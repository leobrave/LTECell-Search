// UHD Includes
#include <uhd/exception.hpp>
#include <uhd/types/tune_request.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <uhd/utils/thread_priority.hpp>

// Boost Includes
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>

// User-defined Includes
#include "MutexFIFO.h"
#include "cfilter.h"
#include "csync.h"

//#include <fstream>

// Use the boost program options (input) namespace
namespace po = boost::program_options;

/***********************************************************************
 * Signal handlers
 **********************************************************************/
static bool stop_signal_called = false;
void sig_int_handler(int){stop_signal_called = true;}

// Setup Structure for RX Block
struct rxBlock {
    unsigned int number;                         //number
    std::vector<std::complex<float>> data;       //data
};

// Prototype - Collect block of RX samples
void collectRXBlock(
        uhd::rx_streamer::sptr &rx_stream,
        size_t total_num_samps, double rate,
        size_t &processedBlocks,
        MutexFIFO<rxBlock> *blockFIFO);

// Prototype - processingFunction
void processingFunction(MutexFIFO<rxBlock> *blockFIFO); //unsigned int pfNum

// Main Routine
int UHD_SAFE_MAIN( int argc, char *argv[] ) {

    // Setup Boost Threads
    uhd::set_thread_priority_safe();

    // Display Starting Message
    std::cout << "Project Name:\tLTE Cell Search"
              << "\nProject Date:\tMarch 16th, 2020"
              << std::endl;

    //variables to be set by po
    std::string args, subdev;
    double seconds_in_future;
    size_t total_num_samps;     // 1.92e6*0.01 = 19200
    size_t threadCount;         // processingFunction  working threads
    double rate, freq, rx_gain; //bw;

    //setup the program options
    po::options_description desc("Allowed options");
    desc.add_options()
            ("args", po::value<std::string>(&args)->default_value(""), "single uhd device address args")
            ("secs", po::value<double>(&seconds_in_future)->default_value(0.01), "number of seconds in the future to receive")
            ("nsamps", po::value<size_t>(&total_num_samps)->default_value(1.92e4), "total number of samples to receive")
            ("rate", po::value<double>(&rate)->default_value(1.92e6), "rate of incoming samples")
            ("freq", po::value<double>(&freq)->default_value(1895.0e6), "RF center frequency in Hz")
            ("rx_gain", po::value<double>(&rx_gain)->default_value(38.0), "RF rx gain in dB")
	        //("bw", po::value<double>(&bw)->default_value(20e6), "analog frontend filter bandwidth in Hz")  //filter bandwidth
            ("threads", po::value<size_t>(&threadCount)->default_value(2), "total worker threads running (excluding main)")
            ("subdev", po::value<std::string>(&subdev), "subdev spec (homogeneous across motherboards)")
            ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    //create a usrp device
    uhd::usrp::multi_usrp::sptr usrp = uhd::usrp::multi_usrp::make(args);

    //always select the subdevice first, the channel mapping affects the other settings
    if (vm.count("subdev")) usrp->set_rx_subdev_spec(subdev);

    //set the rx sample rate (sets across all channels)
    std::cout << boost::format("Desired RX Rate: %04.2f Msps, ") % (rate / 1e6);
    usrp->set_rx_rate(rate);
    std::cout << boost::format("\tActual: %04.2f Msps") % (usrp->get_rx_rate() / 1e6) << std::endl;

    //set the center frequency
    std::cout << boost::format("Desired RX Freq: %04.2f MHz, ") % (freq / 1e6);
    usrp->set_rx_freq(freq);
    std::cout << boost::format("\tActual: %04.2f MHz") % (usrp->get_rx_freq() / 1e6) << std::endl;
    
    //set the rx gain
    std::cout << boost::format("Desired RX Gain: %04.2f dB, ") % (rx_gain);
    usrp->set_rx_gain(rx_gain);
    std::cout << boost::format("\tActual: %04.2f dB") % (usrp->get_rx_gain()) << std::endl;

    //set the IF filter bandwidth
    //std::cout << boost::format("Desired RX Bandwidth: %f MHz...") % (bw/1e6) << std::endl;
    //usrp->set_rx_bandwidth(bw);
    //std::cout << boost::format("Actual RX Bandwidth: %f MHz...") % (usrp->get_rx_bandwidth()/1e6) << std::endl << std::endl;

    // Create a receive streamer (complex floats)
    uhd::stream_args_t stream_args("fc32");
    uhd::rx_streamer::sptr rx_stream = usrp->get_rx_stream(stream_args);

    // Activate CTRL-C Detection
    std::signal(SIGINT, &sig_int_handler);
    std::cout << std::endl << "Press Ctrl-C to stop collecting..." << std::endl;


    // Setup the Stream Command
    uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    stream_cmd.num_samps = total_num_samps;     //接收采样点数
    stream_cmd.stream_now = true;               //现在开始接收
	
	// Send the Stream Command to the USRP
    rx_stream->issue_stream_cmd(stream_cmd);    //配置rx_stream参数

    // Setup FIFO Queue and keep track of processed blocks
    MutexFIFO<rxBlock> *blockFIFO = new MutexFIFO<rxBlock>();
    size_t processedBlocks = 0;

    // Implement 2 copies of processingFunction
    boost::thread_group processingFunctions;
    processingFunctions.create_thread(boost::bind(&processingFunction, blockFIFO));  //add thread
    //processingFunctions.create_thread(boost::bind(&processingFunction, 2, threadCount, total_num_samps, blockFIFO));

    // Keep running until CTRL-C issued
    while (not stop_signal_called)

        // Collect Block of Samples
        collectRXBlock(rx_stream, total_num_samps, rate, processedBlocks, blockFIFO);

    // Requesting the USRP to stop continuous streaming
    std::cout << "\n\nRequesting the USRP to stop continuous mode..." << std::endl;
    rx_stream->issue_stream_cmd(uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS);

    // Join all threads to main
    std::cout << "Requesting worker threads to end..." << std::endl;
    processingFunctions.join_all();

    // Clean our mess
    std::cout << "Flushing remaining blocks in queue..." << std::endl;
    delete blockFIFO;

    // Let the user see any final messages, then exit
    std::cout << std::endl;
    return EXIT_SUCCESS;

}

/* Function:	Collect RX Block
 * Description:	Collect a block of RX samples and append
 *              them to a FIFO queue.
 */

void collectRXBlock(
        uhd::rx_streamer::sptr &rx_stream,
        size_t total_num_samps, double rate,
        size_t &processedBlocks,                 //unsigned long
        MutexFIFO<rxBlock> *blockFIFO) {

    // Setup metadata to collect errors
    uhd::rx_metadata_t md;

    // Setup a buffer to receive the block of samples
    size_t samps_per_buff = rx_stream->get_max_num_samps();

    rxBlock receivedPacket;                       //接受数据包
    receivedPacket.data.resize(total_num_samps);  //resize设置大小

    // Setup collection timeout (factor of sample rate, samples and transmit)
    double timeout = (total_num_samps / rate) + 0.1;

    // Keep track of the accumulated samples
    size_t num_acc_samps = 0;

    while (num_acc_samps < total_num_samps) {

        // Grab a packet from the USRP
        size_t samps_remaining = (num_acc_samps + samps_per_buff < total_num_samps) ? samps_per_buff : (total_num_samps - num_acc_samps);
        size_t num_rx_samps = rx_stream->recv(&receivedPacket.data.at((unsigned int)num_acc_samps), samps_remaining, md, timeout);

        // Did we get our samples?
        switch (md.error_code) {

            case uhd::rx_metadata_t::ERROR_CODE_NONE:
                break;

            case uhd::rx_metadata_t::ERROR_CODE_TIMEOUT:
                std::cout << boost::format("Got timeout before all samples received, possible packet loss...") << std::endl;
                break;

            default:
                std::cout << boost::format("Got error code 0x%x, exiting loop...") % md.error_code << std::endl;
                break;
        }

        // Increment the count of accumulated samples
        num_acc_samps += num_rx_samps;
    }

    // Check if we collected enough samples
    if (num_acc_samps < total_num_samps)

        // We didn't get enough samples for a block, skip it.
        std::cerr << "Receive timeout before all samples received..." << std::endl;  //timeout, 暂停

    else {
        // Set the blockID and push it!
        receivedPacket.number = processedBlocks;

        blockFIFO->push(receivedPacket);
        processedBlocks++;
    }
}

/* Function:	Processing Function (Get the LTE Cell_ID)
 * Description:	From rxBlock, get data of 1 frame ,the number = rate * 10ms, rate = 1.92MHz
 *              For PSS Nid_2, we need to process the data, conv_same + Decimation.
 *              For SSS Nid_1, we can use the data of rxBlock directly.
 */

void processingFunction( MutexFIFO<rxBlock> *blockFIFO)
{
    // int count = 0;
    // std::ofstream outfile;
    // outfile.open("usrp_data.txt", std::ios::out);

    CFilter filt_cx2;
    filt_cx2.init();

    CSync sync;
    sync.init();

    // Keep running until CTRL-C issued
    while (not stop_signal_called) 
    {
        // Grab a block from the queue
        rxBlock topBlock;
        if (true == blockFIFO->pop(&topBlock)) 
        {
            std::vector<std::complex<float>> data_cx2(topBlock.data.size()/2);
            filt_cx2.forward(topBlock.data, data_cx2);

            sync.find_pss(data_cx2);
            sync.find_sss(topBlock.data);
            
            //write to usrp_data.txt, for 10ms * 99
            // if(count < 19)
            // {
            //     for (int i = 0; i < N_SAMPS_10MS * 2; i++)
            //     {
            //         outfile << topBlock.data[i].real() << " " << topBlock.data[i].imag() << std::endl;
            //     }
            //     count ++;
            // }
            // else
            // {
            //     outfile.close();
            // }

            // std::vector<std::complex<float>> data_cx1 = topBlock.data;
            // std::cout << "size " << data_cx1.size() << "capacity " << data_cx1.capacity() << std::endl;

            // std::cout << topBlock.data[0].real() << topBlock.data[0].imag() << std::endl;
            // std::cout << topBlock.data[1].real() << topBlock.data[1].imag() << std::endl;

            // float sum = 0.0;

            // for (i = 0; i < NumSamples_total; i++)
            // {
            //     I_total[i] = topBlock.data[i].real();
            //     Q_total[i] = topBlock.data[i].imag();

            //     sum += (I_total[i] * I_total[i] + Q_total[i] * Q_total[i]);
            // }
            // // printf("Ave %f\n", sum/NumSamples_total);
           // std::cout << Num_rxBlocks << " " << "Cell_ID: " << Cell_ID << std::endl;
        }
    }
}
