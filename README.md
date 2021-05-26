# LTECell-Search
LTECell-Search is a program compiled by C++ to search for Physical Cell Identifier (PCI) in 4G LTE network. 

# Hardware Required
Performing LTECell-Search requires the following:
* PC with a Linux OS
* An RF-frontend device capable of both Tx & Rx
* Antennas 

For this implementation the following equipment was used:
* PC running Ubuntu-18.04-LTS
* USRP B210 with antennas 

# Build Instructions
## Install RF front-end driver
For the USRP driver UHD, you can find UHD installation tutorial [here](https://kb.ettus.com/Building_and_Installing_the_USRP_Open-Source_Toolchain_(UHD_and_GNU_Radio)_on_Linux). We strongly recommend to build UHD from source code, which is flexible in development and prototyping.

The following are the driver links for commonly used SDRs.
* UHD: https://github.com/EttusResearch/uhd
* SoapySDR: https://github.com/pothosware/SoapySDR
* BladeRF: https://github.com/Nuand/bladeRF

## Download and build LTECell-Search
Download and compile the source code of LTECell-Search.
```
git clone https://github.com/painting1213/LTECell-Search
cd LTECell-Search
mkdir build
cd build
cmake ../
make
./LTECell-Search
```
We set the center frequency to `1895MHz`, one can change this frequency in `main.cpp`. 

Note: LTECell-Search currently only works in LTE TDD mode. If you need to run in FDD mode, you could modify the function of `find_sss` in `csync.cpp`.
