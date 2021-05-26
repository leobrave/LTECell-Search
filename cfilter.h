#ifndef _CFILTER_H
#define _CFILTER_H

#include <cmath>
#include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>

#include "common.h"

class CFilter
{    
public:
    CFilter(){}
    ~CFilter(){}

    void init() {
        std::vector<std::complex<float>> tmp(N_FILT_HIST_LEN, std::complex<float>(0.0, 0.0));
        filt_hist.insert(filt_hist.end(), tmp.begin(), tmp.end());
    }

    void forward(std::vector<std::complex<float>> &src, std::vector<std::complex<float>> &dst){
        std::vector<std::complex<float>> hist_src;
        hist_src.insert(hist_src.end(), filt_hist.begin(), filt_hist.end());
        hist_src.insert(hist_src.end(), src.begin(), src.end());

        for(int i = 0; i < src.size(); i+=2){
            std::complex<float> acc(0.0, 0.0);
            for(int j = 0; j < N_FILT_LEN; j++){
                acc += std::complex<float>(filter_ds[j], 0.0) * hist_src[i + N_FILT_LEN - j];  //x(i) * h(n-i)
            }
            dst[i/2] = acc;
        }
        
        filt_hist.assign(src.end() - N_FILT_HIST_LEN, src.end());
    }

private:
    std::vector<std::complex<float>> filt_hist;   // filter history
};

#endif