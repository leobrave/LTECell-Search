#ifndef _CSYNC_H
#define _CSYNC_H

// System
#include <cmath>
#include <complex>
#include <csignal>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <string.h>

#include "common.h"

class CSync
{
public:
    CSync() {}
    ~CSync() {}

    void init()
    {
        memset(corr_pss_val, 0, N_PSS * (N_SAMPS_10MS - N_PSS_TIME_LEN) * sizeof(float));
        //memset(corr_sss_val, 0, N_SSS * 2 * sizeof(float));

        corr_pss_update = 0;
        //corr_sss_update = 0;

        sss_start_flag = 0;
        count = 1;

        gen_pss();
        gen_sss();
    }

    void gen_pss();
    void gen_sss();

    void find_pss(std::vector<std::complex<float>> &src);
    void find_sss(std::vector<std::complex<float>> &src);

private:
    // PSS related
    std::complex<float> pss_f[N_PSS][N_PSS_FREQ_LEN];
    std::complex<float> pss_t[N_PSS][N_PSS_TIME_LEN];
    float corr_pss_val[N_PSS][N_SAMPS_10MS - N_PSS_TIME_LEN];
    int pss_pos[N_PSS];
    float pss_corr[N_PSS];

    int Nid_2_det;
    int pss_pos_det;
    float pss_corr_det;

    // SSS related
    float sss_du_0[N_PCI][N_SSS_SEQU_LEN];
    float sss_du_5[N_PCI][N_SSS_SEQU_LEN];
    float corr_sss_val[N_SSS][2];

    int sss0_pos;
    int sss5_pos;
    float sss0_corr;
    float sss5_corr;

    int Nid_1_det;
    int sss_pos_det;
    float sss_corr_det;
    int subframe;

    // CORR update
    int corr_pss_update; // to sum up 10 frames
    //int corr_sss_update;  // to sum up 10 frames

    // Flag
    int sss_start_flag; //After 10 frames sss_start_flag = 1, begin SSS Sync
    int count;          //Display the result count
};

#endif