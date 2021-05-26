#ifndef _COMMON_H
#define _COMMON_H

// Filter related
const int N_FILT_LEN = 11;
const int N_FILT_HIST_LEN = N_FILT_LEN - 1;
const float filter_ds[N_FILT_LEN] = {0.125, 0, -0.25, 0, 0.625, 1.0, 0.625, 0, -0.25, 0, 0.125};

// PSS related
const int N_PSS = 3;
const int N_PSS_FREQ_LEN = 62;
const int N_PSS_TIME_LEN = 64;
const int N_SAMPS_10MS = 19200 / 2;
const int N_PSS_CORR_FRAMES = 10;

// SSS related
const int N_SSS = 168;
const int N_PCI = 504;
const int N_FFT_LEN = 128;
const int N_SSS_SEQU_LEN = 62;
//const int N_SSS_CORR_FRAMES = 10;

#endif