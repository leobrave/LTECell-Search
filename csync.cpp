#include "csync.h"

// User-defined FFT
#define STB_FFT_IMPLEMENTAION
#include "stb_fft.h"

void CSync::gen_pss()
{
    int root_idx[N_PSS] = {25, 29, 34}; // root_idx

    for (int Nid_2 = 0; Nid_2 < N_PSS; Nid_2++)
    {
        int u = root_idx[Nid_2];

        // freq
        for (int i = 0; i < 31; i++)
        {
            float re = cos(-1 * M_PI * u * i * (i + 1) / 63.0);
            float im = sin(-1 * M_PI * u * i * (i + 1) / 63.0);

            pss_f[Nid_2][i] = std::complex<float>(re, im);
        }

        for (int i = 31; i < 62; i++)
        {
            float re = cos(-1 * M_PI * u * (i + 2) * (i + 1) / 63.0);
            float im = sin(-1 * M_PI * u * (i + 2) * (i + 1) / 63.0);

            pss_f[Nid_2][i] = std::complex<float>(re, im);
        }

        // freq --> time
        std::complex<float> tmp[N_PSS_TIME_LEN];
        for (int i = 0; i < N_PSS_TIME_LEN; i++)
        {
            tmp[i] = std::complex<float>(0.0, 0.0);
        }

        for (int i = 1; i < N_PSS_TIME_LEN / 2; i++)
        {
            tmp[i] = pss_f[Nid_2][30 + i];
            tmp[i + 32] = pss_f[Nid_2][i - 1];
        }

        STB_IFFT((cmplx *)tmp, (cmplx *)pss_t[Nid_2], N_PSS_TIME_LEN);
    }
}

void CSync::gen_sss()
{
    int q = 0;
    int m0 = 0;
    int m1 = 0;
    int q_prime = 0;
    int m_prime = 0;

    int x_s_tilda[31] = {0, 0, 0, 0, 1}; // [0]=0, [1]=0, [2]=0, [3]=0, [4]=1
    int s_tilda[31];                     //not initialize

    int x_c_tilda[31] = {0, 0, 0, 0, 1}; // [0]=0, [1]=0, [2]=0, [3]=0, [4]=1
    int c_tilda[31];                     //not initialize

    int x_z_tilda[31] = {0, 0, 0, 0, 1}; // [0]=0, [1]=0, [2]=0, [3]=0, [4]=1
    int z_tilda[31];

    int s0_m0[31];
    int s1_m1[31];
    int s0_m0_idx = 0; //s0_m0[idx]
    int s1_m1_idx = 0; //s1_m1[idx]

    int c0[31];
    int c1[31];
    int c0_idx = 0; //c0[idx]
    int c1_idx = 0; //c1[idx]

    int z1_m0[31];
    int z1_m1[31];
    int z1_m0_idx = 0; //z1_m1[idx]
    int z1_m1_idx = 0; //z1_m1[idx]

    int sss_pci;

    for (int Nid_1 = 0; Nid_1 < N_SSS; Nid_1++)
    {
        for (int Nid_2 = 0; Nid_2 < N_PSS; Nid_2++)
        {
            sss_pci = Nid_1 * 3 + Nid_2;

            // generate m0 and m1
            q_prime = floor((1.0 * Nid_1) / 30); //floor(x)<=x
            q = floor(1.0 * (Nid_1 + q_prime * (q_prime + 1) / 2) / 30);
            m_prime = Nid_1 + q * (q + 1) / 2;
            m0 = int(m_prime) % 31;
            m1 = int(m0 + floor(m_prime / 31) + 1) % 31;

            // generate s_tilda
            for (int i = 0; i < 26; i++)
            {
                x_s_tilda[i + 5] = (x_s_tilda[i + 2] + x_s_tilda[i]) % 2;
            }
            for (int j = 0; j < 31; j++)
            {
                s_tilda[j] = 1 - 2 * x_s_tilda[j];
            }

            // generate c_tilda
            for (int i = 0; i < 26; i++)
            {
                x_c_tilda[i + 5] = (x_c_tilda[i + 3] + x_c_tilda[i]) % 2;
            }
            for (int j = 0; j < 31; j++)
            {
                c_tilda[j] = 1 - 2 * x_c_tilda[j];
            }

            // generate z_tilda
            for (int i = 0; i < 26; i++)
            {
                x_z_tilda[i + 5] = (x_z_tilda[i + 4] + x_z_tilda[i + 2] + x_z_tilda[i + 1] + x_z_tilda[i]) % 2;
            }
            for (int j = 0; j < 31; j++)
            {
                z_tilda[j] = 1 - 2 * x_z_tilda[j];
            }

            // generate s0_m0[31], s1_m1[31]
            for (int n = 0; n < 31; n++)
            {
                s0_m0_idx = (n + m0) % 31;
                s1_m1_idx = (n + m1) % 31;

                s0_m0[n] = s_tilda[s0_m0_idx];
                s1_m1[n] = s_tilda[s1_m1_idx];
            }

            // generate c0,c1
            for (int n = 0; n < 31; n++)
            {
                c0_idx = (n + Nid_2) % 31;
                c1_idx = (n + Nid_2 + 3) % 31;

                c0[n] = c_tilda[c0_idx];
                c1[n] = c_tilda[c1_idx];
            }

            // generate z1_m0, z1_m1
            for (int n = 0; n < 31; n++)
            {
                z1_m0_idx = (n + (m0 % 8)) % 31;
                z1_m1_idx = (n + (m1 % 8)) % 31;

                z1_m0[n] = z_tilda[z1_m0_idx];
                z1_m1[n] = z_tilda[z1_m1_idx];
            }

            // generate SSS
            // SSS_du_0[62], SSS_du_5[62]
            for (int n = 0; n < 31; n++)
            {
                sss_du_0[sss_pci][2 * n] = s0_m0[n] * c0[n];
                sss_du_5[sss_pci][2 * n] = s1_m1[n] * c0[n];

                sss_du_0[sss_pci][2 * n + 1] = s1_m1[n] * c1[n] * z1_m0[n];
                sss_du_5[sss_pci][2 * n + 1] = s0_m0[n] * c1[n] * z1_m1[n];
            }
        }
    }
}

void CSync::find_pss(std::vector<std::complex<float>> &src)
{
    int n_range = src.size() - N_PSS_TIME_LEN;

    for (int Nid_2 = 0; Nid_2 < N_PSS; Nid_2++)
    {
        for (int i = 0; i < n_range; i++)
        {
            std::complex<float> xx(0.0, 0.0);
            std::complex<float> yy(0.0, 0.0);
            std::complex<float> xy(0.0, 0.0);

            for (int j = 0; j < N_PSS_TIME_LEN; j++)
            {
                xx += pss_t[Nid_2][j] * conj(pss_t[Nid_2][j]);
                xy += src[i + j] * conj(pss_t[Nid_2][j]);
                yy += src[i + j] * conj(src[i + j]);
            }
            corr_pss_val[Nid_2][i] += abs(xy); //abs(xy/(yy-xy/xx));
        }
    }

    if (corr_pss_update >= N_PSS_CORR_FRAMES)
    {
        std::cout << "The num of results order is : " << count;
        std::cout << std::endl;
        count++;

        sss_start_flag = 1;

        for (int Nid_2 = 0; Nid_2 < N_PSS; Nid_2++)
        {

            float max_val = -1.0;
            for (int i = 0; i < n_range; i++)
            {
                if (corr_pss_val[Nid_2][i] > max_val)
                {
                    pss_pos[Nid_2] = i;
                    pss_corr[Nid_2] = corr_pss_val[Nid_2][i];

                    max_val = corr_pss_val[Nid_2][i];
                }
            }
        }

        printf("Nid_2\tPos\tCorr\n0\t%d\t%f\n1\t%d\t%f\n2\t%d\t%f\n", pss_pos[0], pss_corr[0], pss_pos[1], pss_corr[1], pss_pos[2], pss_corr[2]);

        // sel
        float max_val = -1.0;
        for (int i = 0; i < N_PSS; i++)
        {
            if (pss_corr[i] > max_val)
            {
                Nid_2_det = i;
                pss_pos_det = pss_pos[i];
                pss_corr_det = pss_corr[i];

                max_val = pss_corr[i];
            }
        }

        if (pss_pos_det > N_SAMPS_10MS / 2)
        {
            pss_pos_det = pss_pos_det - N_SAMPS_10MS / 2;
        }

        printf("Nid_2_det = %d, pss_pos_det = %d, pss_corr_det = %f\n\n", Nid_2_det, pss_pos_det, pss_corr_det);

        // if (count / 10 == 0)
        // {
        //     printf("Nid_2\tPos\tCorr\n0\t%d\t%f\n1\t%d\t%f\n2\t%d\t%f\n", pss_pos[0], pss_corr[0], pss_pos[1], pss_corr[1], pss_pos[2], pss_corr[2]);
        //     printf("Nid_2_det = %d, pss_pos_det = %d, pss_corr_det = %f\n\n", Nid_2_det, pss_pos_det, pss_corr_det);
        // }

        //
        memset(corr_pss_val, 0, N_PSS * (N_SAMPS_10MS - N_PSS_TIME_LEN) * sizeof(float));
        corr_pss_update = 0;
    }
    else
    {
        corr_pss_update++;
        return;
    }
}

void CSync::find_sss(std::vector<std::complex<float>> &src)
{
    int n_pci = 0;

    //ce
    std::complex<float> ce[N_PSS_FREQ_LEN];
    std::complex<float> ce_pss_t[N_FFT_LEN];
    std::complex<float> ce_pss_f[N_PSS_FREQ_LEN];
    std::complex<float> ce_pss_f_tmp[N_FFT_LEN];

    //sss
    std::complex<float> sss_t[N_FFT_LEN]; //std::complex<float> (0.0, 0.0));
    std::complex<float> sss_f_tmp[N_FFT_LEN];

    std::complex<float> sss_f[N_SSS_SEQU_LEN];
    std::complex<float> sss_c[N_SSS_SEQU_LEN];

    if (sss_start_flag == 1)
    {
        int pss_peak = pss_pos_det * 2;
        int sss_peak = pss_peak - N_FFT_LEN * 3 - 10 - 9 * 2; //TDD, CP0 = 10, CP else = 9

        //ce
        for (int i = 0; i < N_FFT_LEN; i++)
        {
            ce_pss_t[i] = src[pss_peak + i - 1];
        }
        STB_FFT((cmplx *)ce_pss_t, (cmplx *)ce_pss_f_tmp, N_FFT_LEN);

        for (int i = 0; i < N_PSS_FREQ_LEN / 2; i++)
        {
            ce_pss_f[i] = ce_pss_f_tmp[127 - 30 + i];
            ce_pss_f[i + 31] = ce_pss_f_tmp[i + 1];
        }
        for (int i = 0; i < N_PSS_FREQ_LEN; i++)
        {
            ce[i] = ce_pss_f[i] * conj(pss_f[Nid_2_det][i]);
        }

        //sss
        for (int i = 0; i < N_FFT_LEN; i++)
        {
            sss_t[i] = src[sss_peak + i - 1];
        }
        STB_FFT((cmplx *)sss_t, (cmplx *)sss_f_tmp, N_FFT_LEN);

        for (int i = 0; i < N_SSS_SEQU_LEN / 2; i++)
        {
            sss_f[i] = sss_f_tmp[127 - 30 + i];
            sss_f[i + 31] = sss_f_tmp[i + 1];
        }

        for (int i = 0; i < N_SSS_SEQU_LEN; i++)
        {
            sss_c[i] = sss_f[i] * conj(ce[i]);
        }

        for (int Nid_1 = 0; Nid_1 < N_SSS; Nid_1++)
        {
            n_pci = Nid_1 * 3 + Nid_2_det;

            std::complex<float> corr_sss_0(0.0, 0.0);
            std::complex<float> corr_sss_5(0.0, 0.0);
            for (int j = 0; j < N_SSS_SEQU_LEN; j++)
            {
                corr_sss_0 += sss_c[j] * std::complex<float>(sss_du_0[n_pci][j], 0.0);
                corr_sss_5 += sss_c[j] * std::complex<float>(sss_du_5[n_pci][j], 0.0);
            }

            corr_sss_val[Nid_1][0] = abs(corr_sss_0);
            corr_sss_val[Nid_1][1] = abs(corr_sss_5);
        }

        //if (corr_sss_update >= N_SSS_CORR_FRAMES)
        //{
        sss0_corr = -1.0;
        sss5_corr = -1.0;

        for (int i = 0; i < N_SSS; i++)
        {
            if (corr_sss_val[i][0] > sss0_corr)
            {
                sss0_corr = corr_sss_val[i][0];
                sss0_pos = i;
            }

            if (corr_sss_val[i][1] > sss5_corr)
            {
                sss5_corr = corr_sss_val[i][1];
                sss5_pos = i;
            }
        }
        printf("SSS\tPos\tCorr\nsss0\t%d\t%f\nsss5\t%d\t%f\n", sss0_pos, sss0_corr, sss5_pos, sss5_corr);

        if (sss0_corr > sss5_corr)
        {
            subframe = 0;
            sss_pos_det = sss0_pos;
            sss_corr_det = sss0_corr;
        }
        else
        {
            subframe = 5;
            sss_pos_det = sss5_pos;
            sss_corr_det = sss5_corr;
        }
        Nid_1_det = sss_pos_det;

        printf("Nid_1_det = %d, sss_pos_det = %d, sss_corr_det = %f\n", Nid_1_det, sss_pos_det, sss_corr_det);
        printf("Cell_ID   = %d, Subframe = %d\n\n", Nid_1_det * 3 + Nid_2_det, subframe);

        sss_start_flag = 0;
        //     corr_sss_update = 0;
        //     memset(corr_sss_val, 0, N_SSS * 2 * sizeof(float));
        // }
        // else
        // {
        //     corr_sss_update++;
        //     return;
        // }
    }
    else
    {
        return;
    }
}