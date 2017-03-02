//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef HelAmps_sm_ttz_t2lvl_tbar2lvl_z2ll_H
#define HelAmps_sm_ttz_t2lvl_tbar2lvl_z2ll_H

#include <cmath>
#include <complex>

namespace MG5_sm_ttz_t2lvl_tbar2lvl_z2ll
{
void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void txxxxx(double p[4], double tmass, int nhel, int nst, std::complex<double>
    fi[18]);

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

double Sgn(double e, double f); 

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]);

void FFV2_3(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> COUP,
    double M3, double W3, std::complex<double> V3[]);
void FFV2_4_3(std::complex<double> F1[], std::complex<double> F2[], std::complex<double>
    COUP1, std::complex<double> COUP2, double M3, double W3, std::complex<double> V3[]);

void FFV1P0_3(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> COUP,
    double M3, double W3, std::complex<double> V3[]);

void FFV3_1(std::complex<double> F2[], std::complex<double> V3[], std::complex<double> COUP,
    double M1, double W1, std::complex<double> F1[]);

void FFV5_2(std::complex<double> F1[], std::complex<double> V3[], std::complex<double> COUP,
    double M2, double W2, std::complex<double> F2[]);

void FFV2_1(std::complex<double> F2[], std::complex<double> V3[], std::complex<double> COUP,
    double M1, double W1, std::complex<double> F1[]);
void FFV2_5_1(std::complex<double> F2[], std::complex<double> V3[], std::complex<double>
    COUP1, std::complex<double> COUP2, double M1, double W1, std::complex<double> F1[]);
void FFV2_3_1(std::complex<double> F2[], std::complex<double> V3[], std::complex<double>
    COUP1, std::complex<double> COUP2, double M1, double W1, std::complex<double> F1[]);

void FFV1_2(std::complex<double> F1[], std::complex<double> V3[], std::complex<double> COUP,
    double M2, double W2, std::complex<double> F2[]);

void FFV2_0(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, std::complex<double> & vertex);
void FFV2_5_0(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP1, std::complex<double> COUP2, std::complex<double> & vertex);

void FFV5_0(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, std::complex<double> & vertex);

void FFV5_1(std::complex<double> F2[], std::complex<double> V3[], std::complex<double> COUP,
    double M1, double W1, std::complex<double> F1[]);

void FFV1_0(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> V3[],
    std::complex<double> COUP, std::complex<double> & vertex);

void FFV3_2(std::complex<double> F1[], std::complex<double> V3[], std::complex<double> COUP,
    double M2, double W2, std::complex<double> F2[]);

void FFV1_1(std::complex<double> F2[], std::complex<double> V3[], std::complex<double> COUP,
    double M1, double W1, std::complex<double> F1[]);

void FFV4_3(std::complex<double> F1[], std::complex<double> F2[], std::complex<double> COUP,
    double M3, double W3, std::complex<double> V3[]);

void VVV1P0_1(std::complex<double> V2[], std::complex<double> V3[], std::complex<double> COUP,
    double M1, double W1, std::complex<double> V1[]);

void FFV2_2(std::complex<double> F1[], std::complex<double> V3[], std::complex<double> COUP,
    double M2, double W2, std::complex<double> F2[]);
void FFV2_3_2(std::complex<double> F1[], std::complex<double> V3[], std::complex<double>
    COUP1, std::complex<double> COUP2, double M2, double W2, std::complex<double> F2[]);
void FFV2_5_2(std::complex<double> F1[], std::complex<double> V3[], std::complex<double>
    COUP1, std::complex<double> COUP2, double M2, double W2, std::complex<double> F2[]);

}  // end namespace MG5_sm_ttz_t2lvl_tbar2lvl_z2ll

#endif  // HelAmps_sm_ttz_t2lvl_tbar2lvl_z2ll_H
