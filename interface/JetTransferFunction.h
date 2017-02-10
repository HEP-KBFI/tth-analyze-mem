#ifndef BJETTRANSFERFUNCTION_H
#define BJETTRANSFERFUNCTION_H

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // pi()

namespace tthMEM
{
  namespace constants
  {
    const double q_eta_m[2] = { 1.0,  1.0  }; // GeV^-1
    const double q_eta_a[2] = { 0.0,  0.0  }; // GeV
    const double q_eta_b[2] = { 1.56, 1.52 }; // GeV^0.5
    const double q_eta_c[2] = { 0.0,  0.13 }; // 1

    const double fb = 0.65;                     // 1
    const double b_eta_m[2]  = {  1.0,  0.98 }; // GeV^-1 (1?)
    const double b_eta_m_[2] = {  0.94, 0.87 }; // GeV^-1 (1?)
    const double b_eta_n[2]  = { -3.6, -4.3  }; // GeV
    const double b_eta_n_[2] = { -3.3, +9.1  }; // GeV
    const double b_eta_a[2]  = {  5.7,  6.0  }; // GeV
    const double b_eta_a_[2] = {  6.6,  0.0  }; // GeV
    const double b_eta_b[2]  = {  0.99, 1.9  }; // GeV^0.5
    const double b_eta_b_[2] = {  1.7,  1.1  }; // GeV^0.5
    const double b_eta_c[2]  = {  0.0,  0.0  }; // 1
    const double b_eta_c_[2] = {  0.16, 0.23 }; // 1

    const double invSqrt2Pi = 1. / std::sqrt(2. * pi());
  }

  namespace functions
  {
    double
    deltaFunction(double bEnergy,
                  double bEnergyReco,
                  double bEta);

    double
    bJetTF(double bEnergy,
           double bEnergyReco,
           double bEta);

    double
    qJetTF(double qEnergy,
           double qEnergyReco,
           double qEta);

    inline double
    gaussianPDF(double x,
                double mu,
                double sigma);
  }
}

#endif // BJETTRANSFERFUNCTION_H
