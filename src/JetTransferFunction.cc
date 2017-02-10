#include "tthAnalysis/tthMEM/interface/JetTransferFunction.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // pow2()

#include <cmath> // std::fabs()

namespace tthMEM
{
  namespace functions
  {
    double
    deltaFunction(double __attribute__((unused)) bEnergy,
                  double __attribute__((unused)) bEnergyReco,
                  double __attribute__((unused)) bEta)
    {
      return 1.;
    }

    double
    bJetTF(double bEnergy,
           double bEnergyReco,
           double bEta)
    {
      const double & fb = constants::fb;
      const unsigned idx = std::fabs(bEta) < 1.0 ? 0 : 1;
      const double m  = constants::b_eta_m [idx];
      const double m_ = constants::b_eta_m_[idx];
      const double n  = constants::b_eta_n [idx];
      const double n_ = constants::b_eta_n_[idx];
      const double a  = constants::b_eta_a [idx];
      const double a_ = constants::b_eta_a_[idx];
      const double b  = constants::b_eta_b [idx];
      const double b_ = constants::b_eta_b_[idx];
      const double c  = constants::b_eta_c [idx];
      const double c_ = constants::b_eta_c_[idx];

      const double mu = m * bEnergy + n;
      const double mu_ = m_ * bEnergy + n_;
      const double sigma = std::sqrt(
        pow2(a) + pow2(b) * bEnergy + pow2(c * bEnergy)
      );
      const double sigma_ = std::sqrt(
        pow2(a_) + pow2(b_) * bEnergy + pow2(c_ * bEnergy)
      );
      const double pdf = fb * gaussianPDF(bEnergyReco, mu, sigma);
      const double pdf_ = (1. - fb) * gaussianPDF(bEnergyReco, mu_, sigma_);

      return pdf + pdf_;
    }

    double
    qJetTF(double qEnergy,
           double qEnergyReco,
           double qEta)
    {
      const unsigned idx = std::fabs(qEta) < 1.0 ? 0 : 1;
      const double m = constants::q_eta_m[idx];
      const double a = constants::q_eta_a[idx];
      const double b = constants::q_eta_b[idx];
      const double c = constants::q_eta_c[idx];

      const double mu = m * qEnergy;
      const double sigma = std::sqrt(
        pow2(a) + pow2(b) * qEnergy + pow2(c * qEnergy)
      );
      return gaussianPDF(qEnergyReco, mu, sigma);
    }

    double
    gaussianPDF(double x,
                double mu,
                double sigma)
    {
      const double var = pow2((x - mu) / sigma);
      return constants::invSqrt2Pi / sigma * std::exp(-var / 2.);
    }
  }
}
