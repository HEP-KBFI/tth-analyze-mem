#include "tthAnalysis/tthMEM/interface/BJetTransferFunction.h"
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
      const double m  = constants::eta_m[idx];
      const double m_ = constants::eta_m_[idx];
      const double n  = constants::eta_n[idx];
      const double n_ = constants::eta_n_[idx];
      const double a  = constants::eta_a[idx];
      const double a_ = constants::eta_a_[idx];
      const double b  = constants::eta_b[idx];
      const double b_ = constants::eta_b_[idx];
      const double c  = constants::eta_c[idx];
      const double c_ = constants::eta_c_[idx];

      const double mu = m * bEnergy + n;
      const double mu_ = m_ * bEnergy + n_;
      const double sigma = std::sqrt(
        a * a + b * b * bEnergy + c * c * pow2(bEnergy)
      );
      const double sigma_ = std::sqrt(
        a_ * a_ + b_ * b_ * bEnergy + c_ * c_ * pow2(bEnergy)
      );
      const double pdf = fb * gaussianPDF(bEnergyReco, mu, sigma);
      const double pdf_ = (1. - fb) * gaussianPDF(bEnergyReco, mu_, sigma_);

      return pdf + pdf_;
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
