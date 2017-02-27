#include "tthAnalysis/tthMEM/interface/JetTransferFunction.h"

#include <algorithm> // std::any_of()

namespace tthMEM
{
  namespace structs
  {
    qJetParams::qJetParams(double energy,
                           double eta)
    {
      classify(energy, eta);
    }

    void
    qJetParams::classify(double energy,
                         double eta)
    {
      const unsigned idx = std::fabs(eta) < 1.0 ? 0 : 1;
      const double m = constants::q_eta_m[idx];
      const double a = constants::q_eta_a[idx];
      const double b = constants::q_eta_b[idx];
      const double c = constants::q_eta_c[idx];

      mu    = m * energy;
      sigma = std::sqrt(pow2(a) + pow2(b) * energy + pow2(c * energy));
    }

    bJetParams::bJetParams(double energy,
                           double eta)
    {
      classify(energy, eta);
    }

    void
    bJetParams::classify(double energy,
                         double eta)
    {
      const unsigned idx = std::fabs(eta) < 1.0 ? 0 : 1;
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

      first.mu  = m  * energy + n;
      second.mu = m_ * energy + n_;
      first.sigma  = std::sqrt(pow2(a)  + pow2(b)  * energy + pow2(c  * energy));
      second.sigma = std::sqrt(pow2(a_) + pow2(b_) * energy + pow2(c_ * energy));
    }
  }

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
      using namespace tthMEM::structs;
      const double & fb = constants::fb;
      const bJetParams b(bEnergy, bEta);

      const double pdf =        fb  * gaussianPDF(bEnergyReco, b.first.mu,  b.first.sigma );
      const double pdf_ = (1. - fb) * gaussianPDF(bEnergyReco, b.second.mu, b.second.sigma);

      return pdf + pdf_;
    }

    double
    qJetTF(double qEnergy,
           double qEnergyReco,
           double qEta)
    {
      using namespace tthMEM::structs;
      const qJetParams q(qEnergy, qEta);
      return gaussianPDF(qEnergyReco, q.mu, q.sigma);
    }

    double
    gaussianPDF(double x,
                double mu,
                double sigma)
    {
      const double var = pow2((x - mu) / sigma);
      return constants::invSqrt2Pi / sigma * std::exp(-var / 2.);
    }

    double
    qJetAcceptance(const MeasuredJet                 & testQuark,
                   const std::vector<MeasuredLepton> & leptons,
                   const std::vector<MeasuredJet>    & acceptedQuarks,
                   int                               & quarkShiftIdx,
                   const AnalysisCuts                & cuts)
    {
      quarkShiftIdx = -1;

      // if the reconstructed quark is outside of the detector acceptance
      if(testQuark.absEta() > cuts.jetEta)
      {
        return 1.;
      }

      // if the reconstructed quark coincides with an accepted quark
      std::vector<double> dRs;
      std::transform(
        acceptedQuarks.begin(), acceptedQuarks.end(), dRs.begin(),
        [&testQuark](const MeasuredJet & quark) -> double
        {
          return quark.dR(testQuark);
        }
      );
      const int minIdx = std::distance(dRs.begin(), std::min_element(dRs.begin(), dRs.end()));
      if(dRs[minIdx] < cuts.jetAlgoRadius)
      {
        // the energy of an accepted quark must be shifted in its TF
        quarkShiftIdx = minIdx;
        return 1.;
      }

      // check maybe the quark is vetoed by a lepton
      if(std::any_of(
           leptons.begin(), leptons.end(),
           [&cuts, testQuark](const MeasuredLepton & lepton) -> bool
           {
             return testQuark.dR(lepton) < cuts.jetToLepton_dR &&
                    testQuark.energy() / lepton.energy() > cuts.jetToLepton_relIso;
           }
        ))
      {
        return 0.;
      }

      // otherwise use probabilistic value for the quark acceptance
      using namespace tthMEM::structs;
      const qJetParams q(testQuark.energy(), testQuark.eta());
      const double denom = 1. / q.sigma * std::sqrt(2.);
      const double alpha = q.mu * denom;
      const double x     = std::exp(-testQuark.eta());
      return (std::erf(alpha) - std::erf(alpha - cuts.jetPt * std::sqrt(pow2(x) + 1.) * denom / x)) / 2.;
    }
  }
}

