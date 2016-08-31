#include "tthAnalysis/tthMEM/interface/tthMEMrecFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*

#include <cmath> // std::fabs(), std::sqrt()

#include <TVectorD.h> // TVectorD

namespace tthMEM
{
  namespace functions
  {
    double
    nuHtauCosTheta(double nuHtau_en,
                   double hTauEnergy,
                   double hTauMassSquared,
                   double hTauMomentum)
    {
      return (nuHtau_en * hTauEnergy -
               (constants::massTauSquared - hTauMassSquared) / 2.) /
               (hTauMomentum * nuHtau_en);
    }

    double
    nuLeptTauCosTheta(double nuLeptTau_en,
                      double mInvSquared,
                      double nuLeptTau_p,
                      double complLeptEnergy,
                      double complLeptMassSquared,
                      double complLeptMomentum)
    {
      return (nuLeptTau_en * complLeptEnergy -
               (constants::massTauSquared - (complLeptMassSquared + mInvSquared)) / 2.) /
               (complLeptMomentum * nuLeptTau_p);
    }

    double
    bQuarkEnergy(const LorentzVector & W,
                 const Vector & bJetP3Unit,
                 double bJetRecoEnergy)
    {
      const Vector Wp = getVector(W);
      const double a = constants::DeltaFactor / W.e() / constants::massB;
      const double b = W.Beta() * bJetP3Unit.Dot(Wp.unit());
      const double a2 = pow2(a);
      const double b2 = pow2(b);
      const double b_abs = std::fabs(b);
      const double disc = a2 + b2 - 1.;
      const double sep = disc - a2 * b2;

      if((b2 - 1.) >= 1.e-5)
      {
        LOGVRB << "b^2 = " << b2 << " >= 1 => Eb = 0";
        return 0.;
      }
      if(disc < 0.)
      {
        LOGVRB << "a^2 + b^2 - 1 = " << disc << " < 1 => Eb = 0";
        return 0.;
      }

      const double Eb[2] = {
        (a + b_abs * std::sqrt(disc)) / (1. - b2),
        (a - b_abs * std::sqrt(disc)) / (1. - b2)
      };

      const double Eb_p = Eb[0] * constants::massB;
      const double Eb_m = (Eb[1] < 1.0 ? Eb_p : Eb[1]) * constants::massB;

      if(Eb[0] <= 0. || Eb[1] <= 0.)
      {
        LOGWARN << "Eb+ = " << Eb[0] << " <= 0 or Eb- = " << Eb[1] << " <= 0 "
                << "=> Eb = 0";
        LOGWARN << "(for reference: " << lvrap(W) << ")";
        return 0.;
      }

//--- choose physical solution
      if(b > 0. && sep < 0.)
//--- choose solution closest to the resonstructed energy
        return std::fabs(bJetRecoEnergy - Eb_p) < std::fabs(bJetRecoEnergy - Eb_m) ?
               Eb_p : Eb_m;
      else if(b > 0. && sep >= 0.)
        return Eb_p;
      else if(b <= 0. && sep > 0.)
        return Eb_m;
      return 0.;
    }

    double
    tDecayJacobiFactor(const LorentzVector & W,
                       double bQuarkEnergy,
                       double bQuarkP,
                       double nuTopEnergy,
                       double leptonEnergy,
                       const Vector & bQuarkPunit)
    {
      const double bQuarkBeta = bQuarkP / bQuarkEnergy;
      const double invAbsFactor = std::fabs(
        bQuarkPunit.Dot(getVector(W).unit()) * W.Beta() / bQuarkBeta - 1.
      );
      if(invAbsFactor == 0.)
      {
        LOGWARN << "Encountered singularities in the calculation of the"
                << "Jacobi factor of the top decay";
        LOGWARN << "(for reference: " << lvrap("W", W) << ")";
        return 0.;
      }
      return bQuarkEnergy * pow2(nuTopEnergy) /
             (leptonEnergy * W.e() * constants::massWSquared * invAbsFactor);
    }

    double
    MeffSquaredTau2hadrons(double hTauMassSquared)
    {
      const double denom = constants::massTauSquared - hTauMassSquared;
      if(denom <= 0.)
      {
        LOGWARN << "Something's off: hadronic tau mass is greater than "
                << "or equal to the mass of tau lepton";
        return 0.;
      }
      return constants::ttHhadTauPSfactor / denom;
    }

    double
    hadTauPSJacobiFactor(double z,
                         double hTauMassSquared,
                         double hTauInvBeta)
    {
      return MeffSquaredTau2hadrons(hTauMassSquared) /
             (4 * pow6(2 * pi() * z) * hTauInvBeta);
    }

    double
    leptTauPSJacobiFactor(double mInvSquared,
                          double z,
                          double complLeptMassSquared,
                          double complLeptMomentum)
    {
      if(z >= complLeptMassSquared / constants::massTauSquared &&
         z < 1. - mInvSquared / constants::massTauSquared)
      {
        const double complLeptMass = std::sqrt(complLeptMassSquared);
        const double mInv = std::sqrt(mInvSquared);
        const double tau_en =
          (constants::massTauSquared + mInvSquared - complLeptMassSquared) / (2 * mInv);
        const double vis_en = tau_en - mInv;
        if(! (tau_en >= constants::massTau && vis_en >= complLeptMass))
        {
          LOGVRB << "tau energy not greater than or equal to the tau mass "
                 << "(" << tau_en << " < " << constants::massTau << "); or "
                 << "visible energy not greater than or equal to the associated "
                 << "lepton mass (" << vis_en << " < " << complLeptMass << ")";
          return 0.;
        }
    //--- evaluate Iinv in the rest frame of the neutrino pair
        const double Iinv = 2 * constants::GFSquared / pow2(pi()) * mInvSquared *
          (tau_en * vis_en -
           std::sqrt((pow2(tau_en) - constants::massTauSquared) *
                     (pow2(vis_en) - complLeptMassSquared)) / 3.);
        LOGTRC   << "di-neutrino rest frame: tau en = " << tau_en;
        LOGTRC   << "di-neutrino rest frame: vis en = " << vis_en;
        LOGTRC_S << "di-neutrino rest frame: I_inv = "  << Iinv;
        return Iinv / (8 * pow6(2 * pi() * z) * complLeptMomentum);
      }
      else
      {
        LOGVRB << "z = " << z << " not in the physical region "
               << "[" << complLeptMassSquared / constants::massTauSquared << ", "
                      << 1. - mInvSquared / constants::massTauSquared << ")";
        return 0.;
      }
    }

    double
    nuTauEnergy(double z,
                double hlEnergy)
    {
      return hlEnergy * (1. - z) / z;
    }

    LorentzVector
    nuP4(double nuTheta,
         double nuPhi,
         double nuEnergy,
         double nuP,
         const Vector & eX,
         const Vector & eY,
         const Vector & eZ)
    {
      const VectorSpherical nuHtau_loc(nuP, nuTheta, nuPhi);
      const double nu_px = nuHtau_loc.Dot(eX);
      const double nu_py = nuHtau_loc.Dot(eY);
      const double nu_pz = nuHtau_loc.Dot(eZ);
      return LorentzVector(nu_px, nu_py, nu_pz, nuEnergy);
    }

    double
    nuTopEnergy(const VectorSpherical & nuTopPunit,
                const Vector & leptonPunit,
                double leptonEnergy)
    {
      return constants::massWSquared /
             ((1 - leptonPunit.Dot(nuTopPunit)) * leptonEnergy * 2.);
    }

    double
    MET_TF(double nuSumX,
           double nuSumY,
           double METx,
           double METy,
           double MET_TF_denom,
           const TMatrixDSym & invCovMET)
    {
      TVectorD hadRecDiff(2);
      hadRecDiff(0) = METx - nuSumX;
      hadRecDiff(1) = METy - nuSumY;
      const double MET_pull = (invCovMET * hadRecDiff) * hadRecDiff;
      const double MET_TF_value = MET_TF_denom * std::exp(-MET_pull / 2.);

      LOGTRC << "MET_x = "   << METx   << "; MET_y = "   << METy << "; "
             << "nuSum_x = " << nuSumX << "; nuSum_y = " << nuSumY;
      LOGTRC_S << "=> MET_pull = " << MET_pull << " => MET_TF = " << MET_TF_value;

      return MET_TF_value;
    }
  }
}