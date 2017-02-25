#ifndef TTHMEMRECFUNCTIONS_H
#define TTHMEMRECFUNCTIONS_H

#include "tthAnalysis/tthMEM/interface/general/auxFunctions.h" // LorentzVector

#include <TMatrixDSym.h> // TMatrixDSym
#include <TMatrixD.h> // TMatrixD

namespace tthMEM
{
  namespace functions
  {
    /**
     * @brief Functions for reconstructing the underlying event
     *        from reconstructed values in the event or reconstruction variables
     *        from generator level values.
     */

      /**
     * @brief Calculates energy fraction of leptonic tau decay products
     * @param measuredVisMassSquared Measured mass squared of visible tau decay products
     * @param massHiggsOrZsquared    Mass of Higgs or Z
     * @param z1                     Energy fraction of hadronic tau decay products
     * @return The energy fraction of leptonic tau decay products
     *
     * Note that we could use the energy fraction of hadronic tau decay products
     * in the same manner, given leptonic energy fraction. But in the code we use it
     * in the former way.
     */
    double
    z2(double z1,
       double measuredVisMassSquared, /* bind */
       double massHiggsOrZsquared);   /* bind */

    /**
     * @brief A more precise computation of energy fraction in leptonic tau decay system
     * @param z1                  Energy fraction in hadronic tau decay system
     * @param mInvSquared         Squared invisible mass of the neutrino pair in leptonic tau decay system
     * @param nuLtau_phi          Rotation angle of the neutrino pair in leptonic tau decay system
     * @param nuHtau              4-momentum of the neutrino from hadronic tau decay system
     * @param vis                 4-momentum of visible decay products of the tau pair
     * @param complLepton         4-momentum of the lepton in leptonic tau decay system
     * @param nuLtauLocalSystem   Coordinate transformation matrix in leptonic tau decay system
     * @param massHiggsOrZsquared Squared mass of the mother of the tau pair (either Higgs or Z)
     * @return A more precise value for the leptonic energy fraction
     *
     * @note Simple mH^2 = mVis^2 / (z1 * z2) doesn't hold very well; the solution is to expand the mH^2
     *       term in terms of the integration variables and measured quantities. The result is way more complex
     *       than the initial guess (which is pretty good but not good enough since Higgs mass width is so
     *       narrow and even a slight deviation from the PDG mass will render the MEM useless). In principle,
     *       a new set of integration variables should be developed, but this is too much work; hence the hack.
     *
     *       The solution basically boils down to solving a quadratic equation in terms of z2; however, there's
     *       an implicit dependency in the opening angle of the leptonic tau decay system, which also depends on
     *       the energy fraction z2. The approximation we employ is that we use first estimation of z2 in
     *       the computation of this opening angle, which we in turn use to recompute z2.
     *
     *       Since we're dealing with a quadratic equation, we must consider both solutions. If both solutions
     *       are physical, closest value to the initial estimate of z2 is considered. Since this is still
     *       an approximation and not a perfect estimate for z2, there still might be bias in the calculation
     *       (and it looks like it happens more often if the opening angle is relatively large).
     *
     *       An imporovement of this function would calculate all everything up to upper-most parent particle
     *       and would choose the one which is closest to the expected mass of the parent.
     */
    double
    z2(double z1,
       double mInvSquared,
       double nuLtau_phi,
       const LorentzVector & nuHtau,
       const LorentzVector & vis,          /* bind */
       const LorentzVector & complLepton,  /* bind */
       const TMatrixD & nuLtauLocalSystem, /* bind */
       double massHiggsOrZsquared);        /* bind */

    /**
     * @brief Reconstructs cosine of the opening angle between the hadronic tau
     *         decay products (a hadronic tau lepton and a neutrino): tau -> htau nu
     * @param nuHtau_en       Energy of the neutrino (equal to its momentum)
     * @param hTauEnergy      Energy of the hadronic tau lepton (i.e. the decay product)
     * @param hTauMassSquared Squared mass of the hadronic tau lepton
     * @param hTauMomentum    Magnitude of the momentum of the hadronic tau lepton
     * @return Cosine of the opening angle between the decay hadronic tau decay products
     *
     * @note is a special case of nuLeptTauCosTheta() if
     *        - nuLeptTau_en = nuLeptTau_p
     *        - complLeptMassSquared = 0
     *        - complLeptMomentum == complLeptEnergy
     */
    double
    nuHtauCosTheta(double nuHtau_en,
                   double hTauEnergy,      /* bind */
                   double hTauMassSquared, /* bind */
                   double hTauMomentum);   /* bind */

    /**
     * @brief Reconstructs cosine of the opening angle between the leptonic tau
     *        decay products (a lepton and a neutrino pair w/ non-zero mass): tau -> ltau nunu
     * @param nuLeptTau_en         Energy of the neutrino pair
     * @param mInvSquared          Invisible squared mass of the neutrino pair
     * @param nuLeptTau_p          Magnitude of the momentum of the neutrino pair
     * @param complLeptEnergy      Energy of the complementary lepton
     * @param complLeptMassSquared Mass of the complementary lepton
     * @param complLeptMomentum    Magnitude of the momentum of the complementary lepton
     * @return Cosine of the opening angle between the decay leptonic tau decay products
     */
    double
    nuLeptTauCosTheta(double nuLeptTau_en,
                      double mInvSquared,
                      double nuLeptTau_p,
                      double complLeptEnergy,      /* bind */
                      double complLeptMassSquared, /* bind */
                      double complLeptMomentum);   /* bind */

    /**
     * @brief Calculates the b-quark energy from the top decay constituents: t -> W b
     * @param W              4-momentum of the W boson
     * @param bJetP3Unit     Unit vector of the b-jet's 3-momentum
     * @param bJetRecoEnergy Energy of the b-jet
     * @return Energy of the b-quark
     *
     * @note The method imposes the following constraints in the calculation
     *         - b-quark mass is constants::bMass
     *         - t-quark mass is constants::tMass
     *         - direction of the b-quark is perfectly measured (same as associated jet)
     */
    double
    bQuarkEnergy(const LorentzVector & W,
                 const Vector & bJetP3Unit, /* bind */
                 double bJetRecoEnergy);    /* bind */

    /**
     * @brief Calculates Jacobi factor associated w/ the top decay: t -> W b -> l nu b
     * @param W            4-momentum of the W-boson
     * @param bQuarkEnergy Energy of the b-quark
     * @param bQuarkP      Magnitude of the momentum of the b-quark
     * @param nuWEnergy    Neutrino energy originating from the W decay
     * @param leptonEnergy Lepton energy originating from the top decay
     * @param bQuarkPunit  Unit vector of b-jet 3-momentum
     * @return Jacobi factor associated w/ the top decay
     */
    double
    tDecayJacobiFactor(const LorentzVector & W,
                       double bQuarkEnergy,
                       double bQuarkP,
                       double nuWEnergy,
                       double leptonEnergy,         /* bind */
                       const Vector & bQuarkPunit); /* bind */

    /**
     * @brief Calculates the squared matrix elements for hadronic tau decay: tau -> htau nu
     * @param hTauMassSquared Squared mass of the hadronic tau lepton
     * @return The squared matrix element
     *
     * @note The matrix element is calculated in such way that it reproduces
     *       the branching ratio of taus decaying into hadronic particles
     */
    double
    MeffSquaredTau2hadrons(double hTauMassSquared);

    /**
     * @brief Calculates the Jacobi and phase space factor which arises from
     *        the variable changes in the delta function associated with
     *        hadronic tau decays
     * @param z               Energy fraction carried by hadronic tau decay product
     * @param hTauMassSquared Squared mass of the hadronic tau lepton
     * @param hTauInvBeta     Inverse velocity of the hadronic tau lepton
     * @return The phase space x Jacobi factor
     */
    double
    hadTauPSJacobiFactor(double z,
                         double hTauMassSquared, /* bind */
                         double hTauInvBeta);    /* bind */

    /**
     * @brief Calculates the Jacobi and phase space factor which arises from
     *        the variable changes in the delta function associated with
     *        leptonic tau decays
     * @param mInvSquared          Squared invisible mass of the neutrino pair
     * @param z                    Energy fraction carried by lepton w.r.t tau
     * @param complLeptMassSquared Squared mass of the complementary lepton
     * @param complLeptMomentum    Magnitude of the momentum of the complementary lepton
     * @return The phase space x Jacobi factor
     */
    double
    leptTauPSJacobiFactor(double mInvSquared,
                          double z,
                          double complLeptMassSquared, /* bind */
                          double complLeptMomentum);   /* bind */

    /**
     * @brief Calculates neutrino or its system energy from from the energy fraction
     *        carried by the tau lepton decay product (hadronic tau or lepton)
     *        to the tau lepton, z
     * @param z            The energy fraction
     * @param hlEnergy     The lepton or hadronic tau lepton energy
     * @return Neutrino or its system energy
     */
    double
    nuTauEnergy(double z,
                double hlEnergy); /* bind */

    /**
     * @brief Finds the neutrino 4-momentum originating from tau decay
     * @param nuTheta          Polar angle of neutrino 3-momentum
     * @param nuPhi            Azimuthal angle of neutrino 3-momentum
     * @param nuEnergy         Energy of the neutrino
     * @param nuP              Magnitude of the neutrino 3-momentum
     * @param nuTauLocalSystem Local tau decay system in lab coordinates
     * @return Tau neutrino 4-momentum
     *
     * @note nuP and nuEnergy are equal for hadronic tau decays; for leptonic tau
     *       decays the di-neutrino system has non-zero mass in general, thus
     *       these values are not necessarily equal. Also, the nuTauLocalSystem
     *       matrix should be constructed so that the 1st/2nd/3rd row corresponds to
     *       the x/y/z-coordinate axis in the local neutrino system.
     */
    LorentzVector
    nuP4(double nuTheta,
         double nuPhi,
         double nuEnergy,
         double nuP,                         /* (bind) */
         const TMatrixD & nuTauLocalSystem); /* bind */

    /**
     * @brief Computes 4-momentum of the di-neutrino system
     * @param z2                Energy fraction in leptonic tau decay system
     * @param mInvSquared       Squared mass of the neutrino pair
     * @param nuLtau_phi        Rotation angle of the neutrino system
     * @param complLepton       Complementary lepton in leptonic tau decay system
     * @param nuLtauLocalsystem Coordinate transformation matrix of the leptonic
     *                          tau decay system
     * @return Lorentz 4-momentum of the neutrino pair
     */
    LorentzVector
    nuLtau(double z2,
           double mInvSquared,
           double nuLtau_phi,
           const LorentzVector & complLepton,   /* bind */
           const TMatrixD & nuLtauLocalsystem); /* bind */

    /**
     * @brief Calculates lepton neutrino energy originating from W decay:
     *        W -> l nu
     * @param nuWPunit   Unit vector of neutrino 3-momentum
     * @param leptonPunit  Unit vector of lepton 3-momentum
     * @param leptonEnergy The lepton energy
     * @return Neutrino energy
     */
    double
    nuWEnergy(const VectorSpherical & nuWPunit,
              const Vector & leptonPunit,       /* bind */
              double leptonEnergy);             /* bind */

    /**
     * @brief Calculates the value of MET/hadronic recoil transfer function (TF)
     * @param METx_        ,,Reconstructed true'' MET x-component
     * @param METy_        ,,Reconstructed true'' MET y-component
     * @param METx         MET x-component
     * @param METy         MET y-component
     * @param MET_TF_denom Denominator of this TF ( = 1 / sqrt(2 * pi * det(invConvMET)))
     * @param invCovMET    inverse of MET covariance matrix
     * @return Value of the MET/hadronic recoil TF
     */
    double
    MET_TF(double METx_,
           double METy_,
           double METx,                    /* bind */
           double METy,                    /* bind */
           double MET_TF_denom,            /* bind */
           const TMatrixDSym & invCovMET); /* bind */

    /**
     * @brief Calculates rotation angle between visible and invisible tau decay products
     *        in the lab frame
     * @param mother The tau 4-momentum which decays into visible and invisible products
     * @param vis      4-momentum of the visible tau decay product
     * @param beamAxis Collision axis
     * @return The rotation angle in the lab frame
     */
    double
    phiFromLabMomenta(const LorentzVector & mother,
                      const LorentzVector & vis,
                      const Vector & beamAxis);

    /**
     * @brief Calculates energy fraction between its mother particle and its visible product
     * @param mother 4-momentum of the mother particle
     * @param vis    4-momentum of its visible decay product
     * @return The energy fraction E_{vis} / E_{mother}
     */
    double
    z(const LorentzVector & mother,
      const LorentzVector & vis);

    /**
     * @brief Calculates cosine of the angle between the two decay products
     * @param vis 4-momentum of the visible decay product
     * @param inv 4-momentum of the invisible decay product
     * @return The cosine of the angle between the decay products
     */
    double
    cosTheta(const LorentzVector & vis,
             const LorentzVector & inv);

    /**
     * @brief Calculate local neutrino coordinate system arising from tau decay
     * @param beamAxis   3-vector defining the beam axis
     * @param leptP3unit Unit 3-momentum of associated lepton in the tau decay
     * @param str        Optional string for printout
     * @return Matrix where 0,1,2 rows correspond to the local x,y,z coordinate axis
     */
    TMatrixD
    nuLocalSystem(const Vector & beamAxis,
                  const Vector & leptP3unit,
                  const std::string & str = "");
  }
}

#endif // TTHMEMRECFUNCTIONS_H
