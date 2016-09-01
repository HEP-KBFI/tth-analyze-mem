#ifndef TTHMEMRECFUNCTIONS_H
#define TTHMEMRECFUNCTIONS_H

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // LorentzVector

#include <TMatrixDSym.h> // TMatrixDSym

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
     * @param W       4-momentum of the W-boson
     * @param b_en    Energy of the b-quark
     * @param b_p     Magnitude of the momentum of the b-quark
     * @param nuT_en  Neutrino energy originating from the top decay
     * @param lept_en Lepton energy originating from the top decay
     * @param b_Punit Unit vector of b-jet 3-momentum
     * @return Jacobi factor associated w/ the top decay
     */
    double
    tDecayJacobiFactor(const LorentzVector & W,
                       double bQuarkEnergy,
                       double bQuarkP,
                       double nuTopEnergy,
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
     * @param nuTheta  Polar angle of neutrino 3-momentum
     * @param nuPhi    Azimuthal angle of neutrino 3-momentum
     * @param nuEnergy Energy of the neutrino
     * @param nuP      Magnitude of the neutrino 3-momentum
     * @param eX       Unit x-vector of local tau decay system in lab coordinates
     * @param eY       Unit y-vector of local tau decay system in lab coordinates
     * @param eZ       Unit z-vector of local tau decay system in lab coordinates
     * @return Tau neutrino 4-momentum
     *
     * @note nuP and nuEnergy are equal for hadronic tau decays; for leptonic tau
     *       decays the di-neutrino system has non-zero mass in general, thus
     *       these values are not necessarily equal
     */
    LorentzVector
    nuP4(double nuTheta,
         double nuPhi,
         double nuEnergy,
         double nuP,         /* (bind) */
         const Vector & eX,  /* bind */
         const Vector & eY,  /* bind */
         const Vector & eZ); /* bind */

    /**
     * @brief Calculates lepton neutrino energy originating from top decay:
     *        t -> W b -> l nu b
     * @param nuTopPunit   Unit vector of neutrino 3-momentum
     * @param leptonPunit  Unit vector of lepton 3-momentum
     * @param leptonEnergy The lepton energy
     * @return Neutrino energy
     */
    double
    nuTopEnergy(const VectorSpherical & nuTopPunit,
                const Vector & leptonPunit,         /* bind */
                double leptonEnergy);               /* bind */

    /**
     * @brief Calculates the value of MET/hadronic recoil transfer function (TF)
     * @param nuSumX       Sum of x-component of all neutrino momenta
     * @param nuSumY       Sum of y-component of all neutrino momenta
     * @param METx         MET x-component
     * @param METy         MET y-component
     * @param MET_TF_denom Denominator of this TF ( = 1 / sqrt(2 * pi * det(invConvMET)))
     * @param invCovMET    inverse of MET covariance matrix
     * @return Value of the MET/hadronic recoil TF
     */
    double
    MET_TF(double nuSumX,
           double nuSumY,
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
     * @return The energy fraction E_{mother} / E_{vis}
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
  }
}

#endif // TTHMEMRECFUNCTIONS_H
