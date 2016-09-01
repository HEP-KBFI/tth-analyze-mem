#ifndef RECOTRUEEVENT_TTHORZ_3L1TAU_H
#define RECOTRUEEVENT_TTHORZ_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // Vector, LorentzVector

#include <vector> // std::vector<>

namespace tthMEM
{
  class RecoTrueEvent_ttHorZ_3l1tau
  {
  public:
    RecoTrueEvent_ttHorZ_3l1tau() = default;

    LorentzVector hTauLepton, lTauLepton; ///< hadronic tau and lepton from tau decay
    LorentzVector nuHtau, nuLtau;         ///< neutrinos from tau decay
    LorentzVector hTau, lTau;             ///< hadronic tau and leptonic tau
    LorentzVector higgsOrZ;               ///< Higgs boson

    LorentzVector lW[2];  ///< leptons coming from W decay
    LorentzVector nuW[2]; ///< neutrinos from W decay
    LorentzVector W[2];   ///< W-bosons from top decay
    LorentzVector b[2];   ///< b-quarks from top decay
    LorentzVector t[2];   ///< top quark itself
    LorentzVector g[2];   ///< gluons

    /**
     * @brief Returns the same event with all 4-momenta but gluon boosted
     * @param boostVector The boost vector in the direction of which the 4-momenta is boosted
     * @return New instance of this class
     */
    RecoTrueEvent_ttHorZ_3l1tau
    boost(Vector boostVector) const;

    /**
     * @brief Returns the sum of all neutrino 4-momenta in the event
     * @return The sum of 4-momenta of neutrinos from tau end top decays
     */
    LorentzVector
    getNuSum() const;

    /**
     * @brief Returns the sum of 4-momenta of Higgs/Z and both top quarks
     * @return The 4-momentum
     */
    LorentzVector
    getTTHorZ() const;

    /**
     * @brief Returns a new vector containing all 4-momenta of the event subject to MadGraph
     * @return The vector of 4-momenta
     */
    std::vector<LorentzVector>
    getForMG() const;
  };
}

#endif // RECOTRUEEVENT_TTHORZ_3L1TAU_H
