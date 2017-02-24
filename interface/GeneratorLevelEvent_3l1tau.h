#ifndef GENERATORLEVELEVENT_3L1TAU_H
#define GENERATORLEVELEVENT_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/objects/MeasuredLepton.h" // MeasuredLepton
#include "tthAnalysis/tthMEM/interface/VariableManager_3l1tau.h" // VariableManager_3l1tau

#include <array> // std::array<>
#include <unordered_map> // std::unordered_map<,,>

namespace tthMEM
{
  typedef MeasuredLepton GeneratorParticle;

  /**
   * @brief The generator level class corresponding to ttH/Z events of 3l1tau channel
   *
   * @todo: provide interface to retrieve the integration variables
   */
  class GeneratorLevelEvent_3l1tau
  {
  public:
    GeneratorLevelEvent_3l1tau();

    GeneratorParticle genNuFromLtau;   ///< tau neutrino from tau decaying leptonically
    GeneratorParticle genNuLepFromTau; ///< lepton neutrino from tau decaying leptonically
    GeneratorParticle genLepFromTau;   ///< letpon from tau decaying leptonically

    GeneratorParticle genNuFromHtau; ///< tau neutrino from tau decaying hadronically
    GeneratorParticle genHtau;       ///< hadronic tau

    std::array<GeneratorParticle, 2> genTau;        ///< taus

    std::array<GeneratorParticle, 2> genBQuarkFromTop; ///< b quarks
    std::array<GeneratorParticle, 2> genNuFromTop;     ///< neutrinos from top (W) decay
    std::array<GeneratorParticle, 2> genLepFromTop;    ///< leptons from top (W) decay

    std::array<GeneratorParticle, 2> genWBoson; ///< W bosons
    std::array<GeneratorParticle, 2> genTop;    ///< top quarks
    GeneratorParticle genHorZ;                  ///< Higgs or Z boson
    GeneratorParticle genDiNuFromLtau;          ///< di-neutrino from tau decaying leptonically

    std::size_t leptonicTauDecayIdx; ///< index of genTau corresponding to leptonic tau decay
    std::size_t hadronicTauDecayIdx; ///< index of genTau corresponding to hadronic tau decay

    void
    initialize();

    void
    setBranches(TTree * t);

    void
    initNewBranches(TTree * t);

    /**
     * @brief Calculates the generator level integration variables
     * @param vm The variable manager
     *
     * @todo Add the computation of missing variables
     */
    void
    setIntegrationVariables(VariableManager_3l1tau & vm) const;

    /**
     * @brief Sets the beam axis (needed in the calculation of the integration variables)
     * @param beamAxis The beam axis
     */
    void
    setBeamAxis(const Vector & beamAxis);

    friend std::ostream &
    operator<<(std::ostream & os,
               const GeneratorLevelEvent_3l1tau & event);

  private:
    Vector beamAxis_; ///< beam axis in the lab

    std::unordered_map<Var_3l1tau, double, EnumClassHash> genIntVariables;
  };
}

#endif // GENERATORLEVELEVENT_3L1TAU_H
