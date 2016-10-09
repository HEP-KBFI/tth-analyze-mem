#ifndef GENERATORLEVELEVENT_3L1TAU_H
#define GENERATORLEVELEVENT_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/MeasuredLepton.h" // MeasuredLepton

#include <array> // std::array<>

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
    GeneratorLevelEvent_3l1tau() = default;

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

    void
    initialize();

    void
    setBranches(TChain * t);

    void
    initNewBranches(TTree * t);

  private:
    std::size_t leptonicTauDecayIdx; ///< index of genTau corresponding to leptonic tau decay
    std::size_t hadronicTauDecayIdx; ///< index of genTau corresponding to hadronic tau decay
  };
}

#endif // GENERATORLEVELEVENT_3L1TAU_H
