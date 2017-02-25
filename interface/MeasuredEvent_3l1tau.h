#ifndef MEASUREDEVENT_H
#define MEASUREDEVENT_H

#include "tthAnalysis/tthMEM/interface/GeneratorLevelEvent_3l1tau.h" // GeneratorLevelEvent_3l1tau
#include "tthAnalysis/tthMEM/interface/objects/MeasuredMET.h" // MeasuredMET
#include "tthAnalysis/tthMEM/interface/objects/MeasuredJet.h" // MeasuredJet
#include "tthAnalysis/tthMEM/interface/objects/MeasuredHadronicTau.h" // MeasuredHadronicTau
#include "tthAnalysis/tthMEM/interface/objects/MVAVariables.h" // MVAVariables
#include "tthAnalysis/tthMEM/interface/IndexWrapper.h" // IndexWrapper<,>
#include "tthAnalysis/tthMEM/interface/RLESelector.h" // RLESelector<>

#define NOF_RECO_JETS_MEM 3

namespace tthMEM
{
  class DebugPlotter_ttHorZ_3l1tau;

  class
  MeasuredEvent_3l1tau
  {
  public:
    MeasuredEvent_3l1tau() = default;

    Int_t njets;
    TBranch * branch_njets = nullptr;

    mutable RLESelector<> rle;
    MeasuredMET met;
    IndexWrapper<MeasuredLepton, 3> leptons;
    std::array<MeasuredJet, NOF_RECO_JETS_MEM> allJets;
    mutable IndexWrapper<MeasuredJet, 2> jets;
    MeasuredHadronicTau htau;

    MVAVariables mvaVariables;

    unsigned complLeptonIdx;
    std::vector<unsigned> bjetLeptonIdxs;

    mutable DebugPlotter_ttHorZ_3l1tau * debugPlotter = 0;
    std::shared_ptr<GeneratorLevelEvent_3l1tau> generatorLevel;

    void
    initialize();

    void
    setBranches(TTree * t);

    void
    initNewBranches(TTree * t);

    void
    includeGeneratorLevel(bool include);

    bool
    hasNextPermutation() const;

    bool
    isCorrectPermutation() const;

    void
    nextPermutation() const;

    void
    resetPermutation() const;

    unsigned
    getPermutationNumber() const;

    void
    nextJetCombination() const;

    bool
    hasNextJetCombination() const;

    void
    getJetCombination() const;

    void
    resetJetCombination() const;

    unsigned
    getJetCombinationNumber() const;

    void
    printPermutation() const;

    friend std::ostream &
    operator<<(std::ostream & os,
               const MeasuredEvent_3l1tau & event);

  private:
    mutable unsigned currentPermutation_;
    mutable unsigned currentJetCombination_;
    std::vector<IndexWrapper<MeasuredJet, 2>> jetCombinations_;
    std::vector<std::vector<unsigned>> leptonPermIdxs;
    std::vector<std::vector<unsigned>> jetPermIdxs;
    double dR;
  };
}

#endif // MEASUREDEVENT_H
