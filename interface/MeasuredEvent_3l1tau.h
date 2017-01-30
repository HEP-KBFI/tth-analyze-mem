#ifndef MEASUREDEVENT_H
#define MEASUREDEVENT_H

#include <Rtypes.h> // UInt_t, ULong64_t, Float_t, Long64_t
#include <TTree.h> // TTree
#include <TBranch.h> // TBranch

#include "tthAnalysis/tthMEM/interface/MeasuredMET.h" // tthMEM::MeasuredMET
#include "tthAnalysis/tthMEM/interface/MeasuredLepton.h" // tthMEM::MeasuredLepton
#include "tthAnalysis/tthMEM/interface/MeasuredJet.h" // tthMEM::MeasuredJet
#include "tthAnalysis/tthMEM/interface/MeasuredHadronicTau.h" // tthMEM::MeasuredHadronicTau
#include "tthAnalysis/tthMEM/interface/MVAVariables.h" // tthMEM::MVAVariables
#include "tthAnalysis/tthMEM/interface/IndexWrapper.h" // tthMEM::IndexWrapper<,>
#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau
#include "tthAnalysis/tthMEM/interface/GeneratorLevelEvent_3l1tau.h" // GeneratorLevelEvent_3l1tau

#include <ostream> // std::ostream
#include <vector> // std::vector<>
#include <memory> // std::shared_ptr<>
#include <array> // std::array<,>

#define NOF_RECO_JETS_MEM 3

namespace tthMEM
{
  class
  MeasuredEvent_3l1tau
  {
  public:
    MeasuredEvent_3l1tau() = default;

    UInt_t run;
    UInt_t lumi;
    ULong64_t evt;
    Int_t njets;

    TBranch * branch_run   = nullptr;
    TBranch * branch_lumi  = nullptr;
    TBranch * branch_evt   = nullptr;
    TBranch * branch_njets = nullptr;

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

    std::string
    str(bool includePermutation = true) const;

    void
    addFilter(const std::string & rleSelectionFileName);

    bool
    isFiltered() const;

    friend std::ostream &
    operator<<(std::ostream & os,
               const MeasuredEvent_3l1tau & event);

  private:
    mutable unsigned currentPermutation_;
    mutable unsigned currentJetCombination_;
    std::vector<IndexWrapper<MeasuredJet, 2>> jetCombinations_;
    std::vector<std::vector<unsigned>> leptonPermIdxs;
    std::vector<std::vector<unsigned>> jetPermIdxs;
    std::vector<std::string> rleSelection;
    double dR;
  };
}

#endif // MEASUREDEVENT_H
