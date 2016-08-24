#ifndef MEASUREDEVENT_H
#define MEASUREDEVENT_H

#include <Rtypes.h> // UInt_t, ULong64_t, Float_t, Long64_t
#include <TTree.h> // TTree
#include <TChain.h> // TChain
#include <TBranch.h> // TBranch

#include "tthAnalysis/tthMEM/interface/MeasuredMET.h" // tthMEM::MeasuredMET
#include "tthAnalysis/tthMEM/interface/MeasuredLepton.h" // tthMEM::MeasuredLepton
#include "tthAnalysis/tthMEM/interface/MeasuredJet.h" // tthMEM::MeasuredJet
#include "tthAnalysis/tthMEM/interface/MeasuredHadronicTau.h" // tthMEM::MeasuredHadronicTau
#include "tthAnalysis/tthMEM/interface/MVAVariables.h" // tthMEM::MVAVariables
#include "tthAnalysis/tthMEM/interface/IndexWrapper.h" // tthMEM::IndexWrapper<,>
#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau

#include <ostream> // std::ostream
#include <vector> // std::vector<>

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

    TBranch * branch_run = 0;
    TBranch * branch_lumi = 0;
    TBranch * branch_evt = 0;

    MeasuredMET met;
    IndexWrapper<MeasuredLepton, 3> leptons;
    IndexWrapper<MeasuredJet, 2> jets;
    MeasuredHadronicTau htau1;

    MVAVariables mvaVariables;

    unsigned complLeptonIdx;
    std::vector<unsigned> bjetLeptonIdxs;

    mutable DebugPlotter_ttHorZ_3l1tau * debugPlotter_ = 0;

    void
    initialize();

    void
    setBranches(TChain * t);

    void
    initNewBranches(TTree * t);

    bool
    hasNextPermutation() const;

    void
    nextPermutation() const;

    void
    resetPermutation() const;

    std::string
    str() const;

    friend std::ostream &
    operator<<(std::ostream & os,
               const MeasuredEvent_3l1tau & event);

  private:
    mutable unsigned currentPermutation_;
    std::vector<std::vector<unsigned>> leptonPermIdxs;
    std::vector<std::vector<unsigned>> jetPermIdxs;
  };
}

#endif // MEASUREDEVENT_H
