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

#include <vector> // std::vector<>
#include <ostream> // std::ostream

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
    MeasuredLepton lepton1, lepton2, lepton3;
    MeasuredJet jet1, jet2;
    MeasuredHadronicTau htau1;

    MVAVariables mvaVariables;

    std::vector<MeasuredLepton> leptons_;
    std::vector<MeasuredJet>    jets_;

    void
    initialize();

    void
    setBranches(TChain * t);

    void
    initNewBranches(TTree * t);

    const std::vector<MeasuredLepton> &
    leptons() const;

    const std::vector<MeasuredJet> &
    jets() const;

    friend std::ostream &
    operator<<(std::ostream & os,
               const MeasuredEvent_3l1tau & event);
  };
}

#endif // MEASUREDEVENT_H
