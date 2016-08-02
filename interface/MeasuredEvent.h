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

namespace tthMEM
{
  class
  MeasuredEvent
  {
  public:
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

    void
    initialize();

    void
    setBranches(TChain * t);

    void
    initNewBranches(TTree * t);
  };
}

#endif // MEASUREDEVENT_H
