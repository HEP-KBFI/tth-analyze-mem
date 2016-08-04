#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h"

using namespace tthMEM;

MeasuredEvent_3l1tau::MeasuredEvent_3l1tau()
  : covMET(TMatrixD(2, 2))
{
  covMET[0][0] = 100.0; // in GeV
  covMET[1][0] =   0.0;
  covMET[0][1] =   0.0;
  covMET[1][1] = 100.0;
}

void
MeasuredEvent_3l1tau::initialize()
{
  met.initialize();

  lepton1.initialize();
  lepton2.initialize();
  lepton3.initialize();

  jet1.initialize();
  jet2.initialize();

  htau1.initialize();
}

void
MeasuredEvent_3l1tau::setBranches(TChain * t)
{
  t -> SetBranchAddress("run",  &run);
  t -> SetBranchAddress("lumi", &lumi);
  t -> SetBranchAddress("evt",  &evt);

  met.setBranches(t);

  lepton1.setBranches(t, "lepton1");
  lepton2.setBranches(t, "lepton2");
  lepton3.setBranches(t, "lepton3");

  jet1.setBranches(t, "jet1");
  jet2.setBranches(t, "jet2");

  htau1.setBranches(t, "htau1");

  mvaVariables.setBranches(t);
}

void
MeasuredEvent_3l1tau::initNewBranches(TTree * t)
{
  branch_run  = t -> Branch("run",  &run,  "run/i");
  branch_lumi = t -> Branch("lumi", &lumi, "lumi/i");
  branch_evt  = t -> Branch("evt",  &evt,  "evt/l");

  met.initNewBranches(t);

  lepton1.initNewBranches(t, "lepton1");
  lepton2.initNewBranches(t, "lepton2");
  lepton3.initNewBranches(t, "lepton3");

  jet1.initNewBranches(t, "jet1");
  jet2.initNewBranches(t, "jet2");

  htau1.initNewBranches(t, "htau1");

  mvaVariables.initNewBranches(t);
}
