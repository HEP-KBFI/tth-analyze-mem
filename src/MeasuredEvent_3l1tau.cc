#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h"

using namespace tthMEM;

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

  leptons_.clear();
  leptons_.reserve(3);
  leptons_ = { lepton1, lepton2, lepton3 };

  jets_.clear();
  jets_.reserve(2);
  jets_ = { jet1, jet2 };
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

const std::vector<MeasuredLepton> &
MeasuredEvent_3l1tau::leptons() const
{
  return leptons_;
}

const std::vector<MeasuredJet> &
MeasuredEvent_3l1tau::jets() const
{
  return jets_;
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredEvent_3l1tau & event)
  {
    os << "The event "
       << "(" << event.run << ":" << event.lumi << ":" << event.evt << "):\n";
    for(std::size_t i = 0; i < 3; ++i)
      os << "\tLepton " << (i + 1) << event.leptons_[i] << "\n";
    for(std::size_t i = 0; i < 2; ++i)
      os << "\tJet "    << (i + 1) << event.jets_[i]    << "\n";
    os << "\tTau: " << event.htau1 << "\n";
    os << "\tMET: " << event.met   << "\n";
    return os;
  }
}
