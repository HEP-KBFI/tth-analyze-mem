#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

#include "TString.h" // Form()

#include <algorithm> // std::swap(), std::accumulate()
#include <cmath> // std::abs()
#include <exception> // std::runtime_error

using namespace tthMEM;

void
MeasuredEvent_3l1tau::initialize()
{
  met.initialize();

  for(std::size_t i = 0; i < 3; ++i)
    leptons[i].initialize();

  for(std::size_t i = 0; i < 2; ++i)
    jets[i].initialize();

  htau1.initialize();

//--- check whether the sum of lepton charges is +-1
  const int leptonChargeSum = std::accumulate(
    leptons.begin(), leptons.end(), 0,
    [](int sum, const MeasuredLepton & lepton) -> int
    {
      return sum + lepton.charge();
    }
  );
  if(std::abs(leptonChargeSum) != 1)
    throw std::runtime_error("Something's off: the abs of sum of lepton charges is not 1");

//--- retrieve the indices of leptons that have the same charge
  int leptonIdx1 = -1,
      leptonIdx2 = -1;
  for(unsigned i = 0; i < 3; ++i)
    for(unsigned j = i + 1; j < 3; ++j)
      if(leptons[i].charge() == leptons[j].charge())
      {
        leptonIdx1 = i;
        leptonIdx2 = j;
      }
  if(leptonIdx1 < 0 || leptonIdx2 < 0)
    throw std::runtime_error("Something's off: there were no leptons with equal signs");

//--- determine the four index permutations the event can possibly have
  currentPermutation_ = 0;
  leptonPermIdxs = { { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 }, { 0, 1, 2 } };
  jetPermIdxs    = {    { 0, 1 },    { 0, 1 },    { 1, 0 },    { 1, 0 } };
  for(unsigned i = 1; i < 4; i += 2)
    std::swap(leptonPermIdxs[i][leptonIdx1], leptonPermIdxs[i][leptonIdx2]);

//--- set the permutations and their pointers
  leptons.setPermutationPtrs(leptonPermIdxs, &currentPermutation_, 4);
  jets.setPermutationPtrs(jetPermIdxs, &currentPermutation_, 4);

//--- use the first lepton index in the default representation that
//--- has opposite sign w.r.t the tau lepton
  if((leptonChargeSum + htau1.charge()) != 0)
    throw std::runtime_error("The sum of charges of leptonic products is not zero");
  complLeptonIdx = leptonIdx1;
  bjetLeptonIdxs = (complLeptonIdx == 0) ? std::vector<unsigned>{{1, 2}} : (
                   (complLeptonIdx == 1) ? std::vector<unsigned>{{0, 2}} :
                                           std::vector<unsigned>{{1, 2}});
}

void
MeasuredEvent_3l1tau::setBranches(TChain * t)
{
  t -> SetBranchAddress("run",  &run);
  t -> SetBranchAddress("lumi", &lumi);
  t -> SetBranchAddress("evt",  &evt);

  met.setBranches(t);

  for(std::size_t i = 0; i < 3; ++i)
    leptons[i].setBranches(t, Form("lepton%lu", (i + 1)));

  for(std::size_t i = 0; i < 2; ++i)
    jets[i].setBranches(t, Form("jet%lu", (i + 1)));

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

  for(std::size_t i = 0; i < 3; ++i)
    leptons[i].initNewBranches(t, Form("lepton%lu", (i + 1)));

  for(std::size_t i = 0; i < 2; ++i)
    jets[i].initNewBranches(t, Form("jet%lu", (i + 1)));

  htau1.initNewBranches(t, "htau1");

  mvaVariables.initNewBranches(t);
}

bool
MeasuredEvent_3l1tau::hasNextPermutation() const
{
  return currentPermutation_ < 4;
}

void
MeasuredEvent_3l1tau::nextPermutation() const
{
  ++currentPermutation_;
}

void
MeasuredEvent_3l1tau::resetPermutation() const
{
  currentPermutation_ = 0;
}

std::string
MeasuredEvent_3l1tau::str() const
{
  return std::string(Form("%u_%u_%llu_%u", run, lumi, evt, currentPermutation_));
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredEvent_3l1tau & event)
  {
    os << "The event (" << event.run << ":"
                        << event.lumi << ":"
                        << event.evt << "), "
       << "permutation #" << (event.currentPermutation_ + 1) << ":\n";
    for(std::size_t i = 0; i < 3; ++i)
      os << "\tLepton " << (i + 1) << ": " << event.leptons[i] << "\n";
    for(std::size_t i = 0; i < 2; ++i)
      os << "\tJet "    << (i + 1) << ": " << event.jets[i]    << "\n";
    os << "\tTau: " << event.htau1 << "\n";
    os << "\tMET: " << event.met   << "\n";
    return os;
  }
}
