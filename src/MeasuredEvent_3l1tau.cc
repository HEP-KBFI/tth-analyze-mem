#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

#include <TString.h> // Form()

#include <algorithm> // std::swap(), std::accumulate(), std::find()
#include <cmath> // std::abs()
#include <fstream> // std::ifstream

using namespace tthMEM;

void
MeasuredEvent_3l1tau::initialize()
{
  met.initialize();

  for(std::size_t i = 0; i < 3; ++i)
    leptons[i].initialize();

  for(std::size_t i = 0; i < 2; ++i)
    jets[i].initialize();

  htau.initialize();

  if(generatorLevel) generatorLevel -> initialize();

//--- check whether the sum of lepton charges is +-1
  const int leptonChargeSum = std::accumulate(
    leptons.begin(), leptons.end(), 0,
    [](int sum, const MeasuredLepton & lepton) -> int
    {
      return sum + lepton.charge();
    }
  );
  if(std::abs(leptonChargeSum) != 1)
    throw_line("runtime error")
      << "Something's off: the abs of sum of lepton charges is not 1";

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
    throw_line("runtime error")
      << "Something's off: there were no leptons with equal signs";

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
  if((leptonChargeSum + htau.charge()) != 0)
    throw_line("runtime error")
      << "The sum of charges of leptonic products is not zero";
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
    jets[i].setBranches(t, Form("jets%lu", (i + 1)));

  htau.setBranches(t, "htau");

  mvaVariables.setBranches(t);

  if(generatorLevel) generatorLevel -> setBranches(t);
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

  htau.initNewBranches(t, "htau");

  mvaVariables.initNewBranches(t);

  if(generatorLevel) generatorLevel -> initNewBranches(t);
}

void
MeasuredEvent_3l1tau::includeGeneratorLevel(bool include)
{
  if(include)
    generatorLevel = std::make_shared<typename decltype(generatorLevel)::element_type>();
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

void
MeasuredEvent_3l1tau::addFilter(const std::string & rleSelectionFileName)
{
  if(! rleSelectionFileName.empty())
  {
    LOGINFO << "Adding filter file '" << rleSelectionFileName << '\'';
    std::ifstream rleSelectionFile(rleSelectionFileName);
    for(std::string line; std::getline(rleSelectionFile, line);)
    {
      LOGINFO << "Selecting event = '" << line << '\'';
      rleSelection.push_back(line);
    }
  }
}

bool
MeasuredEvent_3l1tau::isFiltered() const
{
  const std::string rle(Form("%u:%u:%llu", run, lumi, evt));
  return rleSelection.size() &&
         std::find(rleSelection.begin(), rleSelection.end(), rle) == rleSelection.end();
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredEvent_3l1tau & event)
  {
    os << "The event (" << event.run << ':'
                        << event.lumi << ':'
                        << event.evt << "), "
       << "permutation #" << (event.currentPermutation_ + 1) << ":\n";
    for(std::size_t i = 0; i < 3; ++i)
      os << "\tLepton " << (i + 1) << ": " << event.leptons[i] << '\n';
    for(std::size_t i = 0; i < 2; ++i)
      os << "\tJet "    << (i + 1) << ": " << event.jets[i]    << '\n';
    os << "\tTau: " << event.htau << '\n';
    os << "\tMET: " << event.met   << '\n';
    if(event.generatorLevel)
      os << event.generatorLevel;
    return os;
  }
}
