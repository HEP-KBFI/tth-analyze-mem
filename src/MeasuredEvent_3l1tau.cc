#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau

using namespace tthMEM;

void
MeasuredEvent_3l1tau::initialize()
{
//--- remove any references to the previous permutations
  leptons.reset();
  jets.reset();

  met.initialize();

  for(int i = 0; i < 3; ++i)
    leptons[i].initialize();

//--- let's read in all jets and print them
  LOGTRC << "All jets";
  for(std::size_t i = 0; i < NOF_RECO_JETS_MEM; ++i)
  {
    allJets[i].initialize();
    LOGTRC << "All jet " << (i + 1) << ": " << allJets[i];
  }

//--- let's build all combinations of 2s out of N valid jets
//--- (the multicplicity varies from event to event)
//--- also, print the combinations just for a quick check
  LOGTRC << "Found " << njets << " jets in total";
  const std::vector<std::array<unsigned, 2>> jetCombinationIndexes = tthMEM::combination<2>(njets);
  if(! jetCombinationIndexes.size())
    throw_line_ext("MeasuredEvent_3l1tau", TTHEXCEPTION_ERR_CODE_MISSING_JET_COMBINATIONS)
      << "No combinations built";
  jetCombinations_.clear();
  for(std::size_t i = 0; i < jetCombinationIndexes.size(); ++i)
  {
    const std::array<unsigned, 2> & iPair = jetCombinationIndexes[i];
    IndexWrapper<MeasuredJet, 2> jetCombination;
    jetCombination[0] = allJets[iPair[0]];
    jetCombination[1] = allJets[iPair[1]];
    jetCombinations_.push_back(jetCombination);

    LOGTRC << "Jet combination #" << (i + 1);
    LOGTRC << "jet 1: " << jetCombination[0];
    LOGTRC << "jet 2: " << jetCombination[1];
  }
  LOGTRC << "Got " << jetCombinations_.size() << " jet pair combinations in total";
  currentJetCombination_ = 0;
//--- select the first combination
  jets = jetCombinations_[currentJetCombination_];

  htau.initialize();

  if(generatorLevel)
  {
    LOGTRC << "Redirecting reconstructed particles to their underlying "
           << "generator level particles";
    generatorLevel -> initialize();
    for(int i = 0; i < 2; ++i)
      leptons[i] = generatorLevel -> genLepFromTop[i];
    leptons[2] = generatorLevel -> genLepFromTau;
    for(int i = 0; i < 2; ++i)
    {
      const GeneratorParticle & jet = generatorLevel -> genBQuarkFromTop[i];
//--- the first combination must include a permutation which passes gen lvl matching
      jetCombinations_[0][i] = MeasuredJet(jet.pt(), jet.eta(), jet.phi(), jet.mass());
    }
//--- remove the rest permutations
    if(jetCombinations_.size() > 1)
    {
      LOGTRC << "Erasing rest of the jet combinations";
      jetCombinations_.erase(jetCombinations_.begin() + 1, jetCombinations_.end());
      LOGTRC << "Forcing only " << jetCombinations_.size() << " jet combination";
      LOGTRC << "Forced gen b 1: " << jetCombinations_[0][0];
      LOGTRC << "Forced gen b 2: " << jetCombinations_[0][1];
    }
    const GeneratorParticle & gHtau = generatorLevel -> genHtau;
    const int htauDecayMode = htau.decayMode(); // use reco decayMode for now
    htau = MeasuredHadronicTau(gHtau.pt(), gHtau.eta(), gHtau.phi(), gHtau.mass(),
                               gHtau.charge(), htauDecayMode);
  }

//--- check whether the sum of lepton charges is +-1
  const int leptonChargeSum = std::accumulate(
    leptons.begin(), leptons.end(), 0,
    [](int sum, const MeasuredLepton & lepton) -> int
    {
      return sum + lepton.charge();
    }
  );
  if(std::abs(leptonChargeSum) != 1)
    throw_line_ext("runtime error", TTHEXCEPTION_ERR_CODE_NONZERO_CHARGE_SUM)
      << "Something's off: the abs of sum of lepton charges is not 1";

//--- Determine the four index permutations the event can possibly have.
//--- The idea here is that we want to retain some kind of consistency in
//--- the signs of the leptons: we require that the first lepton coming from
//--- top decay has always positive sign. If we do not require that, we might
//--- never hit the right kinematic configuration at all. For example, if
//--- hadronic tau has a negative sign, the complementary lepton must have
//--- a positive sign. This leaves a positive-negative lepton pair which must
//--- be associated with the jets. However, we want to permute only the jets
//--- and same-sign leptons. This means that we fix one lepton and permute
//--- the other two. However, if the fixed lepton is associated with a jet
//--- that should not actually be together, none of the permutations would
//--- be a valid kinematic configuration.
//---
//--- The solution is to keep the sign of the first lepton coming from top decay
//--- always positive (we permute the jets anyways, so the association will always
//--- be right). The additional advantage here is that we have to check only two
//--- possible cases in filling the MadGraph matrix element: which signs the tau
//--- leptons have (the first tau lepton must have a positive sign).
//---
//--- The scenario is thus the following:
//---   * rotate the lepton signs until all the positive signs are pushed to the left:
//---       + + -  =>  + + -                 - - + => + - -
//---       + - +  =>  + + -       and       - + - => + - -
//---       - + +  =>  + + -                 + - - => + - -
//---   * if the sum of the signs is +1, we need to permute the first two leptons,
//---     from which follows that the complementary lepton index is either 0 or 1,
//---     and that we have to permute the first two leptons
//---   * if the sum of the signs is -1, we need to permute the last two leptons,
//---     from which follows that the complementary lepton index is either 1 or 2,
//---     and that we have to permute the last two leptons
//---   * since the complementary lepton index must be consitent across one integration,
//---     we choose it to be 1 b/c it would be valid for both cases (sum is +1 or -1)
//---   * therefore the lepton indexes associated with the top decay are 0 and 2
//---   * the permutation pattern for the leptons is thus
//---       (AB) (BA) (AB) (BA)
//---     and for the jets
//---       (CD) (CD) (DC) (DC)
//---     which ensures that all possible lepton-to-jet combinations are considered.
  currentPermutation_ = 0;
  complLeptonIdx = 1;
  bjetLeptonIdxs = std::vector<unsigned>{{0, 2}};
  std::vector<unsigned> permutationUnit{0, 1, 2};
  if(leptonChargeSum == -1)
//--- rotate until we have: + - -
  {
    while(leptons[permutationUnit[0]].charge() != +1)
      std::rotate(
        permutationUnit.begin(), permutationUnit.begin() + 1, permutationUnit.end()
      );
  }
  else // sum of lepton charges is +1
//--- rotate until we have: + + -
  {
    while(leptons[permutationUnit[2]].charge() != -1)
      std::rotate(
        permutationUnit.begin(), permutationUnit.begin() + 1, permutationUnit.end()
      );
  }
  leptonPermIdxs = std::vector<std::vector<unsigned>>(4, permutationUnit );
  jetPermIdxs = std::vector<std::vector<unsigned>>(4, { 0, 1 });
  for(unsigned i = 0; i < 4; ++i)
  {
    if(i % 2 == 0)
    {
      if(leptonChargeSum == -1)
        std::swap(leptonPermIdxs[i][1], leptonPermIdxs[i][2]);
      else
        std::swap(leptonPermIdxs[i][0], leptonPermIdxs[i][1]);
    }
    if(i < 2)
      std::swap(jetPermIdxs[i][0], jetPermIdxs[i][1]);
  }

//--- set the permutations and their pointers
  leptons.setPermutationPtrs(leptonPermIdxs, &currentPermutation_, 4);
  for(IndexWrapper<MeasuredJet, 2> & jetCombination: jetCombinations_)
    jetCombination.setPermutationPtrs(jetPermIdxs, &currentPermutation_, 4);

//--- set dR for permutation checking
  dR = 0.5;
}

void
MeasuredEvent_3l1tau::setBranches(TTree * t)
{
  rle.setBranches(t);
  t -> SetBranchAddress("njets", &njets);

  met.setBranches(t);

  for(int i = 0; i < 3; ++i)
    leptons[i].setBranches(t, Form("lepton%i", (i + 1)));

  for(std::size_t i = 0; i < NOF_RECO_JETS_MEM; ++i)
    allJets[i].setBranches(t, Form("jets%lu", (i + 1)));

  htau.setBranches(t, "htau");

  mvaVariables.setBranches(t);

  if(generatorLevel) generatorLevel -> setBranches(t);
}

void
MeasuredEvent_3l1tau::initNewBranches(TTree * t)
{
  rle.initNewBranches(t);
  branch_njets = t -> Branch("njets", &njets, "njets/I");

  met.initNewBranches(t);

  for(std::size_t i = 0; i < 3; ++i)
    leptons[i].initNewBranches(t, Form("lepton%lu", (i + 1)));

  for(std::size_t i = 0; i < NOF_RECO_JETS_MEM; ++i)
    allJets[i].initNewBranches(t, Form("jets%lu", (i + 1)));

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

bool
MeasuredEvent_3l1tau::isCorrectPermutation() const
{
  if(! generatorLevel)
    return true;
  else if((generatorLevel -> genLepFromTau).dR(leptons[complLeptonIdx]) < dR       &&
          (generatorLevel -> genLepFromTop[0]).dR(leptons[bjetLeptonIdxs[0]]) < dR &&
          (generatorLevel -> genLepFromTop[1]).dR(leptons[bjetLeptonIdxs[1]]) < dR &&
          (generatorLevel -> genBQuarkFromTop[0]).dR(jets[0]) < dR                 &&
          (generatorLevel -> genBQuarkFromTop[1]).dR(jets[1]) < dR)
    return true;
  return false;
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

unsigned
MeasuredEvent_3l1tau::getPermutationNumber() const
{
  return currentPermutation_;
}

void
MeasuredEvent_3l1tau::nextJetCombination() const
{
  ++currentJetCombination_;
}

bool
MeasuredEvent_3l1tau::hasNextJetCombination() const
{
  return currentJetCombination_ < jetCombinations_.size();
}

void
MeasuredEvent_3l1tau::getJetCombination() const
{
  jets = jetCombinations_[currentJetCombination_];
}

void
MeasuredEvent_3l1tau::resetJetCombination() const
{
  currentJetCombination_ = 0;
  jets = jetCombinations_[currentJetCombination_];
}

unsigned
MeasuredEvent_3l1tau::getJetCombinationNumber() const
{
  return currentJetCombination_;
}

void
MeasuredEvent_3l1tau::printPermutation() const
{
  if(generatorLevel)
  {
    LOGTRC << "Printing permutation #" << (currentPermutation_ + 1);
    const double dRlepFromTau  = (generatorLevel -> genLepFromTau).dR(leptons[complLeptonIdx]);
    const double dRlepFromTop1 = (generatorLevel -> genLepFromTop[0]).dR(leptons[bjetLeptonIdxs[0]]);
    const double dRlepFromTop2 = (generatorLevel -> genLepFromTop[1]).dR(leptons[bjetLeptonIdxs[1]]);
    const double dRbQuark1     = (generatorLevel -> genBQuarkFromTop[0]).dR(jets[0]);
    const double dRbQuark2     = (generatorLevel -> genBQuarkFromTop[1]).dR(jets[1]);
    LOGTRC << "Testing lepton from tau";
    LOGTRC << "genLepFromTau:        " << (generatorLevel -> genLepFromTau);
    LOGTRC << "complementary lepton: " << leptons[complLeptonIdx];
    LOGTRC << "\tdR = " << dRlepFromTau << (dRlepFromTau < dR ? " < " : " >= ") << dR
           << " => " << (dRlepFromTau < dR ? "PASS" : "FAIL");
    LOGTRC << "Testing lepton from top (0)";
    LOGTRC << "genLepFromTop[0]: " << (generatorLevel -> genLepFromTop[0]);
    LOGTRC << "b-jet lepton[0]:  " << leptons[bjetLeptonIdxs[0]];
    LOGTRC << "\tdR = " << dRlepFromTop1 << (dRlepFromTop1 < dR ? " < " : " >= ") << dR
           << " => " << (dRlepFromTop1 < dR ? "PASS" : "FAIL");
    LOGTRC << "Testing lepton from top (1)";
    LOGTRC << "genLepFromTop[1]: " << (generatorLevel -> genLepFromTop[1]);
    LOGTRC << "b-jet lepton[1]:  " << leptons[bjetLeptonIdxs[1]];
    LOGTRC << "\tdR = " << dRlepFromTop2 << (dRlepFromTop2 < dR ? " < " : " >= ") << dR
           << " => " << (dRlepFromTop2 < dR ? "PASS" : "FAIL");
    LOGTRC << "Testing b-jet (0)";
    LOGTRC << "genBQuarkFromTop[0]: " << (generatorLevel -> genBQuarkFromTop[0]);
    LOGTRC << "b-jet[0]:            " << jets[0];
    LOGTRC << "\tdR = " << dRbQuark1 << (dRbQuark1 < dR ? " < " : " >= ") << dR
           << " => " << (dRbQuark1 < dR ? "PASS" : "FAIL");
    LOGTRC << "Testing b-jet (1)";
    LOGTRC << "genBQuarkFromTop[1]: " << (generatorLevel -> genBQuarkFromTop[1]);
    LOGTRC << "b-jet[1]:            " << jets[1];
    LOGTRC << "\tdR = " << dRbQuark2 << (dRbQuark2 < dR ? " < " : " >= ") << dR
           << " => " << (dRbQuark2 < dR ? "PASS" : "FAIL");
  }
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredEvent_3l1tau & event)
  {
    os << "The event (" << event.rle.str() << "), "
       << "permutation #" << (event.currentPermutation_ + 1) << ":\n";
    for(int i = 0; i < 3; ++i)
      os << "\tLepton " << (i + 1) << ": " << event.leptons[i] << '\n';
    for(int i = 0; i < 2; ++i)
      os << "\tJet "    << (i + 1) << ":    " << event.jets[i]    << '\n';
    os << "\tTau:      " << event.htau << '\n';
    os << "\tMET:      " << event.met   << '\n';
    if(event.generatorLevel)
      os << *event.generatorLevel;
    return os;
  }
}
