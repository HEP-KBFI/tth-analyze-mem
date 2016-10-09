#include "tthAnalysis/tthMEM/interface/GeneratorLevelEvent_3l1tau.h"

#include <TString.h> // Form()

using namespace tthMEM;

void
GeneratorLevelEvent_3l1tau::initialize()
{
  genNuFromLtau.initialize();
  genNuLepFromTau.initialize();
  genLepFromTau.initialize();

  genNuFromHtau.initialize();
  genHtau.initialize();

  // make sure that the neutrinos are actually neutral
  // by default the ,,charge'' here is actually opposite sign of signum(pdgId)
  genNuFromLtau.setNeutral();
  genNuLepFromTau.setNeutral();
  genNuFromHtau.setNeutral();

  for(std::size_t i = 0; i < 2; ++i)
  {
    genTau[i].initialize();

    genBQuarkFromTop[i].initialize();
    genNuFromTop[i].initialize();
    genLepFromTop[i].initialize();

    genNuFromTop[i].setNeutral();

    genWBoson[i] = genNuFromTop[i] + genLepFromTop[i];
    genTop[i] = genWBoson[i] + genBQuarkFromTop[i];
  }

  genHorZ = genTau[0] + genTau[1];
  genDiNuFromLtau = genNuFromLtau + genNuLepFromTau;

  if(genLepFromTau.charge() == genTau[0].charge())
  {
    leptonicTauDecayIdx = 0;
    hadronicTauDecayIdx = 1;
  }
  else
  {
    leptonicTauDecayIdx = 1;
    hadronicTauDecayIdx = 0;
  }
}

void
GeneratorLevelEvent_3l1tau::setBranches(TChain * t)
{
  genNuFromLtau.setBranches(t, "genNuFromLtau");
  genNuLepFromTau.setBranches(t, "genNuLepFromTau");
  genLepFromTau.setBranches(t, "genLepFromTau");

  genNuFromHtau.setBranches(t, "genNuFromHtau");
  genHtau.setBranches(t, "genHtau");

  for(std::size_t i = 0; i < 2; ++i)
  {
    genTau[i].setBranches(t, Form("genTau%lu", (i + 1)));

    genBQuarkFromTop[i].setBranches(t, Form("genBQuarkFromTop%lu", (i + 1)));
    genNuFromTop[i].setBranches(t, Form("genNuFromTop%lu", (i + 1)));
    genLepFromTop[i].setBranches(t, Form("genLepFromTop%lu", (i + 1)));
  }
}

void
GeneratorLevelEvent_3l1tau::initNewBranches(TTree * t)
{
  genNuFromLtau.initNewBranches(t, "genNuFromLtau");
  genNuLepFromTau.initNewBranches(t, "genNuLepFromTau");
  genLepFromTau.initNewBranches(t, "genLepFromTau");

  genNuFromHtau.initNewBranches(t, "genNuFromHtau");
  genHtau.initNewBranches(t, "genHtau");

  for(std::size_t i = 0; i < 2; ++i)
  {
    genTau[i].initNewBranches(t, Form("genTau%lu", (i + 1)));

    genBQuarkFromTop[i].initNewBranches(t, Form("genBQuarkFromTop%lu", (i + 1)));
    genNuFromTop[i].initNewBranches(t, Form("genNuFromTop%lu", (i + 1)));
    genLepFromTop[i].initNewBranches(t, Form("genLepFromTop%lu", (i + 1)));

    genWBoson[i].initNewBranches(t, Form("genWBoson%lu", (i + 1)));
    genTop[i].initNewBranches(t, Form("genTop%lu", (i + 1)));
  }

  genHorZ.initNewBranches(t, "genHorZ");
  genDiNuFromLtau.initNewBranches(t, "genDiNuFromLtau");
}

