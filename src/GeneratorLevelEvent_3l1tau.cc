#include "tthAnalysis/tthMEM/interface/GeneratorLevelEvent_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/tthMEMrecFunctions.h" // functions::
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // pow2()
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGTRC

#include <TString.h> // Form()

using namespace tthMEM;

GeneratorLevelEvent_3l1tau::GeneratorLevelEvent_3l1tau()
  : beamAxis_(0., 0., 1.)
{
  for(auto var: Enum<Var_3l1tau>())
    genIntVariables[var] = 0.;
}

void
GeneratorLevelEvent_3l1tau::initialize()
{
  genNuFromLtau.initialize();
  genNuLepFromTau.initialize();
  genLepFromTau.initialize();

  genNuFromHtau.initialize();
  genHtau.initialize();

//--- make sure that the neutrinos are actually neutral
//--- by default the ,,charge'' here is actually opposite sign of signum(pdgId)
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

//--- since by construction the first genTau has positive sign, we have to find
//--- the association between tau decay mode and genTau ourselves
  if(genLepFromTau.charge() == genTau[0].charge())
  {
    leptonicTauDecayIdx_ = 0;
    hadronicTauDecayIdx_ = 1;
  }
  else
  {
    leptonicTauDecayIdx_ = 1;
    hadronicTauDecayIdx_ = 0;
  }

//--- calculate the integration variables in advance (we need to fill the TTree)
//--- should the top decay related variables be found w.r.t the beam axis?
  genIntVariables[Var_3l1tau::kBcosTheta1]     = std::cos(genNuFromTop[0].p3().theta());
  genIntVariables[Var_3l1tau::kBphi1]          = genNuFromTop[0].phi();
  genIntVariables[Var_3l1tau::kBcosTheta2]     = std::cos(genNuFromTop[1].p3().theta());
  genIntVariables[Var_3l1tau::kBphi2]          = genNuFromTop[1].phi();
  genIntVariables[Var_3l1tau::kZ1]             = functions::z(
    genTau[hadronicTauDecayIdx_].p4(), genHtau.p4()
  );
  genIntVariables[Var_3l1tau::kTauPhi]         = functions::phiFromLabMomenta(
    genTau[hadronicTauDecayIdx_].p4(), genHtau.p4(), beamAxis_
  );
  genIntVariables[Var_3l1tau::kTauPhiInv]      = functions::phiFromLabMomenta(
    genTau[leptonicTauDecayIdx_].p4(), genLepFromTau.p4(), beamAxis_
  );
  genIntVariables[Var_3l1tau::kTauMinvSquared] = pow2(genDiNuFromLtau.mass());
}

void
GeneratorLevelEvent_3l1tau::setBranches(TTree * t)
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

  for(auto var: Enum<Var_3l1tau>())
  {
    const std::string varName = VariableManager_3l1tau::getVarName(var);
    t -> Branch(Form("gen_%s",   varName.c_str()), &(genIntVariables[var]),
                Form("gen_%s/D", varName.c_str()));
  }
}

void
GeneratorLevelEvent_3l1tau::setIntegrationVariables(VariableManager_3l1tau & vm) const
{
  for(const Var_3l1tau & var: vm.generatorLevels)
    vm.set(var, genIntVariables.at(var));
}

void
GeneratorLevelEvent_3l1tau::setBeamAxis(const Vector & beamAxis)
{
  beamAxis_ = beamAxis;
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const GeneratorLevelEvent_3l1tau & event)
  {
    const std::size_t & lIdx = event.leptonicTauDecayIdx_;
    const std::size_t & hIdx = event.hadronicTauDecayIdx_;
    os << "\tgenNuFromLtau:           " << event.genNuFromLtau   << '\n'
       << "\tgenNuLepFromTau:         " << event.genNuLepFromTau << '\n'
       << "\tgenLepFromTau:           " << event.genLepFromTau   << '\n'
       << "\tgenNuFromHtau:           " << event.genNuFromHtau   << '\n'
       << "\tgenHtau:                 " << event.genHtau         << '\n'
       << "\tgenTau (leptonic decay): " << event.genTau[lIdx]    << '\n'
       << "\tgenTau (hadronic decay): " << event.genTau[hIdx]    << '\n'
       << "\tgenHorZ:                 " << event.genHorZ         << '\n'
       << "\tgenDiNuFromLtau:         " << event.genDiNuFromLtau << '\n';
    for(std::size_t i = 0; i < 2; ++i)
      os << "\tgenBQuarkFromTop[" << (i + 1) << "]:     " << event.genBQuarkFromTop[i] << '\n'
         << "\tgenNuFromTop[" << (i + 1) << "]:         " << event.genNuFromTop[i]     << '\n'
         << "\tgenLepFromTop[" << (i + 1) << "]:        " << event.genLepFromTop[i]    << '\n'
         << "\tgenWBoson[" << (i + 1) << "]:            " << event.genWBoson[i]        << '\n'
         << "\tgenTop[" << (i + 1) << "]:               " << event.genTop[i]           << '\n';
    return os;
  }
}

