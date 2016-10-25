#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau

#include <TString.h> // Form()

#define PLACEHOLDER_DEBUGPLOTTER -9999

using namespace tthMEM;

const decltype(DebugPlotter_ttHorZ_3l1tau::hVarMap_)
DebugPlotter_ttHorZ_3l1tau::hVarMap_ =
{
  { hVar_3l1tau::kZ2,         "z2",        },
  { hVar_3l1tau::kMassHorZ,   "massHorZ"   },
  { hVar_3l1tau::kMassHtau,   "massHtau"   },
  { hVar_3l1tau::kMassLtau,   "massLtau"   },
  { hVar_3l1tau::kB1en,       "B1en"       },
  { hVar_3l1tau::kB2en,       "B2en"       },
  { hVar_3l1tau::kB1RecoEn,   "B1RecoEn"   },
  { hVar_3l1tau::kB2RecoEn,   "B2RecoEn"   },
  { hVar_3l1tau::kB1energyTF, "B1energyTF" },
  { hVar_3l1tau::kB2energyTF, "B2energyTF" },
  { hVar_3l1tau::kXa,         "xa"         },
  { hVar_3l1tau::kXb,         "xb"         },
  { hVar_3l1tau::kFlux,       "flux"       },
  { hVar_3l1tau::kHtauPSF,    "hTauPSF"    },
  { hVar_3l1tau::kLtauPSF,    "lTauPSF"    },
  { hVar_3l1tau::kProbPDF,    "probPDF"    },
  { hVar_3l1tau::kTdecayJF1,  "tDecayJF1"  },
  { hVar_3l1tau::kTdecayJF2,  "tDecayJF2"  },
  { hVar_3l1tau::kJacobiF,    "jacobiF"    },
  { hVar_3l1tau::kMETtf,      "METtf"      },
  { hVar_3l1tau::kMsquared,   "Msquared"   },
  { hVar_3l1tau::kProb,       "prob"       }
};

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file,
                                                       unsigned debugFrequency)
  : file_(file)
  , tree_(0)
  , debugFrequency_(debugFrequency)
  , debugRange_(8)
  , logCounter_(0)
  , log_(false)
{
  reset();
}

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file)
  : DebugPlotter_ttHorZ_3l1tau(file, 1)
{ }

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau()
  : DebugPlotter_ttHorZ_3l1tau(0, 1)
{ }

void
DebugPlotter_ttHorZ_3l1tau::initialize(const std::string & dirName,
                                       const VariableManager_3l1tau & vm)
{
  ++logCounter_;
  log_ = (debugFrequency_ == 1) ||
         (debugFrequency_ > 0 &&
          logCounter_ % debugFrequency_ >= 1 &&
          logCounter_ % debugFrequency_ <= debugRange_);
  if(! log_) return;

  if(file_)
  {
    file_ -> mkdir(dirName.c_str());
    file_ -> cd(dirName.c_str());
    tree_ = new TTree("tree", dirName.c_str());

//--- variables declared in this file
    for(auto var: Enum<hVar_3l1tau>())
    {
      const std::string varName(hVarMap_.at(var));
      tree_ -> Branch(varName.c_str(), &(histRecod_[var]), Form("%s/D", varName.c_str()));
    }

//--- input/sampled variables
    for(auto var: Enum<Var_3l1tau>())
    {
      const std::string varName(Form("%s_sampled", vm.getVarName(var).c_str()));
      tree_ -> Branch(varName.c_str(), &(histSampled_[var]), Form("%s/D", varName.c_str()));
    }
  }
}

DebugPlotter_ttHorZ_3l1tau &
DebugPlotter_ttHorZ_3l1tau::fill(hVar_3l1tau var,
                                 double value)
{
  if(log_)
    histRecod_[var] = value;
  return *this;
}

DebugPlotter_ttHorZ_3l1tau &
DebugPlotter_ttHorZ_3l1tau::fill(const VariableManager_3l1tau & vm,
                                 const double * const x)
{
  for(auto var: Enum<Var_3l1tau>())
    if(log_)
      histSampled_[var] = vm.get(var, x);
  return *this;
}

DebugPlotter_ttHorZ_3l1tau &
DebugPlotter_ttHorZ_3l1tau::fill()
{
  if(tree_)
    tree_ -> Fill();
  return *this;
}

void
DebugPlotter_ttHorZ_3l1tau::write()
{
  if(tree_)
  {
    tree_ -> Write();
    delete tree_;
    tree_ = 0;
    reset();
  }
}

void
DebugPlotter_ttHorZ_3l1tau::reset()
{
  for(auto e: Enum<hVar_3l1tau>()) histRecod_[e]   = PLACEHOLDER_DEBUGPLOTTER;
  for(auto e: Enum<Var_3l1tau>())  histSampled_[e] = PLACEHOLDER_DEBUGPLOTTER;
}
