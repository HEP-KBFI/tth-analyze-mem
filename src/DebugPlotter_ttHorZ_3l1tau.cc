#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <TString.h> // Form()

#include <algorithm> // std::replace()

#define PLACEHOLDER_DEBUGPLOTTER -9999

using namespace tthMEM;

const decltype(DebugPlotter_ttHorZ_3l1tau::hVarMap_)
DebugPlotter_ttHorZ_3l1tau::hVarMap_ =
{
  { hVar_3l1tau::kDz1,             "diff_z1",             },
  { hVar_3l1tau::kDz2,             "diff_z2",             },
  { hVar_3l1tau::kDMassHorZ,       "diff_mHorZ",          },
  { hVar_3l1tau::kDnuHtauCosTheta, "diff_nuHtauCosTheta", },
  { hVar_3l1tau::kDnuLtauCosTheta, "diff_nuLtauCosTheta", },
  { hVar_3l1tau::kDnuHtauPhi,      "diff_nuHtauPhi",      },
  { hVar_3l1tau::kDnuLtauPhi,      "diff_nuLtauPhi",      },
  { hVar_3l1tau::kDmetX,           "diff_METx",           },
  { hVar_3l1tau::kDmetY,           "diff_METy",           },
  { hVar_3l1tau::kZ2,              "z2",                  },
  { hVar_3l1tau::kMassHorZ,        "massHorZ"             },
  { hVar_3l1tau::kMassHtau,        "massHtau"             },
  { hVar_3l1tau::kMassLtau,        "massLtau"             },
  { hVar_3l1tau::kB1en,            "B1en"                 },
  { hVar_3l1tau::kB2en,            "B2en"                 },
  { hVar_3l1tau::kB1RecoEn,        "B1RecoEn"             },
  { hVar_3l1tau::kB2RecoEn,        "B2RecoEn"             },
  { hVar_3l1tau::kB1energyTF,      "B1energyTF"           },
  { hVar_3l1tau::kB2energyTF,      "B2energyTF"           },
  { hVar_3l1tau::kXa,              "xa"                   },
  { hVar_3l1tau::kXb,              "xb"                   },
  { hVar_3l1tau::kFlux,            "flux"                 },
  { hVar_3l1tau::kHtauPSF,         "hTauPSF"              },
  { hVar_3l1tau::kLtauPSF,         "lTauPSF"              },
  { hVar_3l1tau::kProbPDF,         "probPDF"              },
  { hVar_3l1tau::kTdecayJF1,       "tDecayJF1"            },
  { hVar_3l1tau::kTdecayJF2,       "tDecayJF2"            },
  { hVar_3l1tau::kJacobiF,         "jacobiF"              },
  { hVar_3l1tau::kMETx,            "METx"                 },
  { hVar_3l1tau::kMETy,            "METy"                 },
  { hVar_3l1tau::kMETtf,           "METtf"                },
  { hVar_3l1tau::kMsquared,        "Msquared"             },
  { hVar_3l1tau::kProb,            "prob"                 }
};

DebugString_ttHorZ_3l1tau::DebugString_ttHorZ_3l1tau(const std::string & rle_,
                                                     const std::string & me_,
                                                     unsigned leptonPermutation_,
                                                     unsigned bJetCombination_)
  : me(me_)
  , leptonPermutation(leptonPermutation_)
  , bJetCombination(bJetCombination_)
{
  rle = rle_;
  std::replace(rle.begin(), rle.end(), ':', '_');
}

DebugString_ttHorZ_3l1tau::DebugString_ttHorZ_3l1tau()
  : DebugString_ttHorZ_3l1tau("", "", 0, 0)
{}

std::string
DebugString_ttHorZ_3l1tau::str() const
{
  return std::string(Form(
    "%s_%s_%u_%u", rle.c_str(), me.c_str(), leptonPermutation, bJetCombination
  ));
}

bool
DebugString_ttHorZ_3l1tau::hasSameRLE(const DebugString_ttHorZ_3l1tau & other) const
{
  return rle == other.rle;
}

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file,
                                                       unsigned debugFrequency)
  : file_(file)
  , tree_(nullptr)
  , debugFrequency_(debugFrequency)
  , logCounter_(0)
  , log_(false)
  , isFilled_(false)
{
  reset();
}

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file)
  : DebugPlotter_ttHorZ_3l1tau(file, 1)
{ }

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau()
  : DebugPlotter_ttHorZ_3l1tau(nullptr, 1)
{ }

DebugPlotter_ttHorZ_3l1tau::~DebugPlotter_ttHorZ_3l1tau()
{
  write();
}

void
DebugPlotter_ttHorZ_3l1tau::initialize(const DebugString_ttHorZ_3l1tau & metaInfo,
                                       const VariableManager_3l1tau & vm)
{
  if(! metaInfo_.hasSameRLE(metaInfo))
    ++logCounter_;
  log_ = (debugFrequency_ > 1 && (logCounter_ == 1 || logCounter_ % debugFrequency_ == 1)) ||
          debugFrequency_ == 1;
  metaInfo_ = metaInfo;
  if(! log_)
    return;

  if(file_)
  {
    file_ -> mkdir(metaInfo_.str().c_str());
    file_ -> cd(metaInfo_.str().c_str());
    tree_ = new TTree("tree", metaInfo_.str().c_str());

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
  {
    tree_ -> Fill();
    isFilled_ = true;
  }
  return *this;
}

void
DebugPlotter_ttHorZ_3l1tau::write()
{
  if(tree_)
  {
    if(file_ && ! isFilled_)
    {
      delete tree_;
      tree_ = nullptr;
      file_ -> rmdir(metaInfo_.str().c_str());
    }
    else
    {
      tree_ -> Write();
      delete tree_;
      tree_ = nullptr;
    }
    reset();
  }
}

void
DebugPlotter_ttHorZ_3l1tau::reset()
{
  for(auto e: Enum<hVar_3l1tau>()) histRecod_[e]   = PLACEHOLDER_DEBUGPLOTTER;
  for(auto e: Enum<Var_3l1tau>())  histSampled_[e] = PLACEHOLDER_DEBUGPLOTTER;
  isFilled_ = false;
}
