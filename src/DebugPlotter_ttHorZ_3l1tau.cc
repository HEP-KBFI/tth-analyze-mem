#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau

#include <TString.h> // Form()

using namespace tthMEM;

const unsigned DebugPlotter_ttHorZ_3l1tau::nofBins_ = 200;

const decltype(DebugPlotter_ttHorZ_3l1tau::hVarMap_)
DebugPlotter_ttHorZ_3l1tau::hVarMap_ =
{
  { hVar_3l1tau::kZ2,       { "z2",       0.,  1.     } },
  { hVar_3l1tau::kMassHorZ, { "massHorZ", 70., 140.   } },
  { hVar_3l1tau::kMassHtau, { "massHtau", 0.,  3.     } },
  { hVar_3l1tau::kMassLtau, { "massLtau", 0.,  3.     } },
  { hVar_3l1tau::kB1en,     { "B1en",     0.,  600.   } },
  { hVar_3l1tau::kB2en,     { "B2en",     0.,  600.   } },
  { hVar_3l1tau::kB1RecoEn, { "B1RecoEn", 0.,  600.   } },
  { hVar_3l1tau::kB2RecoEn, { "B2RecoEn", 0.,  600.   } },
  { hVar_3l1tau::kMsquared, { "Msquared", 0.,  1.e-6  } },
  { hVar_3l1tau::kProb,     { "prob",     0.,  1.e-45 } }
};

DebugPlotter_ttHorZ_3l1tau::HVarHolder::HVarHolder(const std::string & hVarName,
                                                   double begin,
                                                   double end)
  : hVarName_(hVarName)
  , limits_(begin, end)
{}

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file,
                                                       unsigned debugFrequency)
  : file_(file)
  , debugFrequency_(debugFrequency)
  , debugRange_(8)
  , logCounter_(0)
  , log_(false)
{
  for(auto e: Enum<hVar_3l1tau>()) histRecod_[e] = 0;
  for(auto e: Enum<Var_3l1tau>())  histSampled_[e] = 0;
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
  log_ = debugFrequency_ > 0 &&
         logCounter_ % debugFrequency_ >= 1 &&
         logCounter_ % debugFrequency_ <= debugRange_;
  if(! log_) return;

  if(file_)
  {
    file_ -> mkdir(dirName.c_str());
    file_ -> cd(dirName.c_str());

    for(auto var: Enum<hVar_3l1tau>())
    {
      const TString varName(hVarMap_.at(var).hVarName_.c_str());
      const Limits & limits = hVarMap_.at(var).limits_;
      histRecod_[var] = new TH1D(varName, varName, nofBins_,
                                 limits.begin_, limits.end_);
    }

//--- input/sampled variables
    for(auto var: Enum<Var_3l1tau>())
    {
      const TString varName = Form("%s_sampled", vm.getVarName(var).c_str());
      const Limits & limits = vm.getVarLimits(var);
      histSampled_[var] = new TH1D(varName, varName, nofBins_,
                                   limits.begin_, limits.end_);
    }
  }
}

DebugPlotter_ttHorZ_3l1tau &
DebugPlotter_ttHorZ_3l1tau::fill(hVar_3l1tau var,
                                 double value)
{
  if(log_ && histRecod_[var]) histRecod_[var] -> Fill(value);
  return *this;
}

DebugPlotter_ttHorZ_3l1tau&
DebugPlotter_ttHorZ_3l1tau::fill(const VariableManager_3l1tau & vm,
                                 const double * const x)
{
  for(auto var: Enum<Var_3l1tau>())
    if(log_ && histSampled_[var]) histSampled_[var] -> Fill(vm.get(var, x));
  return *this;
}

void
DebugPlotter_ttHorZ_3l1tau::write()
{
  for(auto & kv: histRecod_)
    if(kv.second)
    {
      kv.second -> Write();
      delete kv.second;
      kv.second = 0;
    }
  for(auto & kv: histSampled_)
    if(kv.second)
    {
      kv.second -> Write();
      delete kv.second;
      kv.second = 0;
    }
}
