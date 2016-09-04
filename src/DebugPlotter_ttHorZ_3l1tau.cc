#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau

using namespace tthMEM;

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file,
                                                       unsigned debugFrequency)
  : file_(file)
  , debugFrequency_(debugFrequency)
  , debugRange_(8)
  , logCounter_(0)
  , log_(false)
{
  for(auto e: Enum<hVar_3l1tau>()) histograms_[e] = 0;
}

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file)
  : DebugPlotter_ttHorZ_3l1tau(file, 1)
{ }

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau()
  : DebugPlotter_ttHorZ_3l1tau(0, 1)
{ }

void
DebugPlotter_ttHorZ_3l1tau::initialize(const std::string & dirName)
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

    histograms_[hVar_3l1tau::kZ2]       = new TH1D("z2",       "z2",       100, 0., 1.);
    histograms_[hVar_3l1tau::kMassHorZ] = new TH1D("massHorZ", "massHorZ", 200, 0., 200.);
    histograms_[hVar_3l1tau::kMassHtau] = new TH1D("massHtau", "massHtau", 200, 0., 3.);
    histograms_[hVar_3l1tau::kMassLtau] = new TH1D("massLtau", "massLtau", 200, 0., 3.);
    histograms_[hVar_3l1tau::kB1en]     = new TH1D("B1en",     "B1en",     200, 0., 600.);
    histograms_[hVar_3l1tau::kB2en]     = new TH1D("B2en",     "B2en",     200, 0., 600.);
    histograms_[hVar_3l1tau::kB1RecoEn] = new TH1D("B1RecoEn", "B1RecoEn", 200, 0., 600.);
    histograms_[hVar_3l1tau::kB2RecoEn] = new TH1D("B2RecoEn", "B2RecoEn", 200, 0., 600.);
    histograms_[hVar_3l1tau::kMsquared] = new TH1D("Msquared", "Msquared", 200, 0., 1.e-10);
    histograms_[hVar_3l1tau::kProb]     = new TH1D("prob",     "prob",     200, 0., 1.e-51);
  }
}

DebugPlotter_ttHorZ_3l1tau &
DebugPlotter_ttHorZ_3l1tau::fill(hVar_3l1tau var,
                                 double value)
{
  if(! log_) return *this;
  if(histograms_[var]) histograms_[var] -> Fill(value);
  return *this;
}

void
DebugPlotter_ttHorZ_3l1tau::write()
{
  for(auto & kv: histograms_)
    if(kv.second)
    {
      kv.second -> Write();
      delete kv.second;
      kv.second = 0;
    }
}
