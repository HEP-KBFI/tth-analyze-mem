#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau

using namespace tthMEM;

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau(TFile * file)
  : file_(file)
{
  for(int i = hVar::kZ1; i < hVar::kProb; ++i)
    histograms_[static_cast<hVar>(i)] = 0;
}

DebugPlotter_ttHorZ_3l1tau::DebugPlotter_ttHorZ_3l1tau()
  : DebugPlotter_ttHorZ_3l1tau(0)
{ }

void
DebugPlotter_ttHorZ_3l1tau::setFile(TFile * file)
{
  file_ = file;
}

void
DebugPlotter_ttHorZ_3l1tau::initialize(const std::string & dirName)
{
  dirName_ = dirName;
  if(file_)
  {
    file_ -> mkdir(dirName_.c_str());
    file_ -> cd(dirName_.c_str());
  }

  histograms_[hVar::kZ1]       = new TH1D("z1", "z1", 100, 0., 1.);
  histograms_[hVar::kZ2]       = new TH1D("z2", "z2", 100, 0., 1.);
  histograms_[hVar::kMassHorZ] = new TH1D("massHorZ", "massHorZ", 200, 0., 600.);
  histograms_[hVar::kMassHtau] = new TH1D("massHtau", "massHtau", 200, 0., 600.);
  histograms_[hVar::kMassLtau] = new TH1D("massLtau", "massLtau", 200, 0., 600.);
  histograms_[hVar::kB1en]     = new TH1D("B1en", "B1en", 200, 0., 600.);
  histograms_[hVar::kB2en]     = new TH1D("B2en", "B2en", 200, 0., 600.);
  histograms_[hVar::kB1RecoEn] = new TH1D("B1RecoEn", "B1RecoEn", 200, 0., 1000.);
  histograms_[hVar::kB2RecoEn] = new TH1D("B2RecoEn", "B2RecoEn", 200, 0., 1000.);
  histograms_[hVar::kMETpull]  = new TH1D("METpull", "METpull", 200, 0., 300.);
  histograms_[hVar::kXa]       = new TH1D("xa", "xa", 100, 0., 1.);
  histograms_[hVar::kXb]       = new TH1D("xb", "xb", 100, 0., 1.);
  histograms_[hVar::kHtauJPS]  = new TH1D("hTauJPS", "hTauJPS", 200, 0., 1.e-5);
  histograms_[hVar::kLtauJPS]  = new TH1D("lTauJPS", "lTauJPS", 200, 0., 1.e-5);
  histograms_[hVar::kMsquared] = new TH1D("Msquared", "Msquared", 200, 0., 1.e-5);
  histograms_[hVar::kProb]     = new TH1D("prob", "prob", 200, 0., 1.e-5);
}

void
DebugPlotter_ttHorZ_3l1tau::fill(hVar var,
                                 double value)
{
  if(histograms_[var]) histograms_[var] -> Fill(value);
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
