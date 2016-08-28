#include <cstdlib> // std::getenv(), EXIT_SUCCESS
#include <vector> // std::vector<>
#include <algorithm> // std::generate()
#include <string> // std::string

#include <boost/filesystem/path.hpp> // boost::filesystem::path
#include <boost/filesystem/operations.hpp> // boost::filesystem::exists(), ...
  // ... boost::filesystem::create_directories()

#include <TCanvas.h> // TCanvas
#include <TGraph.h> // TGraph
#include <TMultiGraph.h> // TMultiGraph
#include <TAxis.h> // SetTitle()
#include <TLegend.h> // TLegend

#include "tthAnalysis/tthMEM/interface/BJetTransferFunction.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

/**
 * @brief Reproduces Fig. 2 in AN2013-313
 */
int
main()
{
  namespace fs = boost::filesystem;
  fs::path plot_dir;
  if(const char * cmssw_base = std::getenv("CMSSW_BASE"))
  {
    plot_dir = fs::path(cmssw_base) / fs::path("src/tthAnalysis/tthMEM/plots");
    if(! fs::exists(plot_dir))
    {
      fs::create_directories(plot_dir);
      LOGINFO << "Created director " << plot_dir;
    }
  }
  else LOGERR << "Environment variable CMSSW_BASE not set!";

  const std::vector<double> eta = { 0.5, 1.5 };
  const double jetEmin = 0.,
               jetEmax = 225.;
  const unsigned nofSteps = 400;
  const std::vector<double> bQuarkEs = { 50., 100., 150. };
  const std::vector<unsigned> colors = { 6, 2, 4 };
  const std::vector<unsigned> lineTypes = { 1, 2 };
  std::vector<double> jetEs(nofSteps);
  std::generate(jetEs.begin(), jetEs.end(),
    [=] {
      static unsigned x = 0;
      return (jetEmax - jetEmin) * (++x) / nofSteps;
    }
  );
  const double * const x = jetEs.data();

  TCanvas * c = new TCanvas("c", "c", 200, 10, 700, 500);
  TGraph * grs[bQuarkEs.size() * eta.size()];
  TMultiGraph * mg = new TMultiGraph();

  unsigned grIdx = 0;
  for(unsigned etaIdx = 0; etaIdx < eta.size(); ++ etaIdx)
    for(unsigned i = 0; i < bQuarkEs.size(); ++i)
    {
      double y[jetEs.size()];
      for(unsigned j = 0; j < jetEs.size(); ++j)
        y[j] = tthMEM::functions::bJetTF(bQuarkEs[i], jetEs[j], eta[etaIdx]);
      grs[grIdx] = new TGraph(jetEs.size(), x, y);
      grs[grIdx] -> SetLineColor(colors[i]);
      grs[grIdx] -> SetLineStyle(lineTypes[etaIdx]);
      grs[grIdx] -> SetFillColor(0);
      mg -> Add(grs[grIdx]);
      ++grIdx;
    }

  mg -> Draw("AC");
  mg -> SetTitle("CMS simulation #sqrt{s}=8 TeV");
  mg -> GetXaxis() -> SetTitle("jet energy (GeV)");
  mg -> GetXaxis() -> SetLimits(jetEmin, jetEmax);
  mg -> GetYaxis() -> SetTitle("transfer function (1/GeV)");
  mg -> SetMinimum(0.);
  mg -> SetMaximum(0.07);

  TLegend * leg = new TLegend(0.15, 0.6, 0.45, 0.8);
  leg -> SetHeader("b-quarks");
  leg -> AddEntry(grs[0], "E_{q} = 50 GeV");
  leg -> AddEntry(grs[1], "E_{q} = 100 GeV");
  leg -> AddEntry(grs[2], "E_{q} = 150 GeV");
  leg -> SetBorderSize(0);
  leg -> Draw();

  TLegend * leg2 = new TLegend(0.63, 0.45, 0.87, 0.55);
  leg2 -> AddEntry(static_cast<TObject *>(0), "Solid: |#eta| < 1.0", "");
  leg2 -> AddEntry(static_cast<TObject *>(0), "Dashed: |#eta| #geq 1.0", "");
  leg2 -> SetBorderSize(0);
  leg2 -> Draw();

  const fs::path plot_path = plot_dir / fs::path("bjetTF.pdf");
  c -> SaveAs(plot_path.c_str());

  return EXIT_SUCCESS;
}
