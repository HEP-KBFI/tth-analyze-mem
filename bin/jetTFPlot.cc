#include "tthAnalysis/tthMEM/interface/JetTransferFunction.h" // bJetTF(), qJetTF()
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

#include <boost/filesystem/path.hpp> // boost::filesystem::path
#include <boost/filesystem/operations.hpp> // boost::filesystem::exists(), ...
 // ... boost::filesystem::create_directories()
#include <boost/program_options/options_description.hpp>
 // boost::program_options::options_description
#include <boost/program_options/value_semantic.hpp>
 // boost::program_options::bool_switch
#include <boost/program_options/variables_map.hpp>
 // boost::program_options::variables_map, boost::program_options::store(),
 // boost::program_options::notify()
#include <boost/program_options/parsers.hpp>
 // boost::program_options::parse_command_line()
#include <boost/program_options/errors.hpp> // boost::program_options::error

#include <TCanvas.h> // TCanvas
#include <TGraph.h> // TGraph
#include <TMultiGraph.h> // TMultiGraph
#include <TAxis.h> // SetTitle()
#include <TLegend.h> // TLegend
#include <TString.h> // Form()

#include <functional> // std::function<>

void
plot(const std::function<double(double, double, double)> & jetTF,
     const std::string & type,
     const std::string & fileName)
{
  namespace fs = boost::filesystem;
  fs::path plot_dir;
  if(const char * cmssw_base = std::getenv("CMSSW_BASE"))
  {
    plot_dir = fs::path(cmssw_base) / fs::path("src/tthAnalysis/tthMEM/data/plots");
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
  const std::vector<double> quarkEs = { 50., 100., 150. };
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

  TCanvas * c = new TCanvas(type.c_str(), type.c_str(), 200, 10, 700, 500);
  TGraph * grs[quarkEs.size() * eta.size()];
  TMultiGraph * mg = new TMultiGraph();

  unsigned grIdx = 0;
  for(unsigned etaIdx = 0; etaIdx < eta.size(); ++ etaIdx)
    for(unsigned i = 0; i < quarkEs.size(); ++i)
    {
      double y[jetEs.size()];
      for(unsigned j = 0; j < jetEs.size(); ++j)
        y[j] = jetTF(quarkEs[i], jetEs[j], eta[etaIdx]);
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
  leg -> SetHeader(Form("%s-quarks", type.c_str()));
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

  const fs::path plot_path = plot_dir / fs::path(fileName);
  c -> SaveAs(plot_path.c_str());

  delete mg;
  delete leg;
  delete leg2;
  delete c;
}

/**
 * @brief Reproduces Fig. 2 in AN2013-313
 */
int
main(int argc,
     char * argv[])
{
  bool type;
  namespace po = boost::program_options;
  try
  {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h",   "produce help message")
      ("light,l",  po::bool_switch(&type) -> default_value(false),
                   "Plot light jet TF (default: b-jet TF)")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    if(vm.count("help"))
    {
      std::cout << desc << '\n';
      return EXIT_SUCCESS;
    }
    po::notify(vm); // might throw
  }
  catch(po::error & e)
  {
    std::cerr << "User error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  catch(std::exception & e)
  {
    std::cerr << "Error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
  catch(...)
  {
    std::cerr << "Unknown error!\n";
    return EXIT_FAILURE;
  }

  if(type)
    plot(tthMEM::functions::qJetTF, "u", "ujetTF.pdf");
  else
    plot(tthMEM::functions::bJetTF, "b", "bjetTF.pdf");

  return EXIT_SUCCESS;
}
