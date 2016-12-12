#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE, setenv(), std::atexit()
#include <string> // std::string
#include <vector> // std::vector<>
#include <algorithm> // std::inner_product()
#include <csignal> // std::signal(), SIG*

#include <boost/filesystem/operations.hpp> // boost::filesystem::exists()

#include <Rtypes.h> // Long64_t
#include <TFile.h> // TFile
#include <TChain.h> // TChain
#include <TTree.h> // TTree
#include <TString.h> // Form()

#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h" // tthMEM::MeasuredEvent_3l1tau
#include "tthAnalysis/tthMEM/interface/MEM_ttHorZ_3l1tau.h" // tthMEM::MEM_ttHorZ_3l1tau
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // tthMEM::findFile()
#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // tthMEM::constants
#include "tthAnalysis/tthMEM/interface/tthMEMenums.h" // tthMEM::ME_mg5_3l1tau::
#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau
#include "tthAnalysis/tthMEM/interface/VariableManager_3l1tau.h" // VariableManager_3l1tau
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

#include <FWCore/ParameterSet/interface/ParameterSet.h> // edm::ParameterSet
#include <FWCore/PythonParameterSet/interface/MakeParameterSets.h> // edm::readPSetsFrom()
#include <DataFormats/FWLite/interface/InputSource.h> // fwlite::InputSource
#include <DataFormats/FWLite/interface/OutputFiles.h> // fwlite::OutputFiles

// JFC!!! fwlite:: doesn't include "FWCore/Utilities/interface/Exception.h"

using namespace tthMEM;

typedef edm::ParameterSet PSet;
typedef std::vector<PSet> vPSet;

int
main(int argc,
     char * argv[])
{
// optimization procedures
//--- needed by std::chrono::
  setenv("TZ", "/etc/localtime", 1);
//--- untie std::cout from std::cin
  std::ios_base::sync_with_stdio(false);
  std::cin.tie(nullptr);
//--- flush before exit
  std::atexit(Logger::flush);
//--- point signal handler to std::exit() which also flushes the stdout
  for(int sig: { SIGABRT, SIGILL, SIGINT, SIGSEGV, SIGTERM, SIGQUIT })
    std::signal(sig, [](int sig) { LOGERR << "Encountered signal " << sig; std::exit(sig); });

//--- parse the configuration file
  if(argc != 2)
  {
    LOGERR << "Usage: " << argv[0] << " [parameters.py]";
    return EXIT_FAILURE;
  }

  const std::string ps = "process";
  if(! edm::readPSetsFrom(argv[1]) -> existsAs<PSet>(ps.c_str()))
    throw_line("tthMEM") << "No ParameterSet '" << ps << "' found "
                         << "in configuration file = " << argv[1];
  const PSet cfg = edm::readPSetsFrom(argv[1]) -> getParameter<PSet>(ps.c_str());

  const PSet cfg_log = cfg.getParameter<PSet>("logging");
  const std::string logLevel = cfg_log.getParameter<std::string>("logLevel");
  const bool enableLogging = cfg_log.getParameter<bool>("enableLogging");
  const bool enableTimeStamp = cfg_log.getParameter<bool>("enableTimeStamp");

  Logger::setLogLevel(logLevel);
  Logger::enableLogging(enableLogging);
  Logger::enableTimeStamp(enableTimeStamp);

  const PSet cfg_tthMEM = cfg.getParameter<PSet>("tthMEM");
  const bool isMC = cfg_tthMEM.getParameter<bool>("isMC");
  const std::string treeName = cfg_tthMEM.getParameter<std::string>("treeName");
  const std::string rleSelectionFileName = cfg_tthMEM.getParameter<std::string>("rleSelectionFile");
  const std::string pdfName = cfg_tthMEM.getParameter<std::string>("pdfName");
  const std::string madgraphFileName = cfg_tthMEM.getParameter<std::string>("madgraphFileName");
  const std::string integrationMode = cfg_tthMEM.getParameter<std::string>("integrationMode");
  const unsigned maxObjFunctionCalls = cfg_tthMEM.getParameter<unsigned>("maxObjFunctionCalls");
  const Long64_t startingFromEntry = cfg_tthMEM.getParameter<Long64_t>("startingFromEntry");
  const unsigned debugPlots = cfg_tthMEM.getParameter<unsigned>("debugPlots");
  const double higgsWidth = cfg_tthMEM.getParameter<double>("higgsWidth");
  const bool is2016 = cfg_tthMEM.getParameter<bool>("is2016");
  bool includeGeneratorLevel = [&]() -> bool
  {
    return cfg_tthMEM.getParameter<bool>("forceGenLevel") && is2016 && isMC;
  }();
  if(! rleSelectionFileName.empty() && ! boost::filesystem::exists(rleSelectionFileName))
  {
    LOGERR << "File '" << rleSelectionFileName << "' does not exists";
    return EXIT_FAILURE;
  }

//--- clamp the variables if needed
  VariableManager_3l1tau vm;
  const vPSet clampVariables = cfg_tthMEM.getParameter<vPSet>("clampVariables");
  for(const auto & cfg_clamp: clampVariables)
  {
    const std::string clampStr = cfg_clamp.getParameter<std::string>("var");
    const bool useGen = cfg_clamp.getParameter<bool>("useGen");
    const bool useCfg = cfg_clamp.getParameter<bool>("useCfg");
    const double clampValue = cfg_clamp.getParameter<double>("val");

    if(useGen && ! isMC)
    {
      LOGERR << "Generator level not accessible for variable '" << clampStr << "' "
             << "as the sample is not a Monte Carlo simulation but a data sample";
      return EXIT_FAILURE;
    }
    if(useGen && ! is2016)
    {
      LOGERR << "Generator level not available for MC samples produced earlier than 2016\n";
      return EXIT_FAILURE;
    }
    if(useGen && useCfg)
    {
      LOGERR << "Conflicting boolean values 'useGen' and 'useCfg' for variable '"
             << clampStr << "'  => pick one!";
      return EXIT_FAILURE;
    }

    if(useCfg)
    {
      if(vm.clamp(clampStr, clampValue)) return EXIT_FAILURE;
    }
    else if(useGen)
    {
      if(vm.clamp(clampStr)) return EXIT_FAILURE;
      includeGeneratorLevel |= true;
    }
  }
  LOGINFO << vm;

  LOGINFO << "PDF name: " << pdfName;
  LOGINFO << "MadGraph file name: " << madgraphFileName;
  LOGINFO << "Integation mode: " << integrationMode;
  LOGINFO << "Maximum number of calls per event: " << maxObjFunctionCalls;

//--- initialize the MEM instance
  LOGINFO << "Initializing the tth&z MEM instance";
  MEM_ttHorZ_3l1tau mem_tt_HandZ(pdfName, findFile(madgraphFileName), vm);
  mem_tt_HandZ.setIntegrationMode(integrationMode);
  mem_tt_HandZ.setMaxObjFunctionCalls(maxObjFunctionCalls);
  mem_tt_HandZ.setBJetTransferFunction(true);
  mem_tt_HandZ.useAvgBjetCombo(true);
  if(mem_tt_HandZ.isMarkovChainIntegrator())
  {
//--- retrieve the parameters for Markov Chain integrator
    const PSet cfg_mx = cfg_tthMEM.getParameter<PSet>("markovChainParams");
    const std::string mxMode = cfg_mx.getParameter<std::string>("mode");
    const unsigned nofBatches = cfg_mx.getParameter<unsigned>("nofBatches");
    const unsigned nofChains = cfg_mx.getParameter<unsigned>("nofChains");
    const unsigned maxCallsStartingPos = cfg_mx.getParameter<unsigned>("maxCallsStartingPos");
    const double epsilon0 = cfg_mx.getParameter<double>("epsilon0");
    const double T0 = cfg_mx.getParameter<double>("T0");
    const double nu = cfg_mx.getParameter<double>("nu");
    mem_tt_HandZ.setMarkovChainParams(mxMode, nofBatches, nofChains,
                                      maxCallsStartingPos, epsilon0, T0, nu);
  }
  if(higgsWidth > 0.)
    mem_tt_HandZ.setHiggsWidth(higgsWidth);

  const fwlite::InputSource inputFiles(cfg);
  const int maxEvents = inputFiles.maxEvents();

  const fwlite::OutputFiles outputFile(cfg);
  const std::string outputFileName = outputFile.file();

//--- create I/O TTrees
  TChain * inputTree = new TChain(treeName.c_str());
  for(const std::string & inputFile: inputFiles.files())
  {
    if(! boost::filesystem::exists(inputFile))
    {
      LOGERR << "Input file '" << inputFile << "' does not exist";
      return EXIT_FAILURE;
    }
    inputTree -> AddFile(inputFile.c_str());
    LOGINFO << "Chained file = " << inputFile;
  }
  inputTree -> LoadTree(0);
  LOGINFO << "Loaded tree '" << treeName << '\'';

  TFile * newFile = new TFile(outputFileName.c_str(), "recreate");
  TTree * newTree = new TTree("tree", Form("Tree created by %s", argv[0]));

  MeasuredEvent_3l1tau measuredEvent;
  measuredEvent.addFilter(rleSelectionFileName);
  measuredEvent.includeGeneratorLevel(includeGeneratorLevel);
  measuredEvent.setBranches(inputTree);
  measuredEvent.initNewBranches(newTree);
  if(debugPlots)
    measuredEvent.debugPlotter = new DebugPlotter_ttHorZ_3l1tau(newFile, debugPlots);

//--- set up the probability variables
  double probSignal,         probSignalErr;
//  double probSignal_th2ww,   probSignalErr_th2ww;
  double probBackground_ttz, probBackgroundErr_ttz;
  double lhRatioNP, lhRatioNP_up, lhRatioNP_down;

  TBranch * probSignalBranch             __attribute__((unused)) =
    newTree -> Branch("probSignal",             &probSignal,            "probSignal/D");
  TBranch * probSignalBranch_err         __attribute__((unused)) =
    newTree -> Branch("probSignal_err",         &probSignalErr,         "probSignal_err/D");
  TBranch * probBackgroundBranch_ttz     __attribute__((unused)) =
    newTree -> Branch("probBackground_ttz",     &probBackground_ttz,    "probBackground_ttz/D");
  TBranch * probBackgroundBranch_ttz_err __attribute__((unused)) =
    newTree -> Branch("probBackground_ttz_err", &probBackgroundErr_ttz, "probBackground_ttz_err/D");
//  TBranch * probSignalBranch_th2ww       __attribute__((unused)) =
//    newTree -> Branch("probSignal_tth2ww",      &probSignal_th2ww,      "probSignal_tth2ww/D");
//  TBranch * probSignalBranch_th2ww_err   __attribute__((unused)) =
//    newTree -> Branch("probSignal_tth2ww_err", &probSignalErr_th2ww,    "probSignal_tth2ww_err/D");
  TBranch * lhRatioNPBranch              __attribute__((unused)) =
    newTree -> Branch("lhRatioNP",              &lhRatioNP,             "lhRatioNP/D");
  TBranch * lhRatioNPBranch_up           __attribute__((unused)) =
    newTree -> Branch("lhRatioNP_up",           &lhRatioNP_up,          "lhRatioNP_up/D");
  TBranch * lhRatioNPBranch_down         __attribute__((unused)) =
    newTree -> Branch("lhRatioNP_down",         &lhRatioNP_down,        "lhRatioNP_down/D");

//--- set up the weights in the denominator of Neyman-Pearson likelihood ratio
  const std::vector<double> bkgWeightNumerator{ constants::xSectionTTZ };
  const std::vector<double> sigWeightNumerator{ constants::xSectionTTH /*, constants::xSectionTTH2diW */ };
  const double bkgWeightDenom = 1. / std::accumulate(bkgWeightNumerator.begin(), bkgWeightNumerator.end(), 0.);
  const double sigWeightDenom = 1. / std::accumulate(sigWeightNumerator.begin(), sigWeightNumerator.end(), 0.);
  std::vector<double> bkgWeights, sigWeights;
  std::transform(
    bkgWeightNumerator.begin(), bkgWeightNumerator.end(), std::back_inserter(bkgWeights),
    [bkgWeightDenom](double weightNumerator) -> double { return weightNumerator * bkgWeightDenom; }
  );
  std::transform(
    sigWeightNumerator.begin(), sigWeightNumerator.end(), std::back_inserter(sigWeights),
    [sigWeightDenom](double weightNumerator) -> double { return weightNumerator * sigWeightDenom; }
  );

//--- start looping over the events
  const Long64_t nof_tree_entries = inputTree -> GetEntries();
  const Long64_t nof_max_entries = maxEvents < 0 ? nof_tree_entries : (startingFromEntry + maxEvents);
  if(nof_max_entries > nof_tree_entries)
  {
    LOGERR << "The requested number of entries to be processed (= "
           << maxEvents << ") starting from entry = "
           << startingFromEntry << " is greater than the total number entries (= "
           << nof_tree_entries << " in the input file(s). Aborting";
    return EXIT_FAILURE;
  }
  LOGINFO << "Processing " << (nof_max_entries - startingFromEntry) << " entries "
          << "(out of " << nof_tree_entries << "), starting from "
          << startingFromEntry << " (event range: [" << startingFromEntry << "; "
          << nof_max_entries << ") )";

  for(Long64_t i = startingFromEntry; i < nof_max_entries; ++i)
  {
    LOGINFO << "Processing " << i << "th event";

    inputTree -> GetEntry(i);
    measuredEvent.initialize();
    if(measuredEvent.isFiltered()) continue;
    LOGINFO << "run:lumi:event = " << measuredEvent.str(false);

    probSignal = 0.;
    probBackground_ttz = 0.;
//    probBackground_th2ww = 0.;

    bool err = false;
    std::array<double, 2> probSignalResult, probBackgroundResult_ttz;
    probSignalResult = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTH, err);
    if(err)
    {
      LOGWARN << "Skipping the event because of errors";
      continue;
    }
    probBackgroundResult_ttz = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTZ, err);
    if(err)
    {
      LOGWARN << "Skipping the event because of errors";
      continue;
    }
//    probSignalResult_th2ww = ...

    probSignal            = probSignalResult[0];
    probSignalErr         = probSignalResult[1];
//    probSignal_th2ww      = probSignalResult_th2ww[0];
//    probSignalErr_th2ww   = probSignalResult_th2ww[1];
    probBackground_ttz    = probBackgroundResult_ttz[0];
    probBackgroundErr_ttz = probBackgroundResult_ttz[1];

    const std::vector<double> probs_signal{
      probSignal//, probSignal_th2ww
    };
    const std::vector<double> probs_signal_err{
      probSignalErr//, probSignalErr_th2ww
    };
    const std::vector<double> prob_background{
      probBackground_ttz
    };
    const std::vector<double> prob_background_err{
      probBackgroundErr_ttz
    };

    const double signal_weighted = std::inner_product(
      probs_signal.begin(), probs_signal.end(), sigWeights.begin(), 0.
    );
    const double signal_weighted_err = std::inner_product(
      probs_signal_err.begin(), probs_signal_err.end(), sigWeights.begin(), 0.
    );
    const double background_weighted = std::inner_product(
      prob_background.begin(), prob_background.end(), bkgWeights.begin(), 0.
    );
    const double background_weighted_err = std::inner_product(
      prob_background_err.begin(), prob_background_err.end(), bkgWeights.begin(), 0.
    );

    const double sig      = signal_weighted;
    const double sig_up   = sig + signal_weighted_err;
    const double sig_down = sig - signal_weighted_err;
    const double bkg      = background_weighted;
    const double bkg_up   = bkg + background_weighted_err;
    const double bkg_down = bkg - background_weighted_err;

    lhRatioNP      = (sig      + bkg)      != 0. ? sig / (sig + bkg)              : 0.;
    lhRatioNP_up   = (sig_up   + bkg_down) != 0. ? sig_up / (sig_up + bkg_down)   : 0.;
    lhRatioNP_down = (sig_down + bkg_up)   != 0. ? sig_down / (sig_down + bkg_up) : 0.;

    newTree -> Fill();
  }
  LOGINFO << "Average time spent on tth & ttz MEM per event: "
          << "Real time: " << mem_tt_HandZ.getAverageComputingTime_real() << " s;  "
          << "CPU time: " << mem_tt_HandZ.getAverageComputingTime_cpu() << " s";

  newFile -> Write();

  LOGINFO << "Wrote results to file = " << outputFileName;
  LOGINFO << "Done";

  return EXIT_SUCCESS;
}
