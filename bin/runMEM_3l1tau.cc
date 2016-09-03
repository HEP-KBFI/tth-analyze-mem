#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE, setenv(), std::atexit()
#include <string> // std::string
#include <vector> // std::vector<>
#include <algorithm> // std::inner_product()
#include <csignal> // std::signal(), SIG*

#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h" // edm::readPSetsFrom()
#include "FWCore/Utilities/interface/Exception.h" // cms::Exception
#include "DataFormats/FWLite/interface/InputSource.h" // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h" // fwlite::OutputFiles

#include <Rtypes.h> // Long64_t
#include <TFile.h> // TFile
#include <TChain.h> // TChain
#include <TTree.h> // TTree
#include <TString.h> // Form()

#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h" // tthMEM::MeasuredEvent_3l1tau
#include "tthAnalysis/tthMEM/interface/MEM_ttHorZ_3l1tau.h" // tthMEM::MEM_ttHorZ_3l1tau
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // tthMEM::findFile()
#include "tthAnalysis/tthMEM/interface/DebugPlotter_ttHorZ_3l1tau.h" // DebugPlotter_ttHorZ_3l1tau
#include "tthAnalysis/tthMEM/interface/VariableManager_3l1tau.h" // VariableManager_3l1tau

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
  if(!edm::readPSetsFrom(argv[1]) -> existsAs<PSet>(ps.c_str()))
    throw cms::Exception("tthMEM")
      << "No ParameterSet '" << ps << "' found in configuration file = " << argv[1] << "\n";
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
  const std::string pdfName = cfg_tthMEM.getParameter<std::string>("pdfName");
  const std::string madgraphFileName = cfg_tthMEM.getParameter<std::string>("madgraphFileName");
  const std::string integrationMode = cfg_tthMEM.getParameter<std::string>("integrationMode");
  const unsigned maxObjFunctionCalls = cfg_tthMEM.getParameter<unsigned>("maxObjFunctionCalls");
  const Long64_t startingFromEntry = cfg_tthMEM.getParameter<Long64_t>("startingFromEntry");
  const unsigned debugPlots = cfg_tthMEM.getParameter<unsigned>("debugPlots");

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
    }
  }
  LOGINFO << vm;

  LOGINFO << "PDF name: " << pdfName;
  LOGINFO << "MadGraph file name: " << madgraphFileName;
  LOGINFO << "Integation mode: " << integrationMode;
  LOGINFO << "Maximum number of calls per event: " << maxObjFunctionCalls;

//--- initialize the MEM instance
  LOGINFO << "Initializing the tth&z MEM instance";
  MEM_ttHorZ_3l1tau mem_tt_HandZ(pdfName, findFile(madgraphFileName), std::move(vm));
  mem_tt_HandZ.setIntegrationMode(integrationMode);
  mem_tt_HandZ.setMaxObjFunctionCalls(maxObjFunctionCalls);
  mem_tt_HandZ.setBJetTransferFunction(true);
  if(mem_tt_HandZ.isMarkovChainIntegrator())
  {
//--- retrieve the parameters for Markov Chain integrator
    const PSet cfg_mx = cfg_tthMEM.getParameter<PSet>("markovChainParams");
    const unsigned nofBatches = cfg_mx.getParameter<unsigned>("nofBatches");
    const unsigned nofChains = cfg_mx.getParameter<unsigned>("nofChains");
    const double epsilon0 = cfg_mx.getParameter<double>("epsilon0");
    const double T0 = cfg_mx.getParameter<double>("T0");
    const double nu = cfg_mx.getParameter<double>("nu");
    mem_tt_HandZ.setMarkovChainParams(nofBatches, nofChains, epsilon0, T0, nu);
  }

  const fwlite::InputSource inputFiles(cfg);
  const int maxEvents = inputFiles.maxEvents();

  const fwlite::OutputFiles outputFile(cfg);
  const std::string outputFileName = outputFile.file();

//--- create I/O TTrees
  TChain * inputTree = new TChain(treeName.c_str());
  for(const std::string & inputFile: inputFiles.files())
  {
    inputTree -> AddFile(inputFile.c_str());
    LOGINFO << "Chained file = " << inputFile;
  }
  inputTree -> LoadTree(0);
  LOGINFO << "Loaded tree '" << treeName << "'";

  TFile * newFile = new TFile(outputFileName.c_str(), "recreate");
  TTree * newTree = new TTree("tree", Form("Tree created by %s", argv[0]));

  MeasuredEvent_3l1tau measuredEvent;
  measuredEvent.setBranches(inputTree);
  measuredEvent.initNewBranches(newTree);
  if(debugPlots)
    measuredEvent.debugPlotter = new DebugPlotter_ttHorZ_3l1tau(newFile, debugPlots);

//--- set up the probability variables
  double probSignal;
  double probBackground_ttz;
//  double probBackground_th2ww;
  double lhRatioNP;

  TBranch * probSignalBranch           __attribute__((unused)) =
    newTree -> Branch("probSignal",            &probSignal,           "probSignal/D");
  TBranch * probBackgroundBranch_ttz   __attribute__((unused)) =
    newTree -> Branch("probBackground_ttz",    &probBackground_ttz,   "probBackground_ttz/D");
//  TBranch * probBackgroundBranch_th2ww __attribute__((unused)) =
//    newTree -> Branch("probBackground_tth2ww", &probBackground_th2ww, "probBackground_tth2ww/D");
  TBranch * lhRatioNPBranch            __attribute__((unused)) =
    newTree -> Branch("lhRatioNP",             &lhRatioNP,            "lhRatioNP/D");

//--- set up the weights in the denominator of Neyman-Pearson likelihood ratio
//  const double bkgWeightDenom = 1. / (constants::xSectionTTZ + constants::xSectionTTH2diW);
  const double bkgWeightDenom = 1. / constants::xSectionTTZ;
  const std::vector<double> denomWeights = {
    1., constants::xSectionTTZ * bkgWeightDenom//, constants::xSectionTTH2diW * bkgWeightDenom
  };

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

    probSignal = 0.;
    probBackground_ttz = 0.;
//    probBackground_th2ww = 0.;

    probSignal = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTH);
    probBackground_ttz = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTZ);

    const std::vector<double> probs = {
      probSignal, probBackground_ttz//, probBackground_th2ww
    };
    const double denom = std::inner_product(
          probs.begin(), probs.end(), denomWeights.begin(), 0.
    );
    lhRatioNP = denom != 0. ? probSignal / denom : 0.;

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
