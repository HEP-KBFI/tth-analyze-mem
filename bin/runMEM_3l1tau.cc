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
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()
#include "tthAnalysis/tthMEM/interface/LikelihoodRatio_3l1tau.h" // LikelihoodRatio_3l1tau, MEMOutput_3l1tau

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
    std::signal(sig, [](int sig_) { LOGERR << "Encountered signal " << sig_; std::exit(sig_); });

//--- parse the configuration file
  if(argc != 2)
  {
    LOGERR << "Usage: " << argv[0] << " [parameters.py]";
    return EXIT_FAILURE;
  }

  const std::string ps = "process";
  if(! edm::readPSetsFrom(argv[1]) -> existsAs<PSet>(ps.c_str()))
    throw_line_ext("tthMEM", TTHEXCEPTION_ERR_CODE_INVALID_PARAMETERSET)
      << "No ParameterSet '" << ps << "' found in configuration file = " << argv[1];
  const PSet cfg = edm::readPSetsFrom(argv[1]) -> getParameter<PSet>(ps.c_str());

  const PSet cfg_log = cfg.getParameter<PSet>("logging");
  const std::string logLevel = cfg_log.getParameter<std::string>("logLevel");
  const bool enableLogging = cfg_log.getParameter<bool>("enableLogging");
  const bool enableTimeStamp = cfg_log.getParameter<bool>("enableTimeStamp");

  Logger::setLogLevel(logLevel);
  Logger::enableLogging(enableLogging);
  Logger::enableTimeStamp(enableTimeStamp);
  Logger::setFloatPrecision(5);

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
  double probSignal_th2ww,   probSignalErr_th2ww;
  double probBackground_ttz, probBackgroundErr_ttz;
  double lhRatioNP, lhRatioNP_err, lhRatioNP_up, lhRatioNP_down;

  TBranch * probSignalBranch             __attribute__((unused)) =
    newTree -> Branch("probSignal",             &probSignal,            "probSignal/D");
  TBranch * probSignalBranch_err         __attribute__((unused)) =
    newTree -> Branch("probSignal_err",         &probSignalErr,         "probSignal_err/D");
  TBranch * probBackgroundBranch_ttz     __attribute__((unused)) =
    newTree -> Branch("probBackground_ttz",     &probBackground_ttz,    "probBackground_ttz/D");
  TBranch * probBackgroundBranch_ttz_err __attribute__((unused)) =
    newTree -> Branch("probBackground_ttz_err", &probBackgroundErr_ttz, "probBackground_ttz_err/D");
  TBranch * probSignalBranch_th2ww       __attribute__((unused)) =
    newTree -> Branch("probSignal_tth2ww",      &probSignal_th2ww,      "probSignal_tth2ww/D");
  TBranch * probSignalBranch_th2ww_err   __attribute__((unused)) =
    newTree -> Branch("probSignal_tth2ww_err", &probSignalErr_th2ww,    "probSignal_tth2ww_err/D");
  TBranch * lhRatioNPBranch              __attribute__((unused)) =
    newTree -> Branch("lhRatioNP",              &lhRatioNP,             "lhRatioNP/D");
  TBranch * lhRatioNPBranch_err          __attribute__((unused)) =
    newTree -> Branch("lhRatioNP_err",          &lhRatioNP_err,         "lhRatioNP_err/D");
  TBranch * lhRatioNPBranch_up           __attribute__((unused)) =
    newTree -> Branch("lhRatioNP_up",           &lhRatioNP_up,          "lhRatioNP_up/D");
  TBranch * lhRatioNPBranch_down         __attribute__((unused)) =
    newTree -> Branch("lhRatioNP_down",         &lhRatioNP_down,        "lhRatioNP_down/D");

//--- log also time and Markov Chain calls before integration
  double realTime_tth, cpuTime_tth,
         realTime_ttz, cpuTime_ttz;
  long long nofMXMCTries_tth, nofMXMCTries_ttz;
  TBranch * realTimeBranch_tth           __attribute__((unused)) =
    newTree -> Branch("realTime_tth",           &realTime_tth,          "realTime_tth/D");
  TBranch * cpuTimeBranch_tth            __attribute__((unused)) =
    newTree -> Branch("cpuTime_tth",            &cpuTime_tth,           "cpuTime_tth/D");
  TBranch * realTimeBranch_ttz           __attribute__((unused)) =
    newTree -> Branch("realTime_ttz",           &realTime_ttz,          "realTime_ttz/D");
  TBranch * cpuTimeBranch_ttz            __attribute__((unused)) =
    newTree -> Branch("cpuTime_ttz",            &cpuTime_ttz,           "cpuTime_ttz/D");
  TBranch * nofTriesBranch_tth           __attribute__((unused)) =
    newTree -> Branch("nofMXMCTries_tth",       &nofMXMCTries_tth,      "nofMXMCTries_tth/L");
  TBranch * nofTriesBranch_ttz           __attribute__((unused)) =
    newTree -> Branch("nofMXMCTries_ttz",       &nofMXMCTries_ttz,      "nofMXMCTries_ttz/L");

//--- log error code
  Int_t errCode;
  TBranch * errCodeBranch                __attribute__((unused)) =
    newTree -> Branch("errCode",                &errCode,               "errCode/I");

  LikelihoodRatio_3l1tau lr_computation(
    { LikelihoodRatio_3l1tau::Hypothesis::tth }, // signal hypotheses
    { LikelihoodRatio_3l1tau::Hypothesis::ttz }  // background hypotheses
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
    errCode = 0;

    inputTree -> GetEntry(i);
    measuredEvent.initialize();
    if(measuredEvent.isFiltered()) continue;
    LOGINFO << "run:lumi:event = " << measuredEvent.str(false);

    MEMOutput_3l1tau result;
    try
    {
      std::array<double, 2> probSignalResult;
      probSignalResult = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTH);
      realTime_tth     = mem_tt_HandZ.getComputingTime_real();
      cpuTime_tth      = mem_tt_HandZ.getComputingTime_cpu();
      nofMXMCTries_tth = mem_tt_HandZ.getNofMXMCTries();

      std::array<double, 2> probBackgroundResult_ttz;
      probBackgroundResult_ttz = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTZ);
      realTime_ttz     = mem_tt_HandZ.getComputingTime_real();
      cpuTime_ttz      = mem_tt_HandZ.getComputingTime_cpu();
      nofMXMCTries_ttz = mem_tt_HandZ.getNofMXMCTries();

      result = lr_computation.compute(
        { { LikelihoodRatio_3l1tau::Hypothesis::tth, probSignalResult         } },
        { { LikelihoodRatio_3l1tau::Hypothesis::ttz, probBackgroundResult_ttz } }
      );
    } catch(const tthMEMexception & exception)
    {
      errCode = exception.getErrCode();
    } catch(...)
    {
      errCode = TTHEXCEPTION_ERR_CODE_DEFAULT;
    }

    probSignal            = result.prob_tth;
    probSignalErr         = result.prob_tth_err;
    probSignal_th2ww      = result.prob_tth_h2ww;
    probSignalErr_th2ww   = result.prob_tth_h2ww_err;
    probBackground_ttz    = result.prob_ttz;
    probBackgroundErr_ttz = result.prob_ttz_err;
    lhRatioNP             = result.lr;
    lhRatioNP_err         = result.lr_err;
    lhRatioNP_up          = result.lr_up;
    lhRatioNP_down        = result.lr_down;

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
