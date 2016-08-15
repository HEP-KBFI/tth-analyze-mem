#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <string> // std::string
#include <vector> // std::vector<>

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
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // tthMEM::findFile(), tthMEM::sqrtS

using namespace tthMEM;

int
main(int argc,
     char * argv[])
{
//--- parse the configuration file
  if(argc != 2)
  {
    LOGERR << "Usage: " << argv[0] << " [parameters.py]";
    return EXIT_FAILURE;
  }

  const std::string ps = "process";
  if(!edm::readPSetsFrom(argv[1]) -> existsAs<edm::ParameterSet>(ps.c_str()))
    throw cms::Exception("tthMEM")
      << "No ParameterSet '" << ps << "' found in configuration file = " << argv[1] << "\n";
  const edm::ParameterSet cfg = edm::readPSetsFrom(argv[1]) -> getParameter<edm::ParameterSet>(ps.c_str());

  const edm::ParameterSet cfg_log = cfg.getParameter<edm::ParameterSet>("logging");
  const std::string logLevel = cfg_log.getParameter<std::string>("logLevel");
  const bool enableLogging = cfg_log.getParameter<bool>("enableLogging");
  const bool enableTimeStamp = cfg_log.getParameter<bool>("enableTimeStamp");

  Logger::setLogLevel(logLevel);
  Logger::enableLogging(enableLogging);
  Logger::enableTimeStamp(enableTimeStamp);

  const edm::ParameterSet cfg_tthMEM = cfg.getParameter<edm::ParameterSet>("tthMEM");
  const std::string treeName = cfg_tthMEM.getParameter<std::string>("treeName");
  const std::string pdfName = cfg_tthMEM.getParameter<std::string>("pdfName");
  const std::string madgraphFileName = cfg_tthMEM.getParameter<std::string>("madgraphFileName");
  const std::string integrationMode = cfg_tthMEM.getParameter<std::string>("integrationMode");
  const unsigned maxObjFunctionCalls = cfg_tthMEM.getParameter<unsigned>("maxObjFunctionCalls");
  const Long64_t startingFromEntry = cfg_tthMEM.getParameter<Long64_t>("startingFromEntry");

  LOGINFO << "PDF name: " << pdfName;
  LOGINFO << "MadGraph file name: " << madgraphFileName;
  LOGINFO << "Integation mode: " << integrationMode;
  LOGINFO << "Maximum number of calls per event: " << maxObjFunctionCalls;

  const fwlite::InputSource inputFiles(cfg);
  const int maxEvents = inputFiles.maxEvents();
  const unsigned reportEvery = inputFiles.reportAfter();

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

//--- set up the probability variables
  double probSignal;
  double probBackground_ttz;

  TBranch * probSignalBranch         = newTree -> Branch(
        "probSignal",         &probSignal,         "probSignal/D");
  TBranch * probBackgroundBranch_ttz = newTree -> Branch(
        "probBackground_ttz", &probBackground_ttz, "probBackground_ttz/D");
  (void) probSignalBranch;         // prevents compilation error
  (void) probBackgroundBranch_ttz; // prevents compilation error

//--- initialize the MEM instance and start looping over the events
  LOGINFO << "Initializing the tth&z MEM instance";
  MEM_ttHorZ_3l1tau mem_tt_HandZ(sqrtS, pdfName, findFile(madgraphFileName));
  mem_tt_HandZ.setIntegrationMode(integrationMode);
  mem_tt_HandZ.setMaxObjFunctionCalls(maxObjFunctionCalls);

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
    if(i > 0 && (i % reportEvery) == 0)
      LOGINFO << "Processing " << i << "th event";

    inputTree -> GetEntry(i);
    measuredEvent.initialize();

    probSignal = -1;
    probBackground_ttz = -1;

    probSignal = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTH);
    probBackground_ttz = mem_tt_HandZ.integrate(measuredEvent, ME_mg5_3l1tau::kTTZ);

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
