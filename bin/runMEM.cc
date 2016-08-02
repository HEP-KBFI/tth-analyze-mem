#include <iostream> // std::cerr
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE
#include <string> // std::string
#include <vector> // std::vector<>

#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h" // edm::readPSetsFrom()
#include "FWCore/Utilities/interface/Exception.h" // cms::Exception
#include "DataFormats/FWLite/interface/InputSource.h" // fwlite::InputSource
#include "DataFormats/FWLite/interface/OutputFiles.h" // fwlite::OutputFiles

#include "boost/algorithm/string/predicate.hpp" // boost::iequals()

#include <Rtypes.h> // Long64_t
#include <TFile.h> // TFile
#include <TChain.h> // TChain
#include <TTree.h> // TTree
#include <TString.h> // Form

#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/MeasuredEvent.h" // tthMEM::MeasuredEvent

int
main(int argc,
     char * argv[])
{
  if(argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " [parameters.py]\n";
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

  tthMEM::Logger::setLogLevel(logLevel);
  tthMEM::Logger::enableLogging(enableLogging);
  tthMEM::Logger::enableTimeStamp(enableTimeStamp);

  const edm::ParameterSet cfg_tthMEM = cfg.getParameter<edm::ParameterSet>("tthMEM");
  const std::string treeName = cfg_tthMEM.getParameter<std::string>("treeName");
  const std::string pdfName = cfg_tthMEM.getParameter<std::string>("pdfName");
  const std::string madgraphFileName = cfg_tthMEM.getParameter<std::string>("madgraphFileName");
  const std::string integrationMode = cfg_tthMEM.getParameter<std::string>("integrationMode");
  //const unsigned maxObjFunctionCalls = cfg_tthMEM.getParameter<unsigned>("maxObjFunctionCalls");

  const fwlite::InputSource inputFiles(cfg);
  const int maxEvents = inputFiles.maxEvents();
  const unsigned reportEvery = inputFiles.reportAfter();

  const fwlite::OutputFiles outputFile(cfg);
  const std::string outputFileName = outputFile.file();

  TChain * inputTree = new TChain(treeName.c_str());
  for(const std::string & inputFile: inputFiles.files())
    inputTree -> AddFile(inputFile.c_str());
  inputTree -> LoadTree(0);

  TFile * newFile = new TFile(outputFileName.c_str(), "recreate");
  TTree * newTree = new TTree("tree", Form("Tree created by %s", argv[0]));

  tthMEM_3l_1tau::MeasuredEvent measuredEvent;
  measuredEvent.setBranches(inputTree);
  measuredEvent.initNewBranches(newTree);

  double probSignal;
  double probBackground;

  TBranch * probSignalBranch     = newTree -> Branch(
        "probSignal",     &probSignal,     "probSignal/D");
  TBranch * probBackgroundBranch = newTree -> Branch(
        "probBackground", &probBackground, "probBackground/D");
  (void) probSignalBranch;     // prevents compilation error
  (void) probBackgroundBranch; // prevents compilation error

  const Long64_t nof_tree_entries = inputTree -> GetEntries();
  const Long64_t nof_entries = maxEvents < 0 ? nof_tree_entries : maxEvents;
  LOGINFO << "Processing " << nof_entries << " entries (out of " << nof_tree_entries << ")";

  for(Long64_t i = 0; i < nof_entries; ++i)
  {
    if(i > 0 && ((i + 1) % reportEvery) == 0)
      LOGINFO << "Processing " << i << "th event";

    inputTree -> GetEntry(i);
    measuredEvent.initialize();

    probSignal = -1;
    probBackground = -1;

    newTree -> Fill();
  }

  newFile -> Write();

  LOGINFO << "Done";

  return EXIT_SUCCESS;
}
