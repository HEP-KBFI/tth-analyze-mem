#include <iostream> // std::cerr
#include <cstdlib> // EXIT_SUCCESS, EXIT_FAILURE

#include "FWCore/ParameterSet/interface/ParameterSet.h" // edm::ParameterSet
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h" // edm::readPSetsFrom()
#include "FWCore/Utilities/interface/Exception.h" // cms::Exception

#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*

using namespace tthMEM;

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
  edm::ParameterSet cfg = edm::readPSetsFrom(argv[1]) -> getParameter<edm::ParameterSet>(ps.c_str());

  edm::ParameterSet cfg_log = cfg.getParameter<edm::ParameterSet>("logging");
  const std::string logLevel = cfg_log.getParameter<std::string>("logLevel");
  const bool enableLogging = cfg_log.getParameter<bool>("enableLogging");
  const bool enableTimeStamp = cfg_log.getParameter<bool>("enableTimeStamp");

  Logger::setLogLevel(logLevel);
  Logger::enableLogging(enableLogging);
  Logger::enableTimeStamp(enableTimeStamp);

  LOGERR << "testing error message " << 1;
  LOGWARN << "testing warning message " << 2;
  LOGINFO << "testing info message " << 3;
  LOGDBG << "testing debug message " << 4;
  return EXIT_SUCCESS;
}
