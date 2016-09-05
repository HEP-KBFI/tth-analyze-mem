#include <cstdlib> // EXIT_FAILURE, EXIT_SUCCESS
#include <string> // std::string
#include <vector> // std::vector<>

#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR, LOGINFO
#include "tthAnalysis/tthMEM/interface/ROC.h" // tthMEM::ROC

#include <FWCore/ParameterSet/interface/ParameterSet.h> // edm::ParameterSet
#include <FWCore/PythonParameterSet/interface/MakeParameterSets.h> // edm::readPSetsFrom()
#include <FWCore/Utilities/interface/Exception.h> // cms::Exception

int
main(int argc,
     char * argv[])
{
  if(argc != 2)
  {
    LOGERR << "Usage: " << argv[0] << " [parameters.py]";
    return EXIT_FAILURE;
  }

  const std::string rocStr = "roc";
  if(! edm::readPSetsFrom(argv[1]) -> existsAs<edm::ParameterSet>(rocStr.c_str()))
    throw_line("tthMEM") << "No ParameterSet '" << rocStr << "' found in "
                         << "configuration file = " << argv[1];
  const edm::ParameterSet cfg = edm::readPSetsFrom(argv[1]) -> getParameter<edm::ParameterSet>(rocStr.c_str());

  const edm::ParameterSet cfg_tthMEM = cfg.getParameter<edm::ParameterSet>("tthMEM");
  const std::string treeName = cfg_tthMEM.getParameter<std::string>("treeName");
  const std::string branchName = cfg_tthMEM.getParameter<std::string>("branchName");

  std::vector<double> legPos;
  if(cfg_tthMEM.exists("legendPosition"))
    legPos = cfg_tthMEM.getParameter<std::vector<double>>("legendPosition");
  if(legPos.size() != 4)
  {
    LOGWARN << "Invalid dimension of legendPosition: is " << legPos.size() << ", "
            << "but must be 4. Resetting legend position";
    legPos.clear();
  }

  typedef std::vector<std::string> vstring;
  const std::string signalFileName = cfg_tthMEM.getParameter<std::string>("signalFile");
  const vstring bkgFileNames = cfg_tthMEM.getParameter<vstring>("bkgFiles");
  if(! bkgFileNames.size())
  {
    LOGERR << "No background files given";
    return EXIT_FAILURE;
  }
  const std::string outputFolderName = cfg_tthMEM.getParameter<std::string>("outFolder");
  const vstring labels = cfg_tthMEM.getParameter<vstring>("labels");
  if(labels.size() != bkgFileNames.size() + 1)
  {
    LOGERR << "Number of labels does not match to the number of input files";
    return EXIT_FAILURE;
  }

  tthMEM::ROC roc(
   signalFileName,  bkgFileNames, outputFolderName,
   treeName, branchName, labels
  );
  if(legPos.size()) roc.setLegendPosition(legPos);
  if(!! roc.plotROC())
  {
    LOGERR << "Not enough input given to draw the ROC curve";
    return EXIT_FAILURE;
  }
  for(unsigned i = 1; i < labels.size(); ++i)
    LOGINFO << "MEM ROC AUC (" << labels[0] << " vs "
            << labels[i] << ") = " << roc.getAUC(i - 1);

  return EXIT_SUCCESS;
}
