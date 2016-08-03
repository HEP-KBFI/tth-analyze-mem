#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

#include <cmath> // std::pow(), std::ceil(), std::log10(), std::fabs(), std::round()
#include <cassert> // assert()

#include "FWCore/ParameterSet/interface/FileInPath.h" // edm::FileInPath

namespace tthMEM
{
  double
  roundToNdigits(double x,
                 int n)
  {
    if(x == 0.) return 0.;
    const double p = std::pow(10., n - std::ceil(std::log10(std::fabs(x))));
    return std::round(p * x) / p;
  }

  std::string
  findFile(const std::string & fileName)
  {
    edm::FileInPath inputFile(fileName);
    if(inputFile.fullPath() == "")
    {
      LOGERR << "Error: cannot find file = " << fileName;
      assert(0);
    }
    return inputFile.fullPath();
  }
}
