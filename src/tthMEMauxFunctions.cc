#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

#include <cfenv> // std::fesetround(), FE_TONEAREST
#include <stdexcept> // std::invalid_argument

#include <boost/range/algorithm/lexicographical_compare.hpp> // ...
  // ... boost::range::lexicographical_compare()
#include <boost/algorithm/string/compare.hpp> // boost::is_iless()

#include <FWCore/ParameterSet/interface/FileInPath.h> // edm::FileInPath

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

  unsigned
  roundToNearestUInt(double x)
  {
    std::fesetround(FE_TONEAREST);
    return static_cast<unsigned>(std::rint(x));
  }

  void
  setMGmomentum(const LorentzVector & lv,
                double * mg)
  {
    mg[0] = lv.e();
    mg[1] = lv.px();
    mg[2] = lv.py();
    mg[3] = lv.pz();
  }

  std::string
  findFile(const std::string & fileName)
  {
    edm::FileInPath inputFile(fileName);
    if(inputFile.fullPath() == "")
    {
      LOGERR << "Error: cannot find file = " << fileName;
      throw std::invalid_argument(__PRETTY_FUNCTION__);
    }
    return inputFile.fullPath();
  }

  bool
  iStrComparator::operator ()(const std::string & lhs,
                              const std::string & rhs) const
  {
    return boost::range::lexicographical_compare(lhs, rhs, boost::is_iless());
  }
}
