#include "tthAnalysis/tthMEM/interface/general/auxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()

#include <FWCore/ParameterSet/interface/FileInPath.h> // edm::FileInPath

#include <boost/range/algorithm/lexicographical_compare.hpp> // ...
  // ... boost::range::lexicographical_compare()
#include <boost/algorithm/string/compare.hpp> // boost::is_iless()

#include <cfenv> // std::fesetround(), FE_TONEAREST

namespace tthMEM
{
  double
  roundToNdigits(double x,
                 int n)
  {
    if(std::fpclassify(x) == FP_ZERO) return 0.;
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
      throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_FILE_NOT_FOUND)
        << "Error: cannot find file = " << fileName;
    return inputFile.fullPath();
  }

  bool
  iStrComparator::operator()(const std::string & lhs,
                             const std::string & rhs) const
  {
    return boost::range::lexicographical_compare(lhs, rhs, boost::is_iless());
  }
}
