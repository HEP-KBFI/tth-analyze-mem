#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"

#include <cmath> // std::pow(), std::ceil(), std::log10(), std::fabs(), std::round()

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
}
