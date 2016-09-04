#include "tthAnalysis/tthMEM/interface/tthMEMvecFunctions.h"

#include <cmath> // std::sqrt()

namespace tthMEM
{
  namespace vec
  {
    double
    l2(const std::vector<double> & vec,
       double shiftValue)
    {
      return l2(vec, 0, vec.size(), shiftValue);
    }

    double
    l2(const std::vector<double> & vec,
       unsigned shiftFromBegin,
       double shiftValue)
    {
      if(! (shiftFromBegin <= vec.size()))
        throw std::range_error(__PRETTY_FUNCTION__);
      return l2(vec, 0, shiftFromBegin, shiftValue);
    }

    double
    l2(const std::vector<double> & vec,
       unsigned shiftFromBegin_begin,
       unsigned shiftFromBegin_end,
       double shiftValue)
    {
      if(! (shiftFromBegin_begin <= shiftFromBegin_end &&
            shiftFromBegin_begin <= vec.size() &&
            shiftFromBegin_end <= vec.size()))
        throw std::range_error(__PRETTY_FUNCTION__);

      return std::sqrt(std::accumulate(
        vec.begin() + shiftFromBegin_begin, vec.begin() + shiftFromBegin_end, 0.,
        [shiftValue](double sum,
                     double element) -> double
        {
          const double eShifted = (element - shiftValue);
          return sum + eShifted * eShifted;
        }
      ));
    }

    double
    avg(const std::vector<double> & vec)
    {
      return avg(vec, 0, vec.size());
    }

    double
    avg(const std::vector<double> & vec,
        unsigned shiftFromBegin)
    {
      return avg(vec, 0, shiftFromBegin);
    }

    double
    avg(const std::vector<double> & vec,
        unsigned shiftFromBegin_begin,
        unsigned shiftFromBegin_end)
    {
      if(! (shiftFromBegin_begin <= shiftFromBegin_end &&
            shiftFromBegin_begin <= vec.size() &&
            shiftFromBegin_end <= vec.size()))
        throw std::range_error(__PRETTY_FUNCTION__);

      if(! vec.size() || shiftFromBegin_begin == shiftFromBegin_end) return 0.;
      return std::accumulate(vec.begin() + shiftFromBegin_begin,
                             vec.begin() + shiftFromBegin_end,
                             0.) /
             (shiftFromBegin_end - shiftFromBegin_begin);
    }

    double
    stdev(const std::vector<double> & vec,
          double average)
    {
      return stdev(vec, 0, vec.size(), average);
    }

    double
    stdev(const std::vector<double> & vec,
          unsigned shiftFromBegin,
          double average)
    {
      return stdev(vec, 0, shiftFromBegin, average);
    }

    double
    stdev(const std::vector<double> & vec,
          unsigned shiftFromBegin_begin,
          unsigned shiftFromBegin_end,
          double average)
    {
      const unsigned range = shiftFromBegin_end - shiftFromBegin_begin;
      return l2(vec, shiftFromBegin_begin, shiftFromBegin_end, average) /
             std::sqrt(range >= 2 ? range * (range - 1) : 1.);
    }
  }
}
