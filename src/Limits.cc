#include "tthAnalysis/tthMEM/interface/Limits.h"

using namespace tthMEM;

Limits::Limits(double begin,
               double end)
  : begin_(begin)
  , end_(end)
{}

bool
Limits::isWithin(double val) const
{
  return begin_ <= val && val <= end_;
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const Limits & limits)
  {
    os << "[" << limits.begin_ << "; " << limits.end_ << "]";
    return os;
  }
}
