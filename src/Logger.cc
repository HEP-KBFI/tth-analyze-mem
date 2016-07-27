#include "tthAnalysis/tthMEM/interface/Logger.h"

namespace tthMEM
{
  Logger
  wrap(std::ostream & os,
       const std::string & level)
  {
    return Logger(os, level);
  }
}
