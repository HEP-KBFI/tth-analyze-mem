#include "tthAnalysis/tthMEM/interface/Logger.h"

namespace tthMEM
{
  Logger
  wrap(std::ostream & os,
       const std::string & level,
       bool enable)
  {
    return Logger(os, level, enable);
  }
}
