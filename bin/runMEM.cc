#include <iostream> // std::cout
#include <cstdlib> // EXIT_SUCCESS

#include "tthAnalysis/tthMEM/interface/Logger.h"

int
main()
{
  LOGERR << "testing error message " << 123;
  return EXIT_SUCCESS;
}
