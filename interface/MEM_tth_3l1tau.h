#ifndef MEM_TTH_3L1TAU_H
#define MEM_TTH_3L1TAU_H

#include <string> // std::string

#include "tthAnalysis/tthMEM/interface/MeasuredEvent.h"

namespace tthMEM
{
  /**
   * @brief Class for performing the phase space integration
   *
   * Sets up the integration boundaries, PDF, ME computed with MadGraph.
   *
   * @todo
   *   - add TF options
   */
  class MEM_tth_3l1tau
  {
  public:
    MEM_tth_3l1tau(double sqrtS,
                   int integrationMode,
                   const std::string & pdfName,
                   const std::string & madgraphFileName);
    ~MEM_tth_3l1tau();

  private:
    int integrationMode_;
  };
}

#endif // MEM_TTH_3L1TAU_H
