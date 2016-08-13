#ifndef MEM_TTH_3L1TAU_H
#define MEM_TTH_3L1TAU_H

#include <string> // std::string

#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/integrand_tth_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/MEMIntegratorBase.h"

#include <TBenchmark.h> // TBenchmark

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
    enum IntegrationMode { kUndefined, kVEGAS, kVAMP };

    MEM_tth_3l1tau(double sqrtS,
                   const std::string & pdfName,
                   const std::string & madgraphFileName);
    ~MEM_tth_3l1tau();

    void
    setMaxObjFunctionCalls(unsigned maxObjFunctionCalls);

    void
    setIntegrationMode(IntegrationMode integrationMode);

    void
    setIntegrationMode(const std::string & integrationModeString);

    double
    getComputingTime_cpu() const;

    double
    getComputingTime_real() const;

    double
    getAverageComputingTime_cpu() const;

    double
    getAverageComputingTime_real() const;

    double
    integrate(const MeasuredEvent_3l1tau & ev);

  private:
    integrand_tth_3l1tau * integrand_;
    double sqrtS_;
    MeasuredEvent_3l1tau ev_;

    IntegrationMode integrationMode_;
    MEMIntegratorBase * intAlgo_;
    unsigned maxObjFunctionCalls_;
    unsigned numDimensions_;
    double precision_;

    double * xl_;
    double * xu_;

    TBenchmark * clock_;
    double numSeconds_cpu_;
    double numSeconds_real_;
    double numSecondsAccumul_cpu_;
    double numSecondsAccumul_real_;
    unsigned nof_calls_;
  };
}

#endif // MEM_TTH_3L1TAU_H
