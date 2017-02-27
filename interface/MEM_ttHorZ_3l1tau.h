#ifndef MEM_TTH_3L1TAU_H
#define MEM_TTH_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h" // MeasuredEvent_3l1tau

class TBenchmark;

namespace tthMEM
{
  class Integrand_ttHorZ_3l1tau;
  class MEMIntegratorBase;

  /**
   * @brief Class for performing the phase space integration
   *
   * Sets up the integration boundaries, PDF, ME computed with MadGraph.
   *
   * @todo
   *   - add TF options
   */
  class MEM_ttHorZ_3l1tau
  {
  public:
    enum class IntegrationMode
    {
      kUndefined = 0,
      kVEGAS = 1,
      kVAMP = 2,
      kMarkovChain = 3
    };

    MEM_ttHorZ_3l1tau() = default;

    MEM_ttHorZ_3l1tau(const std::string & pdfName,
                      const std::string & madgraphFileName,
                      const VariableManager_3l1tau & vm);
    ~MEM_ttHorZ_3l1tau();

    MEM_ttHorZ_3l1tau &
    setMaxObjFunctionCalls(unsigned maxObjFunctionCalls);

    MEM_ttHorZ_3l1tau &
    setMarkovChainParams(const std::string & mxMode,
                         unsigned nofBatches,
                         unsigned nofChains,
                         unsigned maxCallsStartingPos,
                         double epsilon0,
                         double T0,
                         double nu);

    MEM_ttHorZ_3l1tau &
    setIntegrationMode(IntegrationMode integrationMode);

    MEM_ttHorZ_3l1tau &
    setIntegrationMode(const std::string & integrationModeString);

    MEM_ttHorZ_3l1tau &
    setBJetTransferFunction(bool setTF);

    MEM_ttHorZ_3l1tau &
    setHiggsWidth(double higgsWidth);

    MEM_ttHorZ_3l1tau &
    useAvgBjetCombo(bool useAvgB);

    bool
    isMarkovChainIntegrator() const;

    double
    getComputingTime_cpu() const;

    double
    getComputingTime_real() const;

    double
    getAverageComputingTime_cpu() const;

    double
    getAverageComputingTime_real() const;

    long long
    getNofMXMCTries() const;

    std::array<double, 2>
    integrate(const MeasuredEvent_3l1tau & ev,
              ME_mg5_3l1tau currentME);

  private:
    void
    initialize(const std::string & pdfName,
               const std::string madgraphFileName);

    VariableManager_3l1tau vm_;
    Integrand_ttHorZ_3l1tau * integrand_;
    MeasuredEvent_3l1tau ev_;

    IntegrationMode integrationMode_;
    MEMIntegratorBase * intAlgo_;
    unsigned maxObjFunctionCalls_;
    unsigned numCallsGridOpt_;
    unsigned numCallsIntEval_;

    /* for Markov Chain integrator */
    std::string mxMode_;
    unsigned nofChains_;
    unsigned nofIterBurnin_;
    unsigned nofIterSampling_;
    unsigned nofIterSimAnnPhase1_;
    unsigned nofIterSimAnnPhase2_;
    unsigned nofBatches_;
    unsigned maxCallsStartingPos_;
    double T0_;
    double alpha_;
    double epsilon0_;
    double nu_;
    long long nofMXMCTries; // counter variable

    double higgsWidth_;
    ///< use custom Higgs width in case the phase space is too narrow

    unsigned numDimensions_;
    bool setTF_;
    bool useAvgB_;

    TBenchmark * clock_;
    double numSeconds_cpu_;
    double numSeconds_real_;
    double numSecondsAccumul_cpu_;
    double numSecondsAccumul_real_;
    unsigned nof_calls_;
  };
}

#endif // MEM_TTH_3L1TAU_H
