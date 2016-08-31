#ifndef MEMINTEGRATORVEGAS_H
#define MEMINTEGRATORVEGAS_H

#include "tthAnalysis/tthMEM/interface/MEMIntegratorBase.h" // MEMIntegratorBase
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

#include <cassert> // assert()

#include <gsl/gsl_rng.h> // gsl_rng
#include <gsl/gsl_monte.h>       // gsl_monte_function
#include <gsl/gsl_monte_vegas.h> // gsl_monte_vegas_state

namespace tthMEM
{
  /** Interface to VEGAS integration algorithm.
   *
   * The VEGAS algorithm is documented in:
   *  [1] "A New Algorithm for Adaptive Multidimensional Integration",
   *      G.P. Lepage, J. Comput. Phys. 27 (1978) 192.
   *
   * Originally created by
   * @author Christian Veelken, NICPB Tallinn
   * (visit https://github.com/veelken/SVfitMEM)
   */
  class MEMIntegratorVEGAS
    : public MEMIntegratorBase
  {
  public:
    MEMIntegratorVEGAS(unsigned numCallsGridOpt,
                       unsigned numCallsIntEval,
                       unsigned maxIntEvalIter,
                       double maxChi2);
    ~MEMIntegratorVEGAS();

    void
    integrate(MEMIntegratorBase::gPtr_C integrand,
              const double * xl,
              const double * xu,
              unsigned dimension,
              double & integral,
              double & integralErr) override;

    void
    integrate(MEMIntegratorBase::gPtr_Fortran integrand,
              const double * xl,
              const double * xu,
              unsigned dimension,
              double & integral,
              double & integralErr) override
    {
      LOGERR << "You must use the other integrate() not this one";
      assert(0);
    }

  private:
    void
    setIntegrand(MEMIntegratorBase::gPtr_C integrand,
                 const double * xl,
                 const double * xu,
                 unsigned dimension);

    MEMIntegratorBase::gPtr_C integrand_;

    gsl_monte_function * vegasIntegrand_;
    gsl_monte_vegas_state * vegasWorkspace_;
    mutable gsl_rng * vegasRnd_;

    unsigned numCallsGridOpt_;
    unsigned numCallsIntEval_;
    unsigned numDimensions_;

    unsigned maxIntEvalIter_;
    double maxChi2_;
    double precision_;

    mutable double * xl_;
    mutable double * xu_;
  };
}

#endif // MEMINTEGRATORVEGAS_H
