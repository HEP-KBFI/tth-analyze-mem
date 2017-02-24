#ifndef MEMINTEGRATORVEGAS_H
#define MEMINTEGRATORVEGAS_H

#include "tthAnalysis/tthMEM/interface/integrators/MEMIntegratorBase.h" // MEMIntegratorBase
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()

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
   * Originally created by Christian Veelken, NICPB Tallinn
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
    integrate(MEMIntegratorBase::gPtr_Fortran integrand __attribute__((unused)),
              const double * xl __attribute__((unused)),
              const double * xu __attribute__((unused)),
              unsigned dimension __attribute__((unused)),
              double & integral __attribute__((unused)),
              double & integralErr __attribute__((unused))) override
    {
      throw_line_ext("invalid usage", TTHEXCEPTION_ERR_CODE_UNDEFINED_FUNCTION)
        << "You must use integrate(gPtr_C, ...) not this one";
    }

    void *
    metadata() override;

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
