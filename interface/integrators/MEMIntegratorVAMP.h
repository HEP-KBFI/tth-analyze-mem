#ifndef MEMINTEGRATORVAMP_H
#define MEMINTEGRATORVAMP_H

#include "tthAnalysis/tthMEM/interface/integrators/MEMIntegratorBase.h" // MEMIntegratorBase
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()

namespace tthMEM
{
  /** Interface to "Vegas AMPlified: Anisotropy, Multi-channel sampling and Parallelization" (VAMP) integration algorithm.
   *
   * The VAMP algorithm is documented in:
   *  [1] "Vegas revisited: Adaptive Monte Carlo integration beyond factorization",
   *      T. Ohl, J. Comput. Phys. 120 (1999) 13.
   *  [2] https://whizard.hepforge.org/vamp.pdf
   *
   * Originally created by Christian Veelken, NICPB Tallinn
   * (visit https://github.com/veelken/SVfitMEM)
   */
  class MEMIntegratorVAMP
    : public MEMIntegratorBase
  {
  public:
    MEMIntegratorVAMP(unsigned numCallsGridOpt,
                      unsigned numCallsIntEval);
    ~MEMIntegratorVAMP() override;

    void
    integrate(MEMIntegratorBase::gPtr_Fortran integrand,
              const double * xl,
              const double * xu,
              unsigned dimension,
              double & integral,
              double & integralErr) override;

    void
    integrate(MEMIntegratorBase::gPtr_C integrand __attribute__((unused)),
              const double * xl __attribute__((unused)),
              const double * xu __attribute__((unused)),
              unsigned dimension __attribute__((unused)),
              double & integral __attribute__((unused)),
              double & integralErr __attribute__((unused))) override
    {
      throw_line_ext("invalid usage", TTHEXCEPTION_ERR_CODE_UNDEFINED_FUNCTION)
        << "You must use integrate(gPtr_Fortran, ...) not this one";
    }

    void *
    metadata() override;

  private:
    void
    setIntegrand(MEMIntegratorBase::gPtr_Fortran integrand,
                 const double * xl,
                 const double * xu,
                 unsigned dimension);

    MEMIntegratorBase::gPtr_Fortran integrand_;

    unsigned numCallsGridOpt_;
    unsigned numCallsIntEval_;
    unsigned numDimensions_;

    mutable double * xl_;
    mutable double * xu_;
  };
}

#endif // MEMINTEGRATORVAMP_H
