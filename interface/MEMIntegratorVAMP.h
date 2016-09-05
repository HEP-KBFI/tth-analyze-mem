#ifndef MEMINTEGRATORVAMP_H
#define MEMINTEGRATORVAMP_H

#include "tthAnalysis/tthMEM/interface/MEMIntegratorBase.h" // MEMIntegratorBase
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

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
    integrate(MEMIntegratorBase::gPtr_C integrand,
              const double * xl,
              const double * xu,
              unsigned dimension,
              double & integral,
              double & integralErr) override
    {
      throw_line("invalid usage")
        << "You must use integrate(gPtr_Fortran, ...) not this one";
    }

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
