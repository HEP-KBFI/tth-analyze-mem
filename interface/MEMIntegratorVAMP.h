#ifndef MEMINTEGRATORVAMP_H
#define MEMINTEGRATORVAMP_H

#include <cassert> // assert()

#include "tthAnalysis/tthMEM/interface/MEMIntegratorBase.h" // MEMIntegratorBase
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

namespace tthMEM
{
  /** Interface to "Vegas AMPlified: Anisotropy, Multi-channel sampling and Parallelization" (VAMP) integration algorithm.
   *
   * The VAMP algorithm is documented in:
   *  [1] "Vegas revisited: Adaptive Monte Carlo integration beyond factorization",
   *      T. Ohl, J. Comput. Phys. 120 (1999) 13.
   *  [2] https://whizard.hepforge.org/vamp.pdf
   *
   * Originally created by
   * @author Christian Veelken, NICPB Tallinn
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
      LOGERR << "You must use the other integrate() not this one";
      assert(0);
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
