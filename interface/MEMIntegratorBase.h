#ifndef MEMINTEGRATORBASE_H
#define MEMINTEGRATORBASE_H

#include <cstddef> // std::size_t

namespace tthMEM
{
  /**
   * @brief Abstract class for VAMP and VEGAS integration classes
   */
  class MEMIntegratorBase
  {
  public:
    typedef double (*gPtr_C)      (double *,  std::size_t, void *);
    typedef double (*gPtr_Fortran)(double **, std::size_t, void **);

    MEMIntegratorBase() {}
    virtual ~MEMIntegratorBase() {}

    /**
     * @brief C wrapper for function that performs the integration (needed by VEGAS)
     * @param integrand   Function to perform integration on
     * @param xl          Lower left corner of the d-dimensional hypercube
     * @param xu          Upper right corner of the d-dimensional hypercube
     * @param dimension   Dimensionality of the integrand (hypercube)
     * @param integral    Integral
     * @param integralErr Error of integration
     */
    virtual void
    integrate(gPtr_C integrand,
              const double * xl,
              const double * xu,
              unsigned dimension,
              double & integral,
              double & integralErr) = 0;

    /**
     * @brief Fortran wrapper for function that performs the integration (needed by VAMP)
     * @param integrand   Function to perform integration on
     * @param xl          Lower left corner of the d-dimensional hypercube
     * @param xu          Upper right corner of the d-dimensional hypercube
     * @param dimension   Dimensionality of the integrand (hypercube)
     * @param integral    Integral
     * @param integralErr Error of integration
     */
    virtual void
    integrate(gPtr_Fortran integrand,
              const double * xl,
              const double * xu,
              unsigned dimension,
              double & integral,
              double & integralErr) = 0;
  };
}

#endif // MEMINTEGRATORBASE_H
