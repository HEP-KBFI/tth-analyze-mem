#ifndef INTEGRAND_TTH_3L1TAU_LO_H
#define INTEGRAND_TTH_3L1TAU_LO_H

#include <string> // std::string

#include "tthAnalysis/tthMEM/interface/me_tth_3l1tau_lo_mg5.h" // me_tth_3l1tau_lo_mg5

#include "LHAPDF/LHAPDF.h" // LHAPDF::PDF

namespace tthMEM
{
  /**
   * @brief MEM integrand of the process tth, h->tautau, 3l1tau channel (LO)
   *
   * @todo
   *   - add methods for setting transfer functions
   *   - add a method that sets the measured event information
   *   - add measured event objects
   *   - implement the integrand
   */
  class integrand_tth_3l1tau_lo
  {
  public:
    /**
     * @brief Simple constructor
     * @param sqrtS            Center of momentum energy
     * @param pdfName          Name of parton distribution function (pdf)
     * @param madgraphFileName Full path to madgraph5 file (param_card.dat)
     */
    integrand_tth_3l1tau_lo(double sqrtS,
                            const std::string & pdfName,
                            const std::string & madgraphFileName);
    ~integrand_tth_3l1tau_lo();

    /**
     * @brief Evaluates the integrand
     * @return Probability
     */
    double
    eval() const;

  protected:
    double sqrtS_;
    double s_;

    LHAPDF::PDF * pdf_;

    mutable me_tth_3l1tau_lo_mg5 me_madgraph_;
    ///< @note mutable members can be modified by a const function (e.g. eval())
    bool me_madgraph_initialized_;
  };
}

#endif // INTEGRAND_TTH_3L1TAU_LO_H
