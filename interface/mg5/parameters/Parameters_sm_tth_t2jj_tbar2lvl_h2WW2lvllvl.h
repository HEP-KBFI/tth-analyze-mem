//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_sm_tth_t2jj_tbar2lvl_h2WW2lvllvl_H
#define Parameters_sm_tth_t2jj_tbar2lvl_h2WW2lvllvl_H

#include "tthAnalysis/tthMEM/interface/mg5/read_slha.h"

#include <complex>

class Parameters_sm_tth_t2jj_tbar2lvl_h2WW2lvllvl
{
  public:

    static Parameters_sm_tth_t2jj_tbar2lvl_h2WW2lvllvl * getInstance();

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double mdl_WH, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau, mdl_ymt, mdl_ymb, aS,
        mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA, mdl_MT, mdl_MB,
        mdl_conjg__CKM3x3, mdl_CKM3x3, mdl_conjg__CKM1x1, mdl_MZ__exp__2,
        mdl_MZ__exp__4, mdl_sqrt__2, mdl_MH__exp__2, mdl_aEW, mdl_MW,
        mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2,
        mdl_sw, mdl_g1, mdl_gw, mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_yb,
        mdl_yt, mdl_ytau, mdl_muH, mdl_ee__exp__2, mdl_sw__exp__2,
        mdl_cw__exp__2;
    std::complex<double> mdl_complexi, mdl_I1x33, mdl_I2x33, mdl_I3x33,
        mdl_I4x33;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_72, GC_94, GC_100; 
    // Model couplings dependent on aS
    std::complex<double> GC_11, GC_10; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_sm_tth_t2jj_tbar2lvl_h2WW2lvllvl * instance;
}; 

#endif  // Parameters_sm_tth_t2jj_tbar2lvl_h2WW2lvllvl_H

