//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.4.2, 2016-06-10
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "tthAnalysis/tthMEM/interface/Parameters_sm_tth_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

// Initialize static instance
Parameters_sm_tth_3l1tau * Parameters_sm_tth_3l1tau::instance = 0;

// Function to get static instance - only one instance per program
Parameters_sm_tth_3l1tau * Parameters_sm_tth_3l1tau::getInstance()
{
  if (instance == 0)
    instance = new Parameters_sm_tth_3l1tau();

  return instance; 
}

void Parameters_sm_tth_3l1tau::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_WH = slha.get_block_entry("decay", 25, 6.382339e-03); 
  mdl_WW = slha.get_block_entry("decay", 24, 2.047600e+00); 
  mdl_WZ = slha.get_block_entry("decay", 23, 2.441404e+00); 
  mdl_WT = slha.get_block_entry("decay", 6, 1.491500e+00); 
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.730000e+02); 
  mdl_ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  aS = slha.get_block_entry("sminputs", 3, 1.180000e-01); 
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166390e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.325070e+02); 
  mdl_MH = slha.get_block_entry("mass", 25, 1.250000e+02); 
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118800e+01); 
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  mdl_MT = slha.get_block_entry("mass", 6, 1.730000e+02); 
  mdl_MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  mdl_conjg__CKM3x3 = 1.; 
  mdl_CKM3x3 = 1.; 
  mdl_conjg__CKM1x1 = 1.; 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_MZ__exp__2 = pow(mdl_MZ, 2.); 
  mdl_MZ__exp__4 = pow(mdl_MZ, 4.); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__2 = pow(mdl_MH, 2.); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = pow(mdl_MW, 2.); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_vev__exp__2 = pow(mdl_vev, 2.); 
  mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
  mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
  mdl_I1x33 = mdl_yb * mdl_conjg__CKM3x3; 
  mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3; 
  mdl_I3x33 = mdl_CKM3x3 * mdl_yt; 
  mdl_I4x33 = mdl_CKM3x3 * mdl_yb; 
  mdl_ee__exp__2 = pow(mdl_ee, 2.); 
  mdl_sw__exp__2 = pow(mdl_sw, 2.); 
  mdl_cw__exp__2 = pow(mdl_cw, 2.); 
}
void Parameters_sm_tth_3l1tau::setIndependentCouplings()
{
  GC_94 = -((mdl_complexi * mdl_yt)/mdl_sqrt__2); 
  GC_99 = -((mdl_complexi * mdl_ytau)/mdl_sqrt__2); 
  GC_100 = (mdl_ee * mdl_complexi * mdl_conjg__CKM1x1)/(mdl_sw * mdl_sqrt__2); 
}
void Parameters_sm_tth_3l1tau::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = pow(G, 2.); 
}
void Parameters_sm_tth_3l1tau::setDependentCouplings()
{
  GC_11 = mdl_complexi * G; 
  GC_10 = -G; 
}

// Routines for printing out parameters
void Parameters_sm_tth_3l1tau::printIndependentParameters()
{
  LOGDBG <<  "sm model parameters independent of event kinematics:";
  LOGDBG << setw(20) <<  "mdl_WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH;
  LOGDBG << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW;
  LOGDBG << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ;
  LOGDBG << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT;
  LOGDBG << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau;
  LOGDBG << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt;
  LOGDBG << setw(20) <<  "mdl_ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymb;
  LOGDBG << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS;
  LOGDBG << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf;
  LOGDBG << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1;
  LOGDBG << setw(20) <<  "mdl_MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MH;
  LOGDBG << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ;
  LOGDBG << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA;
  LOGDBG << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT;
  LOGDBG << setw(20) <<  "mdl_MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MB;
  LOGDBG << setw(20) <<  "mdl_conjg__CKM3x3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM3x3;
  LOGDBG << setw(20) <<  "mdl_CKM3x3 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM3x3;
  LOGDBG << setw(20) <<  "mdl_conjg__CKM1x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x1;
  LOGDBG << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi;
  LOGDBG << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2;
  LOGDBG << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4;
  LOGDBG << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2;
  LOGDBG << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2;
  LOGDBG << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW;
  LOGDBG << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW;
  LOGDBG << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW;
  LOGDBG << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee;
  LOGDBG << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2;
  LOGDBG << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2;
  LOGDBG << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw;
  LOGDBG << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2;
  LOGDBG << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw;
  LOGDBG << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1;
  LOGDBG << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw;
  LOGDBG << setw(20) <<  "mdl_vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_vev;
  LOGDBG << setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_vev__exp__2;
  LOGDBG << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam;
  LOGDBG << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb;
  LOGDBG << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt;
  LOGDBG << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau;
  LOGDBG << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH;
  LOGDBG << setw(20) <<  "mdl_I1x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1x33;
  LOGDBG << setw(20) <<  "mdl_I2x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2x33;
  LOGDBG << setw(20) <<  "mdl_I3x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3x33;
  LOGDBG << setw(20) <<  "mdl_I4x33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4x33;
  LOGDBG << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2;
  LOGDBG << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2;
  LOGDBG << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2;
}
void Parameters_sm_tth_3l1tau::printIndependentCouplings()
{
  LOGDBG <<  "sm model couplings independent of event kinematics:";
  LOGDBG << setw(20) <<  "GC_94 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_94;
  LOGDBG << setw(20) <<  "GC_99 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_99;
  LOGDBG << setw(20) <<  "GC_100 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_100;
}
void Parameters_sm_tth_3l1tau::printDependentParameters()
{
  LOGDBG <<  "sm model parameters dependent on event kinematics:";
  LOGDBG << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS;
  LOGDBG << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G;
  LOGDBG << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2;
}
void Parameters_sm_tth_3l1tau::printDependentCouplings()
{
  LOGDBG <<  "sm model couplings dependent on event kinematics:";
  LOGDBG << setw(20) <<  "GC_11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11;
  LOGDBG << setw(20) <<  "GC_10 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10;
}


