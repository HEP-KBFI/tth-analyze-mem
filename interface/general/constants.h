#ifndef TTHMEMCONSTANTS_H
#define TTHMEMCONSTANTS_H

#include "tthAnalysis/tthMEM/interface/general/auxFunctions.h" // pow*(), pi()

namespace tthMEM
{
  namespace constants
  {
    // if not specified all masses, energies and widths are given in GeV
    const double sqrtS    __attribute__((unused)) = 13.e+3;
    const double s        __attribute__((unused)) = pow2(sqrtS);
    const double invSqrtS __attribute__((unused)) = 1. / sqrtS;
    const double invS     __attribute__((unused)) = 1. / s;

    const double cTimesHbar       __attribute__((unused)) = 0.1973; // GeV x fm
    const double conversionFactor __attribute__((unused)) = pow2(cTimesHbar) * 1.e+10; ///< 1 pb = 10^-40 m^2 = 10^10 fm^2
    const double GF               __attribute__((unused)) = 1.16637e-5; // Fermi constant, in 1 / GeV^2
    const double GFSquared        __attribute__((unused)) = pow2(GF);

    const double brTau2e       __attribute__((unused)) = 0.1783; ///< taken from PDG booklet (2014, p 15)
    const double brTau2mu      __attribute__((unused)) = 0.1741; ///< taken from PDG booklet (2014, p 15)
    const double brTau2hadrons __attribute__((unused)) = 1 - brTau2e - brTau2mu;
    const double brH2diTau     __attribute__((unused)) = 6.21e-2;
    const double brH2diW       __attribute__((unused)) = 2.26e-1;
    ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2014

    const double xSectionTTH             __attribute__((unused)) = 0.5085; // pb
    ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014#ttH_Process
    const double xSectionTTH2diTau       __attribute__((unused)) = xSectionTTH * brH2diTau;
    const double xSectionTTH2diTauInGeV2 __attribute__((unused)) = xSectionTTH2diTau * conversionFactor; // 1 / GeV^2
    const double xSectionTTH2diW         __attribute__((unused)) = xSectionTTH * brH2diW;
    const double xSectionTTH2diWinGeV2   __attribute__((unused)) = xSectionTTH2diW * conversionFactor;
    const double xSectionTTZ             __attribute__((unused)) = 1.e-3 * 839.3; // second factor given in fb = 10^-3 pb
    ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGTTH#Plans_for_YR4
    const double xSectionTTZinGeV2       __attribute__((unused)) = xSectionTTZ * conversionFactor;

    const double massVisTauMin    __attribute__((unused)) = 0.3;
    const double massVisTauMax    __attribute__((unused)) = 1.65;
    const double massChargedPion  __attribute__((unused)) = 0.139570;        ///< taken from PDG booklet (2014, p 25)
    const double massTau          __attribute__((unused)) = 1.777000e+00;    ///< MG value (MTA)
    const double massTauSquared   __attribute__((unused)) = pow2(massTau);
    const double massHiggs        __attribute__((unused)) = 1.250000e+02;    ///< MG value (MH)
    const double massHiggsSquared __attribute__((unused)) = pow2(massHiggs);
    const double massZ            __attribute__((unused)) = 9.118800e+01;    ///< MG value (MZ)
    const double massZSquared     __attribute__((unused)) = pow2(massZ);
    const double massW            __attribute__((unused)) = 8.02673592E+01;  ///< MG value (MW)
    const double massWSquared     __attribute__((unused)) = pow2(massW);
    const double massB            __attribute__((unused)) = 4.700000e+00;    ///< MG value (MB)
    const double massBSquared     __attribute__((unused)) = pow2(massB);
    const double massT            __attribute__((unused)) = 1.74300000E+02;  ///< MG value (MT)
    const double massTSquared     __attribute__((unused)) = pow2(massT);

    const double cTau             __attribute__((unused)) = 87.03e+9; // in originally in um (10^-6), converted to fm (10^-15)
    ///< mean life, taken from PDG booklet (2014, p 15)
    const double gammaTau         __attribute__((unused)) = cTimesHbar / cTau;
    const double gammaTau2hadrons __attribute__((unused)) = gammaTau * brTau2hadrons;
    const double gammaT           __attribute__((unused)) = 1.491500e+00; ///< MG value (WT)
    const double gammaW           __attribute__((unused)) = 2.047600e+00; ///< MG value (WW)
    const double gammaZ           __attribute__((unused)) = 2.441404e+00; ///< MG value (WZ)
    const double gammaHiggs       __attribute__((unused)) = 6.382339e-03; ///< MG value (WH)

    const double DeltaFactor        __attribute__((unused)) = (massTSquared - massWSquared - massBSquared) / 2.;
    const double resolutionScaleTTH __attribute__((unused)) = massT + massHiggs / 2.;
    const double resolutionScaleTTZ __attribute__((unused)) = massT + massZ / 2.;
    const double ttHhadTauPSfactor  __attribute__((unused)) = 16. * pi() * pow3(massTau);
    const double ttHfactor __attribute__((unused)) = pow2(
      massT * gammaT * gammaW * massHiggs * gammaHiggs / (massW * massTau * gammaTau * s)
    ); // missing factor: 2^-25 * pi^-13
    const double ttZfactor __attribute__((unused)) = pow2(
      massT * gammaT * gammaW * massZ * gammaZ / (massW * massTau * gammaTau * s)
    ); // missing factor: 2^-25 * pi^-13
  }
}

#endif // TTHMEMCONSTANTS_H
