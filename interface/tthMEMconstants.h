#ifndef TTHMEMCONSTANTS_H
#define TTHMEMCONSTANTS_H

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // pow*(), pi()

namespace tthMEM
{
  namespace constants
  {
    // if not specified all masses, energies and widths are given in GeV
    const double sqrtS = 13.e+3;
    const double s = pow2(sqrtS);
    const double invSqrtS = 1. / sqrtS;
    const double invS = 1. / s;

    const double cTimesHbar = 0.1973; // GeV x fm
    const double conversionFactor = pow2(cTimesHbar) * 1.e+10; ///< 1 pb = 10^-40 m^2 = 10^10 fm^2
    const double GF = 1.16637e-5; // Fermi constant, in 1 / GeV^2
    const double GFSquared = pow2(GF);

    const double brTau2e = 0.1783; ///< taken from PDG booklet (2014, p 15)
    const double brTau2mu = 0.1741; ///< taken from PDG booklet (2014, p 15)
    const double brTau2hadrons = 1 - brTau2e - brTau2mu;
    const double brH2diTau = 6.21e-2;
    const double brH2diW = 2.26e-1;
    ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR2014

    const double xSectionTTH = 0.5085; // pb
    ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014#ttH_Process
    const double xSectionTTH2diTau = xSectionTTH * brH2diTau;
    const double xSectionTTH2diTauInGeV2 = xSectionTTH2diTau * conversionFactor; // 1 / GeV^2
    const double xSectionTTH2diW = xSectionTTH * brH2diW;
    const double xSectionTTH2diWinGeV2 = xSectionTTH2diW * conversionFactor;
    const double xSectionTTZ = 1.e-3 * 839.3; // second factor given in fb = 10^-3 pb
    ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGTTH#Plans_for_YR4
    const double xSectionTTZinGeV2 = xSectionTTZ * conversionFactor;

    const double massVisTauMin = 0.3;
    const double massVisTauMax = 1.65;
    const double massChargedPion = 0.139570; ///< taken from PDG booklet (2014, p 25)
    const double massTau = 1.777000e+00;   ///< MG value (MTA)
    const double massTauSquared = pow2(massTau);
    const double massHiggs = 1.250000e+02; ///< MG value (MH)
    const double massHiggsSquared = pow2(massHiggs);
    const double massZ = 9.118800e+01;     ///< MG value (MZ)
    const double massZSquared = pow2(massZ);
    const double massW = 8.02673592E+01;   ///< MG value (MW)
    const double massWSquared = pow2(massW);
    const double massB = 4.700000e+00;     ///< MG value (MB)
    const double massBSquared = pow2(massB);
    const double massT = 1.74300000E+02;   ///< MG value (MT)
    const double massTSquared = pow2(massT);

    const double cTau = 87.03e+9; // in originally in um (10^-6), converted to fm (10^-15)
    ///< mean life, taken from PDG booklet (2014, p 15)
    const double gammaTau = cTimesHbar / cTau;
    const double gammaTau2hadrons = gammaTau * brTau2hadrons;
    const double gammaT = 1.491500e+00;     ///< MG value (WT)
    const double gammaW = 2.047600e+00;     ///< MG value (WW)
    const double gammaZ = 2.441404e+00;     ///< MG value (WZ)
    const double gammaHiggs = 6.382339e-03; ///< MG value (WH)

    const double DeltaFactor = (massTSquared - massWSquared - massBSquared) / 2.;
    const double resolutionScaleTTH = massT + massHiggs / 2.;
    const double resolutionScaleTTZ = massT + massZ / 2.;
    const double ttHhadTauPSfactor = 16. * pi() * pow3(massTau);
    const double ttHfactor = pow2(massT * gammaT * gammaW * massHiggs * gammaHiggs /
                                  massW * massTau * gammaTau * s);
    const double ttZfactor = pow2(massT * gammaT * gammaW * massZ * gammaZ /
                                  massW * massTau * gammaTau * s);
  }
}

#endif // TTHMEMCONSTANTS_H
