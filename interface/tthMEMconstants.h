#ifndef TTHMEMCONSTANTS_H
#define TTHMEMCONSTANTS_H

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
    const double massVisTauMax = 1.5;
    const double massChargedPion = 0.139570; ///< taken from PDG booklet (2014, p 25)
    const double massTau = 1.77682; ///< taken from PDG booklet (2014, p 15)
    const double massTauSquared = pow2(massTau);
    const double massHiggs = 125.7; ///< taken from PDG booklet (2014, p 11)
    const double massHiggsSquared = pow2(massHiggs);
    const double massZ = 91.1876; ///< taken from PDG booklet (2014, p 9)
    const double massZSquared = pow2(massZ);
    const double massW = 80.385; ///< taken from PDG booklet (2014, p 8)
    const double massWSquared = pow2(massW);
    const double massB = 4.18; ///< MS scheme; taken from PDG booklet (2014, p 23)
    const double massBSquared = pow2(massB);
    const double massT = 173.21; ///< taken from PDG booklet (2014, p 23)
    const double massTSquared = pow2(massT);

    const double cTau = 87.03e+9; // in originally in um (10^-6), converted to fm (10^-15)
    ///< mean life, taken from PDG booklet (2014, p 15)
    const double gammaTau = cTimesHbar / cTau;
    const double gammaTau2hadrons = gammaTau * brTau2hadrons;
    const double gammaT = 2.0; ///< taken from PDG booklet (2014, p 23)
    const double gammaW = 2.085; ///< taken from PDG booklet (2014, p 8)
    const double gammaZ = 2.4952; ///< taken from PDG booklet (2014, p 9)
    const double gammaHiggs = 1.e-3 * massHiggs; ///< taken from thin air

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
