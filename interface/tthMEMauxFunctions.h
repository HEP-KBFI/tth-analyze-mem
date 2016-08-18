#ifndef TTHMEMAUXFUNCTIONS_H
#define TTHMEMAUXFUNCTIONS_H

#include "DataFormats/Math/interface/LorentzVector.h" // math::XYZTLorentzVectorD
#include "DataFormats/Math/interface/Vector3D.h" // math::XYZVectorD, math::RThetaPhiVectorD

#include <string> // std::string
#include <iostream> // std::cout
#include <cmath> // std::acos()

namespace tthMEM
{
  typedef math::XYZTLorentzVectorD LorentzVector;
  typedef math::XYZVectorD         Vector;
  typedef math::RThetaPhiVectorD   VectorSpherical;

  typedef struct LorentzVectorWrap
  {
    LorentzVectorWrap(const LorentzVector & v);
    const LorentzVector & v_;

    friend std::ostream &
    operator<<(std::ostream & os,
               const LorentzVectorWrap & v);
  } lvrap;

  typedef struct VectorCartesianWrap
  {
    VectorCartesianWrap(const Vector & v);
    const Vector & v_;

    friend std::ostream &
    operator<<(std::ostream & os,
               const VectorCartesianWrap & v);
  } cvrap;

  typedef struct VectorSphericalWrap
  {
    VectorSphericalWrap(const Vector & v);
    const Vector & v_;

    friend std::ostream &
    operator<<(std::ostream & os,
               const VectorSphericalWrap & v);
  } svrap;

  enum ME_mg5_3l1tau
  {
    kTTH = 0,
    kTTZ = 1
  };

  // if not specified all masses and energies are given in GeV
  const double sqrtS = 13.e+3;

  const double minVisTauMass = 0.3;
  const double maxVisTauMass = 1.5;
  const double chargedPionMass = 0.13957;
  const double tauLeptonMass = 1.77685;
  const double tauLeptonMassSquared = std::pow(tauLeptonMass, 2);
  const double cTimesHbar = 0.1973; // GeV x fm
  const double conversionFactor = std::pow(cTimesHbar, 2) * 1.e+10; ///< 1 pb = 10^-40 m^2 = 10^10 fm^2
  const double xSectionTTH = 0.5085; // pb
  ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt1314TeV2014#ttH_Process
  const double xSectionTTHinGeV2 = xSectionTTH * conversionFactor; // 1 / GeV^2
  const double xSectionTTZ = 1.e-3 * 839.3; // second factor given in fb = 10^-3 pb
  ///< taken from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGTTH#Plans_for_YR4
  const double xSectionTTZinGeV2 = xSectionTTZ * conversionFactor;
  const double massHiggs = 125.7; ///< taken from PDG booklet (2014, p 11)
  const double massHiggsSquared = std::pow(massHiggs, 2);
  const double massZ = 91.1876; ///< taken from PDG booklet (2014, p 9)
  const double massZSquared = std::pow(massZ, 2);
  const double massW = 80.385; ///< taken from PDG booklet (2014, p 8)
  const double massWSquared = std::pow(massW, 2);
  const double massB = 4.18; ///< MS scheme; taken from PDG booklet (2014, p 23)
  const double massBSquared = std::pow(massB, 2);
  const double massT = 173.21; ///< taken from PDG booklet (2014, p 23)
  const double massTSquared = std::pow(massT, 2);
  const double DeltaFactor = (massTSquared - massWSquared - massBSquared) / 2.;
  const double resolutionScaleTTH = massT + massHiggs / 2.;
  const double resolutionScaleTTZ = massT + massZ / 2.;

  /**
   * @brief Rounds double floating point number to N significant digits
   * @param x Floating point number to be rounded
   * @param n Number of significant figures
   * @return The rounded result
   *
   * For examples, see "test/test_tthMEMauxFunctions.cc"
   */
  double
  roundToNdigits(double x,
                 int n = 3);

  /**
   * @brief Rounds a given floating point number to the nearest unsigned integer
   * @param x The floating point number
   */
  unsigned
  roundToNearestUInt(double x);

  /**
   * @brief Calculates pi from relation acos(-1)
   * @return PI
   */
  constexpr double
  pi()
  {
    return std::acos(-1.L);
  }

  void
  setMGmomentum(const LorentzVector & lv,
                double * mg);

  /**
   * @brief Finds the full path of a file in CMSSW directories
   * @param fileName A given file
   * @return Full path of the given file
   */
  std::string
  findFile(const std::string & fileName);
}

#endif // TTHMEMAUXFUNCTIONS_H
