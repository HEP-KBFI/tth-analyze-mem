#ifndef TTHMEMAUXFUNCTIONS_H
#define TTHMEMAUXFUNCTIONS_H

#include "DataFormats/Math/interface/LorentzVector.h" // math::XYZTLorentzVector
#include "DataFormats/Math/interface/Vector3D.h" // math::XYZVectorD

#include <string> // std::string
#include <iostream> // std::cout
#include <cmath> // std::acos()

namespace tthMEM
{
  typedef math::XYZTLorentzVector LorentzVector;
  typedef math::XYZVectorD Vector;

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

  /**
   * @brief Finds the full path of a file in CMSSW directories
   * @param fileName A given file
   * @return Full path of the given file
   */
  std::string
  findFile(const std::string & fileName);
}

#endif // TTHMEMAUXFUNCTIONS_H
