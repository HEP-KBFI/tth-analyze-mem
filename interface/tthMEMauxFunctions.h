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

  // all masses and energies are given in GeV
  const double sqrtS = 13.e+3;
  const double chargedPionMass = 0.13957;
  const double tauLeptonMass = 1.77685;
  const double tauLeptonMassSquared = tauLeptonMass * tauLeptonMass;

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
   * @return
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
