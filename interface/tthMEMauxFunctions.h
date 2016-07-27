#ifndef TTHMEMAUXFUNCTIONS_H
#define TTHMEMAUXFUNCTIONS_H

#include "DataFormats/Math/interface/LorentzVector.h" // math::XYZTLorentzVector
#include "DataFormats/Math/interface/Vector3D.h" // math::XYZVectorD

#include <string> // std::string
#include <iostream> // std::cout

namespace tthMEM
{
  typedef math::XYZTLorentzVector LorentzVector;
  typedef math::XYZVectorD Vector;

  const double chargedPionMass = 0.13957; // GeV

  /**
   * @brief Rounds double floating point number to N significant digits
   * @param x Floating point number to be rounded
   * @param n Number of significant figures
   * @return The rounded result
   *
   * For examples, see "test/test_tthMEMauxFunctions.cc"
   */
  double roundToNdigits(double x,
                        int n = 3);
}

#endif // TTHMEMAUXFUNCTIONS_H
