#ifndef TTHMEMAUXFUNCTIONS_H
#define TTHMEMAUXFUNCTIONS_H

#include <DataFormats/Math/interface/LorentzVector.h> // math::XYZTLorentzVectorD
#include <DataFormats/Math/interface/Vector3D.h> // math::XYZVectorD, math::RThetaPhiVectorD

#include <string> // std::string
#include <cmath> // std::acos()
#include <ostream> // std::ostream

namespace tthMEM
{
  typedef math::XYZTLorentzVectorD LorentzVector;
  typedef math::XYZVectorD         Vector;
  typedef math::RThetaPhiVectorD   VectorSpherical;

  typedef struct LorentzVectorWrap
  {
    LorentzVectorWrap(const LorentzVector & v);
    LorentzVectorWrap(const std::string & name,
                      const LorentzVector & v);
    const std::string name_;
    const LorentzVector & v_;
    std::size_t textFieldWidth_ = 15;

    friend std::ostream &
    operator<<(std::ostream & os,
               const LorentzVectorWrap & v);
  } lvrap;

  typedef struct LorentzMinkowskiWrap
  {
    LorentzMinkowskiWrap(const LorentzVector & v);
    LorentzMinkowskiWrap(const std::string & name,
                         const LorentzVector & v);
    const std::string name_;
    const LorentzVector & v_;
    std::size_t textFieldWidth_ = 15;

    friend std::ostream &
    operator<<(std::ostream & os,
               const LorentzMinkowskiWrap & v);
  } lmvrap;

  typedef struct VectorCartesianWrap
  {
    VectorCartesianWrap(const Vector & v);
    VectorCartesianWrap(const std::string & name,
                        const Vector & v);
    const std::string name_;
    const Vector & v_;
    std::size_t textFieldWidth_ = 15;

    friend std::ostream &
    operator<<(std::ostream & os,
               const VectorCartesianWrap & v);
  } cvrap;

  typedef struct VectorSphericalWrap
  {
    VectorSphericalWrap(const std::string & name,
                        const Vector & v);
    VectorSphericalWrap(const Vector & v);
    const std::string name_;
    const Vector & v_;
    std::size_t textFieldWidth_ = 15;

    friend std::ostream &
    operator<<(std::ostream & os,
               const VectorSphericalWrap & v);
  } svrap;

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
   * @brief Saves 4-momentum to a double array
   * @param lv Lorentz 4-vector
   * @param mg double array
   */
  void
  setMGmomentum(const LorentzVector & lv,
                double * mg);

  /**
   * @brief Constructs LorentzVector from Vector and energy
   * @param v The Vector
   * @param e The energy
   * @return The LorentzVector
   */
  LorentzVector
  getLorentzVector(const Vector & v,
                   const double e);

  /**
   * @brief Returns Vector from LorentzVector
   * @param v The LorentzVector
   * @return The Vector
   */
  Vector
  getVector(const LorentzVector & v);

  /**
   * @brief Finds the full path of a file in CMSSW directories
   * @param fileName A given file
   * @return Full path of the given file
   */
  std::string
  findFile(const std::string & fileName);

  /* lambdas replacing std::pow() as the latter is quite inefficient */
  const auto pow2 = [](double x) { return x * x; };
  const auto pow3 = [](double x) { return x * x * x; };
  const auto pow4 = [](double x) { return x * x * x * x; };
  const auto pow5 = [](double x) { return x * x * x * x * x; };
  const auto pow6 = [](double x) { return x * x * x * x * x * x; };

  /**
   * @brief Case-insensitive string comparison functor
   *
   * Needed in building the by-directional map where the value can be
   * fetched by a string, while ignoring its case
   */
  struct iStrComparator
  {
    bool
    operator()(const std::string & lhs,
               const std::string & rhs) const;
  };
}

#endif // TTHMEMAUXFUNCTIONS_H
