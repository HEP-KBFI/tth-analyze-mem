#ifndef TTHMEMAUXFUNCTIONS_H
#define TTHMEMAUXFUNCTIONS_H

#include "DataFormats/Math/interface/LorentzVector.h" // math::XYZTLorentzVectorD
#include "DataFormats/Math/interface/Vector3D.h" // math::XYZVectorD, math::RThetaPhiVectorD

#include <string> // std::string
#include <cmath> // std::acos()
#include <vector> // std::vector<>
#include <ostream> // std::ostream
#include <iterator> // std::ostream_iterator<>
#include <algorithm> // std::copy()

namespace tthMEM
{
  enum ME_mg5_3l1tau
  {
    kTTH = 0,
    kTTZ = 1
  };

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

  namespace functions
  {
    /**
     * @brief Calculates L2 norm aka Euclidean norm of a given vector of doubles
     * @param v          The vector of doubles
     * @param shiftValue The value which is subtracted from every element beforehand
     * @return Magnitude of the vector in Euclidean space
     */
    double
    l2(const std::vector<double> & v,
       double shiftValue = 0.);

    /**
     * @brief Calculates L2 norm aka Euclidean norm of a given vector of doubles
     * @param v              The vector
     * @param shiftFromBegin End point up to which the norm is calculated
     * @param shiftValue     The value which is subtracted from every
     *                       element beforehand
     * @return The norm of the subvector in the element range [0, shiftFromBegin)
     */
    double
    l2(const std::vector<double> & v,
       unsigned shiftFromBegin,
       double shiftValue = 0.);

    /**
     * @brief Calculates L2 norm aka Euclidean norm of a given vector of doubles
     * @param v                    The vector
     * @param shiftFromBegin_begin The starting point from which the norm is
     *                             calculated
     * @param shiftFromBegin_end   The end point up to which the norm is calculated
     * @param shiftValue           The value which is subtracted from every
     *                             element beforehand
     * @return The norm of the subvector in the element range
     *         [shiftFromBegin_begin, shiftFromBegin_end)
     */
    double
    l2(const std::vector<double> & v,
       unsigned shiftFromBegin_begin,
       unsigned shiftFromBegin_end,
       double shiftValue = 0.);

    /**
     * @brief Calculates average of a given vector of doubles
     * @param v The vector
     * @return Average of the vector
     */
    double
    avg(const std::vector<double> & v);

    /**
     * @brief Calculates average of a given vector of doubles
     * @param v              The vector
     * @param shiftFromBegin End point up to which the average is calculated
     * @return Average of the subvector in the element range [0, shiftFromBegin)
     */
    double
    avg(const std::vector<double> & v,
        unsigned shiftFromBegin);

    /**
     * @brief Calculates average of a given vector of doubles
     * @param v                    The vector
     * @param shiftFromBegin_begin The starting point from which the average is
     *                             calculated
     * @param shiftFromBegin_end   The end point up to which the average is
     *                             calculated
     * @return Average of the subvector in the element range
     *         [shiftFromBegin_begin, shiftFromBegin_end)
     */
    double
    avg(const std::vector<double> & v,
        unsigned shiftFromBegin_begin,
        unsigned shiftFromBegin_end);

    /**
     * @brief Calculates standard deviation of a given vector of doubles
     * @param v       The vector
     * @param average Average of the vector
     * @return Standard deviation of the vector
     */
    double
    stdev(const std::vector<double> & v,
          double average);

    /**
     * @brief Calculates standard deviation of a given vector of doubles
     * @param v              The vector
     * @param shiftFromBegin End point up to which the standard deviation is
     *                       calculated
     * @param average        The average of the subvector
     * @return Standard deviation of the subvector in the element range
     *         [0, shiftFromBegin)
     */
    double
    stdev(const std::vector<double> & v,
          unsigned shiftFromBegin,
          double average);

    /**
     * @brief Calculates standard deviation of a given vector of doubles
     * @param v                    The vector
     * @param shiftFromBegin_begin The starting point from which the standard
     *                             deviation is calculated
     * @param shiftFromBegin_end   The end point up to which the standard
     *                             deviation is calculated
     * @param average              The average of the subvector
     * @return Standard deviation of the subvector in the element range
     *         [shiftFromBegin_begin, shiftFromBegin_end)
     */
    double
    stdev(const std::vector<double> & v,
          unsigned shiftFromBegin_begin,
          unsigned shiftFromBegin_end,
          double average);
  }

  /**
   * @brief Prints comma-separated vector elements
   * @param os Stream to print to
   * @param v  The vector
   * @return Modified stream
   *
   * Assumes that the typename T is printable (otherwise expect compilation error)
   */
  template <typename T>
  std::ostream &
  operator<<(std::ostream & os,
             const std::vector<T> & v)
  {
    if(! v.size()) return os;
    if(v.size() > 1)
      std::copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(os, ", "));
    os << v[v.size() - 1];
    return os;
  }

  /**
   * @brief Hash function for enum classes used in std::unordered_map<>
   *
   * Cheers to the author at http://stackoverflow.com/a/24847480
   */
  struct EnumClassHash
  {
    template <typename T>
    std::size_t operator()(T t) const
    {
      return static_cast<std::size_t>(t);
    }
  };

  /**
   * @brief Meta-class for looping over enum classes which have
   *        First and Last entries specified.
   *
   * Cheers to the author at http://stackoverflow.com/a/8498694
   */
  template <typename T>
  struct Enum
  {
    struct Iterator
    {
      Iterator(int val)
        : val_(val)
      {}

      T
      operator*(void) const
      {
        return static_cast<T>(val_);
      }

      void
      operator++(void)
      {
        ++val_;
      }

      bool
      operator!=(Iterator rhs)
      {
        return val_ != rhs.val_;
      }
    private:
      int val_;
    };
  };

  template <typename T>
  typename Enum<T>::Iterator begin(Enum<T>)
  {
    return typename Enum<T>::Iterator(static_cast<int>(T::First));
  }

  template <typename T>
  typename Enum<T>::Iterator end(Enum<T>)
  {
    return typename Enum<T>::Iterator(static_cast<int>(T::Last) + 1);
  }

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

#endif // TTHMEMAUXFUNCTIONS_H
