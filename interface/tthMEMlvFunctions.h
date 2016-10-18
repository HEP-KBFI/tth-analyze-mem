#ifndef TTHMEMLVFUNCTIONS_H
#define TTHMEMLVFUNCTIONS_H

#include <DataFormats/Math/interface/LorentzVector.h> // math::XYZTLorentzVectorD
#include <DataFormats/Math/interface/Vector3D.h> // math::XYZVectorD, math::RThetaPhiVectorD

#include <TVectorD.h> // TVectorD

#include <ostream> // std::ostream
#include <string> // std::string

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
   * @brief Constructs LorentzVector from Vector and energy
   * @param v The Vector
   * @param e The energy
   * @return The LorentzVector
   */
  LorentzVector
  getLorentzVector(const Vector & v,
                   double e);

  /**
   * @brief Constructs Minkowskian LorentzVector from pT, eta, phi & mass
   * @param v The version of the LorentzVector in pT, eta, phi & mass representation
   * @return The Minkowskian LorentzVector
   */
  LorentzVector
  getLorentzVector(const math::PtEtaPhiMLorentzVector & v);

  /**
   * @brief Returns Vector from LorentzVector
   * @param v The LorentzVector
   * @return The Vector
   */
  Vector
  getVector(const LorentzVector & v);

  /**
   * @brief Returns TVectorD from Vector
   * @param v The Vector
   * @return The TVectorD
   */
  TVectorD
  getVector(const Vector & v);

  /**
   * @brief Returns Vector from TVectorD
   * @param v The TVectorD
   * @return The Vector
   */
  Vector
  getVector(const TVectorD & v);

  /**
   * @brief Returns Vector from TMatrixDRow
   * @param row The TMatrixDRow
   * @return The Vector
   */
  Vector
  getVector(const TMatrixDRow & row);
}

#endif // TTHMEMLVFUNCTIONS_H
