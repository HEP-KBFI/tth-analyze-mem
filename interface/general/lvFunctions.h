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

  enum class LVRepresentation { Cartesian, Cylindrical, Spherical, Extensive };

  template<typename VectorObject,
           LVRepresentation Representation>
  struct AnyVectorWrap
  {
    AnyVectorWrap(const VectorObject & v)
      : AnyVectorWrap("", v)
    {}
    AnyVectorWrap(const std::string & name,
                  const VectorObject & v)
      : name_(name)
      , v_(v)
    {
      textFieldWidth_ = std::max(textFieldWidth_, name_.size()) + 2;
    }

    friend std::ostream &
    operator<<(std::ostream & os,
               const AnyVectorWrap<VectorObject, Representation> & v);

    const std::string name_;
    const VectorObject & v_;
    std::size_t textFieldWidth_ = 15;
  };

  typedef AnyVectorWrap<LorentzVector, LVRepresentation::Cylindrical> LorentzVectorWrap;
  typedef LorentzVectorWrap                                           lvrap;
  typedef AnyVectorWrap<LorentzVector, LVRepresentation::Cartesian> LorentzMinkowskiWrap;
  typedef LorentzMinkowskiWrap                                      lmvrap;
  typedef AnyVectorWrap<LorentzVector, LVRepresentation::Extensive> LorentzVectorExtWrap;
  typedef LorentzVectorExtWrap                                      lextvrap;
  typedef AnyVectorWrap<Vector, LVRepresentation::Cartesian> VectorCartesianWrap;
  typedef VectorCartesianWrap                                cvrap;
  typedef AnyVectorWrap<Vector, LVRepresentation::Spherical> VectorSphericalWrap;
  typedef VectorSphericalWrap                                svrap;

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
