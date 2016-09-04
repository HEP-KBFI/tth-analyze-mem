#ifndef TTHMEMLVFUNCTIONS_H
#define TTHMEMLVFUNCTIONS_H

#include <DataFormats/Math/interface/LorentzVector.h> // math::XYZTLorentzVectorD
#include <DataFormats/Math/interface/Vector3D.h> // math::XYZVectorD, math::RThetaPhiVectorD

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
   * @brief Returns Vector from LorentzVector
   * @param v The LorentzVector
   * @return The Vector
   */
  Vector
  getVector(const LorentzVector & v);
}

#endif // TTHMEMLVFUNCTIONS_H
