#include "tthAnalysis/tthMEM/interface/tthMEMlvFunctions.h"

#include <cmath> // std::max()
#include <iomanip> // std::setw()
#include <vector> // std::vector<>

#define _SPACE std::setw(10) << std::right

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const LorentzVectorWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ':';
    os << "pt ="   << _SPACE << v.v_.pt()  << "; "
       << "eta ="  << _SPACE << v.v_.eta() << "; "
       << "phi ="  << _SPACE << v.v_.phi() << "; "
       << "mass =" << _SPACE << v.v_.mass();
    return os;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const LorentzMinkowskiWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ':';
    os << "en =" << _SPACE << v.v_.energy() << "; "
       << "px =" << _SPACE << v.v_.px()     << "; "
       << "py =" << _SPACE << v.v_.py()     << "; "
       << "pz =" << _SPACE << v.v_.pz();
    return os;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const LorentzVectorExtWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ':';
    os << "pt ="   << _SPACE << v.v_.pt()     << "; "
       << "eta ="  << _SPACE << v.v_.eta()    << "; "
       << "phi ="  << _SPACE << v.v_.phi()    << "; "
       << "mass =" << _SPACE << v.v_.mass()   << "; "
       << "en = "  << _SPACE << v.v_.energy() << "; "
       << "|p| = " << _SPACE << v.v_.P()      << "; "
       << "theta = " << _SPACE << v.v_.theta();
    return os;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const VectorCartesianWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ':';
    os << "x =" << _SPACE << v.v_.x() << "; "
       << "y =" << _SPACE << v.v_.y() << "; "
       << "z =" << _SPACE << v.v_.z();
    return os;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const VectorSphericalWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ':';
    os << "norm ="  << _SPACE << v.v_.r()     << "; "
       << "theta =" << _SPACE << v.v_.theta() << "; "
       << "phi ="   << _SPACE << v.v_.phi();
    return os;
  }

  LorentzVector
  getLorentzVector(const Vector & v,
                   const double e)
  {
    return LorentzVector(v.x(), v.y(), v.z(), e);
  }

  LorentzVector
  getLorentzVector(const math::PtEtaPhiMLorentzVector & v)
  {
    return LorentzVector(v.x(), v.y(), v.z(), v.e());
  }

  Vector
  getVector(const LorentzVector & v)
  {
    return Vector(v.x(), v.y(), v.z());
  }

  TVectorD
  getVector(const Vector & v)
  {
    const std::vector<double> v_({ v.x(), v.y(), v.z() });
    return TVectorD(3, v_.data());
  }

  Vector
  getVector(const TVectorD & v)
  {
    return Vector(v(0), v(1), v(2));
  }

  Vector
  getVector(const TMatrixDRow & row)
  {
    return Vector(row(0), row(1), row(2));
  }
}
