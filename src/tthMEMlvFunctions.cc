#include "tthAnalysis/tthMEM/interface/tthMEMlvFunctions.h"

#include <cmath> // std::max()
#include <iomanip> // std::setw()

namespace tthMEM
{
  LorentzVectorWrap::LorentzVectorWrap(const LorentzVector & v)
    : LorentzVectorWrap("", v)
  {}

  LorentzVectorWrap::LorentzVectorWrap(const std::string & name,
                                       const LorentzVector & v)
    : name_(name)
    , v_(v)
  {
    textFieldWidth_ = std::max(textFieldWidth_, name_.size()) + 2;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const LorentzVectorWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ":";
    os << "En ="   << std::setw(10) << std::right << v.v_.energy()
       << ";" << std::string(5, ' ')
       << "pT ="   << std::setw(10) << std::right << v.v_.pt()
       << ";" << std::string(5, ' ')
       << "mass =" << std::setw(10) << std::right << v.v_.mass();
    return os;
  }

  VectorCartesianWrap::VectorCartesianWrap(const Vector & v)
    : VectorCartesianWrap("", v)
  {}

  LorentzMinkowskiWrap::LorentzMinkowskiWrap(const LorentzVector & v)
    : LorentzMinkowskiWrap("", v)
  {}

  LorentzMinkowskiWrap::LorentzMinkowskiWrap(const std::string & name,
                                             const LorentzVector & v)
    : name_(name)
    , v_(v)
  {
    textFieldWidth_ = std::max(textFieldWidth_, name_.size()) + 2;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const LorentzMinkowskiWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ":";
    os << "En =" << std::setw(10) << std::right << v.v_.energy()
       << ";" << std::string(5, ' ')
       << "px =" << std::setw(10) << std::right << v.v_.px()
       << ";" << std::string(5, ' ')
       << "py =" << std::setw(10) << std::right << v.v_.py()
       << ";" << std::string(5, ' ')
       << "pz =" << std::setw(10) << std::right << v.v_.pz();
    return os;
  }

  VectorCartesianWrap::VectorCartesianWrap(const std::string & name,
                                           const Vector & v)
    : name_(name)
    , v_(v)
  {
    textFieldWidth_ = std::max(textFieldWidth_, name_.size()) + 2;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const VectorCartesianWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ":";
    os << "x =" << std::setw(10) << std::right << v.v_.x()
       << ";" << std::string(5, ' ')
       << "y =" << std::setw(10) << std::right << v.v_.y()
       << ";" << std::string(5, ' ')
       << "z =" << std::setw(10) << std::right << v.v_.z();
    return os;
  }

  VectorSphericalWrap::VectorSphericalWrap(const Vector & v)
    : VectorSphericalWrap("", v)
  {}

  VectorSphericalWrap::VectorSphericalWrap(const std::string & name,
                                           const Vector & v)
    : name_(name)
    , v_(v)
  {
    textFieldWidth_ = std::max(textFieldWidth_, name_.size()) + 2;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const VectorSphericalWrap & v)
  {
    if(v.name_.size() != 0)
      os << std::setw(v.textFieldWidth_) << std::left << v.name_ + ":";
    os << "norm ="  << std::setw(10) << std::right << v.v_.r()
       << ";" << std::string(5, ' ')
       << "theta =" << std::setw(10) << std::right << v.v_.theta()
       << ";" << std::string(5, ' ')
       << "phi ="   << std::setw(10) << std::right << v.v_.phi();
    return os;
  }

  LorentzVector
  getLorentzVector(const Vector & v,
                   const double e)
  {
    return LorentzVector(v.x(), v.y(), v.z(), e);
  }

  Vector
  getVector(const LorentzVector & v)
  {
    return Vector(v.x(), v.y(), v.z());
  }
}
