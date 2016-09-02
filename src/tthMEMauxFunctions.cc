#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGERR

#include <cfenv> // std::fesetround(), FE_TONEAREST
#include <cassert> // assert()

#include "FWCore/ParameterSet/interface/FileInPath.h" // edm::FileInPath

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

  double
  roundToNdigits(double x,
                 int n)
  {
    if(x == 0.) return 0.;
    const double p = std::pow(10., n - std::ceil(std::log10(std::fabs(x))));
    return std::round(p * x) / p;
  }

  unsigned
  roundToNearestUInt(double x)
  {
    std::fesetround(FE_TONEAREST);
    return static_cast<unsigned>(std::rint(x));
  }

  void
  setMGmomentum(const LorentzVector & lv,
                double * mg)
  {
    mg[0] = lv.e();
    mg[1] = lv.px();
    mg[2] = lv.py();
    mg[3] = lv.pz();
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

  std::string
  findFile(const std::string & fileName)
  {
    edm::FileInPath inputFile(fileName);
    if(inputFile.fullPath() == "")
    {
      LOGERR << "Error: cannot find file = " << fileName;
      assert(0);
    }
    return inputFile.fullPath();
  }

  namespace functions
  {
    double
    l2(const std::vector<double> & v,
       double shiftValue)
    {
      return l2(v, 0, v.size(), shiftValue);
    }

    double
    l2(const std::vector<double> & v,
       unsigned shiftFromBegin,
       double shiftValue)
    {
      assert(shiftFromBegin <= v.size());
      return l2(v, 0, shiftFromBegin, shiftValue);
    }

    double
    l2(const std::vector<double> & v,
       unsigned shiftFromBegin_begin,
       unsigned shiftFromBegin_end,
       double shiftValue)
    {
      assert(shiftFromBegin_begin <= shiftFromBegin_end &&
             shiftFromBegin_begin <= v.size() &&
             shiftFromBegin_end <= v.size());
      return std::sqrt(std::accumulate(
        v.begin() + shiftFromBegin_begin, v.end() + shiftFromBegin_end, 0.,
        [shiftValue](double sum,
                     double element) -> double
        {
          return sum + pow2(element - shiftValue);
        }
      ));
    }

    double
    avg(const std::vector<double> & v)
    {
      return avg(v, 0, v.size());
    }

    double
    avg(const std::vector<double> & v,
        unsigned shiftFromBegin)
    {
      return avg(v, 0, shiftFromBegin);
    }

    double
    avg(const std::vector<double> & v,
        unsigned shiftFromBegin_begin,
        unsigned shiftFromBegin_end)
    {
      assert(shiftFromBegin_begin <= shiftFromBegin_end &&
             shiftFromBegin_begin <= v.size() &&
             shiftFromBegin_end <= v.size());
      if(! v.size() || shiftFromBegin_begin == shiftFromBegin_end) return 0.;
      return std::accumulate(v.begin() + shiftFromBegin_begin,
                             v.begin() + shiftFromBegin_end,
                             0.) /
             (shiftFromBegin_end - shiftFromBegin_begin);
    }

    double
    stdev(const std::vector<double> & v, double average)
    {
      return stdev(v, 0, v.size(), average);
    }

    double
    stdev(const std::vector<double> & v,
          unsigned shiftFromBegin,
          double average)
    {
      return stdev(v, 0, shiftFromBegin, average);
    }

    double
    stdev(const std::vector<double> & v,
          unsigned shiftFromBegin_begin,
          unsigned shiftFromBegin_end,
          double average)
    {
      const unsigned range = shiftFromBegin_end - shiftFromBegin_begin;
      return l2(v, shiftFromBegin_begin, shiftFromBegin_end, average) /
             std::sqrt(range >= 2 ? range * (range - 1) : 1.);
    }
  }
}
