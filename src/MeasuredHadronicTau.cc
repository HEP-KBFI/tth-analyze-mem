#include "tthAnalysis/tthMEM/interface/MeasuredHadronicTau.h" // MeasuredHadronicTau, chargedPionMass
                                                              // std::sqrt(), LorentzVector

//#include "tthAnalysis/tthMEM/interface/Logger.h"

using namespace tthMEM;

MeasuredHadronicTau::MeasuredHadronicTau()
  : MeasuredObject()
  , decayMode_(-1)
{}

MeasuredHadronicTau::MeasuredHadronicTau(double pt,
                                         double eta,
                                         double phi,
                                         double mass,
                                         int decayMode)
  : MeasuredObject(pt, eta, phi, mass)
  , decayMode_(decayMode)
{
  initialize();
}

MeasuredHadronicTau::MeasuredHadronicTau(const MeasuredHadronicTau & other)
  : MeasuredObject(other)
  , decayMode_(other.decayMode_)
{
  initialize();
}

MeasuredHadronicTau &
MeasuredHadronicTau::operator=(const MeasuredHadronicTau & other)
{
  MeasuredObject::operator=(other);
  decayMode_ = other.decayMode_;

  initialize();

  return *this;
}

MeasuredHadronicTau::~MeasuredHadronicTau()
{}

int
MeasuredHadronicTau::decayMode() const
{
  return decayMode_;
}

void
MeasuredHadronicTau::initialize()
{
  preciseVisMass_ = mass_;

//  double minVisMass = 0.3; // GeV
//  double maxVisMass = 1.5; // GeV

//  if(decayMode_ == -1)
//    minVisMass = chargedPionMass;
//  else if(decayMode_ == 0)
//  {
//    minVisMass = chargedPionMass;
//    maxVisMass = chargedPionMass;
//  }

//  if(preciseVisMass_ < 0.9 * minVisMass || preciseVisMass_ > 1.1 * maxVisMass)
//    LOGINFO << "Hadronic tau ("
//      << "pt = " << pt_ << "; eta = " << eta_ << "; phi = " << phi_ << "; mass = " << mass_
//      << ") expected in the mass range [" << minVisMass << "; " << maxVisMass << "]";

//  if(preciseVisMass_ < minVisMass) preciseVisMass_ = minVisMass;
//  if(preciseVisMass_ > maxVisMass) preciseVisMass_ = maxVisMass;

  energy_ = std::sqrt(p_ * p_ + preciseVisMass_ * preciseVisMass_);
  p4_ = LorentzVector(px_, py_, pz_, energy_);

}

