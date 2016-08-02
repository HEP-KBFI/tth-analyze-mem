#include "tthAnalysis/tthMEM/interface/MeasuredObject.h" // MeasuredObject, roundToNdigits()

#include <cmath> // std::cos(), std::sin(), std::cosh(), std::sinh(), std::sqrt()

#include <TString.h> // Form()

using namespace tthMEM;

MeasuredObject::MeasuredObject()
  : pt_(0.)
  , eta_(0.)
  , phi_(0.)
  , mass_(0.)
{
  initialize();
}

MeasuredObject::MeasuredObject(double pt,
                               double eta,
                               double phi,
                               double mass)
  : pt_(pt)
  , eta_(eta)
  , phi_(phi)
  , mass_(mass)
{
  initialize();
}

MeasuredObject::MeasuredObject(const MeasuredObject & other)
  : pt_(other.pt_)
  , eta_(other.eta_)
  , phi_(other.phi_)
  , mass_(other.mass_)
{
  initialize();
}

MeasuredObject &
MeasuredObject::operator=(const MeasuredObject & other)
{
  pt_ = other.pt_;
  eta_ = other.eta_;
  phi_ = other.phi_;
  mass_ = other.mass_;

  initialize();

  return *this;
}

MeasuredObject::~MeasuredObject()
{}

double
MeasuredObject::pt() const
{
  return pt_;
}

double
MeasuredObject::eta() const
{
  return eta_;
}

double
MeasuredObject::phi() const
{
  return phi_;
}

double
MeasuredObject::mass() const
{
  return mass_;
}

double
MeasuredObject::energy() const
{
  return energy_;
}

double
MeasuredObject::px() const
{
  return px_;
}

double
MeasuredObject::py() const
{
  return py_;
}

double
MeasuredObject::pz() const
{
  return pz_;
}

double
MeasuredObject::p() const
{
  return p_;
}

double
MeasuredObject::theta() const
{
  return theta_;
}

double
MeasuredObject::cosPhi_sinTheta() const
{
  return cosPhi_sinTheta_;
}

double
MeasuredObject::sinPhi_sinTheta() const
{
  return sinPhi_sinTheta_;
}

double
MeasuredObject::cosTheta() const
{
  return cosTheta_;
}

const LorentzVector &
MeasuredObject::p4() const
{
  return p4_;
}

const Vector &
MeasuredObject::p3() const
{
  return p3_;
}

void
MeasuredObject::setBranches(TChain * t,
                            const std::string & branchName)
{
  t -> SetBranchAddress(Form("%s_pt",   branchName.c_str()), &pt_);
  t -> SetBranchAddress(Form("%s_eta",  branchName.c_str()), &eta_);
  t -> SetBranchAddress(Form("%s_phi",  branchName.c_str()), &phi_);
  t -> SetBranchAddress(Form("%s_mass", branchName.c_str()), &mass_);
}

void
MeasuredObject::initNewBranches(TTree * t,
                                const std::string & branchName)
{
  t -> Branch(Form("%s_pt",   branchName.c_str()), &pt_,   Form("%s_pt/D",   branchName.c_str()));
  t -> Branch(Form("%s_eta",  branchName.c_str()), &eta_,  Form("%s_eta/D",  branchName.c_str()));
  t -> Branch(Form("%s_phi",  branchName.c_str()), &phi_,  Form("%s_phi/D",  branchName.c_str()));
  t -> Branch(Form("%s_mass", branchName.c_str()), &mass_, Form("%s_mass/D", branchName.c_str()));
}

void
MeasuredObject::initialize()
{
  pt_   = roundToNdigits(pt_);
  eta_  = roundToNdigits(eta_);
  phi_  = roundToNdigits(phi_);
  mass_ = roundToNdigits(mass_);

  p_ = pt_ * std::cosh(eta_);
  px_ = pt_ * std::cos(phi_);
  py_ = pt_ * std::sin(phi_);
  pz_ = pt_ * std::sinh(eta_);
  p3_ = Vector(px_, py_, pz_);

  theta_ = p3_.Theta();
  cosPhi_sinTheta_ = std::cos(phi_) * std::sin(theta_);
  sinPhi_sinTheta_ = std::sin(phi_) * std::sin(theta_);
  cosTheta_ = std::cos(theta_);

  energy_ = std::sqrt(p_ * p_ + mass_ * mass_);
  p4_ = LorentzVector(px_, py_, pz_, energy_);

  return;
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
                     const MeasuredObject & o)
  {
    os << "pt = " << o.pt_
       << "; eta = " << o.eta_
       << "; phi = " << o.phi_
       << "; mass = " << o.mass_;
    return os;
  }
}

