#include "tthAnalysis/tthMEM/interface/MeasuredHadronicTau.h" // MeasuredHadronicTau, ...
  // ..., chargedPionMass, std::sqrt(), LorentzVector, minVisTauMass, maxVisTauMass
#include "tthAnalysis/tthMEM/interface/Logger.h"

using namespace tthMEM;

MeasuredHadronicTau::MeasuredHadronicTau()
  : MeasuredLepton()
  , decayMode_(-1)
{}

MeasuredHadronicTau::MeasuredHadronicTau(double pt,
                                         double eta,
                                         double phi,
                                         double mass,
                                         int charge,
                                         int decayMode)
  : MeasuredLepton(pt, eta, phi, mass, charge)
  , decayMode_(decayMode)
{
  initialize();
}

MeasuredHadronicTau::MeasuredHadronicTau(const MeasuredHadronicTau & other)
  : MeasuredLepton(other)
  , decayMode_(other.decayMode_)
{
  initialize();
}

MeasuredHadronicTau &
MeasuredHadronicTau::operator=(const MeasuredHadronicTau & other)
{
  MeasuredLepton::operator=(other);
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
  MeasuredLepton::initialize();

  preciseVisMass_ = mass_;

  double minVisMass = massVisTauMin;
  double maxVisMass = massVisTauMax;

  if(decayMode_ == -1)
    minVisMass = massChargedPion;
  else if(decayMode_ == 0)
  {
    minVisMass = massChargedPion;
    maxVisMass = massChargedPion;
  }

  if(preciseVisMass_ < 0.9 * minVisMass || preciseVisMass_ > 1.1 * maxVisMass)
    LOGWARN << "Hadronic tau (" << *this << ") expected in the mass range "
            << "[" << minVisMass << "; " << maxVisMass << "]";

  if(preciseVisMass_ < minVisMass) preciseVisMass_ = minVisMass;
  if(preciseVisMass_ > maxVisMass) preciseVisMass_ = maxVisMass;

  energy_ = std::sqrt(p_ * p_ + preciseVisMass_ * preciseVisMass_);
  p4_ = LorentzVector(px_, py_, pz_, energy_);
}

void
MeasuredHadronicTau::setBranches(TChain * t,
                                 const std::string & branchName)
{
  MeasuredLepton::setBranches(t, branchName);
  t -> SetBranchAddress(Form("%s_decayMode", branchName.c_str()), &decayMode_);
}

void
MeasuredHadronicTau::initNewBranches(TTree * t,
                                     const std::string & branchName)
{
  MeasuredLepton::initNewBranches(t, branchName);
  t -> Branch(Form("%s_decayMode", branchName.c_str()), &decayMode_,
              Form("%s_decayMode/I", branchName.c_str()));
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredHadronicTau & o)
  {
    os << static_cast<const MeasuredLepton &>(o)
       << "; decayMode = " << o.decayMode_;
    return os;
  }
}

