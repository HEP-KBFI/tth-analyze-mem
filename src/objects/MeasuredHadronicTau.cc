#include "tthAnalysis/tthMEM/interface/objects/MeasuredHadronicTau.h" // MeasuredHadronicTau
#include "tthAnalysis/tthMEM/interface/general/constants.h" // constants::
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGWARN

#include <TTree.h> // TTree
#include <TBranch.h> // TBranch

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

MeasuredHadronicTau::MeasuredHadronicTau(const LorentzVector & lv,
                                         int charge,
                                         int decayMode)
  : MeasuredLepton(lv, charge)
  , decayMode_(decayMode)
{}

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
  preciseVisMass_ = mass_;

  double minVisMass = constants::massVisTauMin;
  double maxVisMass = constants::massVisTauMax;

  if(decayMode_ == -1)
    minVisMass = constants::massChargedPion;
  else if(decayMode_ == 0)
  {
    minVisMass = constants::massChargedPion;
    maxVisMass = constants::massChargedPion;
  }

  if(preciseVisMass_ < 0.9 * minVisMass || preciseVisMass_ > 1.1 * maxVisMass)
    LOGWARN << "Hadronic tau (" << *this << ") expected in the mass range "
            << '[' << minVisMass << "; " << maxVisMass << ']';

  if(preciseVisMass_ < minVisMass) preciseVisMass_ = minVisMass;
  if(preciseVisMass_ > maxVisMass) preciseVisMass_ = maxVisMass;

  // recompute the momenta
  mass_ = preciseVisMass_;
  MeasuredLepton::initialize();
}

void
MeasuredHadronicTau::setBranches(TTree * t,
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

