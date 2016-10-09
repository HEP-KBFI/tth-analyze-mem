#include "tthAnalysis/tthMEM/interface/MeasuredLepton.h"

using namespace tthMEM;

MeasuredLepton::MeasuredLepton()
  : MeasuredObject()
  , charge_(0)
{}

MeasuredLepton::MeasuredLepton(const MeasuredObject & object,
                               int charge)
  : MeasuredObject(object)
  , charge_(charge)
{}

MeasuredLepton::MeasuredLepton(double pt,
                               double eta,
                               double phi,
                               double mass,
                               int charge)
  : MeasuredObject(pt, eta, phi, mass)
  , charge_(charge)
{}

MeasuredLepton::MeasuredLepton(const MeasuredLepton & other)
  : MeasuredObject(other)
  , charge_(other.charge_)
{}

MeasuredLepton &
MeasuredLepton::operator=(const MeasuredLepton & other)
{
  MeasuredObject::operator=(other);
  charge_ = other.charge_;

  return *this;
}

MeasuredLepton::~MeasuredLepton()
{}

int
MeasuredLepton::charge() const
{
  return charge_;
}

void
MeasuredLepton::setBranches(TChain * t,
                            const std::string & branchName)
{
  MeasuredObject::setBranches(t, branchName);
  t -> SetBranchAddress(Form("%s_charge", branchName.c_str()), &charge_);
}

void
MeasuredLepton::initNewBranches(TTree * t,
                                const std::string & branchName)
{
  MeasuredObject::initNewBranches(t, branchName);
  t -> Branch(Form("%s_charge", branchName.c_str()), &charge_,
              Form("%s_charge/I", branchName.c_str()));
}

void
MeasuredLepton::setNeutral()
{
  charge_ = 0;
}

void
MeasuredLepton::flipCharge()
{
  charge_ *= -1;
}

namespace tthMEM
{
  MeasuredLepton
  operator+(const MeasuredLepton & lhs,
            const MeasuredLepton & rhs)
  {
    return MeasuredLepton(static_cast<MeasuredObject>(lhs) + static_cast<MeasuredObject>(rhs),
                          lhs.charge_ + rhs.charge_);
  }

  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredLepton & o)
  {
    os << static_cast<const MeasuredObject &>(o)
       << "; charge = " << o.charge_;
    return os;
  }
}
