#include "tthAnalysis/tthMEM/interface/MeasuredLepton.h"

using namespace tthMEM;

MeasuredLepton::MeasuredLepton()
  : MeasuredObject()
  , charge_(0)
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
