#include "tthAnalysis/tthMEM/interface/MeasuredJet.h"

using namespace tthMEM;

MeasuredJet::MeasuredJet()
  : MeasuredObject()
{}

MeasuredJet::MeasuredJet(double pt,
                         double eta,
                         double phi,
                         double mass)
  : MeasuredObject(pt, eta, phi, mass)
{}

MeasuredJet::MeasuredJet(const MeasuredJet & other)
  : MeasuredObject(other)
{}

MeasuredJet &
MeasuredJet::operator=(const MeasuredJet & other)
{
  MeasuredObject::operator=(other);
  return *this;
}

MeasuredJet::~MeasuredJet()
{}
