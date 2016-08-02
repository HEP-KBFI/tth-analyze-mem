#include "tthAnalysis/tthMEM/interface/MeasuredMET.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // roundToNdigits()

#include <cmath> // std::cos(), std::sin()

using namespace tthMEM;

MeasuredMET::MeasuredMET()
  : pt_(0)
  , phi_(0)
{
  initialize();
}

MeasuredMET::MeasuredMET(double pt,
                         double phi)
  : pt_(pt)
  , phi_(phi)
{
  initialize();
}

double
MeasuredMET::pt() const
{
  return pt_;
}

double
MeasuredMET::phi() const
{
  return phi_;
}

double
MeasuredMET::px() const
{
  return px_;
}

double
MeasuredMET::py() const
{
  return py_;
}

void
MeasuredMET::initialize()
{
  pt_   = roundToNdigits(pt_);
  phi_  = roundToNdigits(phi_);

  px_ = pt_ * std::cos(phi_);
  py_ = pt_ * std::sin(phi_);
}

void
MeasuredMET::setBranches(TChain * t)
{
  t -> SetBranchAddress("met_pt",  &pt_);
  t -> SetBranchAddress("met_phi", &phi_);
}

void
MeasuredMET::initNewBranches(TTree * t)
{
  branch_pt  = t -> Branch("met_pt",  &pt_,  "met_pt/D");
  branch_phi = t -> Branch("met_phi", &phi_, "met_phi/D");
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredMET & o)
  {
    os << "pt = " << o.pt_
       << "; phi = " << o.phi_;
    return os;
  }
}
