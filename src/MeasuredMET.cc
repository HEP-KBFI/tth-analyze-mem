#include "tthAnalysis/tthMEM/interface/MeasuredMET.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // roundToNdigits()

#include <cmath> // std::cos(), std::sin()

using namespace tthMEM;

MeasuredMET::MeasuredMET()
  : pt_(0)
  , phi_(0)
  , covMET_(TMatrixD(2, 2))
  , covMET_eigenVectors_(TMatrixD(2, 2))
  , covMET_eigenValues_(TVectorD(2))
{
  initialize();
}

MeasuredMET::MeasuredMET(double pt,
                         double phi)
  : pt_(pt)
  , phi_(phi)
  , covMET_(TMatrixD(2, 2))
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

  covMET_(0,0) = roundToNdigits(100.0); // in GeV
  covMET_(0,1) = roundToNdigits(0.0);
  covMET_(1,0) = roundToNdigits(0.0);
  covMET_(1,1) = roundToNdigits(100.0);
  calculateEigenVectorsValues();

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

void
MeasuredMET::calculateEigenVectorsValues()
{
//--- the covariance matrix is symmetric and positive semi-definite by construction
  TMatrixDSym covMET_sym(2);
  covMET_sym.SetMatrixArray(covMET_.GetMatrixArray());
  const TMatrixDSymEigen covMET_eigen(covMET_sym);
  covMET_eigenVectors_.SetMatrixArray(covMET_eigen.GetEigenValues().GetMatrixArray());
  covMET_eigenValues_.SetElements(covMET_eigen.GetEigenValues().GetMatrixArray());
//--- eigenvalues of a symmetric positive (semi-)definite matrix are always real and positive
  covMET_eigenValues_(0) = std::sqrt(covMET_eigenValues_(0));
  covMET_eigenValues_(1) = std::sqrt(covMET_eigenValues_(1));
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MeasuredMET & o)
  {
    os << "pt = " << o.pt_ << "; phi = " << o.phi_ << "; "
       << "covMET = {[[" << o.covMET_(0, 0) << "][" << o.covMET_(0, 1) << "]]"
                 << "[[" << o.covMET_(1, 0) << "][" << o.covMET_(1, 1) << "]]};"
       << "covMET eigenVector 1 = { " << o.covMET_eigenVectors_(0, 0) << ", "
                                      << o.covMET_eigenVectors_(1, 0) << " }; "
       << "covMET eigenVector 2 = { " << o.covMET_eigenVectors_(1, 0) << ", "
                                      << o.covMET_eigenVectors_(1, 1) << " }; "
       << "covMET eigenvalues = { " << o.covMET_eigenValues_(0) << ", "
                                    << o.covMET_eigenValues_(0) << " }";
    return os;
  }
}
