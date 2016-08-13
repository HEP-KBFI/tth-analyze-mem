#include "tthAnalysis/tthMEM/interface/MeasuredMET.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // roundToNdigits()
#include "tthAnalysis/tthMEM/interface/Logger.h" // Logger::getFloatPrecision()

#include <cmath> // std::cos(), std::sin()
#include <sstream> // std::stringstream
#include <string> // std::string
#include <iomanip> // std::setw(), std::setprecision()

#include <TMatrixDSymEigen.h> // TMatrixDSymEigen

using namespace tthMEM;

MeasuredMET::MeasuredMET()
  : pt_(0)
  , phi_(0)
  , covMET_(TMatrixDSym(2))
  , covMET_eigenVectors_(TMatrixD(2, 2))
  , covMET_eigenValues_(TVectorD(2))
{
  initialize();
}

MeasuredMET::MeasuredMET(double pt,
                         double phi)
  : pt_(pt)
  , phi_(phi)
  , covMET_(TMatrixDSym(2, 2))
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

const TMatrixDSym &
MeasuredMET::covMET() const
{
  return covMET_;
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
  const TMatrixDSymEigen covMET_eigen(covMET_);
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
    std::stringstream ss;
    ss << std::fixed << std::setprecision(Logger::getFloatPrecision());
    ss << "pt = " << o.pt_ << "; phi = " << o.phi_ << "; ";
    for(unsigned i = 0; i < 2; ++i)
      ss << "covMET eigenVector " << (i + 1) << " = { "
         << o.covMET_eigenVectors_(0, i) << ", "
         << o.covMET_eigenVectors_(1, i) << " }; ";

//--- print the covariance matrix
    const unsigned minIndent = 27;
    const unsigned actualIndent = std::max(minIndent, 16u);
    const unsigned nameIndent = actualIndent - 7;
    const unsigned numberFieldWidth = 12;
    ss << "\n" << std::string(actualIndent, ' ') << "|";
    for(unsigned i = 0; i < 2; ++i)
      ss << std::setw(numberFieldWidth) << i << " |";
    ss << "\n" << std::setw(nameIndent) << "covMET =    "
       << std::string(2 * (numberFieldWidth + 2) + 10, '-') << "\n";
    for(unsigned i = 0; i < 2; ++i)
    {
      ss << std::setw(actualIndent - 1) << i << " |";
      for(unsigned j = 0; j < 2; ++j)
        ss << std::setw(numberFieldWidth) << o.covMET_(i, j) << " |";
      ss << "\n";
    }

    os << ss.str();
    return os;
  }
}
