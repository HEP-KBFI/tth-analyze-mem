#ifndef MEASUREDMET_H
#define MEASUREDMET_H

#include <Rtypes.h> // UInt_t, ULong64_t, Float_t, Long64_t
#include <TTree.h> // TTree
#include <TChain.h> // TChain
#include <TBranch.h> // TBranch
#include <TMatrixD.h> // TMatrixD
#include <TMatrixDSym.h> // TMatrixDSym
#include <TVectorD.h> // TVectorD

namespace tthMEM
{
  /**
   * @brief Class for storing measured MET (missing energy)
   */
  class MeasuredMET
  {
  public:
    /**
     * @brief Default constructor
     * Initializes pt_ and phi_ to zero and calls initialize()
     */
    MeasuredMET();

    /**
     * @brief Default constructor which initializes pt and phi to their actual value
     * @param pt  pT of MET
     * @param phi phi of MET
     */
    MeasuredMET(double pt,
                double phi);

    /* simple getters */
    double pt() const;
    double phi() const;

    double px() const;
    double py() const;

    const TMatrixDSym & covMET() const;

    void
    initialize(); ///< truncates trailing numbers; sets px_ and py_

    void
    setBranches(TChain * t);
    ///< associates the pt and phi with an old input tree

    void
    initNewBranches(TTree * t);
    ///< associates the pt and phi with a new output tree

    friend std::ostream &
    operator<<(std::ostream & os,
               const MeasuredMET & o);
    ///< prints the pt and phi to ostream

  private:
    double pt_;  ///< pT of measured momentum in the lab frame
    double phi_; ///< azimuthal angle of measured momentum in the lab frame

    double px_; ///< x-component of measured momentum in the lab frame
    double py_; ///< y-component of measured momentum in the lab frame

    TMatrixDSym covMET_; ///< soon we'll have a separate branch for it in our Ntuples
                         ///< the covariance matrix is symmetric by construction
    TMatrixD    covMET_eigenVectors_; ///< eigenvectors of covMET_
    TVectorD    covMET_eigenValues_;  ///< eigenvalues of covMET_

    TBranch * branch_pt = 0;  ///< output branch for pt
    TBranch * branch_phi = 0; ///< output branch for phi

    void
    calculateEigenVectorsValues();
  };
}

#endif // MEASUREDMET_H
