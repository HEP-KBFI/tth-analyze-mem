#ifndef MEASUREDMET_H
#define MEASUREDMET_H

#include <TMatrixDSym.h> // TMatrixDSym
#include <TVectorD.h> // TVectorD

class TTree;
class TBranch;

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
     * @brief Default constructor which initializes pt and phi to their actual values
     * @param pt  pT of MET
     * @param phi phi of MET
     */
    MeasuredMET(double pt,
                double phi);

    /**
     * @brief Default constructor which initializes MET pt, phi and covariance matrix
     *        to their actual values
     * @param pt        MET pT
     * @param phi       MET phi
     * @param covMET_XX (0, 0) component of the MET covariance matrix
     * @param covMET_XY (0, 1) and (1, 0) component of the MET covariance matrix
     * @param covMET_YY (1, 1) component of the MET covariance matrix
     */
    MeasuredMET(double pt,
                double phi,
                double covMET_XX,
                double covMET_XY,
                double covMET_YY);

    /* simple getters */
    double pt() const;
    double phi() const;

    double px() const;
    double py() const;

    const TMatrixDSym & covMET() const;

    void
    initialize(); ///< truncates trailing numbers; sets px_ and py_

    void
    setBranches(TTree * t);
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

    double covMET_XX_; ///< (0, 0) component of MET covariance matrix
    double covMET_XY_; ///< (0, 1) and (1, 0) component of MET covariance matrix
    double covMET_YY_; ///< (1, 1) component of MET covariance matrix

    TMatrixDSym covMET_; ///< the covariance matrix is symmetric by construction
    TMatrixD    covMET_eigenVectors_; ///< eigenvectors of covMET_
    TVectorD    covMET_eigenValues_;  ///< eigenvalues of covMET_

    TBranch * branch_pt        = nullptr; ///< output branch for pt
    TBranch * branch_phi       = nullptr; ///< output branch for phi
    TBranch * branch_covMET_XX = nullptr; ///< output branch for (0, 0) component of MET covariance matrix
    TBranch * branch_covMET_XY = nullptr; ///< output branch for (0, 1) and (1, 0) component of MET covariance matrix
    TBranch * branch_covMET_YY = nullptr; ///< output branch for (1, 1) component of MET covariance matrix

    void
    calculateEigenVectorsValues();

    bool
    branchExists(TTree * tree,
                 const std::string & branchName) const;
  };
}

#endif // MEASUREDMET_H
