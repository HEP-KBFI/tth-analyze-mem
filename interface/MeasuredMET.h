#ifndef MEASUREDMET_H
#define MEASUREDMET_H

#include <Rtypes.h> // UInt_t, ULong64_t, Float_t, Long64_t
#include <TTree.h> // TTree
#include <TChain.h> // TChain
#include <TBranch.h> // TBranch

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

    TBranch * branch_pt = 0;  ///< output branch for pt
    TBranch * branch_phi = 0; ///< output branch for phi
  };
}

#endif // MEASUREDMET_H
