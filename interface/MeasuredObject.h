#ifndef MEASUREDOBJECT_H
#define MEASUREDOBJECT_H

#include <string> // std::string
#include <ostream> // std::ostream

#include <TChain.h> // TChain
#include <TTree.h> // TTree
#include <TBranch.h> // TBranch

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // tthMEM::LorentzVector(), tthMEM::Vector()

namespace tthMEM
{
  /**
   * @brief Basic class for storing measured data needed in MEM
   */
  class MeasuredObject
  {
  public:
    /**
     * @brief Default constructor
     * Initializes pt_, eta_, phi_, mass_ to zero and calls initialize()
     */
    MeasuredObject();
    /**
     * @brief Default constructor which initializes the variables to their actual value
     * @param pt   pT of measured momentum in the lab frame
     * @param eta  pseudo-rapidity of measured momentum in the lab frame
     * @param phi  azimuthal angle of measured momentum in the lab frame
     * @param mass measured mass in the lab frame
     * @see initialize()
     */
    MeasuredObject(double pt,
                   double eta,
                   double phi,
                   double mass);
    /**
     * @brief Copy constructor
     * @param other Other object to copy pt_, eta_, phi_, mass_ from
     * @see initialize()
     */
    MeasuredObject(const MeasuredObject & other);
    /**
     * @brief Copy assignment operator
     * @param other Other object to copy pt_, eta_, phi_, mass_ from
     * @return Reference to the current object
     * @see initialize()
     */
    MeasuredObject & operator=(const MeasuredObject & other);
    ~MeasuredObject();

    /* simple getters */
    double pt() const;
    double eta() const;
    double phi() const;
    double mass() const;

    double energy() const;
    double px() const;
    double py() const;
    double pz() const;

    double p() const;
    double theta() const;
    double cosPhi_sinTheta() const;
    double sinPhi_sinTheta() const;
    double cosTheta() const;

    const LorentzVector & p4() const;
    const Vector & p3() const;

    virtual void
    initialize(); ///< sets all momentum components but pt_, eta_, phi_, mass_

    virtual void
    setBranches(TChain * t,
                const std::string & branchName);
    ///< associates the pt, eta, phi and mass with an old input tree

    virtual void
    initNewBranches(TTree * t,
                    const std::string & branchName);
    ///< associates the pt, eta, phi and mass with a new output tree

    friend std::ostream &
    operator<<(std::ostream & os,
               const MeasuredObject & o);
    ///< prints the pt, eta, phi and mass to ostream

  protected:
    double pt_;     ///< pT of measured momentum in the lab frame
    double eta_;    ///< pseudo-rapidity of measured momentum in the lab frame
    double phi_;    ///< azimuthal angle of measured momentum in the lab frame
    double mass_;   ///< measured (visible) mass in the lab frame

    double energy_; ///< (visible) energy in the lab frame
    double px_;     ///< x-component of measured momentum in the lab frame
    double py_;     ///< y-component of measured momentum in the lab frame
    double pz_;     ///< z-component of measured momentum in the lab frame

    double p_;      ///< module of measured momentum in the lab frame
    double theta_;  ///< polar angle of measured momentum in the lab frame
    double cosPhi_sinTheta_; ///< spherical x-component of measured momentum in the lab frame
    double sinPhi_sinTheta_; ///< spherical y-component of measured momentum in the lab frame
    double cosTheta_;        ///< spherical z-component of measured momentum in the lab frame

    LorentzVector p4_; ///< measured 4-momentum in the lab frame
    Vector p3_;        ///< measured 3-momentum in the lab frame

    TBranch * branch_pt   = 0; ///< output branch for pt
    TBranch * branch_eta  = 0; ///< output branch for eta
    TBranch * branch_phi  = 0; ///< output branch for phi
    TBranch * branch_mass = 0; ///< output branch for mass
  };
}

#endif // MEASUREDOBJECT_H
