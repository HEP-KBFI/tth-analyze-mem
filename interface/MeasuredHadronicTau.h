#ifndef MEASUREDHADRONICTAU_H
#define MEASUREDHADRONICTAU_H

#include "tthAnalysis/tthMEM/interface/MeasuredLepton.h"

namespace tthMEM
{
  /**
   * @brief Specialized class for storing measurements of hadronic taus
   *
   * Most of the functionality is inherited from MeasuredObject.
   * However, mass a hadronic tau might need a special treatment.
   *
   * @see MeasuredObject
   */
  class MeasuredHadronicTau
    : public MeasuredLepton
  {
  public:
    MeasuredHadronicTau();
    MeasuredHadronicTau(double pt,
                        double eta,
                        double phi,
                        double mass,
                        int charge,
                        int decayMode);
    MeasuredHadronicTau(const LorentzVector & lv,
                        int charge,
                        int decayMode);
    MeasuredHadronicTau(const MeasuredHadronicTau & other);
    MeasuredHadronicTau & operator=(const MeasuredHadronicTau & other);
    ~MeasuredHadronicTau();

    int decayMode() const;

    virtual void
    initialize() final override;

    virtual void
    setBranches(TTree * t,
                const std::string & branchName) final override;

    virtual void
    initNewBranches(TTree * t,
                    const std::string & branchName) final override;

    friend std::ostream &
    operator<<(std::ostream & os,
               const MeasuredHadronicTau & o);
    ///< prints the pt, eta, phi, mass and decayMode to ostream

  protected:
    int decayMode_;         ///< decay mode of hadronic tau lepton
    double preciseVisMass_; ///< precise visible mass of tau decay products in the lab frame
  };
}

#endif // MEASUREDHADRONICTAU_H
