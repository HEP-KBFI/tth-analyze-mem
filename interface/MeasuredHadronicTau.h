#ifndef MEASUREDHADRONICTAU_H
#define MEASUREDHADRONICTAU_H

#include "tthAnalysis/tthMEM/interface/MeasuredObject.h"

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
    : public MeasuredObject
  {
  public:
    MeasuredHadronicTau();
    MeasuredHadronicTau(double pt,
                        double eta,
                        double phi,
                        double mass,
                        int decayMode);
    MeasuredHadronicTau(const MeasuredHadronicTau & other);
    MeasuredHadronicTau & operator=(const MeasuredHadronicTau & other);
    ~MeasuredHadronicTau();

    int decayMode() const;

    virtual void initialize() override;

  protected:
    int decayMode_;         ///< decay mode of hadronic tau lepton
    double preciseVisMass_; ///< precise visible mass of tau decay products in the lab frame
  };
}

#endif // MEASUREDHADRONICTAU_H
