#ifndef MEASUREDLEPTON_H
#define MEASUREDLEPTON_H

#include "tthAnalysis/tthMEM/interface/MeasuredObject.h"

namespace tthMEM
{
  /**
   * @brief Specialized class for storing measurements of electrons and muons
   *
   * Most of the functionality is inherited from MeasuredObject
   *
   * @see MeasuredObject
   */
  class MeasuredLepton
    : public MeasuredObject
  {
  public:
    MeasuredLepton();
    MeasuredLepton(double pt,
                   double eta,
                   double phi,
                   double mass,
                   int charge);
    MeasuredLepton(const MeasuredLepton & other);
    MeasuredLepton & operator=(const MeasuredLepton & other);
    ~MeasuredLepton();

    int charge() const;
  private:
    int charge_; ///< charge of electron or muon
  };
}

#endif // MEASUREDLEPTON_H
