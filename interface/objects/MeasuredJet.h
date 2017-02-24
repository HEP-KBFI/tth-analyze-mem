#ifndef MEASUREDJET_H
#define MEASUREDJET_H

#include "tthAnalysis/tthMEM/interface/objects/MeasuredObject.h"

namespace tthMEM
{
  /**
   * @brief Specialized class for storing measurements of jets
   *
   * Basically identical to MeasuredObject; this class created for semantic sugarcoating
   *
   * @see MeasuredObject
   */
  class MeasuredJet
    : public MeasuredObject
  {
  public:
    MeasuredJet();
    MeasuredJet(double pt,
                double eta,
                double phi,
                double mass);
    MeasuredJet(const LorentzVector & lv);
    MeasuredJet(const MeasuredJet & other);
    MeasuredJet & operator=(const MeasuredJet & other);
    ~MeasuredJet();
  };
}

#endif // MEASUREDJET_H
