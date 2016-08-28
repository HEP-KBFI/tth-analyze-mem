#ifndef DEBUGPLOTTER_TTHORZ_3L1TAU_H
#define DEBUGPLOTTER_TTHORZ_3L1TAU_H

#include <string> // std::string
#include <unordered_map> // std::unordered_map<,>
#include <functional> // std::hash<>

#include <TFile.h> // TFile
#include <TH1D.h> // TH1D

namespace tthMEM
{
  /**
   * @brief Histogram variable
   */
  enum hVar
  {
    kZ1,       // energy fraction of hadronic tau decay products
    kZ2,       // energy fraction of leptonic tau decay products
    kMassHtau, // mass of hadronic tau lepton
    kMassLtau, // mass of leptonic tau lepton
    kMassHorZ, // mass of Higgs or Z
    kB1en,     // energy of the 1st b-quark
    kB2en,     // energy of the 2nd b-quark
    kB1RecoEn, // energy of the 1st b-jet
    kB2RecoEn, // energy of the 2nd b-jet
    kMETpull,  // exponent in hadronic recoil TF
    kXa,       // 1st Bjorken variable
    kXb,       // 2nd Bjorken variable
    kDrecX,    // x component of the hadronic recoil difference
    kDrecY,    // y component of the hadronic recoil difference
    kMsquared, // squared matrix element amplitude
    kProb      // overall probability in the integrand
  };
  ///< histogram variable

  /**
   * @brief Class for plotting calculated or sampled values in eval() function
   *        of class Integrand_ttHorZ_3l1tau
   */
  class DebugPlotter_ttHorZ_3l1tau
  {
  public:
    DebugPlotter_ttHorZ_3l1tau();
    DebugPlotter_ttHorZ_3l1tau(TFile * file);
    DebugPlotter_ttHorZ_3l1tau(TFile * file,
                               unsigned debugFrequency);

    void
    initialize(const std::string & dirName);

    DebugPlotter_ttHorZ_3l1tau &
    fill(hVar var,
         double value);

    void
    write();

  private:
    std::unordered_map<hVar, TH1D *, std::hash<int>> histograms_;
    TFile * const file_;
    const unsigned debugFrequency_;
    const unsigned debugRange_;
    unsigned logCounter_;
    bool log_;
  };
}

#endif // DEBUGPLOTTER_TTHORZ_3L1TAU_H
