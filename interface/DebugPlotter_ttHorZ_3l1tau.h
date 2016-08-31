#ifndef DEBUGPLOTTER_TTHORZ_3L1TAU_H
#define DEBUGPLOTTER_TTHORZ_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // EnumClassHash

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
  enum class hVar
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
    kMsquared, // squared matrix element amplitude
    kProb,     // overall probability in the integrand
//--- for iteration
    First = kZ1,
    Last = kProb
  };
  ///< histogram variable

  /**
   * @brief Class for plotting calculated or sampled values in eval() function
   *        of class Integrand_ttHorZ_3l1tau
   */
  class DebugPlotter_ttHorZ_3l1tau
  {
  public:
    /**
     * @brief Default constructor; doesn't write the histograms anywhere
     */
    DebugPlotter_ttHorZ_3l1tau();
    /**
     * @brief Default constructor; creates a file and logs all events
     *        and its permutations
     * @param file The file where the histograms will be written
     */
    DebugPlotter_ttHorZ_3l1tau(TFile * file);
    /**
     * @brief Default constructor; create a file and logs every n-th
     *        event and its permutations
     * @param file           The file where the histograms will be written
     * @param debugFrequency Specifies how often every event is logged
     *
     * Let's say that the debugFrequency is n. In this case all events
     * and permutations that are in [n, n + debugRange_] with period n
     * are dumped to the provided file. The debugRange_ variable is hard-
     * coded to 8 in the respective implementation file.
     */
    DebugPlotter_ttHorZ_3l1tau(TFile * file,
                               unsigned debugFrequency);

    /**
     * @brief Creates a new subdirectory for the set of histograms
     *        and initializes the histograms under the directory
     * @param dirName Directory name
     *
     * Also explicitly points the histograms to readily created subdirectory.
     *
     * The user is expected to specify a unique dirName which e.g.
     * holds information about the run, lumi, event, permutation and
     * matrix element type (separated by anything but colon,
     * e.g. 1_123_12345678_3_tth)
     */
    void
    initialize(const std::string & dirName);

    /**
     * @brief Fills a given variable
     * @param var   The variable pointing to the corresponding histogram
     * @param value The value used to filled the histogram (nominal weight)
     * @return Reference to this instance
     */
    DebugPlotter_ttHorZ_3l1tau &
    fill(hVar var,
         double value);

    /**
     * @brief Writes the histograms to the file and resets them to 0
     */
    void
    write();

  private:
    std::unordered_map<hVar, TH1D *, EnumClassHash> histograms_;
    ///< map of all histograms
    TFile * const file_; ///< pointer to the file
    const unsigned debugFrequency_; ///< debugging frequency (see constructor)
    const unsigned debugRange_; ///< debugging range (== 8; see constructor)
    unsigned logCounter_; ///< counter incremented every time initialise() is called
    bool log_; ///< shorthand boolean variable which tells whether to fill or not
  };
}

#endif // DEBUGPLOTTER_TTHORZ_3L1TAU_H
