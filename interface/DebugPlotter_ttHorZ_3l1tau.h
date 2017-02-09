#ifndef DEBUGPLOTTER_TTHORZ_3L1TAU_H
#define DEBUGPLOTTER_TTHORZ_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/tthMEMenums.h" // EnumClassHash, hVar_3l1tau::
#include "tthAnalysis/tthMEM/interface/VariableManager_3l1tau.h" // VariableManager_3l1tau

#include <string> // std::string
#include <unordered_map> // std::unordered_map<,>
#include <functional> // std::hash<>

#include <Rtypes.h> // Double_t
#include <TFile.h> // TFile
#include <TTree.h> // TTree

namespace tthMEM
{
  /**
   * @brief Tree variables for debugging purposes
   */
  enum class hVar_3l1tau
  {
    kDz1,             // difference in energy fraction of hadronic tau decay products
    kDz2,             // difference in energy fraction of leptonic tau decay products
    kDMassHorZ,       // difference in mass of Higgs or Z
    kDnuHtauCosTheta, // difference in cosine of opening angle in hadronic tau system
    kDnuLtauCosTheta, // difference in cosine of opening angle in leptonic tau system
    kDnuHtauPhi,      // difference in revolution angle in hadronic tau system
    kDnuLtauPhi,      // difference in revolution angle in leptonic tau system
    kDmetX,           // difference in x-component of MET
    kDmetY,           // difference in y-component of MET
    kZ2,              // energy fraction of leptonic tau decay products
    kMassHtau,        // mass of hadronic tau lepton
    kMassLtau,        // mass of leptonic tau lepton
    kMassHorZ,        // mass of Higgs or Z
    kB1en,            // energy of the 1st b-quark
    kB2en,            // energy of the 2nd b-quark
    kB1RecoEn,        // energy of the 1st b-jet
    kB2RecoEn,        // energy of the 2nd b-jet
    kB1energyTF,      // TF value for the 1st b-jet energy
    kB2energyTF,      // TF value for the 2nd b-jet energy
    kTdecayJF1,       // top decay Jacobi factor for the 1st top decay branch
    kTdecayJF2,       // top decay Jacobi factor for the 2nd top decay branch
    kHtauPSF,         // phase space factor of the tau decaying hadronically
    kLtauPSF,         // phase space factor of the tau decaying leptonically
    kJacobiF,         // overall Jacobi factor
    kMETx,            // x-component of the MET
    kMETy,            // y-component of the MET
    kMETtf,           // value returned by the MET transfer function
    kMsquared,        // squared matrix element amplitude
    kXa,              // 1st Bjorken fraction
    kXb,              // 2nd Bjorken fraction
    kFlux,            // flux factor
    kProbPDF,         // probability returned by parton distribution fraction
    kProb,            // overall probability in the integrand
//--- for iteration
    First = kDz1,
    Last = kProb
  };

  class DebugString_ttHorZ_3l1tau
  {
  public:
    DebugString_ttHorZ_3l1tau();

    DebugString_ttHorZ_3l1tau(const std::string & rle_,
                              const std::string & me_,
                              unsigned leptonPermutation_,
                              unsigned bJetCombination_);

    std::string
    str() const;

    bool
    hasSameRLE(const DebugString_ttHorZ_3l1tau & other) const;

  private:
    std::string rle;
    std::string me;
    unsigned    leptonPermutation;
    unsigned    bJetCombination;
  };

  /**
   * @brief Class for plotting calculated or sampled values in eval() function
   *        of class Integrand_ttHorZ_3l1tau
   */
  class DebugPlotter_ttHorZ_3l1tau
  {
  public:
    /**
     * @brief Default constructor; doesn't write the trees anywhere
     */
    DebugPlotter_ttHorZ_3l1tau();
    /**
     * @brief Default constructor; creates a file and logs all events
     *        and its permutations
     * @param file The file where the trees will be written
     */
    DebugPlotter_ttHorZ_3l1tau(TFile * file);
    /**
     * @brief Default constructor; create a file and logs every n-th
     *        event and its permutations
     * @param file           The file where the trees will be written
     * @param debugFrequency Specifies how often every event is logged
     */
    DebugPlotter_ttHorZ_3l1tau(TFile * file,
                               unsigned debugFrequency);

    /**
      * @brief Simple destructor
      *
      * Saves unsaved information to the TTree.
      */
    ~DebugPlotter_ttHorZ_3l1tau();

    /**
     * @brief Creates a new subdirectory for the current tree
     *        and initializes the branches under the tree
     * @param dirName Directory name
     * @param vm      Variable manager (holds sampled values)
     *
     * Also explicitly points the current tree to readily created subdirectory.
     *
     * The user is expected to specify a unique dirName which e.g.
     * holds information about the run, lumi, event, permutation and
     * matrix element type (separated by anything but colon,
     * e.g. 1_123_12345678_3_tth)
     */
    void
    initialize(const DebugString_ttHorZ_3l1tau & metaInfo,
               const VariableManager_3l1tau & vm);

    /**
     * @brief Fills a given variable
     * @param var   The variable pointing to the corresponding branch
     * @param value The value used to fill the branch
     * @return Reference to this instance
     */
    DebugPlotter_ttHorZ_3l1tau &
    fill(hVar_3l1tau var,
         double value);

    /**
     * @brief Sets a tree branches to their sampled values
     * @param x The sampled values
     * @return Reference to this instance
     */
    DebugPlotter_ttHorZ_3l1tau &
    fill(const VariableManager_3l1tau & vm,
         const double * const x);

    /**
     * @brief Fill the underlying tree
     * @return Reference to this instance
     */
    DebugPlotter_ttHorZ_3l1tau &
    fill();

    /**
     * @brief Writes the current tree to the file and resets it to 0
     */
    void
    write();

  private:
    static const std::unordered_map<hVar_3l1tau, std::string, EnumClassHash> hVarMap_;
    ///< table for additional variables
    std::unordered_map<Var_3l1tau, Double_t, EnumClassHash> histSampled_;
    ///< map of sampled values
    std::unordered_map<hVar_3l1tau, Double_t, EnumClassHash> histRecod_;
    ///< map of additional variables (listed in the implementation file)
    TFile * const file_;                 ///< pointer to the file
    TTree * tree_;                       ///< pointer to the current tree
    DebugString_ttHorZ_3l1tau metaInfo_; ///< directory name under which the debug tree is written
    const unsigned debugFrequency_;      ///< debugging frequency (see constructor)
    unsigned logCounter_;                ///< counter incremented every time initialise() is called
    bool log_;                           ///< shorthand bool variable which tells whether to fill or not
    bool isFilled_;                      ///< for bookkeeping: true if at least one variable is filled

    /**
     * @brief Reset the maps holding current values to the placeholder value
     */
    void
    reset();
  };
}

#endif // DEBUGPLOTTER_TTHORZ_3L1TAU_H
