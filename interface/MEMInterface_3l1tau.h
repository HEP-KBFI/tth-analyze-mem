#ifndef MEMINTERFACE_3L1TAU_H
#define MEMINTERFACE_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/MeasuredLepton.h" // tthMEM::MeasuredLepton
#include "tthAnalysis/tthMEM/interface/MeasuredJet.h" // tthMEM::MeasuredJet
#include "tthAnalysis/tthMEM/interface/MeasuredHadronicTau.h" // tthMEM::MeasuredHadronicTau
#include "tthAnalysis/tthMEM/interface/MeasuredMET.h" // tthMEM::MeasuredMET
#include "tthAnalysis/tthMEM/interface/MEMOutput_3l1tau.h" // MEMOutput_3l1tau
#include "tthAnalysis/tthMEM/interface/MEM_ttHorZ_3l1tau.h" // tthMEM::MEM_ttHorZ_3l1tau
#include "tthAnalysis/tthMEM/interface/LikelihoodRatio_3l1tau.h" // LikelihoodRatio_3l1tau
#include "tthAnalysis/tthMEM/interface/VariableManager_3l1tau.h" // VariableManager_3l1tau

#define MIN_NOF_RECO_JETS 2
#define MAX_NOF_RECO_JETS 3

using namespace tthMEM;

/**
 * @file The interface for MEM on 3l1tau channel
 *
 * @todo Create a test example for this file
 */

struct MEMInterface_3l1tau
{
  /**
   * @brief Default constructor
   * Sets up the MEM environment
   */
  MEMInterface_3l1tau();

  /**
   * @brief Initializes private data members
   */
  void
  initialize();

  /**
   * @brief Computes MEM score for different hypotheses relevant to 3l1tau channel,
   *        and also likelihood ratio w/ systematics
   * @param selectedJets        All selected jets that passed the analysis cuts (must be >= 2)
   *                            MeasuredJet construct takes 4 floats as an argument:
   *                              { pt, eta, phi, mass }
   * @param leadingLepton       The leading lepton (by pT)
   *                            MeasuredLeptonc constructor takes 4 floats and 1 integer as an argument:
   *                              { pt, eta, phi, mass, charge }
   * @param subLeadingLepton    The subleading lepton (by pT)
   * @param thirdLepton         The subsubleading lepton (by pT)
   * @param selectedLeptons     Exactly three selected leptons that passed the analysis cuts
   *
   * @param selectedHadronicTau One selected hadronic tau
   *                            MeasuredHadronicTau constructor takes 4 floats and 2 integers as an argument:
   *                              { pt, eta, phi, mass, charge, decayMode }
   * @param measuredMET         Measured MET
   * @param run                 Run number (optional)
   * @param lumi                Lumi section number (optional)
   * @param evt                 Event number (optional)
   * @return                    MEMOutput_3l1tau
   *                            A struct containing MEM score for different hypotheses, likelihood ratio, and
   *                            errors thereof
   */
  MEMOutput_3l1tau
  operator()(const std::vector<MeasuredJet> & selectedJets,
             const MeasuredLepton           & leadingLepton,
             const MeasuredLepton           & subLeadingLepton,
             const MeasuredLepton           & thirdLepton,
             const MeasuredHadronicTau      & selectedHadronicTau,
             const MeasuredMET              & measuredMET,
             UInt_t    run  = 0,
             UInt_t    lumi = 0,
             ULong64_t evt  = 0);

  /* MEM environment parameters (do not change unless you know what you're doing) */
  std::string pdfName;
  std::string madgraphFileName;
  std::string integrationMode;
  unsigned maxObjFunctionCalls;
  double higgsWidth;
  bool useBJetTransferFunction;
  bool useAvgBjetCombo;

  /* Parameters for Markov Chain integrator */
  std::string mxMode;
  unsigned nofBatches;
  unsigned nofChains;
  unsigned maxCallsStartingPos;
  double epsilon0;
  double T0;
  double nu;

private:
  MEM_ttHorZ_3l1tau mem_tt_HandZ;
  LikelihoodRatio_3l1tau lr_computation;
  VariableManager_3l1tau vm;
};

#endif // MEMINTERFACE_3L1TAU_H
