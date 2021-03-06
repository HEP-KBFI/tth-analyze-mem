#ifndef INTEGRAND_TTHORZ_3L1TAU_H
#define INTEGRAND_TTHORZ_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/VariableManager_3l1tau.h" // VariableManager_3l1tau
#include "tthAnalysis/tthMEM/interface/RecoTrueEvent_ttHorZ_3l1tau.h" // RecoTrueEvent_ttHorZ_3l1tau

#include <functional> // std::function<>

namespace LHAPDF
{
  class PDF;
}

class mg5_tthz_t2lvl_tbar2lvl_hz2tata;

namespace tthMEM
{
  class MeasuredEvent_3l1tau;

  /**
   * @brief MEM integrand of the process tth, h->tautau, 3l1tau channel (LO)
   *
   * @todo
   *   - transfer function for tau pT
   */
  class Integrand_ttHorZ_3l1tau
  {
  public:
    /**
     * @brief Simple constructor
     * @param pdfName          Name of parton distribution function (pdf)
     * @param madgraphFileName Full path to madgraph5 file (param_card.dat)
     * @param vm               Variable manager (needed to pull new values
     *                         at every iteration)
     */
    Integrand_ttHorZ_3l1tau(const std::string & pdfName,
                            const std::string & madgraphFileName,
                            const VariableManager_3l1tau & vm);
    ~Integrand_ttHorZ_3l1tau();

    /**
     * @brief Set measured/reconstructed event as an input
     * @param measuredEvent The event
     * @return reference to this
     */
    Integrand_ttHorZ_3l1tau &
    setEvent(const MeasuredEvent_3l1tau & measuredEvent);

    /**
     * @brief Updates input parameters changed in the event
     */
    void
    renewInputs();

    /**
     * @brief Sets the mode in which the MEM integration is carried out
     * @param currentME The mode: either kTTH or kTTZ
     * @return Reference to the this instance
     */
    Integrand_ttHorZ_3l1tau &
    setCurrentME(ME_mg5_3l1tau currentME);

    /**
     * @brief Enables or disables the transfer function (TF) for b-quark energy
     * @param setTF If true, uses Gaussian-like TF; if false, uses the delta-function
     * @return Reference to the this instance
     */
    Integrand_ttHorZ_3l1tau &
    setBJetTransferFunction(bool setTF);

    /**
     * @brief Sets Higgs width for TTH matrix element
     * @param higgsWidth The Higgs width
     * @return Reference to this instance
     */
    void
    setHiggsWidth(double higgsWidth) const;

    /**
     * @brief Evaluates the integrand
     * @param x Array of sampled values (generated by the integration lib)
     * @return Probability
     */
    double
    eval(const double * x) const;

    static const Integrand_ttHorZ_3l1tau * gIntegrand;
    ///< static pointer to this instance

    mutable RecoTrueEvent_ttHorZ_3l1tau recoEvent;
    ///< reconstructed event in the lab frame

    const Vector beamAxis_; ///< defines the z-coordinate in the lab frame

  protected:
    LHAPDF::PDF * pdf_;
    ///< pointer to the parton distribution function (PDF)
    ME_mg5_3l1tau currentME_;
    ///< tells whether current integration is for tth or ttz
    mutable mg5_tthz_t2lvl_tbar2lvl_hz2tata * me_madgraph_[2];
    ///< two sets of MadGraph matrix elements: one for tth, the other for ttz

    double Q_;              ///< the resolution scale needed by the PDF

    const MeasuredEvent_3l1tau * measuredEvent_;
    ///< pointer to the measured event
    const VariableManager_3l1tau & vm_;
    ///< variable metadata used to pull new values

    mutable std::vector<double *> mgMomenta_;
    ///< needed by MadGraph to calculate the matrix elements
    std::vector<unsigned> mgMomentaIdxs_;
    ///< specifies the order in which the 4-momenta are assigned in MadGraph

    /* reconstruction functionals */
    std::function<double(double, double, double,
                         const LorentzVector &)>         z2_;
    std::function<double(double)>                        nuHtauCosTheta_;
    std::function<double(double)>                        nuHtauEnergy_;
    std::function<LorentzVector(double, double, double)> nuHTau_;
    std::function<LorentzVector(double, double, double)> nuLtau_;
    std::function<double(const VectorSpherical &)>       nuWEnergy_[2];
    std::function<double(const LorentzVector &)>         bQuarkEnergy_[2];

    /* transfer functionals */
    std::function<double(double, double, double)>        bJetTF_;
    std::function<double(double)>                        bJetTFBound_[2];
    std::function<double(double, double)>                MET_TF_;

    /* Jacobi x phase space factors */
    std::function<double(const LorentzVector &,
                         double, double, double)>        tDecayJacobiFactor_[2];
    std::function<double(double)>                        hadTauPSJacobiFactor_;
    std::function<double(double, double)>                leptTauPSJacobiFactor_;

    /**
     * @brief Sets madgraph momenta
     * @param memVector_p4 The std::vector<LorentzVector> containing 4-vectors
     */
    void
    setMGmomenta(const std::vector<LorentzVector> & memVector_p4) const;
  };
}

#endif // INTEGRAND_TTHORZ_3L1TAU_H
