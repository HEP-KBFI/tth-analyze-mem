#ifndef INTEGRAND_TTHORZ_3L1TAU_H
#define INTEGRAND_TTHORZ_3L1TAU_H

#include <string> // std::string

#include "tthAnalysis/tthMEM/interface/me_ttHorZ_3l1tau_mg5.h" // me_ttHorZ_3l1tau_mg5
#include "tthAnalysis/tthMEM/interface/MeasuredEvent_3l1tau.h" // tthMEM_3l_1tau::MeasuredEvent

#include "LHAPDF/LHAPDF.h" // LHAPDF::PDF

namespace tthMEM
{
  /**
   * @brief MEM integrand of the process tth, h->tautau, 3l1tau channel (LO)
   *
   * @todo
   *   - add methods for setting transfer functions
   *   - implement the integrand
   */
  class Integrand_ttHorZ_3l1tau
  {
  public:
    /**
     * @brief Simple constructor
     * @param sqrtS            Center of momentum energy
     * @param pdfName          Name of parton distribution function (pdf)
     * @param madgraphFileName Full path to madgraph5 file (param_card.dat)
     */
    Integrand_ttHorZ_3l1tau(const std::string & pdfName,
                            const std::string & madgraphFileName);
    ~Integrand_ttHorZ_3l1tau();

    /**
     * @brief Set measured/reconstructed event as an input
     * @param measuredEvent The event
     */
    void
    setEvent(const MeasuredEvent_3l1tau & measuredEvent);

    /**
     * @brief Updates input parameters changed in the event
     */
    void
    renewInputs();

    /* simple setters */
    void setNumDimensions(unsigned numDimensions);
    void setIdxCosTheta1   (int idx);
    void setIdxVarphi1     (int idx);
    void setIdxCosTheta2   (int idx);
    void setIdxVarphi2     (int idx);
    void setIdxZ1          (int idx);
    void setIdxPhi1        (int idx);
    void setIdxPhiInv      (int idx);
    void setIdxMinvSquared (int idx);
    void setCurrentME(ME_mg5_3l1tau currentME);

    /**
     * @brief Evaluates the integrand
     * @param x Array of sampled values (generated by the integration lib)
     * @return Probability
     */
    double
    eval(const double * x) const;

    static const Integrand_ttHorZ_3l1tau * gIntegrand;
    ///< static pointer to this instance

  protected:
    const Vector beamAxis_;

    LHAPDF::PDF * pdf_;

    ME_mg5_3l1tau currentME_;
    mutable me_ttHorZ_3l1tau_mg5 * me_madgraph_[2];
    ///< @note mutable members can be modified by a const function (e.g. eval())
    ///< @note dynamic dispatch is tad bit slower because of additional vlookup

    unsigned numDimensions_;

    const MeasuredEvent_3l1tau * measuredEvent_;
    double measuredVisMassSquared_; ///< mass of visible tau decay products
    // for the hadronic tau:
    double hTauEnergy_;
    double hTauMomentum_;
    double hTauMass_;
    double hTauMassSquared_;
    LorentzVector hTauP4_;
    Vector eX_htau_, eY_htau_, eZ_htau_;
    // for the lepton coming from tau:
    double complLeptEnergy_;
    double complLeptMomentum_;
    double complLeptMass_;
    double complLeptMassSquared_;
    LorentzVector complLeptP4_;
    Vector eX_lept_, eY_lept_, eZ_lept_;
    // for the leptons coming from t decay:
    double leptEnergy_[2];
    LorentzVector leptP4_[2];
    Vector        leptP3_[2];
    Vector        leptP3Unit_[2];
    // for the b quarks coming from t decay:
    double bJetRecoEnergy_[2];
    double bJetP_[2];
    Vector bJetP3Unit_[2];

    double MET_x_, MET_y_;
    double covDet_;
    TMatrixDSym invCovMET_;
    double MET_TF_denom;
    double Q_;

    int idxCosTheta1_;
    int idxVarphi1_;
    int idxCosTheta2_;
    int idxVarphi2_;
    int idxZ1_;
    int idxPhi1_;
    int idxPhiInv_;
    int idxMinvSquared_;

    mutable std::vector<double *> mgMomenta_;
    std::vector<unsigned> mgMomentaIdxs_;

    /**
     * @brief Sets madgraph momenta
     * @param memVector_p4 std::vector<LorentzVector> containing 4-vectors
     */
    void
    setMGmomenta(const std::vector<LorentzVector> & memVector_p4) const;

    /**
     * @brief Calculates cos theta of the neutrino coming from the hadronic tau
     * @param nuHtau_en Energy of the neutrino coming from hadronic tau decay
     * @return Cosine of the opening angle
     *
     * @todo consider moving the function to a separate class, static it and
     *       bind it with missing arguments inside this class
     */
    inline double
    nuHtauCosTheta(double nuHtau_en) const;

    /**
     * @brief Calculates cos theta of the neutrino pair coming from the leptonic tau
     * @param nuLeptTau_en Energy of the neutrino pair coming from leptonic tau decay
     * @return Cosine of the opening angle
     *
     * @todo consider moving the function to a separate class, static it and
     *       bind it with missing arguments inside this class
     */
    inline double
    nuLeptTauCosTheta(double nuLeptTau_en,
                      double mInvSquared,
                      double nuLeptTau_p) const;

    /**
     * @brief Calculates true b quark energy
     * @param W 4-momentum of W boson
     * @param bIdx b-jet index
     * @return The b quark energy
     *
     * @todo consider moving the function to a separate class, static it and
     *       bind it with missing arguments inside this class
     */
    inline double
    bJetEnergy(const LorentzVector & W,
               unsigned bIdx) const;

    /**
     * @brief Calculates Jacobi factor arising from variable change in the delta-function
     *        associated with top decay
     * @param W      the W-boson 4-vector
     * @param b_en   energy of the b-quark
     * @param nuT_en energy of the neutrino
     * @param blIdx  index of the lepton and b-jet (valid values: 0, 1)
     * @return the Jacobi factor
     *
     * @todo consider moving the function to a separate class, static it and
     *       bind it with missing arguments inside this class
     */
    inline double
    tDecayJacobiFactor(const LorentzVector & W,
                       double b_en,
                       double nuT_en,
                       unsigned blIdx) const;

    /**
     * @brief Calculates effective matrix element squared for the process
     *            tau -> hadronic tau + tau neutrino
     *        which is derived such that the BR of the process is reproduced
     * @return
     */
    inline double
    MeffSquaredTau2hadrons() const;

    /**
     * @brief Calculates the Jacobi and phase space factor which arises from
     *        the variable changes in the delta function associated with
     *        hadronic tau decays
     * @param z Energy fraction carried by hadronic tau decay product (in (0, 1))
     * @return The phase space x jacobi factor
     */
    inline double
    hadTauPSJacobiFactor(const double z) const;

    inline double
    leptTauPSJacobiFactor(double mInvSquared,
                          double z) const;
  };
}

#endif // INTEGRAND_TTHORZ_3L1TAU_H
