#ifndef INTEGRAND_TTHORZ_3L1TAU_LO_H
#define INTEGRAND_TTHORZ_3L1TAU_LO_H

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
    Integrand_ttHorZ_3l1tau(double sqrtS,
                            const std::string & pdfName,
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
     * @brief Calculates cos theta of the neutrino coming from the hadronic tau
     * @param nuHtau_en Energy of the neutrino coming from hadronic tau decay
     * @return Cosine of the opening angle
     *
     * @todo consider moving the function to a separate class, static it and
     *       bind it with missing arguments inside this class
     */
    double
    nuHtauCosTheta(double nuHtau_en) const;

    /**
     * @brief Calculates cos theta of the neutrino pair coming from the leptonic tau
     * @param nuLeptTau_en Energy of the neutrino pair coming from leptonic tau decay
     * @return Cosine of the opening angle
     *
     * @todo consider moving the function to a separate class, static it and
     *       bind it with missing arguments inside this class
     */
    double
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
    double
    bJetEnergy(const LorentzVector & W,
               unsigned bIdx) const;

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
    const double sqrtS_;
    const double s_;
    const double invSqrtS_;
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
    unsigned complLeptIdx_; ///< index of lepton coming from tau
    double complLeptEnergy_;
    double complLeptMomentum_;
    double complLeptMass_;
    double complLeptMassSquared_;
    LorentzVector complLeptP4_;
    Vector eX_lept_, eY_lept_, eZ_lept_;
    // for the leptons coming from t decay:
    unsigned lept1Idx_, lept2Idx_;
    double lept1Energy_;
    double lept2Energy_;
    LorentzVector lept1p4_,     lept2p4_;
    Vector        lept1p3_,     lept2p3_;
    Vector        lept1p3Unit_, lept2p3Unit_;
    // for the b quarks coming from t decay:
    double bJetRecoEnergy_[2];
    Vector bJetp3Unit_[2];

    double MET_x_, MET_y_;
    double covDet_;
    TMatrixDSym invCovMET_;
    double MET_TF_denom;

    int idxCosTheta1_;
    int idxVarphi1_;
    int idxCosTheta2_;
    int idxVarphi2_;
    int idxZ1_;
    int idxPhi1_;
    int idxPhiInv_;
    int idxMinvSquared_;

    double * mgGluon1p4_;
    double * mgGluon2p4_;
    double * mgBjet1p4_;
    double * mgLeptonFromBjet1p4_;
    double * mgNeutrinoFromBjet1p4_;
    double * mgBjet2p4_;
    double * mgLeptonFromBjet2p4_;
    double * mgNeutrinoFromBjet2p4_;
    double * mgTau1p4_;
    double * mgTau2p4_;
    mutable std::vector<double *> mgMomenta_;
  };
}

#endif // INTEGRAND_TTHORZ_3L1TAU_LO_H
