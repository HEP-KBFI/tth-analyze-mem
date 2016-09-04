#include "tthAnalysis/tthMEM/interface/Integrand_ttHorZ_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/me_tth_3l1tau_mg5.h"
#include "tthAnalysis/tthMEM/interface/me_ttz_3l1tau_mg5.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/tthMEMrecFunctions.h"
#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h"
#include "tthAnalysis/tthMEM/interface/BJetTransferFunction.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <cmath> // std::sqrt()
#include <cstring> // std::memset()
#include <algorithm> // std::for_each(), std::copy()
#include <sstream> // std::ostringstream
#include <iomanip> // std::setprecision()
#include <exception> // std::invalid_argument, std::domain_error

#include <TMath.h> // TMath::IsNaN() ...
 // ... (why not use std here: http://stackoverflow.com/a/570694)

using namespace tthMEM;

const Integrand_ttHorZ_3l1tau * Integrand_ttHorZ_3l1tau::gIntegrand = 0;

Integrand_ttHorZ_3l1tau::Integrand_ttHorZ_3l1tau(const std::string & pdfName,
                                                 const std::string & madgraphFilename,
                                                 const VariableManager_3l1tau & vm)
  : pdf_(0)
  , currentME_(ME_mg5_3l1tau::kTTH) // default to tth
  , me_madgraph_{0, 0}
  , beamAxis_(0., 0., 1.)
  , measuredEvent_(0)
  , vm_(vm)
  , bJetTF_(functions::deltaFunction)
{
  LOGTRC;

  if(! pdf_ && pdfName != "") pdf_ = LHAPDF::mkPDF(pdfName.c_str(), 0);
  else
  {
    throw std::invalid_argument("PDF file name empty!");
  }

  me_madgraph_[ME_mg5_3l1tau::kTTH] = new me_tth_3l1tau_mg5();
  me_madgraph_[ME_mg5_3l1tau::kTTZ] = new me_ttz_3l1tau_mg5();

  if(madgraphFilename != "")
  {
    for(unsigned i = 0; i < 2; ++i)
      me_madgraph_[i] -> initProc(madgraphFilename);
  }
  else
  {
    throw std::invalid_argument("Madgraph file name empty!");
  }

  for(unsigned i = 0; i < 10; ++i)
    mgMomenta_.push_back(new double[4]);
  std::for_each(mgMomenta_.begin(), mgMomenta_.end(),
    [](double * & d) { std::memset(d, 0., 4 * sizeof(double)); }
  );

  gIntegrand = this;
}

Integrand_ttHorZ_3l1tau::~Integrand_ttHorZ_3l1tau()
{
  LOGTRC;

  if(pdf_) delete pdf_;

  for(unsigned i = 0; i < 2; ++i)
    if(me_madgraph_[i]) delete me_madgraph_[i];

  measuredEvent_ = 0; // no allocation, just the address

  std::for_each(mgMomenta_.begin(), mgMomenta_.end(),
    [](double * & d) { delete d; d = 0; }
  );
}

Integrand_ttHorZ_3l1tau &
Integrand_ttHorZ_3l1tau::setCurrentME(ME_mg5_3l1tau currentME)
{
  LOGTRC;
  currentME_ = currentME;
  Q_ = currentME_ == ME_mg5_3l1tau::kTTH ?
                     constants::resolutionScaleTTH :
                     constants::resolutionScaleTTZ;
  return *this;
}

Integrand_ttHorZ_3l1tau &
Integrand_ttHorZ_3l1tau::setBJetTransferFunction(bool setTF)
{
  LOGTRC;
  bJetTF_ = setTF ? functions::bJetTF : functions::deltaFunction;
  return *this;
}

Integrand_ttHorZ_3l1tau &
Integrand_ttHorZ_3l1tau::setEvent(const MeasuredEvent_3l1tau & measuredEvent)
{
  LOGTRC;
  measuredEvent_ = &measuredEvent;
  mgMomentaIdxs_.clear();

//--- set the variables related to the hadronic tau
  const MeasuredHadronicTau & htau = measuredEvent_ -> htau1;
  recoEvent.hTauLepton = htau.p4();
  const Vector eZ_htau = htau.p3().unit();
  const Vector eY_htau = eZ_htau.Cross(beamAxis_).unit();
  const Vector eX_htau = eY_htau.Cross(eZ_htau).unit();
  // eX should already be unit vector by construction
  LOGVRB << lvrap("htau p4", recoEvent.hTauLepton);
  LOGVRB << svrap("htau eX", eX_htau);
  LOGVRB << svrap("htau eY", eY_htau);
  LOGVRB << svrap("htau eZ", eZ_htau);
  LOGVRB << "htau: " << "eX x eY = " << eX_htau.Cross(eY_htau).r() << " ; "
                     << "eX x eZ = " << eX_htau.Cross(eZ_htau).r() << " ; "
                     << "eY x eZ = " << eY_htau.Cross(eZ_htau).r();
  const Vector eX_htau_(eX_htau.x(), eY_htau.x(), eZ_htau.x());
  const Vector eY_htau_(eX_htau.y(), eY_htau.y(), eZ_htau.y());
  const Vector eZ_htau_(eX_htau.z(), eY_htau.z(), eZ_htau.z());

//--- bind the functional arguments, so that no explicit storage is needed
//--- for the variables which remain constant during the integration
  const double hTauEnergy = htau.energy();
  const double hTauMassSquared = pow2(htau.mass());
  const double hTauMomentum = htau.p();
  const double hTauInvBeta = 1. / htau.p4().Beta();
  nuHtauCosTheta_ = std::bind(functions::nuHtauCosTheta,
                              std::placeholders::_1,
                              hTauEnergy,
                              hTauMassSquared,
                              hTauMomentum);
  hadTauPSJacobiFactor_ = std::bind(functions::hadTauPSJacobiFactor,
                                    std::placeholders::_1,
                                    hTauMassSquared,
                                    hTauInvBeta);
  nuHtauEnergy_ = std::bind(functions::nuTauEnergy,
                            std::placeholders::_1,
                            hTauEnergy);
  nuHTau_ = std::bind(functions::nuP4,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3,
                      std::placeholders::_3,
                      eX_htau_,
                      eY_htau_,
                      eZ_htau_);

//--- set the variables related to the MET/hadronic recoil TF
  const MeasuredMET & met = measuredEvent_ -> met;
  const double MET_x = met.px();
  const double MET_y = met.py();
  TMatrixDSym invCovMET = met.covMET();
  const double covDet = invCovMET.Determinant();
  if(covDet != 0.)
  {
    invCovMET.Invert();
    const double MET_TF_denom = 1. / (2. * pi() * std::sqrt(invCovMET.Determinant()));
    LOGVRB << "MET TF denominator = " << MET_TF_denom;
    MET_TF_ = std::bind(functions::MET_TF,
                        std::placeholders::_1,
                        std::placeholders::_2,
                        MET_x,
                        MET_y,
                        MET_TF_denom,
                        invCovMET);
  }
  else
  {
    MET_TF_ = 0;
    LOGERR << "Cannot invert MET covariance matrix b/c det = 0";
  }

  return *this;
}

void
Integrand_ttHorZ_3l1tau::renewInputs()
{
  LOGTRC;
//--- set the variables related to the leptonic tau
  const unsigned complLeptIdx = measuredEvent_ -> complLeptonIdx;
  const MeasuredLepton & complLepton = measuredEvent_-> leptons[complLeptIdx];

  recoEvent.lTauLepton = complLepton.p4();
  const int complLeptCharge = complLepton.charge();
  const Vector eZ_lept = complLepton.p3().unit();
  const Vector eY_lept = eZ_lept.Cross(beamAxis_).unit();
  const Vector eX_lept = eY_lept.Cross(eZ_lept).unit();
  // eX should already be unit vector by construction
  LOGVRB << lvrap("lept p4", recoEvent.lTauLepton);
  LOGVRB << cvrap("lept p3", complLepton.p3());
  LOGVRB << svrap("lept eX", eX_lept);
  LOGVRB << svrap("lept eY", eY_lept);
  LOGVRB << svrap("lept eZ", eZ_lept);
  LOGVRB << "lept: " << "eX x eY = " << eX_lept.Cross(eY_lept).R() << " ; "
                     << "eX x eZ = " << eX_lept.Cross(eZ_lept).R() << " ; "
                     << "eY x eZ = " << eY_lept.Cross(eZ_lept).R();
  const Vector eX_lept_(eX_lept.x(), eY_lept.x(), eZ_lept.x());
  const Vector eY_lept_(eX_lept.y(), eY_lept.y(), eZ_lept.y());
  const Vector eZ_lept_(eX_lept.z(), eY_lept.z(), eZ_lept.z());

//--- find the measured mass of both visible tau decay products
  const double measuredVisMassSquared =
    (recoEvent.hTauLepton + recoEvent.lTauLepton).mass2();
  const double massHiggsOrZsquared = currentME_ == ME_mg5_3l1tau::kTTH ?
    constants::massHiggsSquared : constants::massZSquared;
  LOGVRB << "measured mass of visible tau decay products: "
         << std::sqrt(measuredVisMassSquared);

//--- find the kinematic variables of the rest of the leptons (from t decay)
//--- and bind the functional arguments, so that no explicit storage is needed
//--- for the variables which remain constant during the integration
  const double complLeptMassSquared = pow2(complLepton.mass());
  const double complLeptEnergy = complLepton.energy();
  const double complLeptP = complLepton.p();
  nuLeptTauCosTheta_ = std::bind(functions::nuLeptTauCosTheta,
                                 std::placeholders::_1,
                                 std::placeholders::_2,
                                 std::placeholders::_3,
                                 complLeptEnergy,
                                 complLeptMassSquared,
                                 complLeptP);
  leptTauPSJacobiFactor_ = std::bind(functions::leptTauPSJacobiFactor,
                                     std::placeholders::_1,
                                     std::placeholders::_2,
                                     complLeptMassSquared,
                                     complLeptP);
  nuLTauEnergy_ = std::bind(functions::nuTauEnergy,
                            std::placeholders::_1,
                            complLeptEnergy);
  nuLTau_ = std::bind(functions::nuP4,
                      std::placeholders::_1,
                      std::placeholders::_2,
                      std::placeholders::_3,
                      std::placeholders::_4,
                      eX_lept_,
                      eY_lept_,
                      eZ_lept_);
  z2_ = std::bind(functions::z2,
                  std::placeholders::_1,
                  measuredVisMassSquared,
                  massHiggsOrZsquared);
  int lept1Charge = 0;
  for(unsigned i = 0; i < 2; ++i)
  {
    const MeasuredJet & bJet_i = measuredEvent_ -> jets[i];
    const double bJetEnergy_i = bJet_i.energy();
    const double bJetEta_i = bJet_i.eta();
    const Vector bJetP3unit_i = bJet_i.p3().unit();

    const unsigned leptIdx_i = measuredEvent_ -> bjetLeptonIdxs[i];
    const MeasuredLepton & lept_i = measuredEvent_ -> leptons[leptIdx_i];
//--- set the lepton energy equal to its momentum, thus setting it massless
    const double leptEnergy_i = lept_i.p();
    const Vector leptP3_i = lept_i.p3();
    const Vector leptP3unit_i = leptP3_i.unit();
    recoEvent.lW[i] = getLorentzVector(leptP3_i, leptEnergy_i);
    LOGVRB << lvrap("W lept " + std::to_string(i + 1) + " p4", recoEvent.lW[i]);
    if(i == 0) lept1Charge = lept_i.charge();

    bQuarkEnergy_[i] = std::bind(functions::bQuarkEnergy,
                                 std::placeholders::_1,
                                 bJetP3unit_i,
                                 bJetEnergy_i);
    bJetTFBound_[i] = std::bind(functions::bJetTF,
                                std::placeholders::_1,
                                bJetEnergy_i,
                                bJetEta_i);
    tDecayJacobiFactor_[i] = std::bind(functions::tDecayJacobiFactor,
                                       std::placeholders::_1,
                                       std::placeholders::_2,
                                       std::placeholders::_3,
                                       std::placeholders::_4,
                                       leptEnergy_i,
                                       bJetP3unit_i);
    nuTopEnergy_[i] = std::bind(functions::nuTopEnergy,
                                std::placeholders::_1,
                                leptP3unit_i,
                                leptEnergy_i);
  }

//--- MadGraph momenta legend:
// 0, g
// 1, g
// 2, b, associated W+ => associated lepton +
// 3, e+, associated W+ => associated lepton +
// 4, ve, associated W+ => associated lepton +
// 5, b~, associated W- => associated lepton -
// 6, e-, associated W- => associated lepton -
// 7, ve~, associated W- => associated lepton -
// 8, ta+ => associated lepton +
// 9, ta- => associated lepton -
  mgMomentaIdxs_.clear();
  if(lept1Charge == +1 && complLeptCharge == +1)
    mgMomentaIdxs_ = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  else
  if(lept1Charge == +1 && complLeptCharge == -1)
    mgMomentaIdxs_ = { 0, 1, 2, 3, 4, 5, 6, 7, 9, 8 };
  else
  if(lept1Charge == -1 && complLeptCharge == +1)
    mgMomentaIdxs_ = { 0, 1, 5, 6, 7, 2, 3, 4, 8, 9 };
  else
  if(lept1Charge == -1 && complLeptCharge == -1)
    mgMomentaIdxs_ = { 0, 1, 5, 6, 7, 2, 3, 4, 9, 8 };

//--- for debugging purposes plot some variables
  if(DebugPlotter_ttHorZ_3l1tau * dPlotter = measuredEvent_ -> debugPlotter)
  {
    dPlotter -> write();
    std::string measuredEventStr = measuredEvent_ -> str();
    measuredEventStr += std::string("_") +
      (currentME_ == ME_mg5_3l1tau::kTTH ? "tth" : "ttz");
    dPlotter -> initialize(measuredEventStr);
  }
}

void
Integrand_ttHorZ_3l1tau::setMGmomenta(const std::vector<LorentzVector> & memVector_p4) const
{
  for(unsigned i = 0; i < mgMomentaIdxs_.size(); ++i)
    setMGmomentum(memVector_p4[i], mgMomenta_[mgMomentaIdxs_[i]]);
}

double
Integrand_ttHorZ_3l1tau::eval(const double * x) const
{
  LOGTRC;
  if(! pdf_)                      LOGERR << "PDF not initialized!";
  if(! me_madgraph_[currentME_])  LOGERR << "Madgraph's ME not initialized!";
  if(! measuredEvent_)            LOGERR << "Measured event not specified!";
  if(mgMomentaIdxs_.size() != 10) LOGERR << "Number of MG momenta indexes not equal to 10";
  if(! pdf_ || ! me_madgraph_[currentME_] || ! measuredEvent_ || mgMomentaIdxs_.size() != 10)
    throw std::domain_error("Insufficient data for eval()");

  LOGVRB << "Current MG5 ME: " << me_madgraph_[currentME_] -> name();
  if(! MET_TF_) return 0.;
  const bool isTTH = currentME_ == ME_mg5_3l1tau::kTTH;

//--- read the sampled values
  LOGVRB << "x = { " << vm_.getArrayString(x) << " }";

  const double cosTheta[2]   = { vm_.get(Var_3l1tau::kBcosTheta1, x),
                                 vm_.get(Var_3l1tau::kBcosTheta2, x) };
  const double varphi[2]     = { vm_.get(Var_3l1tau::kBphi1, x),
                                 vm_.get(Var_3l1tau::kBphi1, x) };
  const double z1            = vm_.get(Var_3l1tau::kZ1, x);
  const double nuHtau_phi    = vm_.get(Var_3l1tau::kTauPhi, x);
  const double nuLTau_phi    = vm_.get(Var_3l1tau::kTauPhiInv, x);
  const double mInvSquared   = vm_.get(Var_3l1tau::kTauMinvSquared, x);

//--- confirm that the energy fraction carried by the tau is indeed in (0,1)
  const double z2 = z2_(z1);
  if(! (z2 >= 1.e-5 && z2 <= 1.))
  {
    LOGVRB << "z2 = " << z2 << " not in (0, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "z2 = " << z2;

//--- compute the neutrino and tau lepton 4-vector from hadronic tau
  const double nuHtau_en = nuHtauEnergy_(z1);
  const double nuHtau_cosTheta = nuHtauCosTheta_(nuHtau_en);
  if(! (nuHtau_cosTheta >= -1. && nuHtau_cosTheta <= +1.))
  {
    LOGVRB << "nuHtau_cosTheta = " << nuHtau_cosTheta << " not in (-1, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "nuHtau_cosTheta = " << nuHtau_cosTheta;
  recoEvent.nuHtau = nuHTau_(std::acos(nuHtau_cosTheta), nuHtau_phi, nuHtau_en);
  recoEvent.hTau = recoEvent.hTauLepton + recoEvent.nuHtau;
  LOGTRC << lvrap("htau nu",      recoEvent.nuHtau);
  LOGTRC << lvrap("hadronic tau", recoEvent.hTau);

//--- compute the neutrino and tau lepton 4-vector from leptonic tau
  const double nuLTau_en = nuLTauEnergy_(z2);
  const double nuLTau_p = std::sqrt(std::max(0., pow2(nuLTau_en) - mInvSquared));
  const double nuLTau_cosTheta = nuLeptTauCosTheta_(nuLTau_en, mInvSquared, nuLTau_p);
  if(! (nuLTau_cosTheta >= -1. && nuLTau_cosTheta <= +1.))
  {
    LOGVRB << "nuLTau_cosTheta = " << nuLTau_cosTheta << " not in (-1, 1) "
           << "=> p = 0";
    return 0.;
  }
  recoEvent.nuLtau = nuLTau_(std::acos(nuLTau_cosTheta), nuLTau_phi, nuLTau_en, nuLTau_p);
  recoEvent.lTau = recoEvent.lTauLepton + recoEvent.nuLtau;
  LOGTRC << lvrap("lept tau nu",  recoEvent.nuLtau);
  LOGTRC << lvrap("leptonic tau", recoEvent.lTau);

  recoEvent.higgsOrZ = recoEvent.hTau + recoEvent.lTau;
  if(isTTH) LOGTRC << lvrap("higgs", recoEvent.higgsOrZ);
  else      LOGTRC << lvrap("Z",     recoEvent.higgsOrZ);

//--- get 4-momenta of top decay products (and top itself)
  double bEnergyTF[2];  // value of b-quark energy transfer function
  for(unsigned i = 0; i < 2; ++i)
  {
    const std::string i_str = std::to_string(i + 1);
    const VectorSpherical nuTP3unit_i(1., std::acos(cosTheta[i]), varphi[i]);
    const double nuTenergy_i = nuTopEnergy_[i](nuTP3unit_i);
    recoEvent.nuW[i] = getLorentzVector(Vector(nuTenergy_i * nuTP3unit_i), nuTenergy_i);
    recoEvent.W[i] = recoEvent.nuW[i] + recoEvent.lW[i];

    LOGTRC << lvrap("W nu " + i_str, recoEvent.nuW[i]);
    LOGTRC << lvrap("W " + i_str,    recoEvent.W[i]);

    const double bEnergy_i = bQuarkEnergy_[i](recoEvent.W[i]);
    if(bEnergy_i == 0.) return 0.;
    const MeasuredJet & bJet_i = measuredEvent_ -> jets[i];
    const Vector bP3_i = std::sqrt(pow2(bEnergy_i) - constants::massBSquared) *
                         bJet_i.p3().unit();
    recoEvent.b[i] = getLorentzVector(bP3_i, bEnergy_i);
    recoEvent.t[i] = recoEvent.b[i] + recoEvent.W[i];

    LOGTRC << lvrap("b " + i_str,           recoEvent.b[i]);
    LOGTRC << lvrap("b " + i_str + " reco", bJet_i.p4());
    LOGTRC << lvrap("t " + i_str,           recoEvent.t[i]);

//--- b-jet energy transfer function
    bEnergyTF[i] = bJetTFBound_[i](bEnergy_i);
  }

  for(unsigned i = 0; i < 2; ++i)
    LOGTRC_S << "TF for #" << i << " b-quark energy = " << bEnergyTF[i];

//--- hadronic recoil transfer function; simplification: use only neutrinos
//--- in the difference of ,,measured'' and ,,true'' hadronic recoil, because
//--- leptons are assumed to be measured perfectly whereas the pT of jets/quarks
//--- introduces uncertainties difficult to handle here
  const LorentzVector nuSum = recoEvent.getNuSum();
  const double MET_TF = MET_TF_(nuSum.x(), nuSum.y());

//--- compute Bjorken x variables
//--- assume that hadronic recoil has only transverse component
  const double hadRecE = 0.;
  const double hadRecPz = 0.;
  const LorentzVector tthOrZ = recoEvent.getTTHorZ();
  const double xa = (hadRecE + tthOrZ.e() + hadRecPz + tthOrZ.pz()) *
                    constants::invSqrtS;
  const double xb = (hadRecE + tthOrZ.e() - hadRecPz - tthOrZ.pz()) *
                    constants::invSqrtS;
  LOGTRC << "xa = " << xa << "; xb = " << xb;
  if(xa <= 0. || xa >= 1. || xb <= 0. || xb >= 1.)
  {
    LOGVRB << "xa or xb have unphysical values";
    return 0.;
  }

  const double fa = pdf_ -> xfxQ(21, xa, Q_) / xa;
  const double fb = pdf_ -> xfxQ(21, xb, Q_) / xb;
  const double probPDF = fa * fb;
  const double flux = constants::invS / (xa * xb);

//--- boost all MG momenta into frame where pT(tth) = 0
//--- note: the boost vector must be a 3-vector of velocities, hence
//---       the division by energy
  const Vector boost(-tthOrZ.px() / tthOrZ.e(), -tthOrZ.py() / tthOrZ.e(), 0.);
  LOGTRC << cvrap("boost vector", boost);
  recoEvent.g[0] = LorentzVector(0., 0., +0.5 * xa * constants::sqrtS,
                                           0.5 * xa * constants::sqrtS);
  recoEvent.g[1] = LorentzVector(0., 0., -0.5 * xb * constants::sqrtS,
                                           0.5 * xb * constants::sqrtS);
  const RecoTrueEvent_ttHorZ_3l1tau recoEventmem = recoEvent.boost(boost);
  const LorentzVector tthOrZ_mem = recoEventmem.getTTHorZ();
  if(isTTH) LOGTRC << lvrap("tth mem", tthOrZ_mem);
  else      LOGTRC << lvrap("ttz mem", tthOrZ_mem);

  const double z1_mem = recoEventmem.hTauLepton.e() / recoEventmem.hTau.e();
  const double z2_mem = recoEventmem.lTauLepton.e() / recoEventmem.lTau.e();
  LOGTRC << "z1_mem = " << z1_mem << "; z2_mem = " << z2_mem;
  if(! (z1_mem >= 1.e-5 && z1_mem <= 1.) ||
     ! (z2_mem >= 1.e-5 && z2_mem <= 1.))
  {
    LOGVRB << "The tau energy fraction z1 and z2 have unphysical values "
           << "when the lab frame is boosted such that pT(tth) = 0";
    return 0.;
  }

//--- set MG momenta & calculate the matrix element
  setMGmomenta(recoEventmem.getForMG());
  me_madgraph_[currentME_] -> setMomenta(mgMomenta_);
  me_madgraph_[currentME_] -> sigmaKin();
  const double prob_ME_mg = me_madgraph_[currentME_] -> getMatrixElements()[0];
  if(TMath::IsNaN(prob_ME_mg) || prob_ME_mg < 0.)
  {
    LOGWARN_S << "MadGraph5 returned NaN or is zero: "
              << "|M|^2 = " << prob_ME_mg << " => skipping event";
    return 0.;
  }
  LOGTRC_S << "prob(PDF) = " << probPDF << "; flux factor = " << flux;
  LOGVRB_S << "|M|^2 = " << prob_ME_mg;

//--- jacobi factors
  double tDecayJacobiFactors[2];
  for(unsigned i = 0; i < 2; ++i)
  {
    tDecayJacobiFactors[i] = tDecayJacobiFactor_[i](
      recoEvent.W[i], recoEvent.b[i].e(), recoEvent.b[i].P(), recoEvent.nuW[i].e()
    );
    LOGTRC_S << "Jacobi factors arising from the top decay: #" << (i + 1) << " = "
             << tDecayJacobiFactors[i];
    if(tDecayJacobiFactors[i] == 0.) return 0.;
  }

  const double hTauPSJacobiFactor = hadTauPSJacobiFactor_(z1);
  LOGTRC_S << "PS x Jacobi factor for hadronic tau decay = "
           << hTauPSJacobiFactor;
  const double lTauPSJacobiFactor = leptTauPSJacobiFactor_(mInvSquared, z2);
  if(lTauPSJacobiFactor == 0.) return 0.;
  LOGTRC_S << "PS x Jacobi factor for leptonic tau decay = "
           << lTauPSJacobiFactor;

  const double jacobiFactor = tDecayJacobiFactors[0] * tDecayJacobiFactors[1] *
    hTauPSJacobiFactor * lTauPSJacobiFactor * z2 *
    (isTTH ? constants::ttHfactor : constants::ttZfactor);
  LOGTRC_S << "Product of all Jacobi factors = " << jacobiFactor;

//--- assemble the integrand
  const double p = prob_ME_mg * probPDF * flux * MET_TF *
                   bEnergyTF[0] * bEnergyTF[1] * jacobiFactor;
  LOGVRB_S << "p = " << p;

//--- for debugging purposes plot some variables
  if(DebugPlotter_ttHorZ_3l1tau * dPlotter = measuredEvent_ -> debugPlotter)
    (*dPlotter).fill(hVar_3l1tau::kZ2,       z2)
               .fill(hVar_3l1tau::kMassHorZ, recoEvent.higgsOrZ.mass())
               .fill(hVar_3l1tau::kMassHtau, recoEvent.hTau.mass())
               .fill(hVar_3l1tau::kMassLtau, recoEvent.lTau.mass())
               .fill(hVar_3l1tau::kB1en,     recoEvent.b[0].e())
               .fill(hVar_3l1tau::kB2en,     recoEvent.b[1].e())
               .fill(hVar_3l1tau::kB1RecoEn, measuredEvent_ -> jets[0].energy())
               .fill(hVar_3l1tau::kB2RecoEn, measuredEvent_ -> jets[1].energy())
               .fill(hVar_3l1tau::kMsquared, prob_ME_mg)
               .fill(hVar_3l1tau::kProb,     p);

  return p;
}
