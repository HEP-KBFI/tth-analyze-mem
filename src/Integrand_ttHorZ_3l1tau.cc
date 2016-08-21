#include "tthAnalysis/tthMEM/interface/Integrand_ttHorZ_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/me_tth_3l1tau_mg5.h"
#include "tthAnalysis/tthMEM/interface/me_ttz_3l1tau_mg5.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <cmath> // std::pow(), std::sqrt()
#include <cstdlib> // std::exit(), EXIT_FAILURE
#include <cstring> // std::memset()
#include <algorithm> // std::for_each(), std::copy()
#include <sstream> // std::ostringstream
#include <iomanip> // std::setprecision()

#include "TMath.h" // TMath::IsNaN() ...
 // ... (why not use std here: http://stackoverflow.com/a/570694)
#include "Math/VectorUtil.h" // ROOT::Math::VectorUtil::boost()

using namespace tthMEM;
namespace VectorUtil = ROOT::Math::VectorUtil;

const Integrand_ttHorZ_3l1tau * Integrand_ttHorZ_3l1tau::gIntegrand = 0;

Integrand_ttHorZ_3l1tau::Integrand_ttHorZ_3l1tau(double sqrtS,
                                                 const std::string & pdfName,
                                                 const std::string & madgraphFilename)
  : sqrtS_(sqrtS)
  , s_(std::pow(sqrtS_, 2))
  , invSqrtS_(1. / sqrtS_)
  , beamAxis_(0., 0., 1.)
  , pdf_(0)
  , currentME_(ME_mg5_3l1tau::kTTH) // default to tth
  , me_madgraph_{0, 0}
  , numDimensions_(0)
  , measuredEvent_(0)
  , invCovMET_(TMatrixDSym(2))
  , idxCosTheta1_(-1)
  , idxVarphi1_(-1)
  , idxCosTheta2_(-1)
  , idxVarphi2_(-1)
  , idxZ1_(-1)
  , idxPhi1_(-1)
  , idxPhiInv_(-1)
  , idxMinvSquared_(-1)
{
  LOGTRC;

  if(! pdf_ && pdfName != "") pdf_ = LHAPDF::mkPDF(pdfName.c_str(), 0);
  else
  {
    LOGERR << "PDF file name empty!";
    std::exit(EXIT_FAILURE);
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
    LOGERR << "Madgraph file name empty!";
    std::exit(EXIT_FAILURE);
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

void
Integrand_ttHorZ_3l1tau::setNumDimensions(unsigned numDimensions)
{
  numDimensions_ = numDimensions;
}

void
Integrand_ttHorZ_3l1tau::setIdxCosTheta1(int idx)
{
  idxCosTheta1_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setIdxVarphi1(int idx)
{
  idxVarphi1_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setIdxCosTheta2(int idx)
{
  idxCosTheta2_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setIdxVarphi2(int idx)
{
  idxVarphi2_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setIdxZ1(int idx)
{
  idxZ1_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setIdxPhi1(int idx)
{
  idxPhi1_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setIdxPhiInv(int idx)
{
  idxPhiInv_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setIdxMinvSquared(int idx)
{
  idxMinvSquared_ = idx;
}

void
Integrand_ttHorZ_3l1tau::setCurrentME(ME_mg5_3l1tau currentME)
{
  LOGTRC;
  currentME_ = currentME;
  Q_ = currentME_ == ME_mg5_3l1tau::kTTH ? resolutionScaleTTH :
                                           resolutionScaleTTZ;
}

void
Integrand_ttHorZ_3l1tau::setEvent(const MeasuredEvent_3l1tau & measuredEvent)
{
  LOGTRC;
  measuredEvent_ = &measuredEvent;
  mgMomentaIdxs_.clear();

//--- set the variables related to the hadronic tau
  const MeasuredHadronicTau & htau = measuredEvent_ -> htau1;
  hTauEnergy_ = htau.energy();
  hTauMomentum_ = htau.p();
  hTauMass_ = htau.mass();
  hTauMassSquared_ = std::pow(hTauMass_, 2);
  hTauP4_ = htau.p4();
  eZ_htau_ = htau.p3().unit();
  eY_htau_ = eZ_htau_.Cross(beamAxis_).unit();
  eX_htau_ = eY_htau_.Cross(eZ_htau_).unit();
  // eX should already be unit vector by construction
  LOGVRB << lvrap("htau p4", hTauP4_);
  LOGVRB << svrap("htau eX", eX_htau_);
  LOGVRB << svrap("htau eY", eY_htau_);
  LOGVRB << svrap("htau eZ", eZ_htau_);
  LOGVRB << "htau: " << "eX x eY = " << eX_htau_.Cross(eY_htau_).r() << " ; "
                     << "eX x eZ = " << eX_htau_.Cross(eZ_htau_).r() << " ; "
                     << "eY x eZ = " << eY_htau_.Cross(eZ_htau_).r();

//--- set the variables related to the MET/hadronic recoil TF
  const MeasuredMET & met = measuredEvent_ -> met;
  MET_x_ = met.px();
  MET_y_ = met.py();
  invCovMET_ = met.covMET();
  covDet_ = invCovMET_.Determinant();
  if(covDet_ != 0.)
  {
    invCovMET_.Invert();
    MET_TF_denom = 1. / (2. * pi() * std::sqrt(invCovMET_.Determinant()));
    LOGVRB << "MET TF denominator = " << MET_TF_denom;
  }
  else
    LOGERR << "Cannot invert MET covariance matrix b/c det = 0";
}

void
Integrand_ttHorZ_3l1tau::renewInputs()
{
  LOGTRC;
//--- set the variables related to the leptonic tau
  const unsigned complLeptIdx = measuredEvent_ -> complLeptonIdx;
  const MeasuredLepton & complLepton = measuredEvent_-> leptons[complLeptIdx];
  complLeptEnergy_ = complLepton.energy();
  complLeptMomentum_ = complLepton.p();
  complLeptMass_ = complLepton.mass();
  complLeptMassSquared_ = std::pow(complLeptMass_, 2);
  complLeptP4_ = complLepton.p4();
  const int complLeptCharge = complLepton.charge();
  eZ_lept_ = complLepton.p3().unit();
  eY_lept_ = eZ_lept_.Cross(beamAxis_).unit();
  eX_lept_ = eY_lept_.Cross(eZ_lept_).unit();
  // eX should already be unit vector by construction
  LOGVRB << lvrap("lept p4", complLeptP4_);
  LOGVRB << cvrap("lept p3", complLepton.p3());
  LOGVRB << svrap("lept eX", eX_lept_);
  LOGVRB << svrap("lept eY", eY_lept_);
  LOGVRB << svrap("lept eZ", eZ_lept_);
  LOGVRB << "lept: " << "eX x eY = " << eX_lept_.Cross(eY_lept_).R() << " ; "
                     << "eX x eZ = " << eX_lept_.Cross(eZ_lept_).R() << " ; "
                     << "eY x eZ = " << eY_lept_.Cross(eZ_lept_).R();

//--- find the measured mass of both visible tau decay products
  measuredVisMassSquared_ = (hTauP4_ + complLeptP4_).mass2();
  LOGVRB << "measured mass of visible tau decay products: "
         << std::sqrt(measuredVisMassSquared_);

//--- find the kinematic variables of the rest of the leptons (from t decay)
  const unsigned lept1Idx = measuredEvent_ -> bjetLeptonIdxs[0];
  const unsigned lept2Idx = measuredEvent_ -> bjetLeptonIdxs[1];
  const MeasuredLepton & lept1 = measuredEvent_ -> leptons[lept1Idx];
  const MeasuredLepton & lept2 = measuredEvent_ -> leptons[lept2Idx];
  lept1p3_ = lept1.p3();
  lept2p3_ = lept2.p3();
  lept1Energy_ = lept1.p();
  lept2Energy_ = lept2.p();
  lept1p3Unit_ = lept1p3_.unit();
  lept2p3Unit_ = lept2p3_.unit();
  lept1p4_ = LorentzVector(lept1p3_.x(), lept1p3_.y(), lept1p3_.z(), lept1Energy_);
  lept2p4_ = LorentzVector(lept2p3_.x(), lept2p3_.y(), lept2p3_.z(), lept2Energy_);
  const int lept1Charge = lept1.charge();
  LOGVRB << lvrap("t lept 1 p4", lept1p4_);
  LOGVRB << lvrap("t lept 2 p4", lept2p4_);

//--- find the kinematic variables of the b quarks
  for(unsigned i = 0; i < 2; ++i)
  {
    const MeasuredJet & bjet_i = measuredEvent_ -> jets[i];
    bJetRecoEnergy_[i] = bjet_i.energy();
    bJetp3Unit_[i] = bjet_i.p3().unit();
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
}

void
Integrand_ttHorZ_3l1tau::setMGmomenta(const std::vector<LorentzVector> & memVector_p4) const
{
  for(unsigned i = 0; i < mgMomentaIdxs_.size(); ++i)
    setMGmomentum(memVector_p4[i], mgMomenta_[mgMomentaIdxs_[i]]);
}

double
Integrand_ttHorZ_3l1tau::nuHtauCosTheta(double nuHtau_en) const
{
  return (nuHtau_en * hTauEnergy_ -
             (tauLeptonMassSquared - hTauMassSquared_) / 2.)
           / (hTauMomentum_ * nuHtau_en);
}

double
Integrand_ttHorZ_3l1tau::nuLeptTauCosTheta(double nuLeptTau_en,
                                           double mInvSquared,
                                           double nuLeptTau_p) const
{
  return (nuLeptTau_en * complLeptEnergy_ -
             (tauLeptonMassSquared - (complLeptMassSquared_ + mInvSquared)) / 2.)
           / (complLeptMomentum_ * nuLeptTau_p);
}

double
Integrand_ttHorZ_3l1tau::bJetEnergy(const LorentzVector & W,
                                    unsigned bIdx) const
{
  const Vector Wp(W.x(), W.y(), W.z());
  const double a_ = DeltaFactor / W.energy() / massB;
  const double b_ = W.Beta() * bJetp3Unit_[bIdx].Dot(Wp.unit());
  const double a_2 = std::pow(a_, 2);
  const double b_2 = std::pow(b_, 2);
  const double b_abs = std::fabs(b_);
  const double a_2b_2_1 = a_2 + b_2 - 1.;

  if((b_2 - 1.) >= 1.e-5)
  {
    LOGVRB << "(b_)^2 = " << b_2 << " >= 1 => Eb[" << bIdx << "] = 0";
    return 0.;
  }
  if(a_2b_2_1 < 0.)
  {
    LOGVRB << "(a_)^2 + (b_)^2 - 1 = " << a_2b_2_1 << " < 1 => Eb[" << bIdx << "] = 0";
    return 0.;
  }

  const double Eb[2] = { massB * (a_ + b_abs * std::sqrt(a_2b_2_1)) / (1 - b_2),
                         massB * (a_ - b_abs * std::sqrt(a_2b_2_1)) / (1 - b_2) };
  if(Eb[0] <= 0. && Eb[1] <= 0.)
  {
    LOGVRB << "Eb+[" << bIdx << "] = " << Eb[0] << " <= 0 and "
           << "Eb-[" << bIdx << "] = " << Eb[1] << " <= 0 "
           << "=> Eb[" << bIdx << "] = 0";
    return 0.;
  }

  const double Eb_a_[2] = { Eb[0] - a_, Eb[1] - a_ };
  std::vector<unsigned> validSolIdx;
  for(unsigned i = 0; i < 2; ++i)
    if(b_ * Eb_a_[i] > 0. && Eb[i] > 0.) validSolIdx.push_back(i);
  if(! validSolIdx.size())
  {
    LOGVRB << "Neither of Eb[" << bIdx << "] values satisfied the solution conditions";
    LOGVRB << "(a' = " << a_ << "; b' = " << b_ << "; "
           << "Eb+ = " << Eb[0] << "; Eb- = " << Eb[1] << ")";
    return 0.;
  }
  else if(validSolIdx.size() == 2)
    return std::max(Eb[0], Eb[1]);

  return Eb[validSolIdx[0]];
}

double
Integrand_ttHorZ_3l1tau::eval(const double * x) const
{
  LOGTRC;
  if(! pdf_)                      LOGERR << "PDF not initialized!";
  if(! me_madgraph_[currentME_])  LOGERR << "Madgraph's ME not initialized!";
  if(! measuredEvent_)            LOGERR << "Measured event not specified!";
  if(! numDimensions_)            LOGERR << "Number of dimensions unspecified!";
  if(mgMomentaIdxs_.size() != 10) LOGERR << "Number of MG momenta indexes not equal to 10";
  if(idxCosTheta1_ < 0)           LOGERR << "Index idxCosTheta1 not set";
  if(idxVarphi1_ < 0)             LOGERR << "Index idxVarphi1 not set";
  if(idxCosTheta2_ < 0)           LOGERR << "Index idxCosTheta2 not set";
  if(idxVarphi2_ < 0)             LOGERR << "Index idxVarphi2 not set";
  if(idxZ1_ < 0)                  LOGERR << "Index idxZ1 not set";
  if(idxPhi1_ < 0)                LOGERR << "Index idxPhi1 not set";
  if(idxPhiInv_ < 0)              LOGERR << "Index idxPhiInv not set";
  if(idxMinvSquared_ < 0)         LOGERR << "Index idxMinvSquared not set";
  if(! pdf_ || ! me_madgraph_[currentME_] || ! measuredEvent_ ||
     idxCosTheta1_ < 0 || idxVarphi1_ < 0 || idxCosTheta2_ < 0 ||
     idxVarphi2_ < 0 || idxZ1_ < 0 || idxPhi1_ < 0 || idxPhiInv_ < 0 ||
     idxMinvSquared_ < 0 || mgMomentaIdxs_.size() != 10)
    std::exit(EXIT_FAILURE);
  LOGVRB << "Current MG5 ME: " << me_madgraph_[currentME_] -> name();
  if(covDet_ == 0.) return 0.;

//--- read the sampled values
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(Logger::getFloatPrecision());
  std::copy(x, x + numDimensions_ - 1, std::ostream_iterator<double>(ss, ", "));
  ss << x[numDimensions_ - 1];
  LOGVRB << "x = { " << ss.str() << " }";

  const double cosTheta1     = x[idxCosTheta1_];
  const double varphi1       = x[idxVarphi1_];
  const double cosTheta2     = x[idxCosTheta2_];
  const double varphi2       = x[idxVarphi2_];
  const double z1            = x[idxZ1_];
  const double nuHtau_phi    = x[idxPhi1_];
  const double nuLTau_phi    = x[idxPhiInv_];
  const double mInvSquared   = x[idxMinvSquared_];

//--- confirm that the energy fraction carried by the tau is indeed in (0,1)
  const double z2 = measuredVisMassSquared_ /
      (z1 * (currentME_ == ME_mg5_3l1tau::kTTH ? massHiggsSquared : massZSquared));
  if(! (z2 >= 1.e-5 && z2 <= 1.))
  {
    LOGVRB << "z2 = " << z2 << " not in (0, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "z2 = " << z2;

//--- compute the neutrino and tau lepton 4-vector from hadronic tau
  const double nuHtau_en = hTauEnergy_ * (1. - z1) / z1;
  const double nuHtau_cosTheta = nuHtauCosTheta(nuHtau_en);
  if(! (nuHtau_cosTheta >= -1. && nuHtau_cosTheta <= +1.))
  {
    LOGVRB << "nuHtau_cosTheta = " << nuHtau_cosTheta << " not in (-1, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "nuHtau_cosTheta = " << nuHtau_cosTheta;
  const double nuHtau_theta = std::acos(nuHtau_cosTheta);

  const VectorSpherical nuHtau_loc(nuHtau_en, nuHtau_theta, nuHtau_phi);
  const double nuHtau_px = nuHtau_loc.Dot(eX_htau_);
  const double nuHtau_py = nuHtau_loc.Dot(eY_htau_);
  const double nuHtau_pz = nuHtau_loc.Dot(eZ_htau_);
  const LorentzVector nuHtau(nuHtau_px, nuHtau_py, nuHtau_pz, nuHtau_en);
  LOGTRC << lvrap("htau nu", nuHtau);

  const LorentzVector hTau = hTauP4_ + nuHtau;
  LOGTRC << lvrap("hadronic tau", hTau);

//--- compute the neutrino and tau lepton 4-vector from leptonic tau
  const double nuLTau_en = complLeptEnergy_ * (1. - z2) / z2;
  const double nuLTau_mass = std::sqrt(mInvSquared);
  const double nuLTau_p = std::sqrt(std::max(0., std::pow(nuLTau_en, 2) -
                                                 std::pow(nuLTau_mass, 2)));
  const double nuLTau_cosTheta = nuLeptTauCosTheta(nuLTau_en, mInvSquared,
                                                   nuLTau_p);
  if(! (nuLTau_cosTheta >= -1. && nuLTau_cosTheta <= +1.))
  {
    LOGVRB << "nuLTau_cosTheta = " << nuLTau_cosTheta << " not in (-1, 1) "
           << "=> p = 0";
    return 0.;
  }
  const double nuLTau_theta = std::acos(nuLTau_cosTheta);
  const VectorSpherical nuLTau_loc(nuLTau_en, nuLTau_theta, nuLTau_phi);
  const double nuLTau_px = nuLTau_loc.Dot(eX_lept_);
  const double nuLTau_py = nuLTau_loc.Dot(eY_lept_);
  const double nuLTau_pz = nuLTau_loc.Dot(eZ_lept_);
  const LorentzVector nuLTau(nuLTau_px, nuLTau_py, nuLTau_pz, nuLTau_en);
  LOGTRC << lvrap("lept tau nu", nuLTau);

  const LorentzVector lTau = complLeptP4_ + nuLTau;
  LOGTRC << lvrap("leptonic tau", lTau);

  const LorentzVector higgs = hTau + lTau;
  LOGTRC << lvrap("higgs", higgs);

//--- get W boson and associated neutrino 4-vectors
  const VectorSpherical nuT1_p3unit(1., std::acos(cosTheta1), varphi1);
  const VectorSpherical nuT2_p3unit(1., std::acos(cosTheta2), varphi2);
  const double nuT1_en = massWSquared / ((1 - lept1p3Unit_.Dot(nuT1_p3unit))
                                         * lept1Energy_ * 2.);
  const double nuT2_en = massWSquared / ((1 - lept2p3Unit_.Dot(nuT2_p3unit))
                                         * lept2Energy_ * 2.);
  const Vector nuT1_p3(nuT1_en * nuT1_p3unit);
  const Vector nuT2_p3(nuT2_en * nuT2_p3unit);
  const LorentzVector nuT1(nuT1_p3.x(), nuT1_p3.y(), nuT1_p3.z(), nuT1_en);
  const LorentzVector nuT2(nuT2_p3.x(), nuT2_p3.y(), nuT2_p3.z(), nuT2_en);
  const LorentzVector W1 = nuT1 + lept1p4_;
  const LorentzVector W2 = nuT2 + lept2p4_;
  LOGTRC << lvrap("t nu 1", nuT1);
  LOGTRC << lvrap("t nu 2", nuT2);
  LOGTRC << lvrap("W 1", W1);
  LOGTRC << lvrap("W 2", W2);

  const double b1_en = bJetEnergy(W1, 0);
  const double b2_en = bJetEnergy(W2, 1);
  if(b1_en == 0. || b2_en == 0.) return 0.;
  const Vector b1_p3 = std::sqrt(std::pow(b1_en, 2) - massBSquared) * bJetp3Unit_[0];
  const Vector b2_p3 = std::sqrt(std::pow(b2_en, 2) - massBSquared) * bJetp3Unit_[1];
  const LorentzVector b1 = LorentzVector(b1_p3.x(), b1_p3.y(), b1_p3.z(), b1_en);
  const LorentzVector b2 = LorentzVector(b2_p3.x(), b2_p3.y(), b2_p3.z(), b2_en);
  LOGTRC << lvrap("b 1", b1);
  LOGTRC << lvrap("b 2", b2);

  const LorentzVector t1 = b1 + W1;
  const LorentzVector t2 = b2 + W2;
  LOGTRC << lvrap("t 1", t1);
  LOGTRC << lvrap("t 2", t2);

  LOGTRC << "b1_en = " << b1_en << " (reco = " << bJetRecoEnergy_[0] << "); "
         << "b2_en = " << b2_en << " (reco = " << bJetRecoEnergy_[1] << ")";

//--- hadronic recoil transfer function; simplification: use only neutrinos
//--- in the difference of ,,measured'' and ,,true'' hadronic recoil, because
//--- leptons are assumed to be measured perfectly whereas the pT of jets/quarks
//--- introduces uncertainties difficult to handle here
  const LorentzVector nuSum = nuHtau + nuLTau + nuT1 + nuT2;
  TVectorD hadRecDiff(2);
  hadRecDiff(0) = MET_x_ - nuSum.x();
  hadRecDiff(1) = MET_y_ - nuSum.y();
  const double MET_pull = (invCovMET_ * hadRecDiff) * hadRecDiff;
  const double MET_TF = MET_TF_denom * std::exp(-MET_pull / 2.);
  LOGTRC << "MET_x = "   << MET_x_    << "; MET_y = "   << MET_y_ << "; "
         << "nuSum_x = " << nuSum.x() << "; nuSum_y = " << nuSum.y();
  LOGTRC_S << "=> MET_pull = " << MET_pull << " => MET_TF = " << MET_TF;

//--- compute Bjorken x variables
//--- assume that hadronic recoil has only transverse component
  const double hadRecE = 0.;
  const double hadRecPz = 0.;
  const LorentzVector tth = higgs + t1 + t2;
  const double xa = invSqrtS_ * (hadRecE + tth.e() + hadRecPz + tth.pz());
  const double xb = invSqrtS_ * (hadRecE + tth.e() - hadRecPz - tth.pz());
  LOGTRC << "xa = " << xa << "; xb = " << xb;

  const double fa = pdf_ -> xfxQ(21, xa, Q_) / xa;
  const double fb = pdf_ -> xfxQ(21, xb, Q_) / xb;
  const double probPDF = fa * fb;
  const double flux = 1. / (s_ * xa * xb);

//--- boost all MG momenta into frame where pT(tth) = 0
  const Vector boost(-tth.px() / tth.e(), -tth.py() / tth.e(), 0.);
  LOGTRC << cvrap("boost vector", boost);
  const LorentzVector gluon1(0., 0., +0.5 * xa * sqrtS_, 0.5 * xa * sqrtS_);
  const LorentzVector gluon2(0., 0., -0.5 * xb * sqrtS_, 0.5 * xb * sqrtS_);
  const LorentzVector b1_mem    = VectorUtil::boost(b1, boost);
  const LorentzVector b2_mem    = VectorUtil::boost(b2, boost);
  const LorentzVector hTau_mem  = VectorUtil::boost(hTau, boost);
  const LorentzVector lTau_mem  = VectorUtil::boost(lTau, boost);
  const LorentzVector nuT1_mem  = VectorUtil::boost(nuT1, boost);
  const LorentzVector nuT2_mem  = VectorUtil::boost(nuT2, boost);
  const LorentzVector lept1_mem = VectorUtil::boost(lept1p4_, boost);
  const LorentzVector lept2_mem = VectorUtil::boost(lept2p4_, boost);

  const LorentzVector t1_mem    = VectorUtil::boost(t1, boost);
  const LorentzVector t2_mem    = VectorUtil::boost(t2, boost);
  const LorentzVector higgs_mem = VectorUtil::boost(higgs, boost);
  const LorentzVector tth_mem = t1_mem + t2_mem + higgs_mem;
  LOGTRC << lvrap("tth mem", tth_mem);

//--- set MG momenta
  const std::vector<LorentzVector> memVector_p4 = {
    gluon1, gluon2, b1_mem, lept1_mem, nuT1_mem,
                    b2_mem, lept2_mem, nuT2_mem,
    lTau_mem, hTau_mem};
  setMGmomenta(memVector_p4);
  me_madgraph_[currentME_] -> setMomenta(mgMomenta_);

  LOGTRC_S << "prob(PDF) = " << probPDF << "; flux factor = " << flux;

//--- calculate the matrix element
  me_madgraph_[currentME_] -> sigmaKin();
  const double prob_ME_mg = me_madgraph_[currentME_] -> getMatrixElements()[0];
  if(TMath::IsNaN(prob_ME_mg) || prob_ME_mg < 0.)
  {
    LOGERR_S << "Warning: MadGraph5 returned NaN or is zero: "
             << "|M|^2 = " << prob_ME_mg << " => skipping event";
    return 0.;
  }
  LOGVRB_S << "|M|^2 = " << prob_ME_mg;

//--- assemble the integrand
  const double p = prob_ME_mg * probPDF * flux * MET_TF;
  LOGVRB_S << "p = " << p;

  return p;
}

