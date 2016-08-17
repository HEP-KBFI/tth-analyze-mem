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

using namespace tthMEM;

const Integrand_ttHorZ_3l1tau * Integrand_ttHorZ_3l1tau::gIntegrand = 0;

Integrand_ttHorZ_3l1tau::Integrand_ttHorZ_3l1tau(double sqrtS,
                                                 const std::string & pdfName,
                                                 const std::string & madgraphFilename)
  : sqrtS_(sqrtS)
  , s_(std::pow(sqrtS_, 2))
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

  mgGluon1p4_            = new double[4];
  mgGluon2p4_            = new double[4];
  mgBjet1p4_             = new double[4];
  mgLeptonFromBjet1p4_   = new double[4];
  mgNeutrinoFromBjet1p4_ = new double[4];
  mgBjet2p4_             = new double[4];
  mgLeptonFromBjet2p4_   = new double[4];
  mgNeutrinoFromBjet2p4_ = new double[4];
  mgTau1p4_              = new double[4];
  mgTau2p4_              = new double[4];

  mgMomenta_.push_back(mgGluon1p4_);
  mgMomenta_.push_back(mgGluon2p4_);
  mgMomenta_.push_back(mgBjet1p4_);
  mgMomenta_.push_back(mgLeptonFromBjet1p4_);
  mgMomenta_.push_back(mgNeutrinoFromBjet1p4_);
  mgMomenta_.push_back(mgBjet2p4_);
  mgMomenta_.push_back(mgLeptonFromBjet2p4_);
  mgMomenta_.push_back(mgNeutrinoFromBjet2p4_);
  mgMomenta_.push_back(mgTau1p4_);
  mgMomenta_.push_back(mgTau2p4_);

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
}

void
Integrand_ttHorZ_3l1tau::setEvent(const MeasuredEvent_3l1tau & measuredEvent)
{
  LOGTRC;
  measuredEvent_ = &measuredEvent;

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
  LOGVRB << "htau p4: " << lvrap(hTauP4_);
  LOGVRB << "htau eX: " << svrap(eX_htau_);
  LOGVRB << "htau eY: " << svrap(eY_htau_);
  LOGVRB << "htau eZ: " << svrap(eZ_htau_);
  LOGVRB << "htau: " << "eX x eY = " << eX_htau_.Cross(eY_htau_).r() << " ; "
                     << "eX x eZ = " << eX_htau_.Cross(eZ_htau_).r() << " ; "
                     << "eY x eZ = " << eY_htau_.Cross(eZ_htau_).r();

//--- set the variables related to the MET TF
  invCovMET_ = measuredEvent_ -> met.covMET();
  const double covDet = invCovMET_.Determinant();
  if(covDet != 0)
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
  complLeptIdx_ = measuredEvent_ -> complLeptonIdx;
  const MeasuredLepton & complLepton = measuredEvent_-> leptons[complLeptIdx_];
  complLeptEnergy_ = complLepton.energy();
  complLeptMomentum_ = complLepton.p();
  complLeptMass_ = complLepton.mass();
  complLeptMassSquared_ = std::pow(complLeptMass_, 2);
  complLeptP4_ = complLepton.p4();
  eZ_lept_ = complLepton.p3().unit();
  eY_lept_ = eZ_lept_.Cross(beamAxis_).unit();
  eX_lept_ = eY_lept_.Cross(eZ_lept_).unit();
  // eX should already be unit vector by construction
  LOGVRB << "lept p4: " << lvrap(complLeptP4_);
  LOGVRB << "lept p3: " << cvrap(complLepton.p3());
  LOGVRB << "lept eX: " << svrap(eX_lept_);
  LOGVRB << "lept eY: " << svrap(eY_lept_);
  LOGVRB << "lept eZ: " << svrap(eZ_lept_);
  LOGVRB << "lept: " << "eX x eY = " << eX_lept_.Cross(eY_lept_).R() << " ; "
                     << "eX x eZ = " << eX_lept_.Cross(eZ_lept_).R() << " ; "
                     << "eY x eZ = " << eY_lept_.Cross(eZ_lept_).R();

//--- find the measured mass of both visible tau decay products
  measuredVisMassSquared_ = (hTauP4_ + complLeptP4_).mass2();
  LOGVRB << "measured mass of visible tau decay products: "
         << std::sqrt(measuredVisMassSquared_);

//--- find the kinematic variables of the rest of the leptons (from t decay)
  lept1Idx_ = measuredEvent_ -> bjetLeptonIdxs[0];
  lept2Idx_ = measuredEvent_ -> bjetLeptonIdxs[1];
  const MeasuredLepton & lept1 = measuredEvent_ -> leptons[lept1Idx_];
  const MeasuredLepton & lept2 = measuredEvent_ -> leptons[lept2Idx_];
  lept1p3_ = lept1.p3();
  lept2p3_ = lept2.p3();
  lept1Energy_ = lept1.p();
  lept2Energy_ = lept2.p();
  lept1p3Unit_ = lept1p3_.unit();
  lept2p3Unit_ = lept2p3_.unit();
  lept1p4_ = LorentzVector(lept1p3_.x(), lept1p3_.y(), lept1p3_.z(), lept1Energy_);
  lept2p4_ = LorentzVector(lept2p3_.x(), lept2p3_.y(), lept2p3_.z(), lept2Energy_);
  LOGVRB << "t lept 1 p4: " << lvrap(lept1p4_);
  LOGVRB << "t lept 2 p4: " << lvrap(lept2p4_);

//--- find the kinematic variables of the b quarks
  for(unsigned i = 0; i < 2; ++i)
  {
    const MeasuredJet & bjet_i = measuredEvent_ -> jets[i];
    bJetRecoEnergy_[i] = bjet_i.energy();
    bJetp3Unit_[i] = bjet_i.p3().unit();
  }
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
    LOGWARN << "(b_)^2 = " << b_2 << " >= 1 => Eb[" << bIdx << "] = 0";
    return 0.;
  }
  if(a_2b_2_1 < 0.)
  {
    LOGWARN << "(a_)^2 + (b_)^2 - 1 = " << a_2b_2_1 << " < 1 => Eb[" << bIdx << "] = 0";
    return 0.;
  }

  const double Eb[2] = { massB * (a_ + b_abs * std::sqrt(a_2b_2_1)) / (1 - b_2),
                         massB * (a_ - b_abs * std::sqrt(a_2b_2_1)) / (1 - b_2) };
  if(Eb[0] <= 0. && Eb[1] <= 0.)
  {
    LOGWARN << "Eb+[" << bIdx << "] = " << Eb[0] << " <= 0 and "
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
    LOGWARN << "Neither of Eb[" << bIdx << "] values satisfied the solution conditions";
    LOGWARN << "(a' = " << a_ << "; b' = " << b_ << "; "
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
  if(! pdf_)                     LOGERR << "PDF not initialized!";
  if(! me_madgraph_[currentME_]) LOGERR << "Madgraph's ME not initialized!";
  if(! measuredEvent_)           LOGERR << "Measured event not specified!";
  if(! numDimensions_)           LOGERR << "Number of dimensions unspecified!";
  if(idxCosTheta1_ < 0)          LOGERR << "Index idxCosTheta1 not set";
  if(idxVarphi1_ < 0)            LOGERR << "Index idxVarphi1 not set";
  if(idxCosTheta2_ < 0)          LOGERR << "Index idxCosTheta2 not set";
  if(idxVarphi2_ < 0)            LOGERR << "Index idxVarphi2 not set";
  if(idxZ1_ < 0)                 LOGERR << "Index idxZ1 not set";
  if(idxPhi1_ < 0)               LOGERR << "Index idxPhi1 not set";
  if(idxPhiInv_ < 0)             LOGERR << "Index idxPhiInv not set";
  if(idxMinvSquared_ < 0)        LOGERR << "Index idxMinvSquared not set";
  if(! pdf_ || ! me_madgraph_[currentME_] || ! measuredEvent_ ||
     idxCosTheta1_ < 0 || idxVarphi1_ < 0 || idxCosTheta2_ < 0 ||
     idxVarphi2_ < 0 || idxZ1_ < 0 || idxPhi1_ < 0 || idxPhiInv_ < 0 ||
     idxMinvSquared_ < 0)
    std::exit(EXIT_FAILURE);
  LOGVRB << "Current MG5 ME: " << me_madgraph_[currentME_] -> name();

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
  const double nuLeptTau_phi = x[idxPhiInv_];
  const double mInvSquared   = x[idxMinvSquared_];

//--- confirm that the energy fraction carried by the tau is indeed in (0,1)
  const double z2 = measuredVisMassSquared_ /
      (z1 * (currentME_ == ME_mg5_3l1tau::kTTH ? massHiggsSquared : massZSquared));
  if(! (z2 >= 1.e-5 && z2 <= 1.))
  {
    LOGDBG << "z2 = " << z2 << " not in (0, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "z2 = " << z2;

//--- compute the neutrino and tau lepton 4-vector from hadronic tau
  const double nuHtau_en = hTauEnergy_ * (1. - z1) / z1;
  const double nuHtau_cosTheta = nuHtauCosTheta(nuHtau_en);
  if(! (nuHtau_cosTheta >= -1. && nuHtau_cosTheta <= +1.))
  {
    LOGDBG << "nuHtau_cosTheta = " << nuHtau_cosTheta << " not in (-1, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "nuHtau_cosTheta = " << nuHtau_cosTheta;
  const double nuHtau_theta = std::acos(nuHtau_cosTheta);

  const VectorSpherical nuHtau_loc(nuHtau_en, nuHtau_theta, nuHtau_phi);
  const double nuHtau_px = nuHtau_loc.Dot(eX_htau_);
  const double nuHtau_py = nuHtau_loc.Dot(eY_htau_);
  const double nuHtau_pz = nuHtau_loc.Dot(eZ_htau_);
  const LorentzVector nuHtau(nuHtau_px, nuHtau_py, nuHtau_pz, nuHtau_en);
  LOGTRC << "htau nu: " << lvrap(nuHtau);

  const LorentzVector hTau = hTauP4_ + nuHtau;
  LOGTRC << "hadronic tau: " << lvrap(hTau);

//--- compute the neutrino and tau lepton 4-vector from leptonic tau
  const double nuLeptTau_en = complLeptEnergy_ * (1. - z2) / z2;
  const double nuLeptTau_mass = std::sqrt(mInvSquared);
  const double nuLeptTau_p = std::sqrt(std::max(0., std::pow(nuLeptTau_en, 2) -
                                                    std::pow(nuLeptTau_mass, 2)));
  const double nuLeptTau_cosTheta = nuLeptTauCosTheta(nuLeptTau_en, mInvSquared,
                                                      nuLeptTau_p);
  if(! (nuLeptTau_cosTheta >= -1. && nuLeptTau_cosTheta <= +1.))
  {
    LOGDBG << "nuLeptTau_cosTheta = " << nuLeptTau_cosTheta << " not in (-1, 1) "
           << "=> p = 0";
    return 0.;
  }
  const double nuLeptTau_theta = std::acos(nuLeptTau_cosTheta);
  const VectorSpherical nuLeptTau_loc(nuLeptTau_en, nuLeptTau_theta, nuLeptTau_phi);
  const double nuLeptTau_px = nuLeptTau_loc.Dot(eX_lept_);
  const double nuLeptTau_py = nuLeptTau_loc.Dot(eY_lept_);
  const double nuLeptTau_pz = nuLeptTau_loc.Dot(eZ_lept_);
  const LorentzVector nuLeptTau(nuLeptTau_px, nuLeptTau_py, nuLeptTau_pz, nuLeptTau_en);
  LOGTRC << "lept tau nu: " << lvrap(nuLeptTau);

  const LorentzVector leptTau = complLeptP4_ + nuLeptTau;
  LOGTRC << "leptonic tau: " << lvrap(leptTau);

  const LorentzVector higgs = hTau + leptTau;
  LOGTRC << "higgs: " << lvrap(higgs);

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
  LOGTRC << "t nu 1: " << lvrap(nuT1);
  LOGTRC << "t nu 2: " << lvrap(nuT2);
  LOGTRC << "W 1: " << lvrap(W1);
  LOGTRC << "W 2: " << lvrap(W2);

  const double b1_en = bJetEnergy(W1, 0);
  const double b2_en = bJetEnergy(W2, 1);
  if(b1_en == 0. || b2_en == 0.) return 0.;
  LOGTRC << "b1_en = " << b1_en << " (reco = " << bJetRecoEnergy_[0] << "); "
         << "b2_en = " << b2_en << " (reco = " << bJetRecoEnergy_[1] << ")";
  const Vector b1_p3 = std::sqrt(std::pow(b1_en, 2) - massBSquared) * bJetp3Unit_[0];
  const Vector b2_p3 = std::sqrt(std::pow(b2_en, 2) - massBSquared) * bJetp3Unit_[1];
  const LorentzVector b1 = LorentzVector(b1_p3.x(), b1_p3.y(), b1_p3.z(), b1_en);
  const LorentzVector b2 = LorentzVector(b2_p3.x(), b2_p3.y(), b2_p3.z(), b2_en);
  LOGTRC << "b 1: " << lvrap(b1);
  LOGTRC << "b 2: " << lvrap(b2);

  const LorentzVector t1 = b1 + W1;
  const LorentzVector t2 = b2 + W2;
  LOGTRC << "t 1: " << lvrap(t1);
  LOGTRC << "t 2: " << lvrap(t2);

//--- get b-quark 4-vectors

//  me_madgraph_[currentME_] -> setMomenta(mgMomenta_);
//  me_madgraph_[currentME_] -> sigmaKin();
//  const double prob_ME_mg = me_madgraph_[currentME_] -> getMatrixElements()[0];
//  if(TMath::IsNaN(prob_ME_mg) || prob_ME_mg < 0.)
//  {
//    LOGERR << "Warning: MadGraph5 returned NaN or is zero: "
//           << "|M|^2 = " << prob_ME_mg << " => skipping event";
//    return 0.;
//  }

  return 0.;
}

