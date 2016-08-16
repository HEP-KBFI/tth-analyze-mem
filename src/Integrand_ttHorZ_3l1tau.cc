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
 // ... (why not use std: http://stackoverflow.com/a/570694)

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
  eZ_htau_ = htau.p3().Unit();
  eY_htau_ = eZ_htau_.Cross(beamAxis_).Unit();
  eX_htau_ = eY_htau_.Cross(eZ_htau_).Unit();
  // eX should already be unit vector by construction
  LOGVRB << "htau p3 = (" << measuredEvent_ -> htau1.p3().x() << ", "
                          << measuredEvent_ -> htau1.p3().y() << ", "
                          << measuredEvent_ -> htau1.p3().z() << ")";
  LOGVRB << "htau eX: " << "theta = " << eX_htau_.theta() << ", "
                        << "phi = "   << eX_htau_.phi() << ", "
                        << "norm = "  << eX_htau_.R();
  LOGVRB << "htau eY: " << "theta = " << eY_htau_.theta() << ", "
                        << "phi = "   << eY_htau_.phi() << ", "
                        << "norm = "  << eY_htau_.R();
  LOGVRB << "htau eZ: " << "theta = " << eZ_htau_.theta() << ", "
                        << "phi = "   << eZ_htau_.phi() << ", "
                        << "norm = "  << eZ_htau_.R();
  LOGVRB << "htau " << "eX x eY = " << eX_htau_.Cross(eY_htau_).R() << " ; "
                    << "eX x eZ = " << eX_htau_.Cross(eZ_htau_).R() << " ; "
                    << "eY x eZ = " << eY_htau_.Cross(eZ_htau_).R();

//--- set the variables related to the MET TF
  invCovMET_ = measuredEvent_ -> met.covMET();
  const double covDet = invCovMET_.Determinant();
  if(covDet != 0)
  {
    invCovMET_.Invert();
    MET_TF_denom = 1. / (2. * pi() * std::sqrt(invCovMET_.Determinant()));
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
  eZ_lept_ = complLepton.p3().Unit();
  eY_lept_ = eZ_lept_.Cross(beamAxis_).Unit();
  eX_lept_ = eY_lept_.Cross(eZ_lept_).Unit();
  // eX should already be unit vector by construction
  LOGVRB << "lept p3 = (" << complLepton.p3().x() << ", "
                          << complLepton.p3().y() << ", "
                          << complLepton.p3().z() << ")";
  LOGVRB << "lept eX: " << "theta = " << eX_lept_.theta() << ", "
                        << "phi = "   << eX_lept_.phi() << ", "
                        << "norm = "  << eX_lept_.R();
  LOGVRB << "lept eY: " << "theta = " << eY_lept_.theta() << ", "
                        << "phi = "   << eY_lept_.phi() << ", "
                        << "norm = "  << eY_lept_.R();
  LOGVRB << "lept eZ: " << "theta = " << eZ_lept_.theta() << ", "
                        << "phi = "   << eZ_lept_.phi() << ", "
                        << "norm = "  << eZ_lept_.R();
  LOGVRB << "lept " << "eX x eY = " << eX_lept_.Cross(eY_lept_).R() << " ; "
                    << "eX x eZ = " << eX_lept_.Cross(eZ_lept_).R() << " ; "
                    << "eY x eZ = " << eY_lept_.Cross(eZ_lept_).R();

//--- find the meaasured mass of both visible tau decay products
  measuredVisMassSquared_ = (hTauP4_ + complLeptP4_).mass2();
  LOGVRB << "measured mass of visible tau decay products: "
         << std::sqrt(measuredVisMassSquared_);
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

//  const double cosTheta1     = x[idxCosTheta1_];
//  const double varphi1       = x[idxVarphi1_];
//  const double cosTheta2     = x[idxCosTheta2_];
//  const double varphi2       = x[idxVarphi2_];
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
  LOGTRC << "htau nu: En = " << nuHtau_en << "; pT = " << nuHtau.pt();

  const LorentzVector hTau = hTauP4_ + nuHtau;
  LOGTRC << "hadronic tau: En = " << hTau.energy() << "; pT = " << hTau.pt();

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
  LOGTRC << "lept tau nu: En = " << nuLeptTau_en << "; pT = " << nuLeptTau.pt();

  const LorentzVector leptTau = complLeptP4_ + nuLeptTau;
  LOGTRC << "leptonic tau: En = " << leptTau.energy() << "; pT = " << leptTau.pt();

//--- steps:
//--- 1) reconstruct missing momenta for MG
//--- 2) assemble the integrand and evaluate it

//

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

