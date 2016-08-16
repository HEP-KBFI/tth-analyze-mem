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
}

void
Integrand_ttHorZ_3l1tau::renewInputs()
{
  LOGTRC;
//--- set the variables related to the hadronic tau
  const MeasuredHadronicTau & htau = measuredEvent_ -> htau1;
  hTauMass_ = htau.mass();
  hTauMassSquared_ = std::pow(hTauMass_, 2);
  const Vector eZ_htau = htau.p3().Unit();
  const Vector eY_htau = eZ_htau.Cross(beamAxis_).Unit();
  const Vector eX_htau = eY_htau.Cross(eZ_htau).Unit();
  // eX should already be unit vector by construction
  LOGVRB << "htau p3 = (" << measuredEvent_ -> htau1.p3().x() << ", "
                          << measuredEvent_ -> htau1.p3().y() << ", "
                          << measuredEvent_ -> htau1.p3().z() << ")";
  LOGVRB << "htau eX: " << "theta = " << eX_htau.theta() << ", "
                        << "phi = "   << eX_htau.phi() << ", "
                        << "norm = "  << eX_htau.R();
  LOGVRB << "htau eY: " << "theta = " << eY_htau.theta() << ", "
                        << "phi = "   << eY_htau.phi() << ", "
                        << "norm = "  << eY_htau.R();
  LOGVRB << "htau eZ: " << "theta = " << eZ_htau.theta() << ", "
                        << "phi = "   << eZ_htau.phi() << ", "
                        << "norm = "  << eZ_htau.R();
  LOGVRB << "htau " << "eX x eY = " << eX_htau.Cross(eY_htau).R() << " ; "
                    << "eX x eZ = " << eX_htau.Cross(eZ_htau).R() << " ; "
                    << "eY x eZ = " << eY_htau.Cross(eZ_htau).R();
  eX_x_htau_ = eX_htau.x();
  eX_y_htau_ = eX_htau.y();
  eX_z_htau_ = eX_htau.z();
  eY_x_htau_ = eY_htau.x();
  eY_y_htau_ = eY_htau.y();
  eY_z_htau_ = eY_htau.z();
  eZ_x_htau_ = eZ_htau.x();
  eZ_y_htau_ = eZ_htau.y();
  eZ_z_htau_ = eZ_htau.z();

  complLeptIdx_ = measuredEvent_ -> complLeptonIdx;
  const MeasuredLepton & complLepton = measuredEvent_-> leptons[complLeptIdx_];
  complLeptMass_ = complLepton.mass();
  complLeptMassSquared_ = std::pow(complLeptMass_, 2);
  const Vector eZ_lept = complLepton.p3().Unit();
  const Vector eY_lept = eZ_lept.Cross(beamAxis_).Unit();
  const Vector eX_lept = eY_lept.Cross(eZ_lept).Unit();
  // eX should already be unit vector by construction
  LOGVRB << "lept p3 = (" << complLepton.p3().x() << ", "
                          << complLepton.p3().y() << ", "
                          << complLepton.p3().z() << ")";
  LOGVRB << "lept eX: " << "theta = " << eX_lept.theta() << ", "
                        << "phi = "   << eX_lept.phi() << ", "
                        << "norm = "  << eX_lept.R();
  LOGVRB << "lept eY: " << "theta = " << eY_lept.theta() << ", "
                        << "phi = "   << eY_lept.phi() << ", "
                        << "norm = "  << eY_lept.R();
  LOGVRB << "lept eZ: " << "theta = " << eZ_lept.theta() << ", "
                        << "phi = "   << eZ_lept.phi() << ", "
                        << "norm = "  << eZ_lept.R();
  LOGVRB << "lept " << "eX x eY = " << eX_lept.Cross(eY_lept).R() << " ; "
                    << "eX x eZ = " << eX_lept.Cross(eZ_lept).R() << " ; "
                    << "eY x eZ = " << eY_lept.Cross(eZ_lept).R();
  eX_x_lept_ = eX_lept.x();
  eX_y_lept_ = eX_lept.y();
  eX_z_lept_ = eX_lept.z();
  eY_x_lept_ = eY_lept.x();
  eY_y_lept_ = eY_lept.y();
  eY_z_lept_ = eY_lept.z();
  eZ_x_lept_ = eZ_lept.x();
  eZ_y_lept_ = eZ_lept.y();
  eZ_z_lept_ = eZ_lept.z();

//--- find the meaasured mass of both visible tau decay products
  measuredVisMassSquared_ = (htau.p4() + complLepton.p4()).mass2();
  LOGVRB << "measured mass of visible tau decay products: "
         << std::sqrt(measuredVisMassSquared_);

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

//  const double cosTheta1      = x[idxCosTheta1_];
//  const double varphi1     = x[idxVarphi1_];
//  const double cosTheta2   = x[idxCosTheta2_];
//  const double varphi2     = x[idxVarphi2_];
  const double z1          = x[idxZ1_];
//  const double phi1        = x[idxPhi1_];
//  const double phiInv      = x[idxPhiInv_];
//  const double mInvSquared = x[idxMinvSquared_];

//--- confirm that the energy fraction carried by the tau is indeed in (0,1)
  const double z2 = measuredVisMassSquared_ / (massHiggsSquared * z1);
  if(! (z2 >= 1.e-5 && z2 <= 1.))
  {
    LOGDBG << "z2 = " << z2 << " is not in (0, 1) => p = 0";
    return 0.;
  }
  else LOGTRC << "z2 = " << z2;

//--- steps:
//--- 1) reconstruct missing momenta for MG
//--- 2) assemble the integrand and evaluate it
//--- 3) later: permutations over leptons and bjets

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

