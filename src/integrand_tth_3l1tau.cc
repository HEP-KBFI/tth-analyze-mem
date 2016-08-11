#include "tthAnalysis/tthMEM/interface/integrand_tth_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <iostream> // std::cerr, std::cout
#include <cmath> // std::pow()
#include <cstdlib> // std::exit(), EXIT_FAILURE
#include <cstring> // std::memset()
#include <algorithm> // std::for_each(), std::copy()
#include <sstream> // std::ostringstream

#include "TMath.h" // TMath::IsNaN() (http://stackoverflow.com/a/570694)

using namespace tthMEM;

const integrand_tth_3l1tau * integrand_tth_3l1tau::gIntegrand = 0;

integrand_tth_3l1tau::integrand_tth_3l1tau(double sqrtS,
                                           const std::string & pdfName,
                                           const std::string & madgraphFilename)
  : sqrtS_(sqrtS)
  , s_(std::pow(sqrtS_, 2))
  , beamAxis_(0., 0., 1.)
  , pdf_(0)
  , me_madgraph_initialized_(false)
  , numDimensions_(0)
  , measuredEvent_(0)
  , invCovMET_(TMatrixD(2, 2))
  , idxCosTheta1_(-1)
  , idxVarphi1_(-1)
  , idxCosTheta2_(-1)
  , idxVarphi2_(-1)
  , idxZ1_(-1)
  , idxTh_(-1)
  , idxPhi1_(-1)
  , idxPhiInv_(-1)
  , idxMinvSquared_(-1)
{
  LOGDBG;

  if(! pdf_ && pdfName != "") pdf_ = LHAPDF::mkPDF(pdfName.c_str(), 0);
  else
  {
    LOGERR << "PDF file name empty!";
    std::exit(EXIT_FAILURE);
  }

  if(madgraphFilename != "")
  {
    me_madgraph_.initProc(madgraphFilename);
    me_madgraph_initialized_ = true;
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
  mgMomenta_.push_back(mgLeptonFromBjet2p4_);
  mgMomenta_.push_back(mgNeutrinoFromBjet2p4_);
  mgMomenta_.push_back(mgTau1p4_);
  mgMomenta_.push_back(mgTau2p4_);

  std::for_each(mgMomenta_.begin(), mgMomenta_.end(),
    [](double * & d) { std::memset(d, 0., 4 * sizeof(double)); }
  );

  gIntegrand = this;
}

integrand_tth_3l1tau::~integrand_tth_3l1tau()
{
  LOGDBG;

  if(pdf_)
  {
    delete pdf_;
    pdf_ = 0;
  }
  measuredEvent_ = 0; // no allocation, just the address
  std::for_each(mgMomenta_.begin(), mgMomenta_.end(),
    [](double * & d) { delete d; d = 0; }
  );
}

void
integrand_tth_3l1tau::setNumDimensions(unsigned numDimensions)
{
  numDimensions_ = numDimensions;
}

void
integrand_tth_3l1tau::setIdxCosTheta1(int idx)
{
  idxCosTheta1_ = idx;
}

void
integrand_tth_3l1tau::setIdxVarphi1(int idx)
{
  idxVarphi1_ = idx;
}

void
integrand_tth_3l1tau::setIdxCosTheta2(int idx)
{
  idxCosTheta2_ = idx;
}

void
integrand_tth_3l1tau::setIdxVarphi2(int idx)
{
  idxVarphi2_ = idx;
}

void
integrand_tth_3l1tau::setIdxZ1(int idx)
{
  idxZ1_ = idx;
}

void
integrand_tth_3l1tau::setIdxTh(int idx)
{
  idxTh_ = idx;
}

void
integrand_tth_3l1tau::setIdxPhi1(int idx)
{
  idxPhi1_ = idx;
}

void
integrand_tth_3l1tau::setIdxPhiInv(int idx)
{
  idxPhiInv_ = idx;
}

void
integrand_tth_3l1tau::setIdxMinvSquared(int idx)
{
  idxMinvSquared_ = idx;
}

void
integrand_tth_3l1tau::setInputs(const MeasuredEvent_3l1tau & measuredEvent)
{
  measuredEvent_ = &measuredEvent; // todo: implement copy constructor?

//--- set the variables related to the hadronic tau
  hTauMassSquared_ = std::pow(measuredEvent_ -> htau1.mass(), 2);
  const Vector eZ = measuredEvent_ -> htau1.p3().Unit();
  const Vector eY = eZ.Cross(beamAxis_).Unit();
  const Vector eX = eY.Cross(eZ).Unit();
  LOGDBG << "eX: theta = " << eX.theta() << ", phi = " << eX.phi() << ", norm = " << eX.R();
  LOGDBG << "eY: theta = " << eY.theta() << ", phi = " << eY.phi() << ", norm = " << eY.R();
  LOGDBG << "eZ: theta = " << eZ.theta() << ", phi = " << eZ.phi() << ", norm = " << eZ.R();
  LOGDBG << "eX x eY = " << eX.Cross(eY).R() << " ; "
         << "eX x eZ = " << eX.Cross(eZ).R() << " ; "
         << "eY x eZ = " << eY.Cross(eZ).R();
  eX_x_ = eX.x();
  eX_y_ = eX.y();
  eX_z_ = eX.z();
  eY_x_ = eY.x();
  eY_y_ = eY.y();
  eY_z_ = eY.z();
  eZ_x_ = eZ.x();
  eZ_y_ = eZ.y();
  eZ_z_ = eZ.z();

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
integrand_tth_3l1tau::eval(const double * x) const
{
  LOGDBG;
  if(! pdf_)                     LOGERR << "PDF not initialized!";
  if(! me_madgraph_initialized_) LOGERR << "Madgraph's ME not initialized!";
  if(! measuredEvent_)           LOGERR << "Measured event not specified!";
  if(! numDimensions_)           LOGERR << "Number of dimensions unspecified!";
  if(idxCosTheta1_ < 0)          LOGERR << "Index idxCosTheta1 not set";
  if(idxVarphi1_ < 0)            LOGERR << "Index idxVarphi1 not set";
  if(idxCosTheta2_ < 0)          LOGERR << "Index idxCosTheta2 not set";
  if(idxVarphi2_ < 0)            LOGERR << "Index idxVarphi2 not set";
  if(idxZ1_ < 0)                 LOGERR << "Index idxZ1 not set";
  if(idxTh_ < 0)                 LOGERR << "Index idxTh not set";
  if(idxPhi1_ < 0)               LOGERR << "Index idxPhi1 not set";
  if(idxPhiInv_ < 0)             LOGERR << "Index idxPhiInv not set";
  if(idxMinvSquared_ < 0)        LOGERR << "Index idxMinvSquared not set";
  if(! pdf_ || ! me_madgraph_initialized_ || ! measuredEvent_ ||
     idxCosTheta1_ < 0 || idxVarphi1_ < 0 || idxCosTheta2_ < 0 ||
     idxVarphi2_ < 0 || idxZ1_ < 0 || idxTh_ < 0 || idxPhi1_ < 0 ||
     idxPhiInv_ < 0 || idxMinvSquared_ < 0)
    std::exit(EXIT_FAILURE);

//--- read the sampled values
  std::ostringstream ss;
  std::copy(x, x + numDimensions_ - 1, std::ostream_iterator<double>(ss, ", "));
  ss << x[numDimensions_ - 1];
  LOGDBG << "x = { " << ss.str() << " }";

//  const double cosTheta1      = x[idxCosTheta1_];
//  const double idxVarphi1     = x[idxVarphi1_];
//  const double idxCosTheta2   = x[idxCosTheta2_];
//  const double idxVarphi2     = x[idxVarphi2_];
//  const double idxZ1          = x[idxZ1_];
//  const double idxTh          = x[idxTh_];
//  const double idxPhi1        = x[idxPhi1_];
//  const double idxPhiInv      = x[idxPhiInv_];
//  const double idxMinvSquared = x[idxMinvSquared_];

//--- steps:
//--- 1) reconstruct missing momenta for MG
//--- 2) assemble the integrand and evaluate it
//--- 3) later: permutations over leptons and bjets

  me_madgraph_.setMomenta(mgMomenta_);
  me_madgraph_.sigmaKin();
  const double prob_ME_mg = me_madgraph_.getMatrixElements()[0];
  if(TMath::IsNaN(prob_ME_mg))
  {
    LOGERR << "Warning: MadGraph5 returned NaN => skipping event";
    return 0.;
  }

  return 0.;
}

