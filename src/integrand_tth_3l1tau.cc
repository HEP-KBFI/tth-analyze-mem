#include "tthAnalysis/tthMEM/interface/integrand_tth_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <iostream> // std::cerr, std::cout
#include <cmath> // std::pow()
#include <cstdlib> // std::exit(), EXIT_FAILURE

using namespace tthMEM;

const integrand_tth_3l1tau * integrand_tth_3l1tau::gIntegrand = 0;

integrand_tth_3l1tau::integrand_tth_3l1tau(double sqrtS,
                                           const std::string & pdfName,
                                           const std::string & madgraphFilename)
  : sqrtS_(sqrtS)
  , s_(std::pow(sqrtS_, 2))
  , pdf_(0)
  , me_madgraph_initialized_(false)
  , measuredEvent_(0)
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

  gIntegrand = this;
}

integrand_tth_3l1tau::~integrand_tth_3l1tau()
{
  LOGDBG;

  delete pdf_;
  pdf_ = 0;
}

void
integrand_tth_3l1tau::setInputs(const MeasuredEvent_3l1tau & measuredEvent)
{
  measuredEvent_ = &measuredEvent;
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

double
integrand_tth_3l1tau::eval(const double * x) const
{
  if(! pdf_)
  {
    LOGERR << "PDF not initialized!";
    std::exit(EXIT_FAILURE);
  }

  if(! me_madgraph_initialized_)
  {
    LOGERR << "Madgraph's ME not initialized!";
    std::exit(EXIT_FAILURE);
  }

  if(! measuredEvent_)
  {
    LOGERR << "Measured event not specified!";
    std::exit(EXIT_FAILURE);
  }

  return 0.;
}

