#include "tthAnalysis/tthMEM/interface/MEM_tth_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOGDBG

#include <cmath> // std::round()

#include <boost/algorithm/string/predicate.hpp> // boost::iequals()

using namespace tthMEM;

MEM_tth_3l1tau::MEM_tth_3l1tau(double sqrtS,
                               const string & pdfName,
                               const string & madgraphFileName)
  : integrand_(new integrand_tth_3l1tau(sqrtS, pdfName, madgraphFileName))
  , sqrtS_(sqrtS)
  , integrationMode_(IntegrationMode::kUndefined)
  , intAlgo_(0)
  , maxObjFunctionCalls_(20000)
  , numDimensions_(0)
  , precision_(1.e-3)
  , xl_(0)
  , xu_(0)
  , clock_(new TBenchmark())
  , numSeconds_cpu_(0)
  , numSeconds_real_(0)
{}

MEM_tth_3l1tau::~MEM_tth_3l1tau()
{
  if(integrand_) delete integrand_;
  if(clock_)     delete clock_;
}

void
MEM_tth_3l1tau::setIntegrationMode(IntegrationMode integrationMode)
{
  integrationMode_ = integrationMode;
}

void
MEM_tth_3l1tau::setIntegrationMode(const string & integrationModeString)
{
  if(boost::iequals(integrationModeString, "vegas"))
    integrationMode_ = IntegrationMode::kVEGAS;
  else if(boost::iequals(integrationModeString, "vamp"))
    integrationMode_ = IntegrationMode::kVAMP;
  else
    integrationMode_ = IntegrationMode::kUndefined;
}

void
MEM_tth_3l1tau::setMaxObjFunctionCalls(unsigned maxObjFunctionCalls)
{
  maxObjFunctionCalls_ = maxObjFunctionCalls;
}

double
MEM_tth_3l1tau::getComputingTime_cpu() const
{
  return numSeconds_cpu_;
}

double
MEM_tth_3l1tau::getComputingTime_real() const
{
  return numSeconds_real_;
}

double
MEM_tth_3l1tau::integrate(const MeasuredEvent_3l1tau & ev)
{
  LOGDBG;
  clock_ -> Reset();
  clock_ -> Start(__PRETTY_FUNCTION__);

  if(integrationMode_ == IntegrationMode::kUndefined)
  {
    LOGERR << "Integration mode not set";
    assert(0);
  }

  /* compute me */

  if(xl_) delete [] xl_;
  if(xu_) delete [] xu_;
  if(intAlgo_) delete intAlgo_;

  clock_ -> Stop(__PRETTY_FUNCTION__);
  numSeconds_cpu_ = clock_ -> GetCpuTime(__PRETTY_FUNCTION__);
  numSeconds_real_ = clock_ -> GetRealTime(__PRETTY_FUNCTION__);
  LOGINFO << "Real time:\t" << std::round(numSeconds_real_ * 100) / 100 << " s;\t"
          << "CPU time:\t" << std::round(numSeconds_cpu_ * 100) / 100 << " s";

  return 0.;
}

