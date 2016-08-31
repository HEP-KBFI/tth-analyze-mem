#include "tthAnalysis/tthMEM/interface/MEM_ttHorZ_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVEGAS.h"
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVAMP.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // ...
  // ... roundToNearestUInt(), pi(), massTauSquared, xSectionTTHinGeV2

#include <cmath> // std::round()

#include <boost/algorithm/string/predicate.hpp> // boost::iequals()

using namespace tthMEM;

double
g_C(double * x,
    std::size_t dimension,
    void * additionalParameters)
{
  const double returnValue = Integrand_ttHorZ_3l1tau::gIntegrand -> eval(x);
  return returnValue;
}

double
g_Fortran(double ** x,
          std::size_t dimension,
          void ** additionalParameters)
{
  const double returnValue = Integrand_ttHorZ_3l1tau::gIntegrand -> eval(*x);
  return returnValue;
}

MEM_ttHorZ_3l1tau::MEM_ttHorZ_3l1tau(const std::string & pdfName,
                                     const std::string & madgraphFileName,
                                     const VariableManager_3l1tau & vm)
  : vm_(vm)
  , integrand_(new Integrand_ttHorZ_3l1tau(pdfName, madgraphFileName, vm_))
  , integrationMode_(IntegrationMode::kUndefined)
  , intAlgo_(0)
  , maxObjFunctionCalls_(20000)
  , numDimensions_(vm_.getCurrentDim())
  , precision_(1.e-3)
  , setTF_(false)
  , clock_(new TBenchmark())
  , numSeconds_cpu_(0.)
  , numSeconds_real_(0.)
  , numSecondsAccumul_cpu_(0.)
  , numSecondsAccumul_real_(0.)
  , nof_calls_(0)
{
  LOGTRC;
}

MEM_ttHorZ_3l1tau::MEM_ttHorZ_3l1tau(const std::string & pdfName,
                                     const std::string & madgraphFileName,
                                     VariableManager_3l1tau && vm)
  : vm_(std::move(vm))
  , integrand_(new Integrand_ttHorZ_3l1tau(pdfName, madgraphFileName, vm_))
  , integrationMode_(IntegrationMode::kUndefined)
  , intAlgo_(0)
  , maxObjFunctionCalls_(20000)
  , numDimensions_(vm_.getCurrentDim())
  , precision_(1.e-3)
  , setTF_(false)
  , clock_(new TBenchmark())
  , numSeconds_cpu_(0.)
  , numSeconds_real_(0.)
  , numSecondsAccumul_cpu_(0.)
  , numSecondsAccumul_real_(0.)
  , nof_calls_(0)
{
  LOGTRC;
}

MEM_ttHorZ_3l1tau::~MEM_ttHorZ_3l1tau()
{
  LOGTRC;
  if(integrand_)
  {
    delete integrand_;
    integrand_ = 0;
  }
  if(clock_)
  {
    delete clock_;
    clock_ = 0;
  }
}

void
MEM_ttHorZ_3l1tau::setIntegrationMode(IntegrationMode integrationMode)
{
  LOGTRC;
  integrationMode_ = integrationMode;
}

void
MEM_ttHorZ_3l1tau::setIntegrationMode(const std::string & integrationModeString)
{
  LOGTRC;
  if(boost::iequals(integrationModeString, "vegas"))
    integrationMode_ = IntegrationMode::kVEGAS;
  else if(boost::iequals(integrationModeString, "vamp"))
    integrationMode_ = IntegrationMode::kVAMP;
  else
    integrationMode_ = IntegrationMode::kUndefined;
}

void
MEM_ttHorZ_3l1tau::setBJetTransferFunction(bool setTF)
{
  setTF_ = setTF;
}

void
MEM_ttHorZ_3l1tau::setMaxObjFunctionCalls(unsigned maxObjFunctionCalls)
{
  LOGTRC;
  maxObjFunctionCalls_ = maxObjFunctionCalls;
}

double
MEM_ttHorZ_3l1tau::getComputingTime_cpu() const
{
  return numSeconds_cpu_;
}

double
MEM_ttHorZ_3l1tau::getComputingTime_real() const
{
  return numSeconds_real_;
}

double
MEM_ttHorZ_3l1tau::getAverageComputingTime_cpu() const
{
  return nof_calls_ != 0 ? numSecondsAccumul_cpu_ / nof_calls_ : 0.;
}

double
MEM_ttHorZ_3l1tau::getAverageComputingTime_real() const
{
  return nof_calls_ != 0 ? numSecondsAccumul_real_ / nof_calls_ : 0.;
}

double
MEM_ttHorZ_3l1tau::integrate(const MeasuredEvent_3l1tau & ev,
                             ME_mg5_3l1tau currentME)
{
  LOGTRC;
  if(integrationMode_ == IntegrationMode::kUndefined)
  {
    LOGERR << "Integration mode not set";
    assert(0);
  }
  const bool isTTH = currentME == ME_mg5_3l1tau::kTTH;

  LOGINFO << (isTTH ? "[TTH]" : "[TTZ]");
  clock_ -> Reset();
  clock_ -> Start(__PRETTY_FUNCTION__);

//--- set the integrand
  (*integrand_).setEvent               (ev)
               .setCurrentME           (currentME)
               .setBJetTransferFunction(setTF_);
  Integrand_ttHorZ_3l1tau::gIntegrand = integrand_;

//--- get the lower and upper corners of the integration hyperspace
  const double * const xl = vm_.getXL();
  const double * const xu = vm_.getXU();

//--- create probability and corresponding error (uncertainty) variable
  double pSum = 0.;
  double pSumErr = 0;

//--- set integration settings
  const unsigned numCallsGridOpt = roundToNearestUInt(0.20 * maxObjFunctionCalls_);
  const unsigned numCallsIntEval = roundToNearestUInt(0.80 * maxObjFunctionCalls_);
  if(numDimensions_)
  {
    if(integrationMode_ == IntegrationMode::kVEGAS)
      intAlgo_ = new MEMIntegratorVEGAS(numCallsGridOpt, numCallsIntEval, 1, 2.);
    else if(integrationMode_ == IntegrationMode::kVAMP)
      intAlgo_ = new MEMIntegratorVAMP(numCallsGridOpt, numCallsIntEval);
  }

//--- loop over different permutations of same-sign leptons and b-jets
///< @note consider moving the permutation part inside integrate()
  do
  {
    LOGDBG << ev;
    integrand_ -> renewInputs();

    double p = 0.;
    double pErr = 0.;

    if(numDimensions_)
    {
      if(integrationMode_ == IntegrationMode::kVEGAS)
        intAlgo_ -> integrate(&g_C, xl, xu, numDimensions_, p, pErr);
      else if(integrationMode_ == IntegrationMode::kVAMP)
        intAlgo_ -> integrate(&g_Fortran, xl, xu, numDimensions_, p, pErr);
    }
    else
      p = integrand_ -> eval(0);

    LOGINFO_S << "p = " << p << "; pErr = " << pErr;
    pSum += p;
    pSumErr += pErr;

    ev.nextPermutation();
  }
  while(ev.hasNextPermutation());
  ev.resetPermutation();

//--- divide by the process cross section (in GeV, i.e. natural units)
//--- and adjust the uncertainty accordingly
  pSum /= isTTH ? constants::xSectionTTH2diTauInGeV2 :
                  constants::xSectionTTZinGeV2;
  pSumErr /= isTTH ? constants::xSectionTTH2diTauInGeV2 :
                     constants::xSectionTTZinGeV2;
  LOGINFO_S << "Summed over permutations: p = " << pSum << "; pErr = " << pSumErr;

//--- cleanup
  if(intAlgo_)
  {
    delete intAlgo_;
    intAlgo_ = 0;
  }

  clock_ -> Stop(__PRETTY_FUNCTION__);
  numSeconds_cpu_ = clock_ -> GetCpuTime(__PRETTY_FUNCTION__);
  numSeconds_real_ = clock_ -> GetRealTime(__PRETTY_FUNCTION__);
  LOGINFO << "Real time:\t" << std::round(numSeconds_real_ * 100) / 100 << " s;\t"
          << "CPU time:\t" << std::round(numSeconds_cpu_ * 100) / 100 << " s";

  numSecondsAccumul_cpu_ += numSeconds_cpu_;
  numSecondsAccumul_real_ += numSeconds_real_;
  ++nof_calls_;

  return pSum;
}

