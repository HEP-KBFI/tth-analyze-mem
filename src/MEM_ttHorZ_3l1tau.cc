#include "tthAnalysis/tthMEM/interface/MEM_ttHorZ_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVEGAS.h"
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVAMP.h"
#include "tthAnalysis/tthMEM/interface/MEMIntegratorMarkovChain.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // ...
  // ... roundToNearestUInt(), pi()
#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // constants::

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
{
  LOGTRC;
  initialize(pdfName, madgraphFileName);
}

MEM_ttHorZ_3l1tau::MEM_ttHorZ_3l1tau(const std::string & pdfName,
                                     const std::string & madgraphFileName,
                                     VariableManager_3l1tau && vm)
  : vm_(std::move(vm))
{
  LOGTRC;
  initialize(pdfName, madgraphFileName);
}

void
MEM_ttHorZ_3l1tau::initialize(const std::string & pdfName,
                              const std::string madgraphFileName)
{
  LOGTRC;
  integrand_ = new Integrand_ttHorZ_3l1tau(pdfName, madgraphFileName, vm_);
  integrationMode_ = IntegrationMode::kUndefined;
  intAlgo_ = 0;
  numDimensions_ = vm_.getCurrentDim();
  setTF_ = false;
  clock_ = new TBenchmark();
  numSeconds_cpu_ = 0.;
  numSeconds_real_ = 0.;
  numSecondsAccumul_cpu_ = 0.;
  numSecondsAccumul_real_ = 0.;
  nof_calls_ = 0;

  mxMode_ = "uniform";
  nofBatches_ = 100;
  nofChains_ = 1;
  maxCallsStartingPos_ = 1000000;
  epsilon0_ = 1.e-2;
  T0_ = 15.;
  nu_ = 0.71;

  higgsWidth_ = -1.;

  setMaxObjFunctionCalls(20000);
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

MEM_ttHorZ_3l1tau &
MEM_ttHorZ_3l1tau::setIntegrationMode(IntegrationMode integrationMode)
{
  LOGTRC;
  integrationMode_ = integrationMode;
  return *this;
}

MEM_ttHorZ_3l1tau &
MEM_ttHorZ_3l1tau::setIntegrationMode(const std::string & integrationModeString)
{
  LOGTRC;
  if(boost::iequals(integrationModeString, "vegas"))
    integrationMode_ = IntegrationMode::kVEGAS;
  else if(boost::iequals(integrationModeString, "vamp"))
    integrationMode_ = IntegrationMode::kVAMP;
  else if(boost::iequals(integrationModeString, "markovchain"))
    integrationMode_ = IntegrationMode::kMarkovChain;
  else
    integrationMode_ = IntegrationMode::kUndefined;
  return *this;
}

MEM_ttHorZ_3l1tau &
MEM_ttHorZ_3l1tau::setBJetTransferFunction(bool setTF)
{
  setTF_ = setTF;
  return *this;
}

MEM_ttHorZ_3l1tau &
MEM_ttHorZ_3l1tau::setHiggsWidth(double higgsWidth)
{
  if(higgsWidth < 0.)
    throw_line("invalid argument")
      << "Provided 'higgsWidth' = " << higgsWidth << ' '
      << "cannot be a negative number";
  higgsWidth_ = higgsWidth;
  return *this;
}

MEM_ttHorZ_3l1tau &
MEM_ttHorZ_3l1tau::setMaxObjFunctionCalls(unsigned maxObjFunctionCalls)
{
  LOGTRC;
  maxObjFunctionCalls_ = maxObjFunctionCalls;
  numCallsGridOpt_     = roundToNearestUInt(0.20 * maxObjFunctionCalls_);
  numCallsIntEval_     = roundToNearestUInt(0.80 * maxObjFunctionCalls_);
//--- Markov Chain stuff
  nofIterBurnin_       = roundToNearestUInt(0.10 * maxObjFunctionCalls_ / nofChains_);
  nofIterSampling_     = roundToNearestUInt(0.90 * maxObjFunctionCalls_ / nofChains_);
  nofIterSimAnnPhase1_ = roundToNearestUInt(0.20 * nofIterBurnin_);
  nofIterSimAnnPhase2_ = roundToNearestUInt(0.60 * nofIterBurnin_);
  alpha_ = 1. - 10. / nofIterBurnin_;

  return *this;
}

MEM_ttHorZ_3l1tau &
MEM_ttHorZ_3l1tau::setMarkovChainParams(const std::string & mxMode,
                                        unsigned nofBatches,
                                        unsigned nofChains,
                                        unsigned maxCallsStartingPos,
                                        double epsilon0,
                                        double T0,
                                        double nu)
{
  mxMode_ = mxMode;
  nofBatches_ = nofBatches;
  nofChains_ = nofChains;
  maxCallsStartingPos_ = maxCallsStartingPos;
  epsilon0_ = epsilon0;
  T0_ = T0;
  nu_ = nu;
//--- update the nofIter* parameters
  return setMaxObjFunctionCalls(maxObjFunctionCalls_);
}

bool
MEM_ttHorZ_3l1tau::isMarkovChainIntegrator() const
{
  return integrationMode_ == IntegrationMode::kMarkovChain;
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
    throw_line("runtime error") << "Integration mode not set";
  const bool isTTH = currentME == ME_mg5_3l1tau::kTTH;

  LOGINFO << (isTTH ? "[TTH]" : "[TTZ]");
  clock_ -> Reset();
  clock_ -> Start(__PRETTY_FUNCTION__);

//--- calculate the integration variables
  if(ev.generatorLevel)
  {
    ev.generatorLevel -> setBeamAxis(integrand_ -> beamAxis_);
    ev.generatorLevel -> setIntegrationVariables(vm_);
  }

//--- set the integrand
  (*integrand_).setEvent               (ev)
               .setCurrentME           (currentME)
               .setBJetTransferFunction(setTF_);
  if(higgsWidth_ > 0.)
    (*integrand_).setHiggsWidth(higgsWidth_);
  Integrand_ttHorZ_3l1tau::gIntegrand = integrand_;

//--- get the lower and upper corners of the integration hyperspace
  const double * const xl = vm_.getXL();
  const double * const xu = vm_.getXU();

//--- create probability and corresponding error (uncertainty) variable
  double pSum = 0.;
  double pSumErr = 0;

//--- integration settings
  if(numDimensions_)
  {
    if(integrationMode_ == IntegrationMode::kVEGAS)
      intAlgo_ = new MEMIntegratorVEGAS(numCallsGridOpt_, numCallsIntEval_, 1, 2.);
    else if(integrationMode_ == IntegrationMode::kVAMP)
      intAlgo_ = new MEMIntegratorVAMP(numCallsGridOpt_, numCallsIntEval_);
    else if(integrationMode_ == IntegrationMode::kMarkovChain)
    {
      intAlgo_ = new MEMIntegratorMarkovChain(
        mxMode_,
        nofIterBurnin_,       nofIterSampling_,
        nofIterSimAnnPhase1_, nofIterSimAnnPhase2_,
        maxCallsStartingPos_,
        T0_,                  alpha_,
        nofChains_,           nofBatches_,
        epsilon0_,            nu_
      );
    }
  }

//--- loop over different permutations of same-sign leptons and b-jets
///< @note consider moving the permutation part inside integrate()
  for(; ev.hasNextPermutation(); ev.nextPermutation())
  {
    if(! ev.isCorrectPermutation())
      continue;

    LOGDBG << ev;
    integrand_ -> renewInputs();

    double p = 0.;
    double pErr = 0.;

    if(numDimensions_)
    {
      if(integrationMode_ == IntegrationMode::kVEGAS ||
         integrationMode_ == IntegrationMode::kMarkovChain)
        intAlgo_ -> integrate(&g_C, xl, xu, numDimensions_, p, pErr);
      else if(integrationMode_ == IntegrationMode::kVAMP)
        intAlgo_ -> integrate(&g_Fortran, xl, xu, numDimensions_, p, pErr);
    }
    else
      p = integrand_ -> eval(0);

    LOGINFO_S << "p = " << p << "; pErr = " << pErr;
    pSum += p;
    pSumErr += pErr;
  }
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

