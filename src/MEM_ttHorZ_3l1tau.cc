#include "tthAnalysis/tthMEM/interface/MEM_ttHorZ_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVEGAS.h"
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVAMP.h"
#include "tthAnalysis/tthMEM/interface/MEMIntegratorMarkovChain.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // ...
  // ... roundToNearestUInt(), pi()
#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // constants::
#include "tthAnalysis/tthMEM/interface/tthMEMvecFunctions.h" // vec::

#include <cmath> // std::round()

#include <boost/algorithm/string/predicate.hpp> // boost::iequals()

using namespace tthMEM;

double
g_C(double * x,
    std::size_t dimension,
    void * additionalParameters)
{
  return Integrand_ttHorZ_3l1tau::gIntegrand -> eval(x);
}

double
g_Fortran(double ** x,
          std::size_t dimension,
          void ** additionalParameters)
{
  return Integrand_ttHorZ_3l1tau::gIntegrand -> eval(*x);
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
  useAvgB_ = false;
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
MEM_ttHorZ_3l1tau::useAvgBjetCombo(bool useAvgB)
{
  useAvgB_ = useAvgB;
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

unsigned
MEM_ttHorZ_3l1tau::getNofMXMCTries() const
{
  return nofMXMCTries;
}

std::array<double, 2>
MEM_ttHorZ_3l1tau::integrate(const MeasuredEvent_3l1tau & ev,
                             ME_mg5_3l1tau currentME,
                             bool & err)
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
  bool hasFoundPermutation = false;
//--- collect the results for different jet combinations and aggregate them (use max or average)
  std::vector<double> pSum_allJc, pSumErr_allJc;
  for(; ev.hasNextJetCombination(); ev.nextJetCombination())
  {
    LOGTRC << "Picking jet combination #" << (ev.getJetCombinationNumber() + 1);
    ev.getJetCombination();
    bool hasFoundLocalPermutation = false;

    double pSum_perJc = 0.;
    double pSumErr_perJc = 0;
    for(; ev.hasNextPermutation(); ev.nextPermutation())
    {
//--- do not evaluate the same permutation twice
      if(hasFoundLocalPermutation && ev.generatorLevel) break;

      ev.printPermutation();
      if(! ev.isCorrectPermutation())
      {
        LOGTRC << "Skipping permutation #" << (ev.getPermutationNumber() + 1);
        continue;
      }

      LOGTRC << "Found correct permutation";
      hasFoundLocalPermutation = true;
      hasFoundPermutation = true;

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
      pSum_perJc += p;
      pSumErr_perJc += pErr;

      if(integrationMode_ == IntegrationMode::kMarkovChain)
        nofMXMCTries = *static_cast<unsigned *>(intAlgo_ -> metadata());
      else
        nofMXMCTries = -1;
    }
    ev.resetPermutation();
    pSum_allJc.push_back(pSum_perJc);
    pSumErr_allJc.push_back(pSumErr_perJc);
  }
  ev.resetJetCombination();

  if(! hasFoundPermutation)
  {
    LOGWARN << "NB! Event was skipped altogether b/c there was no correct permutation found";
    err = true;
    return {{ 0., 0. }};
  }

//--- aggregate probability and corresponding error (uncertainty) variable
  LOGTRC_S << "p    (jet combinations) = " << pSum_allJc;
  LOGTRC_S << "pErr (jet combinations) = " << pSumErr_allJc;
  const unsigned pSum_maxIdx = vec::maxIdx(pSum_allJc);
  const double pSum_max    = pSum_allJc[pSum_maxIdx];
  const double pSumErr_max = pSumErr_allJc[pSum_maxIdx];
  const double pSum_avg    = vec::avg(pSum_allJc);
  const double pSumErr_avg = vec::avg(pSumErr_allJc);
  LOGTRC_S << "avg p = " << pSum_avg << " (" << pSumErr_avg << ')';
  LOGTRC_S << "max p = " << pSum_max << " (" << pSumErr_max << ')';

  double pSum    = useAvgB_ ? pSum_avg    : pSum_max;
  double pSumErr = useAvgB_ ? pSumErr_avg : pSumErr_max;

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

  return {{ pSum, pSumErr }};
}

