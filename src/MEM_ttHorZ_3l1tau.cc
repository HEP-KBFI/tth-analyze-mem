#include "tthAnalysis/tthMEM/interface/MEM_ttHorZ_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVEGAS.h"
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVAMP.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // ...
  // ... roundToNearestUInt(), pi(), tauLeptonMassSquared, xSectionTTHinGeV2

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

MEM_ttHorZ_3l1tau::MEM_ttHorZ_3l1tau(double sqrtS,
                                     const std::string & pdfName,
                                     const std::string & madgraphFileName)
  : integrand_(new Integrand_ttHorZ_3l1tau(sqrtS, pdfName, madgraphFileName))
  , sqrtS_(sqrtS)
  , integrationMode_(IntegrationMode::kUndefined)
  , intAlgo_(0)
  , maxObjFunctionCalls_(20000)
  , numDimensions_(0)
  , precision_(1.e-3)
  , xl_(0)
  , xu_(0)
  , clock_(new TBenchmark())
  , numSeconds_cpu_(0.)
  , numSeconds_real_(0.)
  , numSecondsAccumul_cpu_(0.)
  , numSecondsAccumul_real_(0.)
  , nof_calls_(0)
{}

MEM_ttHorZ_3l1tau::~MEM_ttHorZ_3l1tau()
{
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
  integrationMode_ = integrationMode;
}

void
MEM_ttHorZ_3l1tau::setIntegrationMode(const std::string & integrationModeString)
{
  if(boost::iequals(integrationModeString, "vegas"))
    integrationMode_ = IntegrationMode::kVEGAS;
  else if(boost::iequals(integrationModeString, "vamp"))
    integrationMode_ = IntegrationMode::kVAMP;
  else
    integrationMode_ = IntegrationMode::kUndefined;
}

void
MEM_ttHorZ_3l1tau::setMaxObjFunctionCalls(unsigned maxObjFunctionCalls)
{
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
  if(integrationMode_ == IntegrationMode::kUndefined)
  {
    LOGERR << "Integration mode not set";
    assert(0);
  }

  LOGVRB;
  clock_ -> Reset();
  clock_ -> Start(__PRETTY_FUNCTION__);
  LOGDBG << ev;

//---  it is (rightly) assumed that leptons and jets
//---  have already been sorted by their pt in decreasing order

//--- there are nine variables we have to sample over:
  numDimensions_ = 9;
  const int idxCosTheta1 = 0;   // cosine of polar angle of the 1st neutrino
  const int idxVarphi1 = 1;     // azimuthal angle of the 1st neutrino
  const int idxCosTheta2 = 2;   // cosine of polar angle of the 2nd neutrino
  const int idxVarphi2 = 3;     // azimuthal angle of the 2nd neutrino
  const int idxZ1 = 4;          // energy fraction of hadronic tau system
  const int idxTh = 5;          // aux variable derived from z2
  const int idxPhi1 = 6;        // rotation angle of hadronic tau neutrino
  const int idxPhiInv = 7;      // (invisible) rotation angle of leptonic tau nu
  const int idxMinvSquared = 8; // (invisible) mass of leptonic neutrino pair

  integrand_ -> setNumDimensions(numDimensions_);
  integrand_ -> setInputs(ev);
  integrand_ -> setIdxCosTheta1  (idxCosTheta1);
  integrand_ -> setIdxVarphi1    (idxVarphi1);
  integrand_ -> setIdxCosTheta2  (idxCosTheta2);
  integrand_ -> setIdxVarphi2    (idxVarphi2);
  integrand_ -> setIdxZ1         (idxZ1);
  integrand_ -> setIdxTh         (idxTh);
  integrand_ -> setIdxPhi1       (idxPhi1);
  integrand_ -> setIdxPhiInv     (idxPhiInv);
  integrand_ -> setIdxMinvSquared(idxMinvSquared);
  integrand_ -> setCurrentME(currentME);
  Integrand_ttHorZ_3l1tau::gIntegrand = integrand_;

//--- set integration boundaries
  double xl[9] = { -1., -pi(), -1., -pi(),  0., -pi() / 2, -pi(), -pi(), 0. };
  double xu[9] = { +1., +pi(), +1., +pi(), +1., +pi() / 2, +pi(), +pi(), tauLeptonMassSquared };
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
  std::copy(xl, xl + numDimensions_, xl_);
  std::copy(xu, xu + numDimensions_, xu_);

//--- create probability and corresponding error (uncertainty) variable
  double p = 0.;
  double pErr = 0;

//--- set integration settings
  const unsigned numCallsGridOpt = roundToNearestUInt(0.20 * maxObjFunctionCalls_);
  const unsigned numCallsIntEval = roundToNearestUInt(0.80 * maxObjFunctionCalls_);
  if(integrationMode_ == IntegrationMode::kVEGAS)
  {
    intAlgo_ = new MEMIntegratorVEGAS(numCallsGridOpt, numCallsIntEval, 1, 2.);
    intAlgo_ -> integrate(&g_C, xl_, xu_, numDimensions_, p, pErr);
    delete intAlgo_;
    intAlgo_ = 0;
  }
  else if(integrationMode_ == IntegrationMode::kVAMP)
  {
    intAlgo_ = new MEMIntegratorVAMP(numCallsGridOpt, numCallsIntEval);
    intAlgo_ -> integrate(&g_Fortran, xl_, xu_, numDimensions_, p, pErr);
    delete intAlgo_;
    intAlgo_ = 0;
  }

//--- divide by the process cross section (in GeV, i.e. natural units)
//--- and adjust the uncertainty accordingly
  p /= xSectionTTHinGeV2;
  pErr /= xSectionTTHinGeV2;
  LOGDBG << "p = " << p << "; pErr = " << pErr;

//--- cleanup
  if(xl_)
  {
    delete [] xl_;
    xl_ = 0;
  }
  if(xu_)
  {
    delete [] xu_;
    xu_ = 0;
  }

  clock_ -> Stop(__PRETTY_FUNCTION__);
  numSeconds_cpu_ = clock_ -> GetCpuTime(__PRETTY_FUNCTION__);
  numSeconds_real_ = clock_ -> GetRealTime(__PRETTY_FUNCTION__);
  LOGINFO << "Real time:\t" << std::round(numSeconds_real_ * 100) / 100 << " s;\t"
          << "CPU time:\t" << std::round(numSeconds_cpu_ * 100) / 100 << " s";

  numSecondsAccumul_cpu_ += numSeconds_cpu_;
  numSecondsAccumul_real_ += numSecondsAccumul_real_;
  ++nof_calls_;

  return p;
}

