#include "tthAnalysis/tthMEM/interface/MEMIntegratorVEGAS.h"

using namespace tthMEM;

MEMIntegratorVEGAS::MEMIntegratorVEGAS(unsigned numCallsGridOpt,
                                       unsigned numCallsIntEval,
                                       unsigned maxIntEvalIter,
                                       double maxChi2)
  : integrand_(0)
  , numCallsGridOpt_(numCallsGridOpt)
  , numCallsIntEval_(numCallsIntEval)
  , maxIntEvalIter_(maxIntEvalIter)
  , maxChi2_(maxChi2)
  , precision_(1.e-5)
  , xl_(0)
  , xu_(0)
{}

MEMIntegratorVEGAS::~MEMIntegratorVEGAS()
{}

void
MEMIntegratorVEGAS::setIntegrand(MEMIntegratorBase::gPtr_C integrand,
                                 const double * xl,
                                 const double * xu,
                                 unsigned dimension)
{
  numDimensions_ = dimension;
  xl_ = new double[numDimensions_];
  xu_ = new double[numDimensions_];
  for(unsigned iDim = 0; iDim < numDimensions_; ++iDim)
  {
    xl_[iDim] = xl[iDim];
    xu_[iDim] = xu[iDim];
  }

  integrand_ = integrand;

  vegasIntegrand_ = new gsl_monte_function;
  vegasIntegrand_ -> f = integrand_;
  vegasIntegrand_ -> dim = numDimensions_;
  vegasIntegrand_ -> params = new double[1];

  vegasWorkspace_ = gsl_monte_vegas_alloc(numDimensions_);
  gsl_rng_env_setup();
  vegasRnd_ = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(vegasRnd_, 12345);

  gsl_monte_vegas_init(vegasWorkspace_);
  vegasWorkspace_ -> stage = 0;
  double integral = 0;
  double integralErr = 0;
  gsl_monte_vegas_integrate(vegasIntegrand_, xl_, xu_, numDimensions_,
                            numCallsGridOpt_ / vegasWorkspace_ -> iterations,
                            vegasRnd_, vegasWorkspace_, &integral, &integralErr);
  vegasWorkspace_ -> stage = 1;
}

void
MEMIntegratorVEGAS::integrate(MEMIntegratorBase::gPtr_C integrand,
                              const double * xl,
                              const double *xu,
                              unsigned dimension,
                              double &integral,
                              double &integralErr)
{
  setIntegrand(integrand, xl, xu, dimension);
  if(! integrand_)
  {
    LOGERR << "No integrand function has been set";
    assert(0);
  }

  integral = 0;
  integralErr = 0;

  unsigned iteration = 0;
  double chi2 = -1;
  do
  {
    gsl_monte_vegas_integrate(
      vegasIntegrand_, xl_, xu_, numDimensions_,
      numCallsIntEval_ / vegasWorkspace_ -> iterations,
      vegasRnd_, vegasWorkspace_, &integral, &integralErr);
    vegasWorkspace_ -> stage = 3;
    ++iteration;
    chi2 = vegasWorkspace_ -> chisq;
  }
  while((chi2 > maxChi2_ || integralErr > (integral * precision_)) &&
        iteration < maxIntEvalIter_);

  if(xl_) delete [] xl_;
  if(xu_) delete [] xu_;
  if(vegasIntegrand_ -> params)
    delete [] static_cast<double *>(vegasIntegrand_ -> params);
  if(vegasIntegrand_) delete vegasIntegrand_;
  gsl_monte_vegas_free(vegasWorkspace_);
  gsl_rng_free(vegasRnd_);
}

