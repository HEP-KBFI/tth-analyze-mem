#include "tthAnalysis/tthMEM/interface/integrators/MEMIntegratorVAMP.h"

using namespace tthMEM;

extern "C"
{
  void vamp_integrate_(MEMIntegratorBase::gPtr_Fortran integrand,
                       const double * xl,
                       const double * xu,
                       int * numDimensions,
                       int * numCallsGridOpt,
                       int * numCallsIntEval,
                       double * integral,
                       double * integralErr);
}

MEMIntegratorVAMP::MEMIntegratorVAMP(unsigned numCallsGridOpt,
                                     unsigned numCallsIntEval)
  : integrand_(nullptr)
  , numCallsGridOpt_(numCallsGridOpt)
  , numCallsIntEval_(numCallsIntEval)
  , numDimensions_(0)
  , xl_(nullptr)
  , xu_(nullptr)
{}

MEMIntegratorVAMP::~MEMIntegratorVAMP()
{}

void
MEMIntegratorVAMP::setIntegrand(MEMIntegratorBase::gPtr_Fortran integrand,
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
}

void
MEMIntegratorVAMP::integrate(MEMIntegratorBase::gPtr_Fortran integrand,
                             const double * xl,
                             const double * xu,
                             unsigned dimension,
                             double & integral,
                             double & integralErr)
{
  setIntegrand(integrand, xl, xu, dimension);
  if(! integrand_)
    throw_line_ext("invalid arument", TTHEXCEPTION_ERR_CODE_MISSING_INTEGRAND)
      << "No integrand function has been set";
  int numDimensions = numDimensions_;
  int numCallsGridOpt = numCallsGridOpt_;
  int numCallsIntEval = numCallsIntEval_;

  vamp_integrate_(integrand_, xl, xu, &numDimensions, &numCallsGridOpt,
                  &numCallsIntEval, &integral, &integralErr);

  if(xl_)
  {
    delete [] xl_;
    xl_ = nullptr;
  }
  if(xu_)
  {
    delete [] xu_;
    xu_ = nullptr;
  }
}

void *
MEMIntegratorVAMP::metadata()
{
  return nullptr;
}
