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
setInputs(const tthMEM::MeasuredEvent_3l1tau & measuredEvent)
{
  //
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

  return 0.;
}

