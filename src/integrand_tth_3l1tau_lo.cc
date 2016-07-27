#include "tthAnalysis/tthMEM/interface/integrand_tth_3l1tau_lo.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h"
#include "tthAnalysis/tthMEM/interface/Logger.h"

#include <iostream> // std::cerr, std::cout
#include <cmath> // std::sqrt()
#include <cstdlib> // std::exit(), EXIT_FAILURE

using namespace tthMEM;

integrand_tth_3l1tau_lo::integrand_tth_3l1tau_lo(double sqrtS,
                                                 const std::string & pdfName,
                                                 const std::string & madgraphFilename)
  : sqrtS_(sqrtS)
  , s_(std::sqrt(sqrtS_))
  , pdf_(0)
  , me_madgraph_initialized_(false)
{
  if(! pdf_ && pdfName != "") pdf_ = LHAPDF::mkPDF(pdfName.data(), 0);
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
}

integrand_tth_3l1tau_lo::~integrand_tth_3l1tau_lo()
{
  LOGINFO << "Calling the dtor";

  delete pdf_;
  pdf_ = 0;
}

double
integrand_tth_3l1tau_lo::eval() const
{
  if(! pdf_)
  {
    LOGERR << "PDF not initialized!";
    std::exit(EXIT_FAILURE);
  }

  if(! me_madgraph_initialized_)
  {
    LOGERR << "Madgraph's ME not initialized!" << "aa";
    std::exit(EXIT_FAILURE);
  }

  return 0.;
}

