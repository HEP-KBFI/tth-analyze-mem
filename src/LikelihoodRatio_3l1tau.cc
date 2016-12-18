#include "tthAnalysis/tthMEM/interface/LikelihoodRatio_3l1tau.h"

#include <algorithm> // std::accumulate()
#include <numeric> // std::inner_product()

#include <boost/range/adaptor/map.hpp> // boost::adaptros::map_keys
#include <boost/range/algorithm/copy.hpp> // boost::copy()

#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // constants::, pow2()
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

using namespace tthMEM;

LikelihoodRatio_3l1tau::LikelihoodRatio_3l1tau(const std::vector<Hypothesis> & signalHypotheses,
                                               const std::vector<Hypothesis> & backgroundHypotheses)
  : signalHypotheses_(signalHypotheses)
  , backgroundHypotheses_(backgroundHypotheses)
{
  computeWeights();
}

LikelihoodRatio_3l1tau &
LikelihoodRatio_3l1tau::addSignalHypothesis(Hypothesis hypothesis,
                                            bool debug)
{
  if(debug)
  {
    if(std::find(
         signalHypotheses_.begin(),
         signalHypotheses_.end(),
         hypothesis
      ) != signalHypotheses_.end())
      throw_line("LikelihoodRatio_3l1tau")
        << "Signal hypothesis already added";
    if(std::find(
         backgroundHypotheses_.begin(),
         backgroundHypotheses_.end(),
         hypothesis
      ) != backgroundHypotheses_.end())
      throw_line("LikelihoodRatio_3l1tau")
        << "Signal hypothesis has already been added to list of background hypotheses";
  }
  signalHypotheses_.push_back(hypothesis);
  computeWeights();
  return *this;
}

LikelihoodRatio_3l1tau &
LikelihoodRatio_3l1tau::addBackgroundHypothesis(Hypothesis hypothesis,
                                                bool debug)
{
  if(debug)
  {
    if(std::find(
         backgroundHypotheses_.begin(),
         backgroundHypotheses_.end(),
         hypothesis
      ) != backgroundHypotheses_.end())
      throw_line("LikelihoodRatio_3l1tau")
        << "Background hypothesis already added";
    if(std::find(
         signalHypotheses_.begin(),
         signalHypotheses_.end(),
         hypothesis
      ) != signalHypotheses_.end())
      throw_line("LikelihoodRatio_3l1tau")
        << "Background hypothesis has already been added to list of signal hypotheses";
  }
  backgroundHypotheses_.push_back(hypothesis);
  computeWeights();
  return *this;
}

void
LikelihoodRatio_3l1tau::computeWeights()
{
  computeWeights(signalHypotheses_,     signalWeights);
  computeWeights(backgroundHypotheses_, backgroundWeights);
}

void
LikelihoodRatio_3l1tau::computeWeights(const std::vector<Hypothesis> & hypothesisList,
                                       std::vector<double> & weights)
{
  weights.clear();
  std::vector<double> weightNumerators;
  for(Hypothesis hypothesis: hypothesisList)
    switch(hypothesis)
    {
      case Hypothesis::tth:
        weightNumerators.push_back(constants::xSectionTTH); break;
      case Hypothesis::ttz:
        weightNumerators.push_back(constants::xSectionTTZ); break;
      case Hypothesis::tth_h2ww:
        weightNumerators.push_back(constants::xSectionTTH2diW); break;
    }
  const double weightDenominator = 1. / std::accumulate(weightNumerators.begin(), weightNumerators.end(), 0.);
  std::transform(
    weightNumerators.begin(), weightNumerators.end(), std::back_inserter(weights),
    [weightDenominator](double weightNumerator) -> double
    {
      return weightNumerator * weightDenominator;
    }
  );
  return;
}

std::vector<double>
LikelihoodRatio_3l1tau::getVector(const std::map<Hypothesis, std::array<double, 2> > & vectorMap,
                                  unsigned idx)
{
  if(idx > 2)
    throw_line("LikelihoodRatio_3l1tau") << "Invalid index requested: " << idx;
  std::vector<double> tmp;
  for(const auto & kv: vectorMap)
    tmp.push_back(kv.second[idx]);
  return tmp;
}

void
LikelihoodRatio_3l1tau::fillResults(MEMOutput_3l1tau & result,
                                    const std::map<Hypothesis, std::array<double, 2>> & vectorMap)
{
  for(const auto & kv: vectorMap)
    switch(kv.first)
    {
      case Hypothesis::tth:
        result.prob_tth     = kv.second[0];
        result.prob_tth_err = kv.second[1];
        break;
      case Hypothesis::ttz:
        result.prob_ttz = kv.second[0];
        result.prob_ttz = kv.second[1];
        break;
      case Hypothesis::tth_h2ww:
        result.prob_tth_h2ww     = kv.second[0];
        result.prob_tth_h2ww_err = kv.second[1];
        break;
    }
}

MEMOutput_3l1tau
LikelihoodRatio_3l1tau::compute(const std::map<Hypothesis, std::array<double, 2>> & signalResults,
                                const std::map<Hypothesis, std::array<double, 2>> & backgroundResults,
                                bool debug) const
{
  if(debug)
  {
    if(signalResults.size() != signalWeights.size())
      throw_line("LikelihoodRatio_3l1tau")
        << "Number of signal results doesn't match with the number of weights in the LR";
    if(backgroundResults.size() != backgroundWeights.size())
      throw_line("LikelihoodRatio_3l1tau")
        << "Number of background results doesn't match with the number of weights in the LR";

    const std::vector<Hypothesis> signalKeys = [&signalResults]()
    {
      std::vector<Hypothesis> tmp;
      boost::copy(signalResults | boost::adaptors::map_keys, std::back_inserter(tmp));
      return tmp;
    }();
    if(std::any_of(
         signalKeys.begin(), signalKeys.end(),
         [&backgroundResults](Hypothesis hypothesis) -> bool
         {
           return backgroundResults.count(hypothesis) > 0;
         }
      ))
      throw_line("LikelihoodRatio_3l1tau")
        << "Results of the given hypotheses are interleaving";
  }

  std::vector<double> signalProb        = getVector(signalResults, 0);
  std::vector<double> signalProbErr     = getVector(signalResults, 1);
  std::vector<double> backgroundProb    = getVector(backgroundResults, 0);
  std::vector<double> backgroundProbErr = getVector(backgroundResults, 1);

  auto dot = [](const std::vector<double> & lhs,
                const std::vector<double> & rhs) -> double
  {
    return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), 0.);
  };

  const double signalWeighted        = dot(signalProb,        signalWeights);
  const double signalErrWeighted     = dot(signalProbErr,     signalWeights);
  const double backgroundWeighted    = dot(backgroundProb,    backgroundWeights);
  const double backgroundErrWeighted = dot(backgroundProbErr, backgroundWeights);

  const double sig      = signalWeighted;
  const double sig_err  = signalErrWeighted;
  const double sig_up   = signalWeighted + signalErrWeighted;
  const double sig_down = signalWeighted - signalErrWeighted;
  const double bkg      = backgroundWeighted;
  const double bkg_err  = backgroundErrWeighted;
  const double bkg_up   = backgroundWeighted + backgroundErrWeighted;
  const double bkg_down = backgroundWeighted - backgroundErrWeighted;

  const double lr      = (sig      + bkg)      != 0. ? sig      / (sig      + bkg)      : 0.;
  const double lr_up   = (sig_up   + bkg_down) != 0. ? sig_up   / (sig_up   + bkg_down) : 0.;
  const double lr_down = (sig_down + bkg_up)   != 0. ? sig_down / (sig_down + bkg_up)   : 0.;
  const double lr_err  = (sig + bkg)           != 0. ?
                         std::sqrt(pow2(sig * sig_err) + pow2(bkg * bkg_err)) / pow2(sig + bkg) : 0.;

  MEMOutput_3l1tau result;
  fillResults(result, signalResults);
  fillResults(result, backgroundResults);
  result.sig      = sig;
  result.sig_err  = sig_err;
  result.sig_up   = sig_up;
  result.sig_down = sig_down;
  result.bkg      = bkg;
  result.bkg_err  = bkg_err;
  result.bkg_up   = bkg_up;
  result.bkg_down = bkg_down;
  result.lr       = lr;
  result.lr_err   = lr_err;
  result.lr_up    = lr_up;
  result.lr_down  = lr_down;

  return result;
}

