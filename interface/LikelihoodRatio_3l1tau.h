#ifndef LIKELIHOODRATIO_3L1TAU_H
#define LIKELIHOODRATIO_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/MEMOutput_3l1tau.h" // MEMOutput_3l1tau

#include <array> // std::array<,>
#include <vector> // std::vector<>
#include <map> // std::map<>

namespace tthMEM
{
  class LikelihoodRatio_3l1tau
  {
public:
    enum class Hypothesis
    {
      tth      = 0,
      ttz      = 1,
      tth_h2ww = 2
    };

    LikelihoodRatio_3l1tau() = default;
    LikelihoodRatio_3l1tau(const std::vector<Hypothesis> & signalHypotheses,
                           const std::vector<Hypothesis> & backgroundHypotheses);

    LikelihoodRatio_3l1tau &
    addSignalHypothesis(Hypothesis hypothesis,
                        bool debug = true);

    LikelihoodRatio_3l1tau &
    addBackgroundHypothesis(Hypothesis hypothesis,
                            bool debug = true);

    MEMOutput_3l1tau
    compute(const std::map<Hypothesis, std::array<double, 2>> & signalResults,
            const std::map<Hypothesis, std::array<double, 2>> & backgroundResults,
            bool debug = true) const;

private:
    void
    computeWeights();

    void
    computeWeights(const std::vector<Hypothesis> & hypothesisList,
                   std::vector<double> & weights);

    static std::vector<double>
    getVector(const std::map<Hypothesis, std::array<double, 2>> & vectorMap,
              unsigned idx);

    static void
    fillResults(MEMOutput_3l1tau & result,
                const std::map<Hypothesis, std::array<double, 2>> & vectorMap);

    std::vector<Hypothesis> signalHypotheses_;
    std::vector<Hypothesis> backgroundHypotheses_;

    std::vector<double> signalWeights;
    std::vector<double> backgroundWeights;
  };
}

#endif // LIKELIHOODRATIO_3L1TAU_H
