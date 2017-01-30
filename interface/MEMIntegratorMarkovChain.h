#ifndef MEMINTEGRATORMARKOVCHAIN_H
#define MEMINTEGRATORMARKOVCHAIN_H

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // iStrComparator
#include "tthAnalysis/tthMEM/interface/MEMIntegratorBase.h" // MEMIntegratorBase
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*
#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line_ext()

#include <vector> // std::vector<>
#include <ostream> // std::ostream

#include <boost/bimap/bimap.hpp> // boost::bimaps::bimap<,>
#include <boost/bimap/set_of.hpp> // boost::bimaps::set_of<,>

#include <TRandom3.h> // TRandom3

namespace tthMEM
{
  enum class MarkovChainMode
  {
    kUniform = 0,
    kGaussian = 1,
    kNone = 2
  };

  /**
   * @brief Generic class to perform Markov Chain type integration in N-dimensional space.
   *
   * The code is implemented following the description in:
   *  [1] "Probabilistic Inference Using Markov Chain Monte Carlo Methods",
   *      R. Neal, http://www.cs.toronto.edu/pub/radford/review.pdf
   *  [2] "Bayesian Training of Backpropagation Networks by the Hybrid Monte Carlo Method",
   *      R. Neal, http://www.cs.toronto.edu/pub/radford/bbp.ps
   *
   * Originally created by Christian Veelken, NICPB Tallinn
   * (visit https://github.com/veelken/SVfitMEM)
   */
  class MEMIntegratorMarkovChain
    : public MEMIntegratorBase
  {
  public:
    MEMIntegratorMarkovChain(const std::string & modeStr,
                             unsigned nofIterBurnin,
                             unsigned nofIterSampling,
                             unsigned nofIterSimAnnPhase1,
                             unsigned nofIterSimAnnPhase2,
                             unsigned maxCallsStartingPos,
                             double T0,
                             double alpha,
                             unsigned nofChains,
                             unsigned nofBatches,
                             double epsilon0,
                             double nu);
    ~MEMIntegratorMarkovChain();

    void
    integrate(gPtr_C integrand,
              const double * xl,
              const double * xu,
              unsigned dimension,
              double & integral,
              double & integralErr) override;

    void
    integrate(gPtr_Fortran integrand __attribute__((unused)),
              const double * xl __attribute__((unused)),
              const double * xu __attribute__((unused)),
              unsigned dimension __attribute__((unused)),
              double & integral __attribute__((unused)),
              double & integralErr __attribute__((unused))) override
    {
      throw_line_ext("invalid usage", TTHEXCEPTION_ERR_CODE_UNDEFINED_FUNCTION)
        << "You must use integrate(gPtr_C, ...) not this one";
    }

    void *
    metadata() override;

    friend std::ostream &
    operator<<(std::ostream & os,
               const MEMIntegratorMarkovChain & mx);

  private:

    void
    setIntegrand(gPtr_C integrand,
                 const double * xl,
                 const double * xu,
                 unsigned dimension);

    void
    initialize();

    void
    makeStochasticMove(unsigned idxMove,
                       bool & isAccepted);

    void
    update_x(const std::vector<double> & q);

    double
    evalProb(const std::vector<double> & q);

    MarkovChainMode mode_;
    ///< flag indicating how initial position of Markov Chain (MX) is chosen
    const unsigned nofIterBurnin_;
    ///< number of ,,stochastic moves'' performed to reach ,,ergodic'' state of MX
    ///< the ,,burnin'' iterations do not enter the integral computation;
    ///< their purpose is to make the computed integral value somewhat independent
    ///< of the initial position at which the MX is started
    const unsigned nofIterSampling_;
    ///< number of ,,stochastic moves'' used to compute the integral
    const unsigned nofIterSimAnnPhase1_;
    ///< ,,simulated annealing'' phase: number of ,,stochastic moves'' performed
    ///<                                at high temperature during ,,burnin'' stage
    const unsigned nofIterSimAnnPhase2_;
    ///< ,,simulated annealing'' phase: number of ,,stochastic moves'' performed
    ///<                                at high temperature during ,,burnin'' stage
    const unsigned nofIterSimAnnPhaseSum_;
    ///< ,,simulated annealing'' phase: number of ,,stochastic moves'' performed
    ///<                                at high temperature during ,,burnin'' stage
    const unsigned maxCallsStartingPos_;
    ///< max number of attempts to find a valid starting position for the MX
    ///< aka initial position of non-zero probability
    const double T0_;
    ///< ,,simulated annealing'' phase: initial annealing temperature
    const double sqrtT0_;
    ///< ,,simulated annealing'' phase: square root of the initial annealing temperature
    const double alpha_;
    ///< ,,simulated annealing'' phase: speed at which the temperature decreases
    const double alphaSquared_;
    ///< ,,simulated annealing'' phase: squared speed at which the temperature decreses
    const unsigned nofChains_;
    ///< number of MXs run in parallel
    const unsigned nofBatches_;
    ///< number of batches per chain (used for uncertainty of computed integral value)
    const unsigned nofChainsXnofBatches_;
    ///< product of the nofChains_ and nofBatches_ (total number of batches)
    const unsigned nofIterPerBatch_;
    ///< number of iterations per batch
    const double epsilon0_;
    ///< average step size of Metropolis moves
    const double nu_;
    ///< variation of step-size for individual Metropolis moves
    std::vector<double> epsilon0s_;
    ///< vector of step-sizes for individual Metropolis moves

    gPtr_C integrand_;
    ///< integrand functional
    unsigned nofDimensions_;
    ///< dimensionality of the integration space
    std::vector<double> x_;
    ///< current point in the integration space
    std::vector<double> xMin_;
    ///< lower boundaries of integration region
    std::vector<double> xMax_;
    ///< upper boundaries of integration region

    /* internal variables storing current MX state */
    std::vector<double> p_, q_;
    double prob_;

    /* temporary variables used in computations */
    std::vector<double> pProposal_, qProposal_, probSum_, integral_;

    /* counter variables */
    unsigned long nofMoves_accepted_,      nofMoves_rejected_,
                  nofMoves_acceptedTotal_, nofMoves_rejectedTotal_,
                  nofIntegrationCalls_;
    unsigned nofTries_;
    ///< number of calls made until positive probability has been found (summed over chains)

    TRandom3 prng_;
    ///< PRNG (Mersenne Twister)

    static const boost::bimaps::bimap<
      boost::bimaps::set_of<std::string, iStrComparator>, MarkovChainMode
    > mxModeStrings_;
    ///< bimap where the keys in left presentation are case-insensitive strings
    ///< needed to initialize the mode_ variable in the constructor
  };
}

#endif // MEMINTEGRATORMARKOVCHAIN_H
