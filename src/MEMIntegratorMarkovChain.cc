#include "tthAnalysis/tthMEM/interface/MEMIntegratorMarkovChain.h"
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // pow*()
#include "tthAnalysis/tthMEM/interface/tthMEMvecFunctions.h" // tthMEM::vec::

#include <cmath> // std::sqrt(), std::exp()
#include <algorithm> // std::copy(), std::fill_n(), std::generate(), ...
  // ..., std::transform(), std::find_if_not(), std::generate_n()
#include <limits> // std::numeric_limist<>
#include <iterator> // std::back_inserter()

#include <boost/assign/list_of.hpp> // boost::assign::list_of<>

#include <TMath.h> // TMath::IsNaN(), TMath::Finite()

using namespace tthMEM;

const decltype(MEMIntegratorMarkovChain::mxModeStrings_)
  MEMIntegratorMarkovChain::mxModeStrings_ =
  boost::assign::list_of<decltype(mxModeStrings_)::relation>
  ( "uniform",  MarkovChainMode::kUniform  )
  ( "gaussian", MarkovChainMode::kGaussian )
  ( "none",     MarkovChainMode::kNone     );

MEMIntegratorMarkovChain::MEMIntegratorMarkovChain(const std::string & modeStr,
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
                                                   double nu)
  : nofIterBurnin_(nofIterBurnin)
  , nofIterSampling_(nofIterSampling)
  , nofIterSimAnnPhase1_(nofIterSimAnnPhase1)
  , nofIterSimAnnPhase2_(nofIterSimAnnPhase2)
  , nofIterSimAnnPhaseSum_(nofIterSimAnnPhase1_ + nofIterSimAnnPhase2_)
  , maxCallsStartingPos_(maxCallsStartingPos)
  , T0_(T0)
  , sqrtT0_(T0 > 0. ? std::sqrt(T0) : 0.)
  , alpha_(alpha)
  , alphaSquared_(pow2(alpha_))
  , nofChains_(nofChains)
  , nofBatches_(nofBatches)
  , nofChainsXnofBatches_(nofChains_ * nofBatches_)
  , nofIterPerBatch_(nofBatches_ > 0 ? nofIterSampling_ / nofBatches_ : 0)
  , epsilon0_(epsilon0)
  , nu_(nu)
{
  const decltype(mxModeStrings_)::left_const_iterator it =
    mxModeStrings_.left.find(modeStr);
  if(it == mxModeStrings_.left.end())
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_MODE)
      << "No such Markov Chain mode: '" << modeStr << '\'';
  mode_ = it -> second;
  if(nofIterSimAnnPhaseSum_ > nofIterBurnin_)
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_ANNSUM)
      << "Invalid configuration parameters: "
      << "'nofIterSimAnnPhase1' = " << nofIterSimAnnPhase1_ << "; "
      << "'nofIterSimAnnPhase2' = " << nofIterSimAnnPhase2_ << "; "
      << "'nofIterBurnin' = " << nofIterBurnin_ << '\n'
      << " => Annealing and sampling stages must not overlap, i.e. "
      << "the number of ,,burnin'' iterations must be less than "
      << "the sum of annealing and sampling iteration numbers";
  if(! (alpha_ > 0. && alpha_ < 1.))
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_ALPHA)
      << "Invalid configuration parameter: "
      << "'alpha' = " << alpha_ << " must be in (0, 1)";
  if(! nofChains_)
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_NOF_CHAINS)
      << "Invalid configuration parameter: "
      << "'nofChains' = " << nofChains_ << " must be above zero";
  if(! nofBatches_)
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_NOF_BATCHES)
      << "Invalid configuration parameter: "
      << "'nofBatches' = " << nofBatches_ << " must be above zero";
  if(nofIterSampling_ % nofBatches_ != 0)
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_BATCH_DIV)
      << "Invalid configuration parameters: "
      << "'nofIterSampling' = " << nofIterSampling_ << "; "
      << "'nofBatches' = " << nofBatches_ << '\n'
      << " => 'nofIterSampling' must be divisible by 'nofBatches'";
  if(T0 <= 0.)
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_TEMP)
      << "Invalid configuration parameter: "
      << "'T0' = " << T0 << " must be above zero";
  /* not sure though but place the constrain anyways */
  if(! (nu_ > 0. && nu_ < 1.))
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_INVALID_NU)
      << "Invalid configuration parameter: "
      << "'nu' = " << nu_ << " must be in (0, 1)";
}

MEMIntegratorMarkovChain::~MEMIntegratorMarkovChain()
{
  LOGVRB_S << "integration calls = " << nofIntegrationCalls_ << "; "
           << "accepted moves = " << nofMoves_acceptedTotal_ << "; "
           << "rejected moves = " << nofMoves_rejectedTotal_ << ' '
           << "(acceptance rate = "
           << 100. * nofMoves_acceptedTotal_ /
              (nofMoves_acceptedTotal_ + nofMoves_rejectedTotal_) << "%)";
}

void
MEMIntegratorMarkovChain::setIntegrand(gPtr_C integrand,
                                       const double * xl,
                                       const double * xu,
                                       unsigned dimension)
{
//--- reset current dimensionality
  nofDimensions_ = dimension;
//--- clean old stuff (quick n dirty)
  for(auto & vec: { &x_, &xMin_, &xMax_, &epsilon0s_, &p_, &q_,
                    &pProposal_, &qProposal_, &probSum_, &integral_ })
    (*vec).clear();
//--- copy the integration space boundaries over
  std::copy(xl, xl + nofDimensions_, std::back_inserter(xMin_));
  std::copy(xu, xu + nofDimensions_, std::back_inserter(xMax_));
  std::fill_n(std::back_inserter(epsilon0s_), nofDimensions_, epsilon0_);
//--- 1st N entries ,,significant'' components, last N ,,dummy'' components
  std::fill_n(std::back_inserter(p_), 2 * nofDimensions_, 0.);
//--- reset other vectors as well (quick n dirty)
//--- ,,potential energy'' E(q) depends only on the 1st N ,,significant components''
  for(auto & vec: { &x_, &q_, &pProposal_, &qProposal_ })
    std::fill_n(std::back_inserter(*vec), nofDimensions_, 0.);
//--- reset the current probability and integral
  prob_ = 0.;
  std::fill_n(std::back_inserter(probSum_), nofChainsXnofBatches_, 0.);
  std::fill_n(std::back_inserter(integral_), nofChainsXnofBatches_, 0.);
//--- actually set the integrand
  integrand_ = integrand;
//--- printout
  LOGVRB << "Markov Chain mode = "     << mxModeStrings_.right.at(mode_);
  LOGVRB << "nofIterBurnin = "         << nofIterBurnin_;
  LOGVRB << "nofIterSampling = "       << nofIterSampling_;
  LOGVRB << "nofIterSimAnnPhase1 = "   << nofIterSimAnnPhase1_;
  LOGVRB << "nofIterSimAnnPhase2 = "   << nofIterSimAnnPhase2_;
  LOGVRB << "nofIterSimAnnPhaseSum = " << nofIterSimAnnPhaseSum_;
  LOGVRB << "maxCallsStartingPos = "   << maxCallsStartingPos_;
  LOGVRB << "nofChains = "             << nofChains_;
  LOGVRB << "nofBatches = "            << nofBatches_;
  LOGVRB << "nofIterPerBatch = "       << nofIterPerBatch_;
  LOGVRB << "T0 = "                    << T0_;
  LOGVRB << "alpha = "                 << alpha_;
  LOGVRB << "epsilon0 = "              << epsilon0_;
  LOGVRB << "nu = "                    << nu_;
}

void
MEMIntegratorMarkovChain::initialize()
{
//--- randomly choose starting position of MX in the N-dimensional space
  q_.clear();
  std::generate_n(std::back_inserter(q_), nofDimensions_,
    [this]() -> double
    {
      while(true)
      {
        const double q0 = mode_ == MarkovChainMode::kGaussian ?
                          prng_.Gaus(0.5, 0.5) : prng_.Uniform(0., 1.);
        if(q0 > 0. && q0 < 1.) return q0;
      }
    }
  );
//--- printout
  LOGTRC_S << "q = " << q_;
}

void
MEMIntegratorMarkovChain::integrate(gPtr_C integrand,
                                    const double * xl,
                                    const double * xu,
                                    unsigned dimension,
                                    double & integral,
                                    double & integralErr)
{
  setIntegrand(integrand, xl, xu, dimension);
  if(! integrand_)
    throw_line_ext("invalid argument", TTHEXCEPTION_ERR_CODE_MXMC_MISSING_INTEGRAND)
      << "No integrand functional has been set";

//--- set PRNG seed to the same number for each integration so that the integration
//--- retsults do not depend on the pervious integration history
//--- also good for reproducability
  prng_.SetSeed(12345);

//--- initialize counters
  nofMoves_accepted_ = 0;
  nofMoves_rejected_ = 0;

//--- loop over MX
  nofTries_ = 0;
  for(unsigned iChain = 0; iChain < nofChains_; ++iChain)
  {
    bool isValidStartPos = false;
    if(mode_ == MarkovChainMode::kNone)
    {
      prob_ = evalProb(q_);
      if(prob_ > 0.)
      {
        const bool isWithinBounds = std::find_if_not(
          q_.begin(), q_.end(),
          [](double qi) -> bool
          {
            return (qi > 0. && qi < 1.);
          }
        ) != q_.end();

        if(isWithinBounds)
          isValidStartPos = true;
        else
          LOGVRB_S << "Requested start position = " << q_ << ' '
                   << "not within interval (0, 1) => "
                   << "searching for a valid alternative";
      } // prob_ > 0.
      else
        LOGVRB_S << "Requested start position = " << q_ << ' '
                 << "returned zero probability => "
                 << "searching for a valid alternative";
    } // mode_ == MarkovChainMode::kNone

    unsigned nofTries = 0;
    while(! isValidStartPos && nofTries < maxCallsStartingPos_)
    {
      initialize();
      prob_ = evalProb(q_);
      if(prob_ > 0.)
        isValidStartPos = true;
      else if(nofTries > 0 && nofTries % 100000 == 0)
        LOGTRC << "attempt #" << iChain << " did not find a valid "
               << "start position, yet";
      ++nofTries;
    } // ! isValidStartPos && nofTries < maxCallsStartingPos_
    nofTries_ += nofTries;
    if(! isValidStartPos) continue;

//--- propose MX transition to a new, randomly chosen point
    for(unsigned iMove = 0; iMove < nofIterBurnin_; ++iMove)
    {
      bool isAccepted = false; // never used!
      makeStochasticMove(iMove, isAccepted);
    }

//--- propose MX transition to a new, randomly chosen point
    unsigned iBatch = iChain * nofBatches_;
    for(unsigned iMove = 0; iMove < nofIterSampling_; ++iMove)
    {
      bool isAccepted = false;
      makeStochasticMove(nofIterBurnin_ + iMove, isAccepted);
      if(isAccepted) ++nofMoves_accepted_;
      else           ++nofMoves_rejected_;

      update_x(q_);

      if(iMove > 0 && iMove % nofIterPerBatch_ == 0)
        ++iBatch;
      if(iBatch >= probSum_.size())
        throw_line_ext("runtime error", TTHEXCEPTION_ERR_CODE_MXMC_RUNTIME)
          << "Something's off: "
          << "'iBatch' = " << iBatch << " >= 'probSum_' size = " << probSum_.size();
      probSum_[iBatch] += prob_;
    } // nofIterSampling
  } // nofChains

  integral_ = probSum_ / nofIterPerBatch_;
  LOGTRC_S << "integral = " << integral_;

//--- compute integral value and its uncertainty ((6.39) and (6.40) in [1])
  integral = vec::avg(integral_);
  integralErr = vec::stdev(integral_, integral);
  LOGVRB_S << "integral = " << integral << " +/- " << integralErr;

  ++nofIntegrationCalls_;
  nofMoves_acceptedTotal_ += nofMoves_accepted_;
  nofMoves_rejectedTotal_ += nofMoves_rejected_;

  LOGTRC << *this;
}

void
MEMIntegratorMarkovChain::makeStochasticMove(unsigned idxMove,
                                             bool & isAccepted)
{
//--- performs ,,stochastic move'' ((24) in [2])
  LOGTRC << "idxMove = " << idxMove;

//--- perform random updates of momentum components
  const std::function<double()> gaus = [this]() -> double
    {
      return prng_.Gaus(0., 1.);
    };
  if     (idxMove < nofIterSimAnnPhase1_)
    std::generate(p_.begin(), p_.end(), gaus);
  else if(idxMove < nofIterSimAnnPhaseSum_)
  {
//--- sample random numbers spherically (?)
    std::vector<double> u = vec::genv(gaus, p_.size());
    const double uL2 = vec::l2(u);
//--- divide by the magnitude of the vector, uL2
    u /= uL2;
    const double pL2 = vec::l2(p_);
    p_ = alpha_ * pL2 * u + (1. - alphaSquared_) * vec::genv(gaus, p_.size());
  }
  else
    std::generate(p_.begin(), p_.end(), gaus);

//--- choose random stpe size
  double expNuTimesC = 0.;
  do
  {
    expNuTimesC = std::exp(nu_ * prng_.BreitWigner(0., 1.));
  } while(TMath::IsNaN(expNuTimesC) ||
          ! TMath::Finite(expNuTimesC) ||
          expNuTimesC > 1.e+6);
  const std::vector<double> epsilons = expNuTimesC * epsilon0s_;

//--- Metropolis algorithm move ((27) in [2])

//--- update position components by single step of chosen size
//--- in the direction of the momentum components
  qProposal_ = q_ + epsilons * vec::subv(p_, nofDimensions_);

//--- ensure that proposed new point is within integration region
//--- (assume ,,periodic'' integration domain)
  for(unsigned iDim = 0; iDim < qProposal_.size(); ++iDim)
  {
    const double qi = qProposal_[iDim] - std::floor(qProposal_[iDim]);
    if(! (qi >= 0. && qi <= 1.))
      throw_line_ext("runtime error", TTHEXCEPTION_ERR_CODE_MXMC_RUNTIME)
        << "Encountered position component q[" << iDim << "] = " << qi << ", "
        << "which is not in [0, 1] => bailing out";
    qProposal_[iDim] = qi;
  }
  LOGTRC << "epsilons = " << epsilons;
  LOGTRC << "qProposal = " << qProposal_;
  LOGTRC << "q = " << q_;
  LOGTRC << "p = " << p_;

//--- check if proposed move of MX to a new position is accepted or not: compute change
//--- in the phase space volume for ,,dummy'' momentum components ((25) in [2])
  const double probProposal = evalProb(qProposal_);
  double deltaEnergy = 0.;
  if     (probProposal > 0. && prob_ > 0.)
    deltaEnergy = -std::log(probProposal / prob_);
  else if(probProposal > 0.)
    deltaEnergy = -std::numeric_limits<double>::max();
  else if(                     prob_ > 0.)
    deltaEnergy = +std::numeric_limits<double>::max();
  else
    throw_line_ext("runtime error", TTHEXCEPTION_ERR_CODE_MXMC_RUNTIME)
      << "Encountered negative change in the phase space voulme: "
      << "'probProposal' = " << probProposal << "; "
      << "'prob_' = " << prob_;

//--- Metropolis move: accept the move if generated number drawn
//--- from uniform distribution is less than the probability
  const double pAccept = std::exp(-deltaEnergy);
  const double u = prng_.Uniform(0., 1.);
  if(u < pAccept)
  {
    q_ = qProposal_;
    prob_ = probProposal;
    isAccepted = true;
  }
  else
    isAccepted = false;
}

void
MEMIntegratorMarkovChain::update_x(const std::vector<double> & q)
{
  x_ = (1. - q) * xMin_ + q * xMax_;
}

double
MEMIntegratorMarkovChain::evalProb(const std::vector<double> & q)
{
  update_x(q);
  return (*integrand_)(x_.data(), nofDimensions_, 0);
}

void *
MEMIntegratorMarkovChain::metadata()
{
  return &nofTries_;
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MEMIntegratorMarkovChain & MX)
  {
    os << "Moves: accepted = " << MX.nofMoves_accepted_ << "; "
       << "rejected = " << MX.nofMoves_rejected_ << ' '
       << "(acceptance rate " << 100. * MX.nofMoves_accepted_ /
                                 (MX.nofMoves_accepted_ + MX.nofMoves_rejected_)
       << " %)\n";

    for(unsigned i = 0; i < MX.nofChains_; ++i)
    {
//--- find the average per batch
      const double integral = vec::avg(
        MX.integral_, MX.nofBatches_ * i, MX.nofBatches_ * (i + 1)
      );
//--- find the standard deviation per batch
      const double integralErr = vec::stdev(
        MX.integral_, MX.nofBatches_ * i, MX.nofBatches_ * (i + 1), integral
      );
      os << std::string(10, ' ') << "chain #" << i << ": "
         << "integral = " << integral << "; "
         << "error = " << integralErr << '\n';
    }

    return os;
  }
}
