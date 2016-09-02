#include "tthAnalysis/tthMEM/interface/MEMIntegratorMarkovChain.h"
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*, Logger::
#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // pow*(), functions::

#include <cmath> // std::sqrt(), std::exp()
#include <stdexcept> // std::invalid_argument, std::runtime_error
#include <algorithm> // std::copy(), std::fill(), std::generate(), ...
  // ..., std::transform(), find_if_not()
#include <functional> // std::divides<>, std::bind2nd()
#include <limits> // std::numeric_limist<>

#include <TMath.h> // TMath::IsNaN(), TMath::Finite()

using namespace tthMEM;

MEMIntegratorMarkovChain::MEMIntegratorMarkovChain(MarkovChainMode mode,
                                                   unsigned nofIterBurnin,
                                                   unsigned nofIterSampling,
                                                   unsigned nofIterSimAnnPhase1,
                                                   unsigned nofIterSimAnnPhase2,
                                                   double T0,
                                                   double alpha,
                                                   unsigned nofChains,
                                                   unsigned nofBatches,
                                                   double epsilon0,
                                                   double nu)
  : mode_(mode)
  , nofIterBurnin_(nofIterBurnin)
  , nofIterSampling_(nofIterSampling)
  , nofIterSimAnnPhase1_(nofIterSimAnnPhase1)
  , nofIterSimAnnPhase2_(nofIterSimAnnPhase2)
  , nofIterSimAnnPhaseSum_(nofIterSimAnnPhase1_ + nofIterSimAnnPhase2_)
  , T0_(T0)
  , sqrtT0_(std::sqrt(T0))
  , alpha_(alpha)
  , alphaSquared_(pow2(alpha_))
  , nofChains_(nofChains)
  , nofBatches_(nofBatches)
  , nofChainsXnofBatches_(nofChains_ * nofBatches_)
  , epsilon0_(epsilon0)
  , nu_(nu)
  , maxCallsStartingPos_(1000000)
  , x_(0)
{
  if(nofIterSimAnnPhaseSum_ > nofIterBurnin_)
  {
    LOGERR << "Invalid configuration parameters: "
           << "'nofIterSimAnnPhase1' = " << nofIterSimAnnPhase1_ << "; "
           << "'nofIterSimAnnPhase2' = " << nofIterSimAnnPhase2_ << "; "
           << "'nofIterBurnin' = " << nofIterBurnin_;
    LOGERR << " => Annealing and sampling stages must not overlap, i.e. "
           << "the number of ,,burnin'' iterations must be less than "
           << "the sum of annealing and sampling iteration numbers";
    throw std::invalid_argument(__PRETTY_FUNCTION__);
  }
  if(! (alpha_ > 0. && alpha_ < 1.))
  {
    LOGERR << "Invalid configuration parameter: "
           << "'alpha' = " << alpha_ << " "
           << "must be in (0, 1)";
    throw std::invalid_argument(__PRETTY_FUNCTION__);
  }
  if(! nofChains_)
  {
    LOGERR << "Invalid configuration parameter: "
           << "'nofChains' = " << nofChains_ << " "
           << "must be above zero";
    throw std::invalid_argument(__PRETTY_FUNCTION__);
  }
  if(! nofBatches_)
  {
    LOGERR << "Invalid configuration parameter: "
           << "'nofBatches' = " << nofBatches_ << " "
           << "must be above zero";
    throw std::invalid_argument(__PRETTY_FUNCTION__);
  }
  if(nofIterSampling_ % nofBatches_ != 0)
  {
    LOGERR << "Invalid configuration parameters: "
           << "'nofIterSampling' = " << nofIterSampling_ << "; "
           << "'nofBatches' = " << nofBatches_;
    LOGERR << " => 'nofIterSampling' must be divisible by 'nofBatches'";
    throw std::invalid_argument(__PRETTY_FUNCTION__);
  }
}

MEMIntegratorMarkovChain::~MEMIntegratorMarkovChain()
{
  LOGVRB_S << "integration calls = " << nofIntegrationCalls_ << "; "
           << "accepted moves = " << nofMoves_acceptedTotal_ << "; "
           << "rejected moves = " << nofMoves_rejectedTotal_ << " "
           << "(acceptance rate = "
           << static_cast<double>(nofMoves_acceptedTotal_) /
              (nofMoves_acceptedTotal_ + nofMoves_rejectedTotal_) << ")";
  if(x_)
  {
    delete [] x_;
    x_ = 0;
  }
}

void
MEMIntegratorMarkovChain::setIntegrand(gPtr_C integrand,
                                       const double * xl,
                                       const double * xu,
                                       unsigned dimension)
{
//--- reset current dimensionality
  nofDimensions_ = dimension;
//--- clean old stuff
  if(x_) delete [] x_;
  xMin_.clear();
  xMax_.clear();
//--- ...
  x_ = new double[nofDimensions_];
//--- copy the integration space boundaries over
  xMin_.resize(nofDimensions_);
  xMax_.resize(nofDimensions_);
  std::copy(xl, xl + nofDimensions_, std::back_inserter(xMin_));
  std::copy(xu, xu + nofDimensions_, std::back_inserter(xMax_));
  std::fill(epsilon0s_.begin(), epsilon0s_.end(), epsilon0_);
//--- 1st N entries ,,significant'' components, last N ,,dummy'' components
  p_.resize(2 * nofDimensions_);
  u_.resize(2 * nofDimensions_);
//--- ,,potential energy'' E(q) depends on only the 1st N ,,significant components''
  q_.resize(nofDimensions_);
  pProposal_.resize(nofDimensions_);
  qProposal_.resize(nofDimensions_);
//--- reset the current probability
  prob_ = 0.;
  probSum_.resize(nofChainsXnofBatches_);
  std::fill(probSum_.begin(), probSum_.end(), 0.);
//--- ...
  integral_.resize(nofChains_ * nofBatches_);
//--- actually set the integrand
  integrand_ = integrand;
}

void
MEMIntegratorMarkovChain::initialize()
{
//--- randomly choose starting position of MX in the N-dimensional space
  q_.clear();
  q_.reserve(nofDimensions_);
  std::generate(
    q_.begin(), q_.end(),
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
  {
    LOGERR << "No integrand functional has been set";
    throw std::invalid_argument(__PRETTY_FUNCTION__);
  }

//--- set PRNG seed to the same number for each integration so that the integration
//--- retsults do not depend on the pervious integration history
//--- also good for reproducability
  prng_.SetSeed(12345);

//--- initialize counters
  nofMoves_accepted_ = 0;
  nofMoves_rejected_ = 0;
  nofChainsRun_ = 0;

//--- TTree stuff (tbd)

//--- loop over MX
  const unsigned nofIterPerBatch = nofBatches_ / nofIterSampling_;
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
          LOGVRB_S << "Requested start position = " << q_ << " "
                   << "not within interval (0, 1) => "
                   << "searching for a valid alternative";
      } // prob_ > 0.
      else
        LOGVRB_S << "Requested start position = " << q_ << " "
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

//--- TTree stuff (tbd)

      if(iMove > 0 && iMove % nofIterPerBatch == 0) ++iBatch;
      if(iBatch >= probSum_.size())
      {
        LOGERR << "Something's off: 'iBatch' = " << iBatch << " >= "
               << "'probSum_' size = " << probSum_.size();
        throw std::runtime_error(__PRETTY_FUNCTION__);
      }
      probSum_[iBatch] += prob_;
    } // nofIterSampling

    ++nofChainsRun_;
  } // nofChains

  std::transform(probSum_.begin(), probSum_.end(), integral_.begin(),
                 std::bind2nd(std::divides<double>(), nofIterPerBatch));
  LOGTRC_S << "integral = " << integral_;

//--- compute integral value and its uncertainty
  integral = functions::avg(integral_);
  integralErr = functions::stdev(integral_, integral);
  LOGVRB_S << "integral = " << integral << " +/- " << integralErr;

  ++nofIntegrationCalls_;
  nofMoves_acceptedTotal_ += nofMoves_accepted_;
  nofMoves_rejectedTotal_ += nofMoves_rejected_;

//--- TTree stuff (tbd)
}

void
MEMIntegratorMarkovChain::makeStochasticMove(unsigned idxMove,
                                             bool & isAccepted)
{
//--- performs ,,stochastic move''

//--- perform random updates of momentum components
  if     (idxMove < nofIterSimAnnPhase1_)
    std::generate(p_.begin(), p_.end(),
                  [this]() -> double
                  {
                    return sqrtT0_ * prng_.Gaus(0., 1.);
                  });
  else if(idxMove < nofIterSimAnnPhaseSum_)
  {
//--- sample random numbers spherically (?)
    std::generate(u_.begin(), u_.end(),
                  [this]() -> double
                  {
                    return prng_.Gaus(0., 1.);
                  });
    const double uL2 = functions::l2(u_);
//--- divide by the magnitude of the vector, uL2
    std::transform(u_.begin(), u_.end(), u_.begin(),
                   std::bind2nd(std::divides<double>(), uL2));
    const double pL2 = functions::l2(p_);
    for(unsigned i = 0; i < p_.size(); ++i)
      p_[i] = alpha_ * pL2 * u_[i] + (1. - alphaSquared_) * prng_.Gaus(0., 1.);
  }
  else
    std::generate(p_.begin(), p_.end(),
                  [this]() -> double
                  {
                    return prng_.Gaus(0., 1.);
                  });

//--- choose random stpe size
  double expNuTimesC = 0.;
  do
  {
    expNuTimesC = std::exp(nu_ * prng_.BreitWigner(0., 1.));
  } while(! TMath::IsNaN(expNuTimesC) ||
          ! TMath::Finite(expNuTimesC) ||
          expNuTimesC > 1.e+6);
  std::vector<double> epsilons;
  std::transform(epsilon0s_.begin(), epsilon0s_.end(), std::back_inserter(epsilons),
                 [expNuTimesC](double epsilon0) -> double
                 {
                   return expNuTimesC * epsilon0;
                 }
  );

//--- Metropolis algorithm move

//--- update position components by single step of chosen size
//--- in the direction of the momentum components
  for(unsigned i = 0; i < qProposal_.size(); ++i)
    qProposal_[i] = q_[i] + epsilons[i] * p_[i];

//--- ensure that proposed new point is within integration region
//--- (assume ,,periodic'' integration domain)
  for(unsigned i = 0; i < qProposal_.size(); ++i)
  {
    const double qi = qProposal_[i] - std::floor(qProposal_[i]);
    if(! (qi >= 0. && qi <= 1.))
    {
      LOGERR << "Encountered position component q[" << i << "] = " << qi << ", "
             << "which is not in [0, 1] => bailing out";
      throw std::runtime_error(__PRETTY_FUNCTION__);
    }
    qProposal_[i] = qi;
  }

//--- check if proposed move of MX to a new position is accepted or not:
//--- compute change in the phase space volume for ,,dummy'' momentum components
  const double probProposal = evalProb(qProposal_);
  double deltaEnergy = 0.;
  if     (probProposal > 0. && prob_ > 0.)
    deltaEnergy = -std::log(probProposal / prob_);
  else if(probProposal > 0.)
    deltaEnergy = -std::numeric_limits<double>::max();
  else if(                     prob_ > 0.)
    deltaEnergy = +std::numeric_limits<double>::max();
  else
  {
    LOGERR << "Encountered negative change in the phase space voulme: "
           << "'probProposal' = " << probProposal << "; "
           << "'prob_' = " << prob_;
    throw std::runtime_error(__PRETTY_FUNCTION__);
  }

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
  for(unsigned i = 0; i < nofDimensions_; ++i)
    x_[i] = (1. - q[i]) * xMin_[i] + q[i] * xMax_[i];
}

double
MEMIntegratorMarkovChain::evalProb(const std::vector<double> & q)
{
  update_x(q);
  return (*integrand_)(x_, nofDimensions_, 0);
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const MEMIntegratorMarkovChain & MX)
  {
    os << "Moves: accepted = " << MX.nofMoves_accepted_ << "; "
       << "rejected = " << MX.nofMoves_rejected_ << " "
       << "(acceptance rate " << static_cast<double>(MX.nofMoves_accepted_) /
                                 (MX.nofMoves_accepted_ + MX.nofMoves_rejected_)
       << " %)\n";

    for(unsigned i = 0; i < MX.nofChains_; ++i)
    {
//--- find the average per batch
      const double integral = functions::avg(
        MX.integral_, MX.nofBatches_ * i, MX.nofBatches_ * (i + 1)
      );
//--- find the standard deviation per batch
      const double integralErr = functions::stdev(
        MX.integral_, MX.nofBatches_ * i, MX.nofBatches_ * (i + 1), integral
      );
      os << std::string(10, ' ') << "chain #" << i << ": "
         << "integral = " << integral << "; "
         << "error = " << integralErr << "\n";
    }

    return os;
  }
}

