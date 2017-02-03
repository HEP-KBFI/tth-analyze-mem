#include "tthAnalysis/tthMEM/interface/MEMIntegratorMarkovChain.h" // MEMIntegratorMarkovChain
#include "tthAnalysis/tthMEM/interface/MEMIntegratorVEGAS.h" // MEMIntegratorVEGAS

using namespace tthMEM;

double eval_x  (double * x,
                std::size_t dimension       __attribute__((unused)),
                void * additionalParameters __attribute__((unused)))
{ return      x[0];               }
double eval_x2 (double * x,
                std::size_t dimension       __attribute__((unused)),
                void * additionalParameters __attribute__((unused)))
{ return pow2(x[0]);              }
double eval_x3 (double * x,
                std::size_t dimension       __attribute__((unused)),
                void * additionalParameters __attribute__((unused)))
{ return pow2(x[0]);              }
double eval_xy (double * x,
                std::size_t dimension       __attribute__((unused)),
                void * additionalParameters __attribute__((unused)))
{ return      x[0]  *      x[1];  }
double eval_x2y(double * x,
                std::size_t dimension       __attribute__((unused)),
                void * additionalParameters __attribute__((unused)))
{ return pow2(x[0]) *      x[1];  }
double eval_xy2(double * x,
                std::size_t dimension       __attribute__((unused)),
                void * additionalParameters __attribute__((unused)))
{ return      x[0]  * pow2(x[1]); }

struct IntMeta
{
  IntMeta()
    : ptr(nullptr)
    , name("")
    , dim(0)
    , xl({})
    , xu({})
    , value(0.)
  {}
  IntMeta(MEMIntegratorBase::gPtr_C ptr_,
          const std::string & name_,
          unsigned dim_,
          const std::vector<double> xl_,
          const std::vector<double> xu_,
          double value_)
    : ptr(ptr_)
    , name(name_)
    , dim(dim_)
    , xl(xl_)
    , xu(xu_)
    , value(value_)
  {}
  MEMIntegratorBase::gPtr_C ptr;
  std::string               name;
  unsigned                  dim;
  std::vector<double>       xl;
  std::vector<double>       xu;
  double                    value;
};

struct IntResult
{
  IntResult()
    : p(0.)
    , pErr(0.)
    , intAlgo(nullptr)
  {}

  void
  reset()
  {
    if(intAlgo)
    {
      delete intAlgo;
      intAlgo = nullptr;
    }
    p = 0.;
    pErr = 0;
  }

  void
  compute()
  {
    intAlgo -> integrate(info.ptr, info.xl.data(), info.xu.data(), info.dim, p, pErr);
  }

  friend std::ostream &
  operator<<(std::ostream & os,
             const IntResult & result)
  {
    os << result.p << " (" << result.pErr;
    const double absErr = result.pErr > 1.e-6 ? std::fabs(result.info.value - result.p) : 0.;
    const bool isOk = absErr < result.pErr;
    os << (isOk ? " >= " : " <  ") << absErr << " => " << (isOk ? "OK" : "FAIL") << ')';
    return os;
  }

  double p;
  double pErr;
  MEMIntegratorBase * intAlgo;
  IntMeta info;
};

int
main()
{
  Logger::setLogLevel("warning");
  std::cout << std::setprecision(6);
  std::cout.setf(std::ios::fixed, std::ios::floatfield);

  std::vector<IntMeta> integrands = {
    { eval_x,   "f(x)    = x      ", 1, { 0.     }, { std::sqrt(2.)                              }, 1. },
    { eval_x2,  "f(x)    = x^2    ", 1, { 0.     }, { std::cbrt(3.)                              }, 1. },
    { eval_x3,  "f(x)    = x^3    ", 1, { 0.     }, { std::sqrt(2.)                              }, 1. },
    { eval_xy,  "f(x, y) = x * y  ", 2, { 0., 0. }, { std::sqrt(2.),        std::sqrt(2.)        }, 1. },
    { eval_x2y, "f(x, y) = x^2 * y", 2, { 0., 0. }, { std::pow(6., 1. / 5), std::pow(6., 1. / 5) }, 1. }
  };
  IntResult result;

  const unsigned nof_chains   = 1;
  const unsigned nof_batches  = 100;
  const unsigned nof_sweeping = 1000;

  for(const auto & integrand: integrands)
  {
    result.info = integrand;
    for(unsigned nof_calls: {1000, 10000, 100000})
    {
      result.intAlgo = new MEMIntegratorVEGAS(
        roundToNearestUInt(nof_calls * 0.2),
        roundToNearestUInt(nof_calls * 0.8),
        1,
        2.
      );
      result.compute();
      std::cout << "VEGAS: " << integrand.name << " @ " << nof_calls << " calls: "
                << std::string(73, ' ') << result << '\n';
      result.reset();

      for(const std::string & mode: { "uniform", "gaussian", "none" })
      {
        for(double nu: { 0.50, 0.71, 0.99 })
        {
          for(double T0: { 1., 15. })
          {
            for(double epsilon0: { 1.e-2, 1.e-4 })
            {
              const unsigned nof_burnin   = roundToNearestUInt(0.10 * nof_calls / nof_chains);
              const unsigned nof_sampling = roundToNearestUInt(0.90 * nof_calls);
              const unsigned nof_ann1     = roundToNearestUInt(0.20 * nof_burnin);
              const unsigned nof_ann2     = roundToNearestUInt(0.60 * nof_burnin);
              const double alpha = 1. - 10. / nof_burnin;
              result.intAlgo = new MEMIntegratorMarkovChain(
                mode,
                nof_burnin,
                nof_sampling,
                nof_ann1,
                nof_ann2,
                nof_sweeping,
                T0,
                alpha,
                nof_chains,
                nof_batches,
                epsilon0,
                nu
              );
              result.compute();
              std::cout << "MXMC:  " << integrand.name << " @ " << nof_calls << " calls, "
                           " mode = " << std::setw(10) << std::left << mode
                        << "  nu = " << nu << ", T0 = " << std::setw(9) << std::left
                        << T0 << ", epsilon0  = " << epsilon0 << ": " << result << '\n';
              result.reset();
            }
          }
        }
      }
    }
  }
  return EXIT_SUCCESS;
}
