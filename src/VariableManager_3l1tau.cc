#include "tthAnalysis/tthMEM/interface/VariableManager_3l1tau.h"
#include "tthAnalysis/tthMEM/interface/tthMEMconstants.h" // constants::
#include "tthAnalysis/tthMEM/interface/Logger.h" // LOG*

#include <algorithm> // std::count_if(), std:find_if()
#include <utility> // std::move()
#include <iomanip> // std::setw()
#include <sstream> // std::stringstream
#include <cstdlib> // EXIT_FAILURE, EXIT_SUCCESS

#include <boost/assign/list_of.hpp> // boost::assign::list_of<>

using namespace tthMEM;

//--- complete the static declarations
const std::unordered_map<Var_3l1tau, VariableManager_3l1tau::Limits, EnumClassHash>
VariableManager_3l1tau::varLimits_ = {
  { Var_3l1tau::kBcosTheta1,     {   -1., +1.                       } },
  { Var_3l1tau::kBphi1,          { -pi(), +pi()                     } },
  { Var_3l1tau::kBcosTheta2,     {   -1., +1.                       } },
  { Var_3l1tau::kBphi2,          { -pi(), +pi()                     } },
  { Var_3l1tau::kZ1,             {    0., +1.                       } },
  { Var_3l1tau::kTauPhi,         { -pi(), +pi()                     } },
  { Var_3l1tau::kTauPhiInv,      { -pi(), +pi()                     } },
  { Var_3l1tau::kTauMinvSquared, {    0., constants::massTauSquared } }
};

const boost::bimaps::bimap<std::string, Var_3l1tau>
VariableManager_3l1tau::varNames_ =
  boost::assign::list_of<decltype(VariableManager_3l1tau::varNames_)::relation>
  ( "bCosTheta1",     Var_3l1tau::kBcosTheta1     )
  ( "bPhi1",          Var_3l1tau::kBphi1          )
  ( "bCosTheta2",     Var_3l1tau::kBcosTheta2     )
  ( "bPhi2",          Var_3l1tau::kBphi2          )
  ( "z1",             Var_3l1tau::kZ1             )
  ( "tauPhi",         Var_3l1tau::kTauPhi         )
  ( "tauPhiInv",      Var_3l1tau::kTauPhiInv      )
  ( "tauMinvSquared", Var_3l1tau::kTauMinvSquared );

//--- private inner classes
VariableManager_3l1tau::Limits::Limits(double begin,
                                       double end)
  : begin_(begin)
  , end_(end)
{}

bool
VariableManager_3l1tau::Limits::isWithin(double val) const
{
  return begin_ <= val && val <= end_;
}

namespace tthMEM
{
  std::ostream &
  operator<<(std::ostream & os,
             const VariableManager_3l1tau::Limits & limits)
  {
    os << "[" << limits.begin_ << "; " << limits.end_ << "]";
    return os;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const VariableManager_3l1tau::Variable & var)
  {
    const static std::unordered_map<VarMode_3l1tau,
      std::string, EnumClassHash> modeMap =
    {
      { VarMode_3l1tau::kFree,      "[free]"      },
      { VarMode_3l1tau::kFixed,     "[fixed]"     },
      { VarMode_3l1tau::kGenerator, "[generator]" }
    };
    os << modeMap.at(var.mode_) << std::setw(12) << std::left
       << std::string("(" + std::to_string(var.idx_) + ")");
    if(var.mode_ == VarMode_3l1tau::kFixed)
      os << std::setw(5) << std::left << var.value_;
    return os;
  }

  std::ostream &
  operator<<(std::ostream & os,
             const VariableManager_3l1tau & vm)
  {
    os << "Variable index mapping (3l1tau):\n";
    for(const auto & kv: VariableManager_3l1tau::varNames_.left)
      os << std::string(10, ' ')
         << std::setw(15) << std::left << kv.first
         << std::setw(12) << std::left << vm.variables_.at(kv.second)
         << "\n";
    os << "Total integration dimension: " << vm.numDimensions_;
    return os;
  }
}

VariableManager_3l1tau::Variable::Variable(VarMode_3l1tau mode,
                                           double value,
                                           int idx)
  : mode_(mode)
  , value_(value)
  , idx_(idx)
{}

//--- public functions
VariableManager_3l1tau::VariableManager_3l1tau()
  : xl_(0)
  , xu_(0)
{
  for(auto var: Enum<Var_3l1tau>())
    variables_[var] = { VarMode_3l1tau::kFree, 0., static_cast<int>(var) };
  updateIndices();
}

VariableManager_3l1tau::VariableManager_3l1tau(const VariableManager_3l1tau & vm) noexcept
  : variables_(vm.variables_)
  , xl_(vm.xl_)
  , xu_(vm.xu_)
  , numDimensions_(vm.numDimensions_)
{}

VariableManager_3l1tau::VariableManager_3l1tau(VariableManager_3l1tau && vm) noexcept
  : variables_(std::move(vm.variables_))
  , xl_(std::move(vm.xl_))
  , xu_(std::move(vm.xu_))
  , numDimensions_(vm.numDimensions_)
{}

int
VariableManager_3l1tau::clamp(const std::string & varName)
{
  if(varNames_.left.count(varName))
  {
    const Var_3l1tau var = varNames_.left.find(varName) -> second;
    variables_[var].mode_ = VarMode_3l1tau::kGenerator;
    variables_[var].idx_ = -1;
    updateIndices();
  }
  else
  {
    LOGERR << "No such integration variable: '" << varName << "'";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

int
VariableManager_3l1tau::clamp(const std::string & varName,
                              double value)
{
  if(varNames_.left.count(varName))
  {
    const Var_3l1tau var = varNames_.left.find(varName) -> second;
    if(! varLimits_.at(var).isWithin(value))
    {
      LOGERR << "Variable '" << varName << "' value = " << value << " "
             << "is not within expected limits: " << varLimits_.at(var);
      return EXIT_FAILURE;
    }
    variables_[var].mode_ = VarMode_3l1tau::kFixed;
    variables_[var].value_ = value;
    variables_[var].idx_ = -1;
    updateIndices();
  }
  else
  {
    LOGERR << "No such integration variable: '" << varName << "'";
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}

unsigned
VariableManager_3l1tau::getCurrentDim() const
{
  return numDimensions_;
}

double
VariableManager_3l1tau::get(Var_3l1tau var,
                            const double * const x) const
{
  const double val = variables_.at(var).mode_ == VarMode_3l1tau::kFree ?
                     x[variables_.at(var).idx_] : variables_.at(var).value_;
  if(! varLimits_.at(var).isWithin(val))
  {
    const std::string varName = varNames_.right.find(var) -> second;
    LOGERR << "Fetched value '" << varName << "' = " << val << " "
           << "is not within expected limits: " << varLimits_.at(var);
    throw EXIT_FAILURE;
  }
  return val;
}

std::string
VariableManager_3l1tau::getArrayString(const double * const x) const
{
  std::stringstream ss;
  ss << std::fixed << std::setprecision(Logger::getFloatPrecision());
  for(auto var: Enum<Var_3l1tau>())
    ss << get(var, x) << (var != Var_3l1tau::Last ? ", " : "");
  return ss.str();
}

int
VariableManager_3l1tau::set(Var_3l1tau var,
                            double value)
{
  if(variables_[var].mode_ != VarMode_3l1tau::kGenerator)
  {
    const std::string varName = varNames_.right.find(var) -> second;
    LOGERR << "Variable '" << varName << "' not configured "
           << "to read from generator level";
    return EXIT_FAILURE;
  }
  if(! varLimits_.at(var).isWithin(value))
  {
    const std::string varName = varNames_.right.find(var) -> second;
    LOGERR << "Variable '" << varName << "' = " << value << " "
           << "not within expected limits: " << varLimits_.at(var);
    return EXIT_FAILURE;
  }
  variables_[var].value_ = value;
  return EXIT_SUCCESS;
}

const double * const
VariableManager_3l1tau::getXL() const
{
  return xl_.data();
}

const double * const
VariableManager_3l1tau::getXU() const
{
  return xu_.data();
}

//--- private functions
void
VariableManager_3l1tau::updateIndices()
{
  xl_.clear();
  xu_.clear();
  int i = 0;
  for(auto & kv: variables_)
    if(kv.second.mode_ == VarMode_3l1tau::kFree)
    {
      xl_.push_back(varLimits_.at(kv.first).begin_);
      xu_.push_back(varLimits_.at(kv.first).end_);
      kv.second.idx_ = i++;
    }
  numDimensions_ = i;
}
