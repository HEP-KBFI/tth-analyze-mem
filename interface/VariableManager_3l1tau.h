#ifndef VARIABLEMANAGER_3L1TAU_H
#define VARIABLEMANAGER_3L1TAU_H

#include <string> // std::string
#include <unordered_map> // std::unordered_map<,,>
#include <ostream> // std::ostream
#include <vector> // std::vector<>

#include <boost/bimap/bimap.hpp> // boost::bimaps::bimap<,>

#include "tthAnalysis/tthMEM/interface/tthMEMauxFunctions.h" // EnumClassHash

namespace tthMEM
{
  enum class Var_3l1tau
  {
    kBcosTheta1     = 0, // cosine of polar angle of the 1st neutrino
    kBphi1          = 1, // azimuthal angle of the 1st neutrino
    kBcosTheta2     = 2, // cosine of polar angle of the 2nd neutrino
    kBphi2          = 3, // azimuthal angle of the 2nd neutrino
    kZ1             = 4, // energy fraction of hadronic tau system
    kTauPhi         = 5, // rotation angle of hadronic tau neutrino
    kTauPhiInv      = 6, // (invisible) rotation angle of leptonic tau nu
    kTauMinvSquared = 7, // (invisible) mass of leptonic neutrino pair
//--- for iteration
    First = kBcosTheta1,
    Last = kTauMinvSquared
  };

  enum class VarMode_3l1tau
  {
    kFree = 0,
    kFixed = 1,
    kGenerator = 2
  };

  class VariableManager_3l1tau
  {
  public:
    VariableManager_3l1tau();
    VariableManager_3l1tau(const VariableManager_3l1tau & vm);
    VariableManager_3l1tau(VariableManager_3l1tau && vm) noexcept;

    int
    clamp(const std::string & varName);

    int
    clamp(const std::string & varName,
          double value);

    unsigned
    getCurrentDim() const;

    double
    get(Var_3l1tau var,
        const double * const x) const;

    std::string
    getArrayString(const double * const x) const;

    int
    set(Var_3l1tau var,
        double value);

    const double * const
    getXU() const;

    const double * const
    getXL() const;

  private:
    struct Limits
    {
      Limits() = default;
      Limits(double begin,
             double end);

      bool
      isWithin(double val) const;

      friend std::ostream &
      operator<<(std::ostream & os,
                 const Limits & limits);

      double begin_, end_;
    };

    struct Variable
    {
      Variable() = default;
      Variable(VarMode_3l1tau mode,
               double value,
               int idx);

      VarMode_3l1tau mode_;
      double value_;
      int idx_;

      friend std::ostream &
      operator<<(std::ostream & os,
                 const Variable & var);
    };

    void
    updateIndices();

    static const std::unordered_map<Var_3l1tau, Limits, EnumClassHash> varLimits_;
    static const boost::bimaps::bimap<std::string, Var_3l1tau> varNames_;

    std::unordered_map<Var_3l1tau, Variable, EnumClassHash> variables_;
    std::vector<double> xl_;
    std::vector<double> xu_;
    unsigned numDimensions_;

  public:
    friend std::ostream &
    operator<<(std::ostream &,
               const VariableManager_3l1tau::Limits & limits);

    friend std::ostream &
    operator<<(std::ostream &,
               const VariableManager_3l1tau::Variable & var);

    friend std::ostream &
    operator<<(std::ostream &,
               const VariableManager_3l1tau & vm);
  };
}

#endif // VARIABLEMANAGER_3L1TAU_H
