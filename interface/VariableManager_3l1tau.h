#ifndef VARIABLEMANAGER_3L1TAU_H
#define VARIABLEMANAGER_3L1TAU_H

#include "tthAnalysis/tthMEM/interface/tthMEMenums.h" // EnumClassHash
#include "tthAnalysis/tthMEM/interface/Limits.h" // Limits

#include <string> // std::string
#include <unordered_map> // std::unordered_map<,,>
#include <ostream> // std::ostream
#include <vector> // std::vector<>

#include <boost/bimap/bimap.hpp> // boost::bimaps::bimap<,>

namespace tthMEM
{
  /**
   * @brief The enum class used to acces and set the underlying values
   */
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

  /**
   * @brief Enum class that specificies, whether a given variable is
   *          - kFree i.e. expected to be sampled
   *          - kFixed i.e. set to a provided value
   *          - kGenerator i.e. read from the underlying generator level
   */
  enum class VarMode_3l1tau
  {
    kFree = 0,
    kFixed = 1,
    kGenerator = 2
  };

  /**
   * @brief Class for managing variables -- their access and values
   */
  class VariableManager_3l1tau
  {
  public:
    /**
     * @brief Default constructor
     *
     * Sets all variables to kFree by default. Also defines
     * the integration boundaries.
     */
    VariableManager_3l1tau();
    /**
     * @brief Simple copy constructor
     * @param vm Variable manager the members of which are copied over
     *
     * Copies all non-static member variables over to the this class
     */
    VariableManager_3l1tau(const VariableManager_3l1tau & vm) noexcept;

    /**
     * @brief Clamps a given variable to kGenerator mode
     * @param varName The variable name
     * @return EXIT_FAILURE (!= 0) if the variable does not exist
     *         in varNames_, EXIT_SUCCESS (== 0) otherwise
     */
    int
    clamp(const std::string & varName);

    /**
     * @brief Clamps a given variable to kFixed mode and sets its
     *        value to the provided one
     * @param varName The variable name
     * @param value   The value the variable is expected to take
     * @return EXIT_FAILURE (!= 0) if the variable does not exist
     *         in varNames_ or if the provided value is not within
     *         the variable range limits specified by varLimits_,
     *         EXIT_SUCCESS (== 0) otherwise
     */
    int
    clamp(const std::string & varName,
          double value);

    /**
     * @brief Returns integration dimension
     *
     * Basically counts, how many variables are in the mode kFree.
     * If none are, then the function promptly returns a 0.
     */
    unsigned
    getCurrentDim() const;

    /**
     * @brief Returns the variable value
     * @param var Variable type
     * @param x   The C-type array containing the sampled values
     * @return    The true value of the variable
     *
     * If the variable is in kFree mode, its value is read
     * from the provided array x; otherwise, the value is read
     * locally. For mode kGenerator it is assumed that the user
     * updates the generator level value using set() function.
     *
     * If the fetched value is not withing the actual boundaries,
     * the function throws. Don't expect to catch it but fix
     * the problem in your configuration file!
     */
    double
    get(Var_3l1tau var,
        const double * const x) const;

    /**
     * @brief Returns a comma-separated list of values
     *        of all the variables
     * @param x The C-type array containing the sampled values
     * @return The final string
     */
    std::string
    getArrayString(const double * const x) const;

    /**
     * @brief Returns the name of the variable
     * @param var Variable enum
     * @return The variable name
     */
    static std::string
    getVarName(Var_3l1tau var);

    /**
     * @brief Returns the integration limits of a variable
     * @param var The variable
     * @return The limits
     */
    Limits
    getVarLimits(Var_3l1tau var) const;

    /**
     * @brief Sets the variable to a given value at runtime
     * @param var   The variable to which the new value is assigned
     * @param value The new value
     * @return EXIT_FAILURE (!= 0) if the variable is not in kGenerator
     *         mode or if the provided value is not within the expected
     *         limits, EXIT_SUCCESS (== 0) otherwise.
     */
    int
    set(Var_3l1tau var,
        double value);

    /**
     * @brief Returns the upper right corner of the integration space
     * @return The C-array containing the upper limits of the integration
     *         variables
     *
     * Note that the returned value includes the upper limit of variables
     * which are declared as kFree.
     *
     * The function is needed by the integration libraries such as VEGAS
     * and VAMP.
     */
    const double * const
    getXU() const;

    /**
     * @brief Returns the lower left corner of the integration space
     * @return The C-array containing the lower limits of the integration
     *         variables
     *
     * Note that the returned value includes the lower limit of variables
     * which are declared as kFree.
     *
     * The function is needed by the integration libraries such as VEGAS
     * and VAMP.
     */
    const double * const
    getXL() const;

    std::vector<Var_3l1tau> generatorLevels;

  private:
    /**
     * @brief Simple helper struct holding variable mode, its value
     *        and array index. The latter is used to read the actual
     *        variable value from a C-style array whenever the mode
     *        is set to kFree. Note that if the mode is not kFree,
     *        the idx_ value is expected to be -1 (thus rendering
     *        the attempts at using it as an array subscript index
     *        futile). Conversely, the value is ignored if the mode
     *        is set to kFree.
     */
    struct Variable
    {
      /**
        * Default constructor needed by templated containers
        */
      Variable() = default;
      /**
       * @brief The default constructor
       * @param mode  The variable mode (@see VarMode_3l1tau)
       * @param value Its value
       * @param idx   The variable index
       */
      Variable(VarMode_3l1tau mode,
               double value,
               int idx);

      VarMode_3l1tau mode_; ///< Variable mode  (kFree, kFixed or kGenerator)
      double value_;        ///< Variable value (ignored if not kFree)
      int idx_;             ///< Variable index (used to read from provided array)

      /**
       * @brief Simple ostream operator for pertty-printing the variable meta-data
       * @param os  The std::ostream to print to
       * @param var The given Variable instance
       * @return The modified ostream
       */
      friend std::ostream &
      operator<<(std::ostream & os,
                 const Variable & var);
    };

    /**
     * @brief Recalculates the integration dimension and limits by
     *        taking only the variable into account which are in kFree mode.
     */
    void
    updateIndices();

    static const std::unordered_map<Var_3l1tau, Limits, EnumClassHash> varLimits_;
    ///< hardcoded map integration limits
    static const boost::bimaps::bimap<std::string, Var_3l1tau> varNames_;
    ///< hardcoded bi-directional map between variable type and corresponding string

    std::unordered_map<Var_3l1tau, Variable, EnumClassHash> variables_;
    ///< map containing all variables and corresponding meta-data
    std::vector<double> xl_; ///< lower left integration limits
    std::vector<double> xu_; ///< upper right integration limits
    unsigned numDimensions_; ///< the integration dimensionality

  /* just so that the compiler sees the inner friend classes */
  public:
    friend std::ostream &
    operator<<(std::ostream &,
               const Limits & limits);

    friend std::ostream &
    operator<<(std::ostream &,
               const VariableManager_3l1tau::Variable & var);

    friend std::ostream &
    operator<<(std::ostream &,
               const VariableManager_3l1tau & vm);
  };
}

#endif // VARIABLEMANAGER_3L1TAU_H
