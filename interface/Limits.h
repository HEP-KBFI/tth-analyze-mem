#ifndef LIMITS_H
#define LIMITS_H

#include <ostream> // std::ostream

namespace tthMEM
{
  /**
   * @brief Helper POD struct which holds the physical (or equivalently,
   *        the integration) limits of a particular variable.
   */
  struct Limits
  {
    Limits() = default;
    Limits(double begin,
           double end);

    /**
     * @brief Checks whether a given value falls within the limits
     * @param val The given value
     * @return True, if the value is indeed between the limits,
     *         false otherwise
     */
    bool
    isWithin(double val) const;

    /**
     * @brief Simple ostream operator for pretty-printing the limits
     * @param os     The std::ostream to print to
     * @param limits The instance of this class
     * @return The modified ostream
     */
    friend std::ostream &
    operator<<(std::ostream & os,
               const Limits & limits);

    double begin_, end_; ///< end points
  };
}

#endif // LIMITS_H
