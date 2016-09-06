#ifndef TTHMEMVECFUNCTIONS_H
#define TTHMEMVECFUNCTIONS_H

#include <vector> // std::vector<>
#include <ostream> // std::ostream
#include <iterator> // std::ostream_iterator<>, std::back_inserter()
#include <algorithm> // std::copy(), std::transform()
#include <functional> // std::bind1st(), std::bind2nd(), std::multiplies<>, ...
  // ..., std::plus<>, std::minus<>, std::divides<>

#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

namespace tthMEM
{
  namespace vec
  {
    /**
     * @brief Calculates L2 norm aka Euclidean norm of a given vector of doubles
     * @param vec        The vector of doubles
     * @param shiftValue The value which is subtracted from every element beforehand
     * @return Magnitude of the vector in Euclidean space
     */
    double
    l2(const std::vector<double> & vec,
       double shiftValue = 0.);

    /**
     * @brief Calculates L2 norm aka Euclidean norm of a given vector of doubles
     * @param vec            The vector
     * @param shiftFromBegin End point up to which the norm is calculated
     * @param shiftValue     The value which is subtracted from every
     *                       element beforehand
     * @return The norm of the subvector in the element range [0, shiftFromBegin)
     */
    double
    l2(const std::vector<double> & vec,
       unsigned shiftFromBegin,
       double shiftValue = 0.);

    /**
     * @brief Calculates L2 norm aka Euclidean norm of a given vector of doubles
     * @param vec                  The vector
     * @param shiftFromBegin_begin The starting point from which the norm is
     *                             calculated
     * @param shiftFromBegin_end   The end point up to which the norm is calculated
     * @param shiftValue           The value which is subtracted from every
     *                             element beforehand
     * @return The norm of the subvector in the element range
     *         [shiftFromBegin_begin, shiftFromBegin_end)
     */
    double
    l2(const std::vector<double> & vec,
       unsigned shiftFromBegin_begin,
       unsigned shiftFromBegin_end,
       double shiftValue = 0.);

    /**
     * @brief Calculates average of a given vector of doubles
     * @param vec The vector
     * @return Average of the vector
     */
    double
    avg(const std::vector<double> & vec);

    /**
     * @brief Calculates average of a given vector of doubles
     * @param vec            The vector
     * @param shiftFromBegin End point up to which the average is calculated
     * @return Average of the subvector in the element range [0, shiftFromBegin)
     */
    double
    avg(const std::vector<double> & vec,
        unsigned shiftFromBegin);

    /**
     * @brief Calculates average of a given vector of doubles
     * @param vec                  The vector
     * @param shiftFromBegin_begin The starting point from which the average is
     *                             calculated
     * @param shiftFromBegin_end   The end point up to which the average is
     *                             calculated
     * @return Average of the subvector in the element range
     *         [shiftFromBegin_begin, shiftFromBegin_end)
     */
    double
    avg(const std::vector<double> & vec,
        unsigned shiftFromBegin_begin,
        unsigned shiftFromBegin_end);

    /**
     * @brief Calculates standard deviation of a given vector of doubles
     * @param vec     The vector
     * @param average Average of the vector
     * @return Standard deviation of the vector
     */
    double
    stdev(const std::vector<double> & vec,
          double average);

    /**
     * @brief Calculates standard deviation of a given vector of doubles
     * @param vec            The vector
     * @param shiftFromBegin End point up to which the standard deviation is
     *                       calculated
     * @param average        The average of the subvector
     * @return Standard deviation of the subvector in the element range
     *         [0, shiftFromBegin)
     */
    double
    stdev(const std::vector<double> & vec,
          unsigned shiftFromBegin,
          double average);

    /**
     * @brief Calculates standard deviation of a given vector of doubles
     * @param vec                  The vector
     * @param shiftFromBegin_begin The starting point from which the standard
     *                             deviation is calculated
     * @param shiftFromBegin_end   The end point up to which the standard
     *                             deviation is calculated
     * @param average              The average of the subvector
     * @return Standard deviation of the subvector in the element range
     *         [shiftFromBegin_begin, shiftFromBegin_end)
     */
    double
    stdev(const std::vector<double> & vec,
          unsigned shiftFromBegin_begin,
          unsigned shiftFromBegin_end,
          double average);

    /**
     * @brief Returns a subvector
     * @param vec                  The vector the subvector of which is taken
     * @param shiftFromBegin_begin The first element of the subvector
     * @param shiftFromBegin_end   The last element of the subvector
     * @return A newly allocated subvector in the element range
     *         [shiftFromBegin_begin, shiftFromBegin_end)
     */
    template <typename VectorElementType>
    std::vector<VectorElementType>
    subv(const std::vector<VectorElementType> & vec,
         unsigned shiftFromBegin_begin,
         unsigned shiftFromBegin_end)
    {
      if(! (shiftFromBegin_begin <= shiftFromBegin_end &&
            shiftFromBegin_begin <= vec.size() &&
            shiftFromBegin_end <= vec.size()))
        throw_line("range error")
          << "Invalid arguments: "
          << "'shiftFromBegin_begin' = " << shiftFromBegin_begin << "; "
          << "'shiftFromBegin_end' = " << shiftFromBegin_end << "; "
          << "vector size = " << vec.size();

      return std::vector<VectorElementType>(vec.begin() + shiftFromBegin_begin,
                                            vec.begin() + shiftFromBegin_end);
    }

    /**
     * @brief Returns a subvector
     * @param vec            The vector the subvector of which is taken
     * @param shiftFromBegin The number of elements are taken from the beginning
     * @return A newly allocated subvector in the element range [0, shiftFromBegin)
     */
    template <typename VectorElementType>
    std::vector<VectorElementType>
    subv(const std::vector<VectorElementType> & vec,
         unsigned shiftFromBegin)
    {
      return subv(vec, 0, shiftFromBegin);
    }

    /**
     * @brief genv
     * @param functor
     * @param nofElements
     * @return
     */
    template <typename VectorElementType>
    std::vector<VectorElementType>
    genv(std::function<VectorElementType()> functor,
         std::size_t nofElements)
    {
      std::vector<VectorElementType> result;
      std::generate_n(std::back_inserter(result), nofElements, functor);
      return result;
    }
  }

  /**
   * @brief Prints comma-separated vector elements
   * @param os  Stream to print to
   * @param vec The vector
   * @return Modified stream
   *
   * Assumes that the typename T is printable (otherwise expect compilation error)
   */
  template <typename VectorElementType>
  std::ostream &
  operator<<(std::ostream & os,
             const std::vector<VectorElementType> & vec)
  {
    if(! vec.size()) return os << "{}";
    os << "{ ";
    if(vec.size() > 1)
      std::copy(vec.begin(), vec.end() - 1,
                std::ostream_iterator<VectorElementType>(os, ", "));
    os << vec[vec.size() - 1];
    os << " }";
    return os;
  }


  /**
   * @brief Computes Hadamard product between two vectors
   * @param lhs Left-hand side vector
   * @param rhs Right-hand side vector
   * @return A copy of the Hadamard product
   *
   * Check https://en.wikipedia.org/wiki/Hadamard_product_(matrices)
   * It's not dot product as can be seen from the return type!
   * The complementary inverse operation (division) is not defined
   * at the time of writing
   */
  template <typename VectorElementType>
  std::vector<VectorElementType>
  operator*(const std::vector<VectorElementType> & lhs,
            const std::vector<VectorElementType> & rhs)
  {
    if(lhs.size() != rhs.size())
      throw_line("range error") << "LHS size = " << lhs.size() << "; "
                                << "RHS size = " << rhs.size();
    std::vector<VectorElementType> result;
    std::transform(
      lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter(result),
      std::multiplies<VectorElementType>()
    );
    return result;
  }

  /**
   * @brief Computes sum of two vectors
   * @param lhs Left-hand side vector
   * @param rhs Right-hand side vector
   * @return A copy of the sum
   */
  template <typename VectorElementType>
  std::vector<VectorElementType>
  operator+(const std::vector<VectorElementType> & lhs,
            const std::vector<VectorElementType> & rhs)
  {
    if(lhs.size() != rhs.size())
      throw_line("range error") << "LHS size = " << lhs.size() << "; "
                                << "RHS size = " << rhs.size();
    std::vector<VectorElementType> result;
    std::transform(
      lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter(result),
      std::plus<VectorElementType>()
    );
    return result;
  }

  /**
   * @brief Computes difference of two vectors
   * @param lhs Left-hand side vector
   * @param rhs Right-hand side vector
   * @return A copy of the difference
   */
  template <typename VectorElementType>
  std::vector<VectorElementType>
  operator-(const std::vector<VectorElementType> & lhs,
            const std::vector<VectorElementType> & rhs)
  {
    if(lhs.size() != rhs.size())
      throw_line("range error") << "LHS size = " << lhs.size() << "; "
                                << "RHS size = " << rhs.size();
    std::vector<VectorElementType> result;
    std::transform(
      lhs.begin(), lhs.end(), rhs.begin(), std::back_inserter(result),
      std::minus<VectorElementType>()
    );
    return result;
  }

  /**
   * @brief Computes the sum of a vector and a scalar by adding
   *        the scalar to each element in the vector
   * @param lhs The scalar
   * @param rhs The vector
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator+(ScalarType lhs,
            const std::vector<VectorElementType> & rhs)
  {
    std::vector<VectorElementType> result;
    std::transform(
      rhs.begin(), rhs.end(), std::back_inserter(result),
      std::bind1st(std::plus<VectorElementType>(),
                   static_cast<VectorElementType>(lhs))
    );
    return result;
  }

  /**
   * @brief Computes the sum of a vector and a scalar by adding
   *        the scalar to each element in the vector
   * @param lhs The vector
   * @param rhs The scalar
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator+(const std::vector<VectorElementType> & lhs,
            ScalarType rhs)
  {
    return rhs + lhs;
  }

  /**
   * @brief Computes the difference of a scalar and a vector by
   *        subtracting each element in the vector from the scalar
   * @param lhs The scalar
   * @param rhs The vector
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator-(ScalarType lhs,
            const std::vector<VectorElementType> & rhs)
  {
    std::vector<VectorElementType> result;
    std::transform(
      rhs.begin(), rhs.end(), std::back_inserter(result),
      std::bind1st(std::minus<VectorElementType>(),
                   static_cast<VectorElementType>(lhs))
    );
    return result;
  }

  /**
   * @brief Computes the difference of a vector and a scalar by
   *        subtracting the scalar from each element in the vector
   * @param lhs The vector
   * @param rhs The scalar
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator-(const std::vector<VectorElementType> & lhs,
            ScalarType rhs)
  {
    std::vector<VectorElementType> result;
    std::transform(
      lhs.begin(), lhs.end(), std::back_inserter(result),
      std::bind2nd(std::minus<VectorElementType>(),
                   static_cast<ScalarType>(rhs))
    );
    return result;
  }

  /**
   * @brief Computes the product of a vector and a scalar by
   *        multiplying the scalar with the vector element-wise
   * @param lhs The scalar
   * @param rhs The vector
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator*(ScalarType lhs,
            const std::vector<VectorElementType> & rhs)
  {
    std::vector<VectorElementType> result;
    std::transform(
      rhs.begin(), rhs.end(), std::back_inserter(result),
      std::bind1st(std::multiplies<VectorElementType>(),
                   static_cast<VectorElementType>(lhs))
    );
    return result;
  }

  /**
   * @brief Computes the product of a vector and a scalar by
   *        multiplying the scalar with the vector element-wise
   * @param lhs The vector
   * @param rhs The scalar
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator*(const std::vector<VectorElementType> & lhs,
            ScalarType rhs)
  {
    return rhs * lhs;
  }

  /**
   * @brief Computes the quotient of a scalar and a vector by
   *        dividing the scalar with every element in the vector
   * @param lhs The scalar
   * @param rhs The vector
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator/(ScalarType lhs,
            const std::vector<VectorElementType> & rhs)
  {
    std::vector<VectorElementType> result;
    std::transform(
      rhs.begin(), rhs.end(), std::back_inserter(result),
      [&lhs](const VectorElementType element) -> VectorElementType
      {
        if(element == 0)
          throw_line("invalid argument")
            << "Dividing with a vector that contains a zero element";
        return static_cast<VectorElementType>(lhs) / element;
      }
    );
    return result;
  }

  /**
   * @brief Computes the quotient of a scalar and a vector by
   *        dividing every element in the vector with the scalar
   * @param lhs The vector
   * @param rhs The scalar
   * @return A copy of the resulting vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType>
  operator/(const std::vector<VectorElementType> & lhs,
            ScalarType rhs)
  {
    if(rhs == 0)
      throw_line("invalid argument")
        << "Passed zero divisor";
    std::vector<VectorElementType> result;
    std::transform(
      lhs.begin(), lhs.end(), std::back_inserter(result),
      std::bind2nd(std::divides<VectorElementType>(),
                   static_cast<VectorElementType>(rhs))
    );
    return result;
  }

  /**
   * @brief Adds a scalar to each element of a vector
   * @param lhs The vector
   * @param rhs The scalar
   * @return The modified vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType> &
  operator+=(std::vector<VectorElementType> & lhs,
             ScalarType rhs)
  {
    std::transform(
      lhs.begin(), lhs.end(), lhs.begin(),
      std::bind2nd(std::plus<VectorElementType>(),
                   static_cast<VectorElementType>(rhs))
    );
    return lhs;
  }

  /**
   * @brief Subtracts a scalar from each element of a vector
   * @param lhs The vector
   * @param rhs The scalar
   * @return The modified vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType> &
  operator-=(std::vector<VectorElementType> & lhs,
             ScalarType rhs)
  {
    std::transform(
      lhs.begin(), lhs.end(), lhs.begin(),
      std::bind2nd(std::minus<VectorElementType>(),
                   static_cast<VectorElementType>(rhs))
    );
    return lhs;
  }

  /**
   * @brief Multiplies each element in a vector with a scalar
   * @param lhs The vector
   * @param rhs The scalar
   * @return The modified vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType> &
  operator*=(std::vector<VectorElementType> & lhs,
             ScalarType rhs)
  {
    std::transform(
      lhs.begin(), lhs.end(), lhs.begin(),
      std::bind2nd(std::multiplies<VectorElementType>(),
                   static_cast<VectorElementType>(rhs))
    );
    return lhs;
  }

  /**
   * @brief Divides each element in a vector with a scalar
   * @param lhs The vector
   * @param rhs The scalar
   * @return The modified vector
   *
   * Note that the scalar type ScalarType is always casted into
   * vector element type VectorElementType.
   */
  template <typename VectorElementType,
            typename ScalarType>
  std::vector<VectorElementType> &
  operator/=(std::vector<VectorElementType> & lhs,
             ScalarType rhs)
  {
    if(rhs == 0)
      throw_line("invalid argument")
        << "Division by zero";
    std::transform(
      lhs.begin(), lhs.end(), lhs.begin(),
      std::bind2nd(std::divides<VectorElementType>(),
                   static_cast<VectorElementType>(rhs))
    );
    return lhs;
  }

  /**
   * @brief Adds vector together by modifying the left operand
   * @param lhs The left vector
   * @param rhs The right vector
   * @return The modified (left) vector
   */
  template <typename VectorElementType>
  std::vector<VectorElementType> &
  operator+=(std::vector<VectorElementType> & lhs,
             const std::vector<VectorElementType> & rhs)
  {
    if(lhs.size() != rhs.size())
      throw_line("range error") << "LHS size = " << lhs.size() << "; "
                                << "RHS size = " << rhs.size();
    std::transform(
      lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
      std::plus<VectorElementType>()
    );
    return lhs;
  }

  /**
   * @brief Subtracts right vector from the left vector
   * @param lhs The left vwector
   * @param rhs The right vector
   * @return The modified (left) vector
   */
  template <typename VectorElementType>
  std::vector<VectorElementType> &
  operator-=(std::vector<VectorElementType> & lhs,
             const std::vector<VectorElementType> & rhs)
  {
    if(lhs.size() != rhs.size())
      throw_line("range error") << "LHS size = " << lhs.size() << "; "
                                << "RHS size = " << rhs.size();
    std::transform(
      lhs.begin(), lhs.end(), rhs.begin(), lhs.begin(),
      std::minus<VectorElementType>()
    );
    return lhs;
  }
}

#endif // TTHMEMVECFUNCTIONS_H
