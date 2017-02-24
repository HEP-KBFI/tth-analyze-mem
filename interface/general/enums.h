#ifndef TTHMEMENUMS_H
#define TTHMEMENUMS_H

#include <cstddef> // std::size_t

namespace tthMEM
{
  /**
   * @brief For selecting ttH or ttZ matrix element in 3l1tau channel
   */
  enum ME_mg5_3l1tau
  {
    kTTH = 0,
    kTTZ = 1
  };

  /**
   * @brief Hash function for enum classes used in std::unordered_map<>
   *
   * Cheers to the author at http://stackoverflow.com/a/24847480
   */
  struct EnumClassHash
  {
    template <typename T>
    std::size_t operator()(T t) const
    {
      return static_cast<std::size_t>(t);
    }
  };

  /**
   * @brief Meta-class for looping over enum classes which have
   *        First and Last entries specified.
   *
   * Cheers to the author at http://stackoverflow.com/a/8498694
   */
  template <typename T>
  struct Enum
  {
    struct Iterator
    {
      Iterator(int val)
        : val_(val)
      {}

      T
      operator*(void) const
      {
        return static_cast<T>(val_);
      }

      void
      operator++(void)
      {
        ++val_;
      }

      bool
      operator!=(Iterator rhs)
      {
        return val_ != rhs.val_;
      }
    private:
      int val_;
    };
  };

  template <typename T>
  typename Enum<T>::Iterator begin(Enum<T>)
  {
    return typename Enum<T>::Iterator(static_cast<int>(T::First));
  }

  template <typename T>
  typename Enum<T>::Iterator end(Enum<T>)
  {
    return typename Enum<T>::Iterator(static_cast<int>(T::Last) + 1);
  }
}

#endif // TTHMEMENUMS_H
