#ifndef INDEXWRAPPER_H
#define INDEXWRAPPER_H

#include "tthAnalysis/tthMEM/interface/Exception.h" // throw_line()

#include <numeric> // std::iota()

template<typename T,
         unsigned NofObjects>
struct IndexWrapper
{
  IndexWrapper()
    : currentPermutation_(0)
    , maxCurrentPermutation_(0)
  {
    std::iota(defaultPermutation, defaultPermutation + NofObjects, 0);
  }

  void
  setPermutationPtrs(const std::vector<std::vector<unsigned>> & permutations,
                     unsigned * currentPermutation,
                     unsigned maxCurrentPermutation)
  {
    if(! permutations.size())
      throw_line("runtime error") << "Permutation size is zero";
    for(const std::vector<unsigned> & permutation: permutations)
      if(permutation.size() != NofObjects)
        throw_line("runtime error")
          << "Permutation size = " << permutation.size() << " does not equal to "
          << "the number of objects = " << NofObjects;
    if(! currentPermutation)
      throw_line("runtime error") << "Passed nullptr for 'currentPermutation'";
    if(! maxCurrentPermutation)
      throw_line("runtime error") << "Passed 0 for 'maxCurrentPermutation'";

    permutations_ = permutations;
    currentPermutation_ = currentPermutation;
    maxCurrentPermutation_ = maxCurrentPermutation;
  }

  T &
  operator[](const int index)
  {
    if(currentPermutation_)
    {
      if(*currentPermutation_ >= maxCurrentPermutation_)
        throw_line("runtime error")
          << "'currentPermutation' ( = " << *currentPermutation_ << ") >= "
          << "'maxCurrentPermutation' ( = " << maxCurrentPermutation_ << ')';
      return objects[permutations_[*currentPermutation_][index]];
    }
    return objects[defaultPermutation[index]];
  }

  const T &
  operator[](const int index) const
  {
    if(currentPermutation_)
    {
      if(*currentPermutation_ >= maxCurrentPermutation_)
        throw_line("runtime error")
          << "'currentPermutation' ( = " << *currentPermutation_ << ") >= "
          << "'maxCurrentPermutation' ( = " << maxCurrentPermutation_ << ')';
      return objects[permutations_[*currentPermutation_][index]];
    }
    return objects[defaultPermutation[index]];
  }

  const T * const
  begin() const
  {
    return &objects[defaultPermutation[0]];
  }

  const T * const
  end() const
  {
    return &objects[defaultPermutation[NofObjects - 1] + 1];
  }

  void
  reset()
  {
    currentPermutation_ = 0;
    maxCurrentPermutation_ = 0;
  }

  T objects[NofObjects];
  unsigned defaultPermutation[NofObjects];
  unsigned * currentPermutation_;
  unsigned maxCurrentPermutation_;
  std::vector<std::vector<unsigned>> permutations_;
};

#endif // INDEXWRAPPER_H
