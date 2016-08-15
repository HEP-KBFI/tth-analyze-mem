#ifndef INDEXWRAPPER_H
#define INDEXWRAPPER_H

#include <numeric> // std::iota()
#include <cassert> // assert()

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
    if(! permutations.size()) assert(0);
    for(const std::vector<unsigned> & permutation: permutations)
      if(permutation.size() != NofObjects) assert(0);
    assert(currentPermutation);
    assert(maxCurrentPermutation);

    permutations_ = permutations;
    currentPermutation_ = currentPermutation;
    maxCurrentPermutation_ = maxCurrentPermutation;
  }

  T &
  operator[](const int index)
  {
    if(currentPermutation_)
    {
      if(*currentPermutation_ >= maxCurrentPermutation_) assert(0);
      return objects[permutations_[*currentPermutation_][index]];
    }
    return objects[defaultPermutation[index]];
  }

  const T &
  operator[](const int index) const
  {
    if(currentPermutation_)
    {
      if(*currentPermutation_ >= maxCurrentPermutation_) assert(0);
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

  T objects[NofObjects];
  unsigned defaultPermutation[NofObjects];
  unsigned * currentPermutation_;
  unsigned maxCurrentPermutation_;
  std::vector<std::vector<unsigned>> permutations_;
};

#endif // INDEXWRAPPER_H
