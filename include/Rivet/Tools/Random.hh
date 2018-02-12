// -*- C++ -*-
#ifndef RIVET_Random_HH
#define RIVET_Random_HH

#include <random>
// #if defined(_OPENMP)
// #include "omp.h"
// #endif

namespace Rivet {


  /// Return a thread-safe random number generator (mainly for internal use)
  mt19937& rng();

  /// Return a uniformly sampled random number between 0 and 1
  double rand01();

  /// Return a Gaussian/normal sampled random number with the given mean and width
  double randnorm(double loc, double scale);

  /// Return a log-normal sampled random number
  double randlognorm(double loc, double scale);


}

#endif
