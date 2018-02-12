// -*- C++ -*-
#include "Rivet/Config/RivetCommon.hh"
#include <random>
#if defined(_OPENMP)
#include "omp.h"
#endif

namespace Rivet {


  // Return a thread-safe random number generator
  mt19937& rng() {
    #if defined(_OPENMP)
    static map<int,mt19937> gens;
    const int nthread = omp_get_thread_num();
    if (gens.find(nthread) == gens.end()) {
      seed_seq seq{1,2,3,4,5};
      vector<uint32_t> seeds(nthread+1);
      seq.generate(seeds.begin(), seeds.end());
      gens[nthread] = mt19937(seeds[nthread]);
      // cout << "Thread " << nthread+1 << ", seed=" << seeds[nthread] << " (" << gens.size() << " RNGs)" << endl;
    }
    mt19937& g = gens[nthread];
    #else
    static mt19937 g(12345);
    #endif
    return g;
  }


  // Return a uniformly sampled random number between 0 and 1
  double rand01() {
    // return rand() / (double)RAND_MAX;
    return generate_canonical<double, 32>(rng()); ///< @todo What's the "correct" number of bits of randomness?
  }


  // Return a Gaussian/normal sampled random number with the given mean and width
  double randnorm(double loc, double scale) {
    normal_distribution<> d(loc, scale);
    return d(rng());
  }


  // Return a log-normal sampled random number
  double randlognorm(double loc, double scale) {
    lognormal_distribution<> d(loc, scale);
    return d(rng());
  }


}
