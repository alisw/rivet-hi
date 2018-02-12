#include <iostream>
#include <limits>
#include <cassert>

#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"

using namespace std;
using namespace Rivet;

int main() {

  // Angle tests
  assert(fuzzyEquals(angle(FourMomentum(1,0,0,1), FourMomentum(1,0,0,1))/M_PI, 0.0));
  assert(fuzzyEquals(angle(FourMomentum(1,0,0,1), FourMomentum(1,0,1,0))/M_PI, 0.5));
  assert(fuzzyEquals(angle(FourMomentum(1,0,0,1), FourMomentum(1,0,0,-1))/M_PI, 1.0));
  // Test with vectors of different magnitude
  assert(fuzzyEquals(angle(FourMomentum(3,0,0,3), FourMomentum(1,0,0,1))/M_PI, 0.0));
  assert(fuzzyEquals(angle(FourMomentum(5,0,0,5), FourMomentum(1,0,1,0))/M_PI, 0.5));
  assert(fuzzyEquals(angle(FourMomentum(7,0,0,7), FourMomentum(1,0,0,-1))/M_PI, 1.0));

  ////////////

  assert(linspace(50, 0, 10).size() == 51);
  assert(logspace(50, 0.000001, 1.0).back() == 1.0);

  assert(binIndex(3, vector<double>{0, 1, 2, 3, 4, 5}) == 3);
  assert(binIndex(2.99, vector<double>{0, 1, 2, 3, 4, 5}) == 2);
  assert(binIndex(-4, vector<double>{0, 1, 2, 3, 4, 5}) == -1);
  assert(binIndex(5.0, vector<double>{0, 1, 2, 3, 4, 5}) == -1);
  assert(binIndex(5.1, vector<double>{0, 1, 2, 3, 4, 5}, true) == 5);

  ////////////

  inRange(1, 0, 2);
  inRange(1, 0.0, 2);
  inRange(1, 0, 2.0);
  inRange(1, 0.0, 2.0);
  inRange(1.0, 0, 2);
  inRange(1.0, 0.0, 2);
  inRange(1.0, 0, 2.0);
  inRange(1.0, 0.0, 2.0);

  assert(isZero(1e-15));
  assert(fuzzyEquals(1e-15, 0.0));
  assert(fuzzyEquals(2.0, 2.0));
  assert(!fuzzyEquals(1, 2.0));
  assert(fuzzyGtrEquals(1.0, 1.0 - 1e-15));

  return EXIT_SUCCESS;
}
