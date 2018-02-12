// -*- C++ -*-
#ifndef RIVET_MT2_HH
#define RIVET_MT2_HH

#include "Rivet/Math/Vector4.hh"

namespace Rivet {


  /// @brief Compute asymm mT2**2 using the bisection method
  ///
  /// If the second invisible mass is not given, symm mT2**2 will be calculated.
  ///
  /// @note Cheng/Han arXiv:0810.5178, Lester arXiv:1411.4312
  double mT2Sq(const FourMomentum& a, const FourMomentum& b, const Vector3& ptmiss,
               double invisiblesMass, double invisiblesMass2=-1);

  /// Override for mT2Sq with FourMomentum ptmiss
  inline double mT2Sq(const FourMomentum& a, const FourMomentum& b, const FourMomentum& ptmiss,
                      double invisiblesMass, double invisiblesMass2) {
    return mT2Sq(a, b, ptmiss.perpVec(), invisiblesMass, invisiblesMass2=-1);
  }


  /// @brief Compute asymm mT2 using the bisection method
  ///
  /// If the second invisible mass is not given, symm mT2 will be calculated.
  ///
  /// @note Cheng/Han arXiv:0810.5178, Lester arXiv:1411.4312
  inline double mT2(const FourMomentum& a, const FourMomentum& b, const Vector3& ptmiss,
                    double invisiblesMass, double invisiblesMass2) {
    const double mt2sq = mT2Sq(a, b, ptmiss, invisiblesMass, invisiblesMass2=-1);
    return mt2sq >= 0 ? sqrt(mt2sq) : -1;
  }

  /// Override for mT2 with FourMomentum ptmiss
  inline double mT2(const FourMomentum& a, const FourMomentum& b, const FourMomentum& ptmiss,
                    double invisiblesMass, double invisiblesMass2=-1) {
    return mT2(a, b, ptmiss.perpVec(), invisiblesMass, invisiblesMass2);
  }


}

#endif
