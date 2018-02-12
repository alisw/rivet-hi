// -*- C++ -*-
#ifndef RIVET_MomentumSmearingFunctions_HH
#define RIVET_MomentumSmearingFunctions_HH

#include "Rivet/Math/Vector4.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {


  /// @name FourMomentum efficiency and smearing functions
  //@{

  /// @name Typedef for FourMomentum smearing functions/functors
  typedef std::function<FourMomentum(const FourMomentum&)> P4SmearFn;

  /// @name Typedef for FourMomentum efficiency functions/functors
  typedef std::function<double(const FourMomentum&)> P4EffFn;


  /// Take a FourMomentum and return 0
  inline double P4_EFF_ZERO(const FourMomentum& ) { return 0; }
  /// @deprecated Alias for P4_EFF_ZERO
  inline double P4_FN0(const FourMomentum& ) { return 0; }

  /// Take a FourMomentum and return 1
  inline double P4_EFF_ONE(const FourMomentum& ) { return 1; }
  /// @deprecated Alias for P4_EFF_ONE
  inline double P4_FN1(const FourMomentum& ) { return 1; }

  /// Take a FourMomentum and return a constant number
  struct P4_EFF_CONST {
    P4_EFF_CONST(double x) : _x(x) {}
    double operator () (const FourMomentum& )  const { return _x; }
    double _x;
  };


  /// Take a FourMomentum and return it unmodified
  inline FourMomentum P4_SMEAR_IDENTITY(const FourMomentum& p) { return p; }

  /// Smear a FourMomentum's energy using a Gaussian of absolute width @a resolution
  /// @todo Also make jet versions that update/smear constituents?
  inline FourMomentum P4_SMEAR_E_GAUSS(const FourMomentum& p, double resolution) {
    const double mass = p.mass2() > 0 ? p.mass() : 0; //< numerical carefulness...
    const double smeared_E = max(randnorm(p.E(), resolution), mass); //< can't let the energy go below the mass!
    return FourMomentum::mkEtaPhiME(p.eta(), p.phi(), mass, smeared_E);
  }

  /// Smear a FourMomentum's transverse momentum using a Gaussian of absolute width @a resolution
  inline FourMomentum P4_SMEAR_PT_GAUSS(const FourMomentum& p, double resolution) {
    const double smeared_pt = max(randnorm(p.pT(), resolution), 0.);
    const double mass = p.mass2() > 0 ? p.mass() : 0; //< numerical carefulness...
    return FourMomentum::mkEtaPhiMPt(p.eta(), p.phi(), mass, smeared_pt);
  }

  /// Smear a FourMomentum's mass using a Gaussian of absolute width @a resolution
  inline FourMomentum P4_SMEAR_MASS_GAUSS(const FourMomentum& p, double resolution) {
    const double smeared_mass = max(randnorm(p.mass(), resolution), 0.);
    return FourMomentum::mkEtaPhiMPt(p.eta(), p.phi(), smeared_mass, p.pT());
  }

  //@}



  /// @name FourMomentum efficiency and smearing functions
  //@{

  /// Take a Vector3 and return 0
  inline double P3_EFF_ZERO(const Vector3& p) { return 0; }
  /// @deprecated Alias for P3_EFF_ZERO
  inline double P3_FN0(const Vector3& p) { return 0; }

  /// Take a Vector3 and return 1
  inline double P3_EFF_ONE(const Vector3& p) { return 1; }
  /// @deprecated Alias for P3_EFF_ONE
  inline double P3_FN1(const Vector3& p) { return 1; }

  /// Take a Vector3 and return a constant number
  struct P3_EFF_CONST {
    P3_EFF_CONST(double x) : _x(x) {}
    double operator () (const Vector3& )  const { return _x; }
    double _x;
  };


  /// Take a Vector3 and return it unmodified
  inline Vector3 P3_SMEAR_IDENTITY(const Vector3& p) { return p; }

  /// Smear a Vector3's length using a Gaussian of absolute width @a resolution
  inline Vector3 P3_SMEAR_LEN_GAUSS(const Vector3& p, double resolution) {
    const double smeared_mod = max(randnorm(p.mod(), resolution), 0.); //< can't let the energy go below the mass!
    return smeared_mod * p.unit();
  }

  //@}


}

#endif
