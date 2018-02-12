#ifndef RIVET_ParticleBase_HH
#define RIVET_ParticleBase_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Jet.fhh"
#include "Rivet/Tools/Cuts.fhh"
#include "Rivet/Math/Vectors.hh"

namespace Rivet {


  /// @brief Base class for particle-like things like Particle and Jet
  class ParticleBase {
  public:

    /// Default constructor
    ParticleBase() { }

    /// Virtual destructor
    virtual ~ParticleBase() { }


    // /// @name Constituent accessors
    // //@{

    // /// @todo Can't do this because a) ParticleBase is pure-virtual; b) inheritance causality for Particle... urk
    // virtual const vector<ParticleBase>& constituents() const = 0;
    // virtual const vector<ParticleBase>& rawConstituents() const = 0;

    // //@}


    /// @name Effective momentum accessors
    //@{

    /// Get equivalent single momentum four-vector (const).
    virtual const FourMomentum& momentum() const = 0;
    /// Get equivalent single momentum four-vector (const) (alias).
    const FourMomentum& mom() const { return momentum(); };

    /// Cast operator for conversion to FourMomentum
    operator const FourMomentum& () const { return momentum(); }

    //@}


    /// @name Convenience access to the effective 4-vector properties
    //@{

    /// Get the energy directly.
    double E() const { return momentum().E(); }
    /// Get the energy directly (alias).
    double energy() const { return momentum().E(); }

    /// Get the energy-squared.
    double E2() const { return momentum().E2(); }
    /// Get the energy-squared (alias).
    double energy2() const { return momentum().E2(); }

    /// Get the \f$ p_T \f$ directly.
    double pt() const { return momentum().pt(); }
    /// Get the \f$ p_T \f$ directly (alias).
    double pT() const { return pt(); }
    /// Get the \f$ p_T \f$ directly (alias).
    double perp() const { return pt(); }

    /// Get the \f$ p_T^2 \f$ directly.
    double pt2() const { return momentum().pt2(); }
    /// Get the \f$ p_T^2 \f$ directly (alias).
    double pT2() const { return pt2(); }
    /// Get the \f$ p_T^2 \f$ directly (alias).
    double perp2() const { return pt2(); }

    /// Get the \f$ E_T \f$ directly.
    double Et() const { return momentum().Et(); }
    /// Get the \f$ E_T^2 \f$ directly.
    double Et2() const { return momentum().Et2(); }

    /// Get the mass directly.
    double mass() const { return momentum().mass(); }
    /// Get the mass**2 directly.
    double mass2() const { return momentum().mass2(); }

    /// Get the \f$ \eta \f$ directly.
    double pseudorapidity() const { return momentum().eta(); }
    /// Get the \f$ \eta \f$ directly (alias).
    double eta() const { return momentum().eta(); }
    /// Get the \f$ |\eta| \f$ directly.
    double abspseudorapidity() const { return momentum().abspseudorapidity(); }
    /// Get the \f$ |\eta| \f$ directly (alias).
    double abseta() const { return momentum().abseta(); }

    /// Get the \f$ y \f$ directly.
    double rapidity() const { return momentum().rapidity(); }
    /// Get the \f$ y \f$ directly (alias).
    double rap() const { return momentum().rapidity(); }
    /// Get the \f$ |y| \f$ directly.
    double absrapidity() const { return momentum().absrapidity(); }
    /// Get the \f$ |y| \f$ directly (alias).
    double absrap() const { return momentum().absrap(); }

    /// Azimuthal angle \f$ \phi \f$.
    double azimuthalAngle(const PhiMapping mapping=ZERO_2PI) const { return momentum().azimuthalAngle(mapping); }
    /// Get the \f$ \phi \f$ directly.
    double phi(const PhiMapping mapping=ZERO_2PI) const { return momentum().phi(mapping); }

    /// Get the 3-momentum directly.
    Vector3 p3() const { return momentum().vector3(); }
    /// Get the 3-momentum magnitude directly.
    double p() const { return momentum().p(); }
    /// Get the 3-momentum magnitude-squared directly.
    double p2() const { return momentum().p2(); }

    /// Get the transverse 3-momentum directly.
    Vector3 ptvec() const { return momentum().ptvec(); }
    /// Get the transverse 3-momentum directly.
    Vector3 pTvec() const { return momentum().pTvec(); }

    /// x component of momentum.
    double px() const { return momentum().x(); }
    /// y component of momentum.
    double py() const { return momentum().y(); }
    /// z component of momentum.
    double pz() const { return momentum().z(); }

    /// x component of momentum, squared.
    double px2() const { return momentum().x2(); }
    /// y component of momentum, squared.
    double py2() const { return momentum().y2(); }
    /// z component of momentum, squared.
    double pz2() const { return momentum().z2(); }

    /// Angle subtended by the 3-vector and the z-axis.
    double polarAngle() const { return momentum().polarAngle(); }
    /// Synonym for polarAngle.
    double theta() const { return momentum().theta(); }

    /// Angle between this vector and another
    double angle(const ParticleBase& v) const { return momentum().angle(v.momentum()); }
    /// Angle between this vector and another
    double angle(const FourVector& v) const { return momentum().angle(v); }
    /// Angle between this vector and another (3-vector)
    double angle(const Vector3& v3) const { return momentum().angle(v3); }

    //@}

  };


  /// @name deltaR, deltaEta, deltaPhi functions specifically for Particle/Jet arguments
  //@{

  inline double deltaR(const ParticleBase& p1, const ParticleBase& p2,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(p1.momentum(), p2.momentum(), scheme);
  }

  inline double deltaR(const ParticleBase& p, const FourMomentum& v,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(p.momentum(), v, scheme);
  }

  inline double deltaR(const ParticleBase& p, const FourVector& v,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(p.momentum(), v, scheme);
  }

  inline double deltaR(const ParticleBase& p, const Vector3& v) {
    return deltaR(p.momentum(), v);
  }

  inline double deltaR(const ParticleBase& p, double eta, double phi) {
    return deltaR(p.momentum(), eta, phi);
  }

  inline double deltaR(const FourMomentum& v, const ParticleBase& p,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(v, p.momentum(), scheme);
  }

  inline double deltaR(const FourVector& v, const ParticleBase& p,
                       RapScheme scheme = PSEUDORAPIDITY) {
    return deltaR(v, p.momentum(), scheme);
  }

  inline double deltaR(const Vector3& v, const ParticleBase& p) {
    return deltaR(v, p.momentum());
  }

  inline double deltaR(double eta, double phi, const ParticleBase& p) {
    return deltaR(eta, phi, p.momentum());
  }


  inline double deltaPhi(const ParticleBase& p1, const ParticleBase& p2) {
    return deltaPhi(p1.momentum(), p2.momentum());
  }

  inline double deltaPhi(const ParticleBase& p, const FourMomentum& v) {
    return deltaPhi(p.momentum(), v);
  }

  inline double deltaPhi(const ParticleBase& p, const FourVector& v) {
    return deltaPhi(p.momentum(), v);
  }

  inline double deltaPhi(const ParticleBase& p, const Vector3& v) {
    return deltaPhi(p.momentum(), v);
  }

  inline double deltaPhi(const ParticleBase& p, double phi) {
    return deltaPhi(p.momentum(), phi);
  }

  inline double deltaPhi(const FourMomentum& v, const ParticleBase& p) {
    return deltaPhi(v, p.momentum());
  }

  inline double deltaPhi(const FourVector& v, const ParticleBase& p) {
    return deltaPhi(v, p.momentum());
  }

  inline double deltaPhi(const Vector3& v, const ParticleBase& p) {
    return deltaPhi(v, p.momentum());
  }

  inline double deltaPhi(double phi, const ParticleBase& p) {
    return deltaPhi(phi, p.momentum());
  }


  inline double deltaEta(const ParticleBase& p1, const ParticleBase& p2) {
    return deltaEta(p1.momentum(), p2.momentum());
  }

  inline double deltaEta(const ParticleBase& p, const FourMomentum& v) {
    return deltaEta(p.momentum(), v);
  }

  inline double deltaEta(const ParticleBase& p, const FourVector& v) {
    return deltaEta(p.momentum(), v);
  }

  inline double deltaEta(const ParticleBase& p, const Vector3& v) {
    return deltaEta(p.momentum(), v);
  }

  inline double deltaEta(const ParticleBase& p, double eta) {
    return deltaEta(p.momentum(), eta);
  }

  inline double deltaEta(const FourMomentum& v, const ParticleBase& p) {
    return deltaEta(v, p.momentum());
  }

  inline double deltaEta(const FourVector& v, const ParticleBase& p) {
    return deltaEta(v, p.momentum());
  }

  inline double deltaEta(const Vector3& v, const ParticleBase& p) {
    return deltaEta(v, p.momentum());
  }

  inline double deltaEta(double eta, const ParticleBase& p) {
    return deltaEta(eta, p.momentum());
  }


  inline double deltaRap(const ParticleBase& p1, const ParticleBase& p2) {
    return deltaRap(p1.momentum(), p2.momentum());
  }

  inline double deltaRap(const ParticleBase& p, const FourMomentum& v) {
    return deltaRap(p.momentum(), v);
  }

  inline double deltaRap(const ParticleBase& p, double y) {
    return deltaRap(p.momentum(), y);
  }

  inline double deltaRap(const FourMomentum& v, const ParticleBase& p) {
    return deltaRap(v, p.momentum());
  }

  inline double deltaRap(double y, const ParticleBase& p) {
    return deltaRap(y, p.momentum());
  }

  //@}


}

#endif
