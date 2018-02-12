#ifndef RIVET_PARTICLEBASEUTILS_HH
#define RIVET_PARTICLEBASEUTILS_HH

#include "Rivet/ParticleBase.hh"

namespace Rivet {



  /// @name ParticleBase classifier -> bool functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to any() or all() e.g. any(jets, DeltaRLess(electron, 0.4))
  //@{

  /// std::function instantiation for functors taking a ParticleBase and returning a bool
  using ParticleBaseSelector = function<bool(const ParticleBase&)>;
  /// std::function instantiation for functors taking two ParticleBase and returning a bool
  using ParticleBaseSorter = function<bool(const ParticleBase&, const ParticleBase&)>;


  /// Base type for Particle -> bool functors
  struct BoolParticleBaseFunctor {
    virtual bool operator()(const ParticleBase& p) const = 0;
  };


  /// Transverse momentum greater-than functor
  struct PtGtr : public BoolParticleBaseFunctor {
    PtGtr(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() > ptcut; }
    double ptcut;
  };
  using pTGtr = PtGtr;
  using ptGtr = PtGtr;

  /// Transverse momentum less-than functor
  struct PtLess : public BoolParticleBaseFunctor {
    PtLess(double pt) : ptcut(pt) { }
    bool operator()(const ParticleBase& p) const { return p.pT() < ptcut; }
    double ptcut;
  };
  using pTLess = PtLess;
  using ptLess = PtLess;

  /// Transverse momentum in-range functor
  struct PtInRange : public BoolParticleBaseFunctor {
    PtInRange(pair<double, double> ptcuts) : ptcut(ptcuts) { }
    PtInRange(double ptlow, double pthigh) : PtInRange(make_pair(ptlow, pthigh)) { }
    bool operator()(const ParticleBase& p) const { return p.pT() >= ptcut.first && p.pT() < ptcut.second; }
    pair<double,double> ptcut;
  };
  using pTInRange = PtInRange;
  using ptInRange = PtInRange;


  /// Pseudorapidity greater-than functor
  struct EtaGtr : public BoolParticleBaseFunctor {
    EtaGtr(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() > etacut; }
    double etacut;
  };
  using etaGtr = EtaGtr;

  /// Pseudorapidity less-than functor
  struct EtaLess : public BoolParticleBaseFunctor {
    EtaLess(double eta) : etacut(eta) { }
    bool operator()(const ParticleBase& p) const { return p.eta() < etacut; }
    double etacut;
  };
  using etaLess = EtaLess;

  /// Pseudorapidity in-range functor
  struct EtaInRange : public BoolParticleBaseFunctor {
    EtaInRange(pair<double, double> etacuts) : etacut(etacuts) { }
    EtaInRange(double etalow, double etahigh) : EtaInRange(make_pair(etalow, etahigh)) { }
    bool operator()(const ParticleBase& p) const { return p.eta() >= etacut.first && p.eta() < etacut.second; }
    pair<double,double> etacut;
  };
  using etaInRange = EtaInRange;


  /// Abs pseudorapidity greater-than functor
  struct AbsEtaGtr : public BoolParticleBaseFunctor {
    AbsEtaGtr(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() > absetacut; }
    double absetacut;
  };
  using absEtaGtr = AbsEtaGtr;
  using absetaGtr = AbsEtaGtr;

  /// Abs pseudorapidity momentum less-than functor
  struct AbsEtaLess : public BoolParticleBaseFunctor {
    AbsEtaLess(double abseta) : absetacut(abseta) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() < absetacut; }
    double absetacut;
  };
  using absEtaLess = AbsEtaLess;
  using absetaLess = AbsEtaLess;

  /// Abs pseudorapidity in-range functor
  struct AbsEtaInRange : public BoolParticleBaseFunctor {
    AbsEtaInRange(const pair<double, double>& absetacuts) : absetacut(absetacuts) { }
    AbsEtaInRange(double absetalow, double absetahigh) : AbsEtaInRange(make_pair(absetalow, absetahigh)) { }
    bool operator()(const ParticleBase& p) const { return p.abseta() >= absetacut.first && p.abseta() < absetacut.second; }
    pair<double,double> absetacut;
  };
  using absEtaInRange = AbsEtaInRange;
  using absetaInRange = AbsEtaInRange;


  /// Rapidity greater-than functor
  struct RapGtr : public BoolParticleBaseFunctor {
    RapGtr(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() > rapcut; }
    double rapcut;
  };
  using rapGtr = RapGtr;

  /// Rapidity momentum less-than functor
  struct RapLess : public BoolParticleBaseFunctor {
    RapLess(double rap) : rapcut(rap) { }
    bool operator()(const ParticleBase& p) const { return p.rap() < rapcut; }
    double rapcut;
  };
  using rapLess = RapLess;

  /// Rapidity in-range functor
  struct RapInRange : public BoolParticleBaseFunctor {
    RapInRange(const pair<double, double>& rapcuts) : rapcut(rapcuts) { }
    RapInRange(double raplow, double raphigh) : RapInRange(make_pair(raplow, raphigh)) { }
    bool operator()(const ParticleBase& p) const { return p.rap() >= rapcut.first && p.rap() < rapcut.second; }
    pair<double,double> rapcut;
  };
  using rapInRange = RapInRange;


  /// Abs rapidity greater-than functor
  struct AbsRapGtr : public BoolParticleBaseFunctor {
    AbsRapGtr(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() > absrapcut; }
    double absrapcut;
  };
  using absRapGtr = AbsRapGtr;
  using absrapGtr = AbsRapGtr;

  /// Abs rapidity momentum less-than functor
  struct AbsRapLess : public BoolParticleBaseFunctor {
    AbsRapLess(double absrap) : absrapcut(absrap) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() < absrapcut; }
    double absrapcut;
  };
  using absRapLess = AbsRapLess;
  using absrapLess = AbsRapLess;

  /// Abs rapidity in-range functor
  struct AbsRapInRange : public BoolParticleBaseFunctor {
    AbsRapInRange(const pair<double, double>& absrapcuts) : absrapcut(absrapcuts) { }
    AbsRapInRange(double absraplow, double absraphigh) : AbsRapInRange(make_pair(absraplow, absraphigh)) { }
    bool operator()(const ParticleBase& p) const { return p.absrap() >= absrapcut.first && p.absrap() < absrapcut.second; }
    pair<double,double> absrapcut;
  };
  using absRapInRange = AbsRapInRange;
  using absrapInRange = AbsRapInRange;



  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) greater-than functor
  struct DeltaRGtr : public BoolParticleBaseFunctor {
    DeltaRGtr(const ParticleBase& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRGtr(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRGtr(const Vector3& vec, double dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) > drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using deltaRGtr = DeltaRGtr;

  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) less-than functor
  struct DeltaRLess : public BoolParticleBaseFunctor {
    DeltaRLess(const ParticleBase& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRLess(const FourMomentum& vec, double dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRLess(const Vector3& vec, double dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    bool operator()(const ParticleBase& p) const { return deltaR(p, refvec, rapscheme) < drcut; }
    FourMomentum refvec;
    double drcut;
    RapScheme rapscheme;
  };
  using deltaRLess = DeltaRLess;

  /// @f$ \Delta R @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaRInRange : public BoolParticleBaseFunctor {
    DeltaRInRange(const ParticleBase& vec, const pair<double,double>& dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec.mom()), drcut(dr), rapscheme(scheme) { }
    DeltaRInRange(const ParticleBase& vec, double drmin, double drmax, RapScheme scheme=PSEUDORAPIDITY)
      : DeltaRInRange(vec, make_pair(drmin, drmax), scheme) { }
    DeltaRInRange(const FourMomentum& vec, const pair<double,double>& dr, RapScheme scheme=PSEUDORAPIDITY)
      : refvec(vec), drcut(dr), rapscheme(scheme) { }
    DeltaRInRange(const FourMomentum& vec, double drmin, double drmax, RapScheme scheme=PSEUDORAPIDITY)
      : DeltaRInRange(vec, make_pair(drmin, drmax), scheme) { }
    DeltaRInRange(const Vector3& vec, const pair<double,double>& dr)
      : drcut(dr), rapscheme(PSEUDORAPIDITY) { refvec.setPx(vec.x()); refvec.setPy(vec.y()); refvec.setPz(vec.z()); }
    DeltaRInRange(const Vector3& vec, double drmin, double drmax)
      : DeltaRInRange(vec, make_pair(drmin, drmax)) { }
    bool operator()(const ParticleBase& p) const {
      const double dR = deltaR(p, refvec, rapscheme);
      return dR >= drcut.first && dR < drcut.second;
    }
    FourMomentum refvec;
    pair<double,double> drcut;
    RapScheme rapscheme;
  };
  using deltaRInRange = DeltaRInRange;


  /// @f$ |\Delta \phi| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaPhiGtr : public BoolParticleBaseFunctor {
    DeltaPhiGtr(const ParticleBase& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiGtr(const FourMomentum& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiGtr(const Vector3& vec, double dphi)
      : refvec(vec), dphicut(dphi) { }
    bool operator()(const ParticleBase& p) const { return deltaPhi(p, refvec) > dphicut; }
    Vector3 refvec;
    double dphicut;
  };
  using deltaPhiGtr = DeltaPhiGtr;

  /// @f$ |\Delta \phi| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaPhiLess : public BoolParticleBaseFunctor {
    DeltaPhiLess(const ParticleBase& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiLess(const FourMomentum& vec, double dphi)
      : refvec(vec.p3()), dphicut(dphi) { }
    DeltaPhiLess(const Vector3& vec, double dphi)
      : refvec(vec), dphicut(dphi) { }
    bool operator()(const ParticleBase& p) const { return deltaPhi(p, refvec) < dphicut; }
    Vector3 refvec;
    double dphicut;
  };
  using deltaPhiLess = DeltaPhiLess;

  /// @f$ \Delta \phi @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaPhiInRange : public BoolParticleBaseFunctor {
    DeltaPhiInRange(const ParticleBase& vec, const pair<double,double>& dphi)
      : refvec(vec.mom()), dphicut(dphi) { }
    DeltaPhiInRange(const ParticleBase& vec, double dphimin, double dphimax)
      : DeltaPhiInRange(vec, make_pair(dphimin, dphimax)) { }
    DeltaPhiInRange(const FourMomentum& vec, const pair<double,double>& dphi)
      : refvec(vec), dphicut(dphi) { }
    DeltaPhiInRange(const FourMomentum& vec, double dphimin, double dphimax)
      : DeltaPhiInRange(vec, make_pair(dphimin, dphimax)) { }
    DeltaPhiInRange(const Vector3& vec, const pair<double,double>& dphi)
      : refvec(vec), dphicut(dphi) { }
    DeltaPhiInRange(const Vector3& vec, double dphimin, double dphimax)
      : DeltaPhiInRange(vec, make_pair(dphimin, dphimax)) { }
    bool operator()(const ParticleBase& p) const {
      const double dphi = deltaPhi(p, refvec);
      return dphi >= dphicut.first && dphi < dphicut.second;
    }
    Vector3 refvec;
    pair<double,double> dphicut;
  };
  using deltaPhiInRange = DeltaPhiInRange;


  /// @f$ |\Delta \eta| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaEtaGtr : public BoolParticleBaseFunctor {
    DeltaEtaGtr(const ParticleBase& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaGtr(const FourMomentum& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaGtr(const Vector3& vec, double deta)
      : refvec(vec), detacut(deta) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaEta(p, refvec)) > detacut; }
    Vector3 refvec;
    double detacut;
  };
  using deltaEtaGtr = DeltaEtaGtr;

  /// @f$ |\Delta \eta| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaEtaLess : public BoolParticleBaseFunctor {
    DeltaEtaLess(const ParticleBase& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaLess(const FourMomentum& vec, double deta)
      : refvec(vec.p3()), detacut(deta) { }
    DeltaEtaLess(const Vector3& vec, double deta)
      : refvec(vec), detacut(deta) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaEta(p, refvec)) < detacut; }
    Vector3 refvec;
    double detacut;
  };
  using deltaEtaLess = DeltaEtaLess;

  /// @f$ \Delta \eta @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaEtaInRange : public BoolParticleBaseFunctor {
    DeltaEtaInRange(const ParticleBase& vec, const pair<double,double>& deta)
      : refvec(vec.mom()), detacut(deta) { }
    DeltaEtaInRange(const ParticleBase& vec, double detamin, double detamax)
      : DeltaEtaInRange(vec, make_pair(detamin, detamax)) { }
    DeltaEtaInRange(const FourMomentum& vec, const pair<double,double>& deta)
      : refvec(vec), detacut(deta) { }
    DeltaEtaInRange(const FourMomentum& vec, double detamin, double detamax)
      : DeltaEtaInRange(vec, make_pair(detamin, detamax)) { }
    DeltaEtaInRange(const Vector3& vec, const pair<double,double>& deta)
      : refvec(vec), detacut(deta) { }
    DeltaEtaInRange(const Vector3& vec, double detamin, double detamax)
      : DeltaEtaInRange(vec, make_pair(detamin, detamax)) { }
    bool operator()(const ParticleBase& p) const {
      const double deta = deltaEta(p, refvec);
      return deta >= detacut.first && deta < detacut.second;
    }
    Vector3 refvec;
    pair<double,double> detacut;
  };
  using deltaEtaInRange = DeltaEtaInRange;


  /// @f$ |\Delta y| @f$ (with respect to another momentum, @a vec) greater-than functor
  struct DeltaRapGtr : public BoolParticleBaseFunctor {
    DeltaRapGtr(const ParticleBase& vec, double drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapGtr(const FourMomentum& vec, double drap)
      : refvec(vec), drapcut(drap) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaRap(p, refvec)) > drapcut; }
    FourMomentum refvec;
    double drapcut;
  };
  using deltaRapGtr = DeltaRapGtr;

  /// @f$ |\Delta y| @f$ (with respect to another momentum, @a vec) less-than functor
  struct DeltaRapLess : public BoolParticleBaseFunctor {
    DeltaRapLess(const ParticleBase& vec, double drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapLess(const FourMomentum& vec, double drap)
      : refvec(vec), drapcut(drap) { }
    bool operator()(const ParticleBase& p) const { return std::abs(deltaRap(p, refvec)) < drapcut; }
    FourMomentum refvec;
    double drapcut;
  };
  using deltaRapLess = DeltaRapLess;

  /// @f$ \Delta y @f$ (with respect to another 4-momentum, @a vec) in-range functor
  struct DeltaRapInRange : public BoolParticleBaseFunctor {
    DeltaRapInRange(const ParticleBase& vec, const pair<double,double>& drap)
      : refvec(vec.mom()), drapcut(drap) { }
    DeltaRapInRange(const ParticleBase& vec, double drapmin, double drapmax)
      : DeltaRapInRange(vec, make_pair(drapmin, drapmax)) { }
    DeltaRapInRange(const FourMomentum& vec, const pair<double,double>& drap)
      : refvec(vec), drapcut(drap) { }
    DeltaRapInRange(const FourMomentum& vec, double drapmin, double drapmax)
      : DeltaRapInRange(vec, make_pair(drapmin, drapmax)) { }
    bool operator()(const ParticleBase& p) const {
      const double drap = deltaRap(p, refvec);
      return drap >= drapcut.first && drap < drapcut.second;
    }
    FourMomentum refvec;
    pair<double,double> drapcut;
  };
  using deltaRapInRange = DeltaRapInRange;

  //@}


  /// @name ParticleBase comparison -> double functors
  /// @todo Move to FourMomentum functions
  ///
  /// To be passed to transform()any(jets, DeltaRLess(electron, 0.4))
  //@{

  /// Base type for Particle -> double functors
  struct DoubleParticleBaseFunctor {
    virtual double operator()(const ParticleBase& p) const = 0;
  };

  /// Calculator of @f$ \Delta R @f$ with respect to a given momentum
  struct DeltaRWRT : public DoubleParticleBaseFunctor {
    DeltaRWRT(const ParticleBase& pb, RapScheme scheme=PSEUDORAPIDITY) : p(pb.mom()), rapscheme(scheme) {}
    DeltaRWRT(const FourMomentum& p4, RapScheme scheme=PSEUDORAPIDITY) : p(p4), rapscheme(scheme) {}
    DeltaRWRT(const Vector3& p3) : p(p3.mod(), p3.x(), p3.y(), p3.z()), rapscheme(PSEUDORAPIDITY) {}
    double operator()(const ParticleBase& pb) const { return deltaR(p, pb, rapscheme); }
    double operator()(const FourMomentum& p4) const { return deltaR(p, p4, rapscheme); }
    double operator()(const Vector3& p3) const { return deltaR(p, p3); }
    const FourMomentum p;
    RapScheme rapscheme;
  };
  using deltaRWRT = DeltaRWRT;

  /// Calculator of @f$ \Delta \phi @f$ with respect to a given momentum
  struct DeltaPhiWRT : public DoubleParticleBaseFunctor {
    DeltaPhiWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    DeltaPhiWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    DeltaPhiWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaPhi(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaPhi(p, p4); }
    double operator()(const Vector3& p3) const { return deltaPhi(p, p3); }
    const Vector3 p;
  };
  using deltaPhiWRT = DeltaPhiWRT;

  /// Calculator of @f$ \Delta \eta @f$ with respect to a given momentum
  struct DeltaEtaWRT : public DoubleParticleBaseFunctor {
    DeltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    DeltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    DeltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return deltaEta(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaEta(p, p4); }
    double operator()(const Vector3& p3) const { return deltaEta(p, p3); }
    const Vector3 p;
  };
  using deltaEtaWRT = DeltaEtaWRT;

  /// Calculator of @f$ |\Delta \eta| @f$ with respect to a given momentum
  struct AbsDeltaEtaWRT : public DoubleParticleBaseFunctor {
    AbsDeltaEtaWRT(const ParticleBase& pb) : p(pb.mom().vector3()) {}
    AbsDeltaEtaWRT(const FourMomentum& p4) : p(p4.vector3()) {}
    AbsDeltaEtaWRT(const Vector3& p3) : p(p3) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaEta(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaEta(p, p4)); }
    double operator()(const Vector3& p3) const { return fabs(deltaEta(p, p3)); }
    const Vector3 p;
  };
  using absDeltaEtaWRT = AbsDeltaEtaWRT;

  /// Calculator of @f$ \Delta y @f$ with respect to a given momentum
  struct DeltaRapWRT : public DoubleParticleBaseFunctor {
    DeltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    DeltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return deltaRap(p, pb); }
    double operator()(const FourMomentum& p4) const { return deltaRap(p, p4); }
    const FourMomentum p;
  };
  using deltaRapWRT = DeltaRapWRT;

  /// Calculator of @f$ |\Delta y| @f$ with respect to a given momentum
  struct AbsDeltaRapWRT : public DoubleParticleBaseFunctor {
    AbsDeltaRapWRT(const ParticleBase& pb) : p(pb.mom()) {}
    AbsDeltaRapWRT(const FourMomentum& p4) : p(p4) {}
    double operator()(const ParticleBase& pb) const { return fabs(deltaRap(p, pb)); }
    double operator()(const FourMomentum& p4) const { return fabs(deltaRap(p, p4)); }
    const FourMomentum p;
  };
  using absDeltaRapWRT = AbsDeltaRapWRT;

  //@}


  /// @name Non-PID particle properties, via unbound functions
  /// @todo Mostly move to functions on FourMomentum
  /// @note In a sub-namespace (imported by default) for protection
  //@{
  namespace Kin {

    /// Unbound function access to momentum
    inline FourMomentum mom(const ParticleBase& p) { return p.mom(); }
    /// Unbound function access to momentum
    inline FourMomentum p4(const ParticleBase& p) { return p.mom(); }

    /// Unbound function access to p3
    inline Vector3 p3(const ParticleBase& p) { return p.p3(); }

    /// Unbound function access to pTvec
    inline Vector3 pTvec(const ParticleBase& p) { return p.pTvec(); }

    /// Unbound function access to p
    inline double p(const ParticleBase& p) { return p.p(); }

    /// Unbound function access to pT
    inline double pT(const ParticleBase& p) { return p.pT(); }

    /// Unbound function access to ET
    inline double Et(const ParticleBase& p) { return p.Et(); }

    /// Unbound function access to eta
    inline double eta(const ParticleBase& p) { return p.eta(); }

    /// Unbound function access to abseta
    inline double abseta(const ParticleBase& p) { return p.abseta(); }

    /// Unbound function access to rapidity
    inline double rap(const ParticleBase& p) { return p.rap(); }

    /// Unbound function access to abs rapidity
    inline double absrap(const ParticleBase& p) { return p.absrap(); }

  }
  //@}


  // Import Kin namespace into Rivet
  using namespace Kin;


}

#endif
