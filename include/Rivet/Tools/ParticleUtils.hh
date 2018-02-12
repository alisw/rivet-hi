#ifndef RIVET_PARTICLEUTILS_HH
#define RIVET_PARTICLEUTILS_HH

#include "Rivet/Particle.hh"
#include "Rivet/Tools/ParticleBaseUtils.hh"
#include "Rivet/Tools/ParticleIdUtils.hh"

// Macros to map Rivet::Particle functions to PID:: functions of the same name
#define PARTICLE_TO_PID_BOOLFN(fname) inline bool fname (const Particle& p) { return PID:: fname (p.pid()); }
#define PARTICLE_TO_PID_INTFN(fname) inline int fname (const Particle& p) { return PID:: fname (p.pid()); }
#define PARTICLE_TO_PID_DBLFN(fname) inline double fname (const Particle& p) { return PID:: fname (p.pid()); }

namespace Rivet {


  /// @name Particle classifier functions
  //@{

  /// Unbound function access to PID code
  inline int pid(const Particle& p) { return p.pid(); }

  /// Unbound function access to abs PID code
  inline int abspid(const Particle& p) { return p.abspid(); }


  /// Is this particle species charged?
  PARTICLE_TO_PID_BOOLFN(isCharged)

  /// Is this particle species neutral?
  PARTICLE_TO_PID_BOOLFN(isNeutral)


  /// Is this a neutrino?
  PARTICLE_TO_PID_BOOLFN(isNeutrino)

  /// Determine if the PID is that of a charged lepton
  PARTICLE_TO_PID_BOOLFN(isChargedLepton)
  PARTICLE_TO_PID_BOOLFN(isChLepton)

  /// Determine if the PID is that of a lepton (charged or neutral)
  PARTICLE_TO_PID_BOOLFN(isLepton)

  /// Determine if the PID is that of a photon
  PARTICLE_TO_PID_BOOLFN(isPhoton)

  /// Determine if the PID is that of an electron or positron
  PARTICLE_TO_PID_BOOLFN(isElectron)

  /// Determine if the PID is that of an muon or antimuon
  PARTICLE_TO_PID_BOOLFN(isMuon)

  /// Determine if the PID is that of an tau or antitau
  PARTICLE_TO_PID_BOOLFN(isTau)

  /// Determine if the PID is that of a hadron
  PARTICLE_TO_PID_BOOLFN(isHadron)

  /// Determine if the PID is that of a meson
  PARTICLE_TO_PID_BOOLFN(isMeson)

  /// Determine if the PID is that of a baryon
  PARTICLE_TO_PID_BOOLFN(isBaryon)

  /// Determine if the PID is that of a quark
  PARTICLE_TO_PID_BOOLFN(isQuark)

  /// Determine if the PID is that of a parton (quark or gluon)
  PARTICLE_TO_PID_BOOLFN(isParton)



  /// Determine if the PID is that of a W+
  PARTICLE_TO_PID_BOOLFN(isWplus)

  /// Determine if the PID is that of a W-
  PARTICLE_TO_PID_BOOLFN(isWminus)

  /// Determine if the PID is that of a W+-
  PARTICLE_TO_PID_BOOLFN(isW)

  /// Determine if the PID is that of a Z0
  PARTICLE_TO_PID_BOOLFN(isZ)

  /// Determine if the PID is that of an SM/lightest SUSY Higgs
  PARTICLE_TO_PID_BOOLFN(isHiggs)

  /// Determine if the PID is that of an s/sbar
  PARTICLE_TO_PID_BOOLFN(isStrange)

  /// Determine if the PID is that of a c/cbar
  PARTICLE_TO_PID_BOOLFN(isCharm)

  /// Determine if the PID is that of a b/bbar
  PARTICLE_TO_PID_BOOLFN(isBottom)

  /// Determine if the PID is that of a t/tbar
  PARTICLE_TO_PID_BOOLFN(isTop)


  /// Determine if the particle is a heavy flavour hadron or parton
  PARTICLE_TO_PID_BOOLFN(isHeavyFlavour)

  /// Determine if the PID is that of a heavy parton (c,b,t)
  PARTICLE_TO_PID_BOOLFN(isHeavyParton)

  /// Determine if the PID is that of a light parton (u,d,s)
  PARTICLE_TO_PID_BOOLFN(isLightParton)


  /// Determine if the PID is that of a heavy flavour (b or c) meson
  PARTICLE_TO_PID_BOOLFN(isHeavyMeson)

  /// Determine if the PID is that of a heavy flavour (b or c) baryon
  PARTICLE_TO_PID_BOOLFN(isHeavyBaryon)

  /// Determine if the PID is that of a heavy flavour (b or c) hadron
  PARTICLE_TO_PID_BOOLFN(isHeavyHadron)


  /// Determine if the PID is that of a light flavour (not b or c) meson
  PARTICLE_TO_PID_BOOLFN(isLightMeson)

  /// Determine if the PID is that of a light flavour (not b or c) baryon
  PARTICLE_TO_PID_BOOLFN(isLightBaryon)

  /// Determine if the PID is that of a light flavour (not b or c) hadron
  PARTICLE_TO_PID_BOOLFN(isLightHadron)


  /// Determine if the PID is that of a b-meson.
  PARTICLE_TO_PID_BOOLFN(isBottomMeson)

  /// Determine if the PID is that of a b-baryon.
  PARTICLE_TO_PID_BOOLFN(isBottomBaryon)

  /// Determine if the PID is that of a b-hadron.
  PARTICLE_TO_PID_BOOLFN(isBottomHadron)


  /// @brief Determine if the PID is that of a c-meson.
  ///
  /// Specifically, the _heaviest_ quark is a c: a B_c is a b-meson and NOT a c-meson.
  /// Charmonia (closed charm) are counted as c-mesons here.
  PARTICLE_TO_PID_BOOLFN(isCharmMeson)

  /// @brief Determine if the PID is that of a c-baryon.
  ///
  /// Specifically, the _heaviest_ quark is a c: a baryon containing a b & c
  /// is a b-baryon and NOT a c-baryon. To test for the simpler case, just use
  /// a combination of hasCharm() and isBaryon().
  PARTICLE_TO_PID_BOOLFN(isCharmBaryon)

  /// Determine if the PID is that of a c-hadron.
  PARTICLE_TO_PID_BOOLFN(isCharmHadron)


  // /// Determine if the PID is that of a strange meson
  // PARTICLE_TO_PID_BOOLFN(isStrangeMeson)

  // /// Determine if the PID is that of a strange baryon
  // PARTICLE_TO_PID_BOOLFN(isStrangeBaryon)

  // /// Determine if the PID is that of a strange hadron
  // PARTICLE_TO_PID_BOOLFN(isStrangeHadron)



  /// Is this a pomeron, odderon, or generic reggeon?
  PARTICLE_TO_PID_BOOLFN(isReggeon)

  /// Determine if the PID is that of a diquark (used in hadronization models)
  PARTICLE_TO_PID_BOOLFN(isDiquark)

  /// Determine if the PID is that of a pentaquark (hypothetical hadron)
  PARTICLE_TO_PID_BOOLFN(isPentaquark)

  /// Is this a fundamental SUSY particle?
  PARTICLE_TO_PID_BOOLFN(isSUSY)

  /// Is this an R-hadron?
  PARTICLE_TO_PID_BOOLFN(isRhadron)

  /// Is this a technicolor particle?
  PARTICLE_TO_PID_BOOLFN(isTechnicolor)

  /// Is this an excited (composite) quark or lepton?
  PARTICLE_TO_PID_BOOLFN(isExcited)

  /// Is this a Kaluza-Klein excitation?
  PARTICLE_TO_PID_BOOLFN(isKK)

  /// Is this a graviton?
  PARTICLE_TO_PID_BOOLFN(isGraviton)

  /// Is this a BSM particle (including graviton)?
  PARTICLE_TO_PID_BOOLFN(isBSM)



  /// Determine if the PID is in the generator-specific range
  PARTICLE_TO_PID_BOOLFN(isGenSpecific)

  /// Determine if the PID is that of an EW scale resonance
  PARTICLE_TO_PID_BOOLFN(isResonance)

  /// Check the PID for usability in transport codes like Geant4
  PARTICLE_TO_PID_BOOLFN(isTransportable)



  /// Does this particle contain an up quark?
  PARTICLE_TO_PID_BOOLFN(hasUp)

  /// Does this particle contain a down quark?
  PARTICLE_TO_PID_BOOLFN(hasDown)

  /// Does this particle contain a strange quark?
  PARTICLE_TO_PID_BOOLFN(hasStrange)

  /// Does this particle contain a charm quark?
  PARTICLE_TO_PID_BOOLFN(hasCharm)

  /// Does this particle contain a bottom quark?
  PARTICLE_TO_PID_BOOLFN(hasBottom)

  /// Does this particle contain a top quark?
  PARTICLE_TO_PID_BOOLFN(hasTop)



  /// jSpin returns 2J+1, where J is the total spin
  PARTICLE_TO_PID_INTFN(jSpin)

  /// sSpin returns 2S+1, where S is the spin
  PARTICLE_TO_PID_INTFN(sSpin)

  /// lSpin returns 2L+1, where L is the orbital angular momentum
  PARTICLE_TO_PID_INTFN(lSpin)


  /// Return the charge
  PARTICLE_TO_PID_DBLFN(charge)

  /// Return 3 times the charge (3 x quark charge is an int)
  PARTICLE_TO_PID_INTFN(charge3)

  /// Return the absolute charge
  PARTICLE_TO_PID_DBLFN(abscharge)

  /// Return 3 times the abs charge (3 x quark charge is an int)
  PARTICLE_TO_PID_INTFN(abscharge3)

  /// Alias for charge3
  /// @deprecated Use charge3
  PARTICLE_TO_PID_INTFN(threeCharge)


  /// Get the atomic number (number of protons) in a nucleus/ion
  PARTICLE_TO_PID_INTFN(nuclZ)

  /// Get the atomic weight (number of nucleons) in a nucleus/ion
  PARTICLE_TO_PID_INTFN(nuclA)

  /// If this is a nucleus (ion), get nLambda
  PARTICLE_TO_PID_INTFN(nuclNlambda)

  //@}


  /// @name Particle charge/sign comparison functions
  //@{

  /// @brief Return true if Particles @a a and @a b have the opposite charge sign
  /// @note Two neutrals returns false
  inline bool oppSign(const Particle& a, const Particle& b) {
    return sign(a.charge3()) == -sign(b.charge3()) && sign(a.charge3()) != ZERO;
  }

  /// Return true if Particles @a a and @a b have the same charge sign
  /// @note Two neutrals returns true
  inline bool sameSign(const Particle& a, const Particle& b) {
    return sign(a.charge3()) == sign(b.charge3());
  }

  /// Return true if Particles @a a and @a b have the exactly opposite charge
  /// @note Two neutrals returns false
  inline bool oppCharge(const Particle& a, const Particle& b) {
    return a.charge3() == -b.charge3() && a.charge3() != 0;
  }

  /// Return true if Particles @a a and @a b have the same charge (including neutral)
  /// @note Two neutrals returns true
  inline bool sameCharge(const Particle& a, const Particle& b) {
    return a.charge3() == b.charge3();
  }

  /// Return true if Particles @a a and @a b have a different (not necessarily opposite) charge
  inline bool diffCharge(const Particle& a, const Particle& b) {
    return a.charge3() != b.charge3();
  }

  //@}



  //////////////////////////////////////



  /// @name Non-PID particle properties, via unbound functions
  //@{

  /// @brief Determine whether a particle is the first in a decay chain to meet the function requirement
  inline bool isFirstWith(const Particle& p, const ParticleSelector& f) {
    return p.isFirstWith(f);
  }

  /// @brief Determine whether a particle is the first in a decay chain not to meet the function requirement
  inline bool isFirstWithout(const Particle& p, const ParticleSelector& f) {
    return p.isFirstWithout(f);
  }


  /// @brief Determine whether a particle is the last in a decay chain to meet the function requirement
  inline bool isLastWith(const Particle& p, const ParticleSelector& f) {
    return p.isLastWith(f);
  }

  /// @brief Determine whether a particle is the last in a decay chain not to meet the function requirement
  inline bool isLastWithout(const Particle& p, const ParticleSelector& f) {
    return p.isLastWithout(f);
  }



  /// @brief Determine whether a particle has an ancestor which meets the function requirement
  inline bool hasAncestorWith(const Particle& p, const ParticleSelector& f) {
    return p.hasAncestorWith(f);
  }

  /// @brief Determine whether a particle has an ancestor which doesn't meet the function requirement
  inline bool hasAncestorWithout(const Particle& p, const ParticleSelector& f) {
    return p.hasAncestorWithout(f);
  }


  /// @brief Determine whether a particle has a parent which meets the function requirement
  inline bool hasParentWith(const Particle& p, const ParticleSelector& f) {
    return p.hasParentWith(f);
  }

  /// @brief Determine whether a particle has a parent which doesn't meet the function requirement
  inline bool hasParentWithout(const Particle& p, const ParticleSelector& f) {
    return p.hasParentWithout(f);
  }


  /// @brief Determine whether a particle has a child which meets the function requirement
  inline bool hasChildWith(const Particle& p, const ParticleSelector& f) {
    return p.hasChildWith(f);
  }

  /// @brief Determine whether a particle has a child which doesn't meet the function requirement
  inline bool hasChildWithout(const Particle& p, const ParticleSelector& f) {
    return p.hasChildWithout(f);
  }


  /// @brief Determine whether a particle has a descendant which meets the function requirement
  inline bool hasDescendantWith(const Particle& p, const ParticleSelector& f) {
    return p.hasDescendantWith(f);
    // return !p.allDescendants(f).empty();
  }

  /// @brief Determine whether a particle has a descendant which doesn't meet the function requirement
  inline bool hasDescendantWithout(const Particle& p, const ParticleSelector& f) {
    return p.hasDescendantWithout(f);
  }


  /// @brief Determine whether a particle has a stable descendant which meets the function requirement
  inline bool hasStableDescendantWith(const Particle& p, const ParticleSelector& f) {
    return p.hasStableDescendantWith(f);
  }

  /// @brief Determine whether a particle has a stable descendant which doesn't meet the function requirement
  inline bool hasStableDescendantWithout(const Particle& p, const ParticleSelector& f) {
    return p.hasStableDescendantWithout(f);
  }



  /// Is this particle potentially visible in a detector?
  inline bool isVisible(const Particle& p) { return p.isVisible(); }

  /// @brief Decide if a given particle is direct, via Particle::isDirect()
  ///
  /// A "direct" particle is one directly connected to the hard process. It is a
  /// preferred alias for "prompt", since it has no confusing implications about
  /// distinguishability by timing information.
  ///
  /// The boolean arguments allow a decay lepton to be considered direct if
  /// its parent was a "real" direct lepton.
  inline bool isDirect(const Particle& p, bool allow_from_direct_tau=false, bool allow_from_direct_mu=false) {
    return p.isDirect(allow_from_direct_tau, allow_from_direct_mu);
  }

  /// @brief Decide if a given particle is prompt, via Particle::isPrompt()
  ///
  /// The boolean arguments allow a decay lepton to be considered prompt if
  /// its parent was a "real" prompt lepton.
  inline bool isPrompt(const Particle& p, bool allow_from_prompt_tau=false, bool allow_from_prompt_mu=false) {
    return p.isPrompt(allow_from_prompt_tau, allow_from_prompt_mu);
  }


  /// Decide if a given particle is stable, via Particle::isStable()
  inline bool isStable(const Particle& p) { return p.isStable(); }

  /// Decide if a given particle decays hadronically
  inline bool hasHadronicDecay(const Particle& p) {
    if (p.isStable()) return false;
    if (p.hasChildWith(isHadron)) return true;
    return false;
  }

  /// Decide if a given particle decays leptonically (decays, and no hadrons)
  inline bool hasLeptonicDecay(const Particle& p) {
    if (p.isStable()) return false;
    if (p.hasChildWith(isHadron)) return false;
    return true;
  }


  /// Check whether a given PID is found in the particle's ancestor list
  /// @deprecated Prefer hasAncestorWith
  inline bool hasAncestor(const Particle& p, PdgId pid)  { return p.hasAncestor(pid); }

  /// Determine whether the particle is from a b-hadron decay
  inline bool fromBottom(const Particle& p) { return p.fromBottom(); }

  /// @brief Determine whether the particle is from a c-hadron decay
  inline bool fromCharm(const Particle& p) { return p.fromCharm(); }

  /// @brief Determine whether the particle is from a hadron decay
  inline bool fromHadron(const Particle& p) { return p.fromHadron(); }

  /// @brief Determine whether the particle is from a tau decay
  inline bool fromTau(const Particle& p, bool prompt_taus_only=false) {
    return p.fromTau(prompt_taus_only);
  }

  /// @brief Determine whether the particle is from a prompt tau decay
  inline bool fromPromptTau(const Particle& p) { return p.fromPromptTau(); }

  /// @brief Determine whether the particle is from a hadron or tau decay
  /// @deprecated Too vague: use fromHadron or fromHadronicTau
  inline bool fromDecay(const Particle& p) { return p.fromDecay(); }

  //@}


  /// @name Particle classifier -> bool functors
  ///
  /// To be passed to any() or all() e.g. any(p.children(), HasPID(PID::MUON))
  //@{

  /// Base type for Particle -> bool functors
  struct BoolParticleFunctor {
    virtual bool operator()(const Particle& p) const = 0;
  };

  struct BoolParticleAND : public BoolParticleFunctor {
    BoolParticleAND(const std::vector<ParticleSelector>& sels) : selectors(sels) {}
    BoolParticleAND(const ParticleSelector& a, const ParticleSelector& b) : selectors({a,b}) {}
    BoolParticleAND(const ParticleSelector& a, const ParticleSelector& b, const ParticleSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Particle& p) const {
      for (const ParticleSelector& sel : selectors) if (!sel(p)) return false;
      return true;
    }
    std::vector<ParticleSelector> selectors;
  };

  struct BoolParticleOR : public BoolParticleFunctor {
    BoolParticleOR(const std::vector<ParticleSelector>& sels) : selectors(sels) {}
    BoolParticleOR(const ParticleSelector& a, const ParticleSelector& b) : selectors({a,b}) {}
    BoolParticleOR(const ParticleSelector& a, const ParticleSelector& b, const ParticleSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Particle& p) const {
      for (const ParticleSelector& sel : selectors) if (sel(p)) return true;
      return false;
    }
    std::vector<ParticleSelector> selectors;
  };

  struct BoolParticleNOT : public BoolParticleFunctor {
    BoolParticleNOT(const ParticleSelector& sel) : selector(sel) {}
    bool operator()(const Particle& p) const { return !selector(p); }
    ParticleSelector selector;
  };


  /// PID matching functor
  struct HasPID : public BoolParticleFunctor {
    HasPID(PdgId pid) : targetpid(pid) { }
    bool operator()(const Particle& p) const { return p.pid() == targetpid; }
    PdgId targetpid;
  };
  using hasPID = HasPID;

  /// |PID| matching functor
  struct HasAbsPID : public BoolParticleFunctor {
    HasAbsPID(PdgId pid) : targetpid(abs(pid)) { }
    bool operator()(const Particle& p) const { return p.abspid() == abs(targetpid); }
    PdgId targetpid;
  };
  using hasAbsPID = HasAbsPID;


  /// Determine whether a particle is the first in a decay chain to meet the cut/function
  struct FirstParticleWith : public BoolParticleFunctor {
    FirstParticleWith(const ParticleSelector& f) : fn(f) { }
    FirstParticleWith(const Cut& c);
    bool operator()(const Particle& p) const { return isFirstWith(p, fn); }
    ParticleSelector fn;
  };
  using firstParticleWith = FirstParticleWith;

  /// Determine whether a particle is the first in a decay chain not to meet the cut/function
  struct FirstParticleWithout : public BoolParticleFunctor {
    FirstParticleWithout(const ParticleSelector& f) : fn(f) { }
    FirstParticleWithout(const Cut& c);
    bool operator()(const Particle& p) const { return isFirstWithout(p, fn); }
    ParticleSelector fn;
  };
  using firstParticleWithout = FirstParticleWithout;


  /// Determine whether a particle is the last in a decay chain to meet the cut/function
  struct LastParticleWith : public BoolParticleFunctor {
    template <typename FN>
    LastParticleWith(const FN& f) : fn(f) { }
    LastParticleWith(const Cut& c);
    bool operator()(const Particle& p) const { return isLastWith(p, fn); }
    std::function<bool(const Particle&)> fn;
  };
  using lastParticleWith = LastParticleWith;

  /// Determine whether a particle is the last in a decay chain not to meet the cut/function
  struct LastParticleWithout : public BoolParticleFunctor {
    LastParticleWithout(const ParticleSelector& f) : fn(f) { }
    LastParticleWithout(const Cut& c);
    bool operator()(const Particle& p) const { return isLastWithout(p, fn); }
    ParticleSelector fn;
  };
  using lastParticleWithout = LastParticleWithout;


  /// Determine whether a particle has an ancestor which meets the cut/function
  struct HasParticleAncestorWith : public BoolParticleFunctor {
    HasParticleAncestorWith(const ParticleSelector& f) : fn(f) { }
    HasParticleAncestorWith(const Cut& c);
    bool operator()(const Particle& p) const { return hasAncestorWith(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleAncestorWith = HasParticleAncestorWith;

  /// Determine whether a particle has an ancestor which doesn't meet the cut/function
  struct HasParticleAncestorWithout : public BoolParticleFunctor {
    HasParticleAncestorWithout(const ParticleSelector& f) : fn(f) { }
    HasParticleAncestorWithout(const Cut& c);
    bool operator()(const Particle& p) const { return hasAncestorWithout(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleAncestorWithout = HasParticleAncestorWithout;


  /// Determine whether a particle has an parent which meets the cut/function
  struct HasParticleParentWith : public BoolParticleFunctor {
    HasParticleParentWith(const ParticleSelector& f) : fn(f) { }
    HasParticleParentWith(const Cut& c);
    bool operator()(const Particle& p) const { return hasParentWith(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleParentWith = HasParticleParentWith;

  /// Determine whether a particle has an parent which doesn't meet the cut/function
  struct HasParticleParentWithout : public BoolParticleFunctor {
    HasParticleParentWithout(const ParticleSelector& f) : fn(f) { }
    HasParticleParentWithout(const Cut& c);
    bool operator()(const Particle& p) const { return hasParentWithout(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleParentWithout = HasParticleParentWithout;


  /// Determine whether a particle has a child which meets the cut/function
  struct HasParticleChildWith : public BoolParticleFunctor {
    HasParticleChildWith(const ParticleSelector& f) : fn(f) { }
    HasParticleChildWith(const Cut& c);
    bool operator()(const Particle& p) const { return hasChildWith(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleChildWith = HasParticleChildWith;

  /// Determine whether a particle has a child which doesn't meet the cut/function
  struct HasParticleChildWithout : public BoolParticleFunctor {
    HasParticleChildWithout(const ParticleSelector& f) : fn(f) { }
    HasParticleChildWithout(const Cut& c);
    bool operator()(const Particle& p) const { return hasChildWithout(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleChildWithout = HasParticleChildWithout;


  /// Determine whether a particle has a descendant which meets the cut/function
  struct HasParticleDescendantWith : public BoolParticleFunctor {
    HasParticleDescendantWith(const ParticleSelector& f) : fn(f) { }
    HasParticleDescendantWith(const Cut& c);
    bool operator()(const Particle& p) const { return hasDescendantWith(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleDescendantWith = HasParticleDescendantWith;

  /// Determine whether a particle has a descendant which doesn't meet the cut/function
  struct HasParticleDescendantWithout : public BoolParticleFunctor {
    HasParticleDescendantWithout(const ParticleSelector& f) : fn(f) { }
    HasParticleDescendantWithout(const Cut& c);
    bool operator()(const Particle& p) const { return hasDescendantWithout(p, fn); }
    ParticleSelector fn;
  };
  using hasParticleDescendantWithout = HasParticleDescendantWithout;

  //@}


  /// @name Unbound functions for filtering particles
  //@{

  /// Filter a particle collection in-place to the subset that passes the supplied Cut
  Particles& ifilter_select(Particles& particles, const Cut& c);
  /// Alias for ifilter_select
  /// @deprecated Use ifilter_select
  inline Particles& ifilterBy(Particles& particles, const Cut& c) { return ifilter_select(particles, c); }

  /// Filter a particle collection in-place to the subset that passes the supplied Cut
  inline Particles filter_select(const Particles& particles, const Cut& c) {
    Particles rtn = particles;
    return ifilter_select(rtn, c);
  }
  /// Alias for ifilter_select
  /// @deprecated Use filter_select
  inline Particles filterBy(const Particles& particles, const Cut& c) { return filter_select(particles, c); }

  /// Filter a particle collection in-place to the subset that passes the supplied Cut
  inline Particles filter_select(const Particles& particles, const Cut& c, Particles& out) {
    out = filter_select(particles, c);
    return out;
  }
  /// Alias for ifilter_select
  /// @deprecated Use filter_select
  inline Particles filterBy(const Particles& particles, const Cut& c, Particles& out) { return filter_select(particles, c, out); }


  /// Filter a particle collection in-place to the subset that fails the supplied Cut
  Particles& ifilter_discard(Particles& particles, const Cut& c);

  /// Filter a particle collection in-place to the subset that fails the supplied Cut
  inline Particles filter_discard(const Particles& particles, const Cut& c) {
    Particles rtn = particles;
    return ifilter_discard(rtn, c);
  }

  /// Filter a particle collection in-place to the subset that fails the supplied Cut
  inline Particles filter_discard(const Particles& particles, const Cut& c, Particles& out) {
    out = filter_discard(particles, c);
    return out;
  }

  //@}



  /// @name Particle pair functions
  //@{

  /// Get the PDG ID codes of a ParticlePair
  /// @todo Make ParticlePair a custom class instead?
  inline PdgIdPair pids(const ParticlePair& pp) {
    return make_pair(pp.first.pid(), pp.second.pid());
  }

  //@}



  /// @name Operations on collections of Particle
  /// @note This can't be done on generic collections of ParticleBase -- thanks, C++ :-/
  //@{
  namespace Kin {

    inline double sumPt(const Particles& ps) {
      return sum(ps, pT, 0.0);
    }

    inline FourMomentum sumP4(const Particles& ps) {
      return sum(ps, p4, FourMomentum());
    }

    inline Vector3 sumP3(const Particles& ps) {
      return sum(ps, p3, Vector3());
    }

    /// @todo Min dPhi, min dR?
    /// @todo Isolation routines?

  }
  //@}


  // Import Kin namespace into Rivet
  using namespace Kin;


}

#endif
