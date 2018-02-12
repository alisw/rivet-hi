#ifndef RIVET_JETUTILS_HH
#define RIVET_JETUTILS_HH

#include "Rivet/Jet.hh"
#include "Rivet/Tools/ParticleBaseUtils.hh"

namespace Rivet {


  /// @name Unbound functions for converting between Jets, Particles and PseudoJets
  //@{

  inline PseudoJets mkPseudoJets(const Particles& ps) {
    PseudoJets rtn; rtn.reserve(ps.size());
    for (const Particle& p : ps)
      rtn.push_back(p);
    return rtn;
  }

  inline PseudoJets mkPseudoJets(const Jets& js) {
    PseudoJets rtn; rtn.reserve(js.size());
    for (const Jet& j : js)
      rtn.push_back(j);
    return rtn;
  }

  inline Jets mkJets(const PseudoJets& pjs) {
    Jets rtn; rtn.reserve(pjs.size());
    for (const PseudoJet& pj : pjs)
      rtn.push_back(pj);
    return rtn;
  }

  //@}


  /// @name Jet classifier -> bool functors
  //@{

  /// std::function instantiation for functors taking a Jet and returning a bool
  using JetSelector = function<bool(const Jet&)>;
  /// std::function instantiation for functors taking two Jets and returning a bool
  using JetSorter = function<bool(const Jet&, const Jet&)>;


  /// Base type for Jet -> bool functors
  struct BoolJetFunctor {
    virtual bool operator()(const Jet& p) const = 0;
  };

  struct BoolJetAND : public BoolJetFunctor {
    BoolJetAND(const std::vector<ParticleSelector>& sels) : selectors(sels) {}
    BoolJetAND(const ParticleSelector& a, const ParticleSelector& b) : selectors({a,b}) {}
    BoolJetAND(const ParticleSelector& a, const ParticleSelector& b, const ParticleSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Particle& p) const {
      for (const ParticleSelector& sel : selectors) if (!sel(p)) return false;
      return true;
    }
    std::vector<ParticleSelector> selectors;
  };

  struct BoolJetOR : public BoolJetFunctor {
    BoolJetOR(const std::vector<ParticleSelector>& sels) : selectors(sels) {}
    BoolJetOR(const ParticleSelector& a, const ParticleSelector& b) : selectors({a,b}) {}
    BoolJetOR(const ParticleSelector& a, const ParticleSelector& b, const ParticleSelector& c) : selectors({a,b,c}) {}
    bool operator()(const Particle& p) const {
      for (const ParticleSelector& sel : selectors) if (sel(p)) return true;
      return false;
    }
    std::vector<ParticleSelector> selectors;
  };

  struct BoolJetNOT : public BoolJetFunctor {
    BoolJetNOT(const ParticleSelector& sel) : selector(sel) {}
    bool operator()(const Particle& p) const { return !selector(p); }
    ParticleSelector selector;
  };




  /// B-tagging functor, with a tag selection cut as the stored state
  struct HasBTag : BoolJetFunctor {
    HasBTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasBTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.bTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasBTag = HasBTag;

  /// C-tagging functor, with a tag selection cut as the stored state
  struct HasCTag : BoolJetFunctor {
    HasCTag(const Cut& c=Cuts::open()) : cut(c) {}
    // HasCTag(const std::function<bool(const Jet& j)>& f) : selector(f) {}
    bool operator() (const Jet& j) const { return j.cTagged(cut); }
    // const std::function<bool(const Jet& j)> selector;
    const Cut cut;
  };
  using hasCTag = HasCTag;

  //@}


  /// @name Unbound functions for filtering jets
  //@{

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  Jets& ifilter_select(Jets& jets, const Cut& c);
  /// Alias for ifilter_select
  /// @deprecated Use ifilter_select
  inline Jets& ifilterBy(Jets& jets, const Cut& c) { return ifilter_select(jets, c); }

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets filter_select(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return ifilter_select(rtn, c);
  }
  /// Alias for ifilter_select
  /// @deprecated Use filter_select
  inline Jets filterBy(const Jets& jets, const Cut& c) { return filter_select(jets, c); }

  /// Filter a jet collection in-place to the subset that passes the supplied Cut
  inline Jets filter_select(const Jets& jets, const Cut& c, Jets& out) {
    out = filter_select(jets, c);
    return out;
  }
  /// Alias for ifilter_select
  /// @deprecated Use filter_select
  inline Jets filterBy(const Jets& jets, const Cut& c, Jets& out) { return filter_select(jets, c, out); }


  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  Jets& ifilter_discard(Jets& jets, const Cut& c);

  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets filter_discard(const Jets& jets, const Cut& c) {
    Jets rtn = jets;
    return ifilter_discard(rtn, c);
  }

  /// Filter a jet collection in-place to the subset that fails the supplied Cut
  inline Jets filter_discard(const Jets& jets, const Cut& c, Jets& out) {
    out = filter_discard(jets, c);
    return out;
  }

  //@}



  /// @name Operations on collections of Jet
  /// @note This can't be done on generic collections of ParticleBase -- thanks, C++ :-/
  //@{
  namespace Kin {

    inline double sumPt(const Jets& js) {
      return sum(js, pT, 0.0);
    }

    inline FourMomentum sumP4(const Jets& js) {
      return sum(js, p4, FourMomentum());
    }

    inline Vector3 sumP3(const Jets& js) {
      return sum(js, p3, Vector3());
    }

    /// @todo Min dPhi, min dR?
    /// @todo Isolation routines?

  }
  //@}


  // Import Kin namespace into Rivet
  using namespace Kin;


}

#endif
