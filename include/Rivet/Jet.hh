// -*- C++ -*-
#ifndef RIVET_Jet_HH
#define RIVET_Jet_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Jet.fhh"
#include "Rivet/Particle.hh"
#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Tools/Utils.hh"
#include "Rivet/Tools/RivetFastJet.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include <numeric>

namespace Rivet {


  /// @brief Representation of a clustered jet of particles.
  class Jet : public ParticleBase {
  public:

    /// @name Constructors
    //@{

    /// Constructor from a FastJet PseudoJet, with optional full particle constituents information.
    Jet(const fastjet::PseudoJet& pj, const Particles& particles=Particles(), const Particles& tags=Particles()) {
      setState(pj, particles, tags);
    }

    /// Set the jet data, with optional full particle information.
    Jet(const FourMomentum& pjet, const Particles& particles=Particles(), const Particles& tags=Particles()) {
      setState(pjet, particles, tags);
    }

    /// Default constructor -- only for STL storability
    Jet() { clear(); }

    //@}


    /// @name Access jet constituents
    //@{

    /// Number of particles in this jet.
    size_t size() const { return _particles.size(); }

    /// Get the particles in this jet.
    Particles& particles() { return _particles; }
    /// Get the particles in this jet (const version)
    const Particles& particles() const { return _particles; }
    /// Get the particles in this jet which pass a cut (const)
    const Particles particles(const Cut& c) const { return filter_select(_particles, c); }
    /// Get the particles in this jet which pass a filtering functor (const)
    const Particles particles(const ParticleSelector& s) const { return filter_select(_particles, s); }

    /// Get the particles in this jet (FastJet-like alias)
    Particles& constituents() { return particles(); }
    /// Get the particles in this jet (FastJet-like alias, const version)
    const Particles& constituents() const { return particles(); }
    /// Get the particles in this jet which pass a cut (FastJet-like alias, const)
    const Particles constituents(const Cut& c) const { return particles(c); }
    /// Get the particles in this jet which pass a filtering functor (FastJet-like alias, const)
    const Particles constituents(const ParticleSelector& s) const { return particles(s); }

    /// Check whether this jet contains a particular particle.
    bool containsParticle(const Particle& particle) const;
    /// Nicer alias for containsParticleId
    bool containsPID(const Particle& particle) const { return containsParticle(particle); }

    /// Check whether this jet contains a certain particle type.
    bool containsParticleId(PdgId pid) const;
    /// Nicer alias for containsParticleId
    bool containsPID(PdgId pid) const { return containsParticleId(pid); }

    /// Check whether this jet contains at least one of certain particle types.
    bool containsParticleId(const vector<PdgId>& pids) const;
    /// Nicer alias for containsParticleId
    bool containsPID(const vector<PdgId>& pids) const { return containsParticleId(pids); }

    //@}


    /// @name Tagging
    ///
    /// @note General sources of tag particles are planned. The default jet finding
    /// adds b-hadron, c-hadron, and tau tags by ghost association.
    //@{

    /// @brief Particles which have been tag-matched to this jet
    Particles& tags() { return _tags; }
    /// @brief Particles which have been tag-matched to this jet (const version)
    const Particles& tags() const { return _tags; }
    /// @brief Particles which have been tag-matched to this jet _and_ pass a selector function
    ///
    /// @note Note the less efficient return by value, due to the filtering.
    Particles tags(const ParticleSelector& f) const { return filter_select(tags(), f); }
    /// @brief Particles which have been tag-matched to this jet _and_ pass a Cut
    ///
    /// @note Note the less efficient return by value, due to the cut-pass filtering.
    Particles tags(const Cut& c) const;


    /// @brief b particles which have been tag-matched to this jet (and pass an optional Cut)
    ///
    /// The default jet finding adds b-hadron tags by ghost association.
    Particles bTags(const Cut& c=Cuts::open()) const;
    /// @brief b particles which have been tag-matched to this jet _and_ pass a selector function
    Particles bTags(const ParticleSelector& f) const { return filter_select(bTags(), f); }

    /// Does this jet have at least one b-tag (that passes an optional Cut)?
    bool bTagged(const Cut& c=Cuts::open()) const { return !bTags(c).empty(); }
    /// Does this jet have at least one b-tag (that passes the supplied selector function)?
    bool bTagged(const ParticleSelector& f) const { return !bTags(f).empty(); }


    /// @brief c (and not b) particles which have been tag-matched to this jet (and pass an optional Cut)
    ///
    /// The default jet finding adds c-hadron tags by ghost association.
    Particles cTags(const Cut& c=Cuts::open()) const;
    /// @brief c (and not b) particles which have been tag-matched to this jet and pass a selector function
    Particles cTags(const ParticleSelector& f) const { return filter_select(cTags(), f); }

    /// Does this jet have at least one c-tag (that passes an optional Cut)?
    bool cTagged(const Cut& c=Cuts::open()) const { return !cTags(c).empty(); }
    /// Does this jet have at least one c-tag (that passes the supplied selector function)?
    bool cTagged(const ParticleSelector& f) const { return !cTags(f).empty(); }


    /// @brief Tau particles which have been tag-matched to this jet (and pass an optional Cut)
    ///
    /// The default jet finding adds tau tags by ghost association.
    Particles tauTags(const Cut& c=Cuts::open()) const;
    /// @brief Tau particles which have been tag-matched to this jet and pass a selector function
    Particles tauTags(const ParticleSelector& f) const { return filter_select(tauTags(), f); }

    /// Does this jet have at least one tau-tag (that passes an optional Cut)?
    bool tauTagged(const Cut& c=Cuts::open()) const { return !tauTags(c).empty(); }
    /// Does this jet have at least one tau-tag (that passes the supplied selector function)?
    bool tauTagged(const ParticleSelector& f) const { return !tauTags(f).empty(); }


    /// @brief Check whether this jet contains a bottom-flavoured hadron.
    ///
    /// @deprecated The bTags() or bTagged() function is probably what you want
    /// for tagging. This one ignores the tags() list and draws conclusions
    /// based directly on the jet constituents; the other gives a much better match
    /// to typical experimental methods.
    ///
    /// @note The decision is made by first trying to find a bottom-flavoured particle
    /// in the particles list. Most likely this will fail unless bottom hadrons
    /// are set stable. If @a include_decay_products is true (the default), a
    /// fallback is attempted, using the post-hadronization ancestor history of
    /// all constituents.
    DEPRECATED("Prefer the bTags() or bTagged() function")
    bool containsBottom(bool include_decay_products=true) const;

    /// @brief Check whether this jet contains a charm-flavoured hadron.
    ///
    /// @deprecated The cTags() or cTagged() function is probably what you want
    /// for tagging. This one ignores the tags() list and draws conclusions
    /// based directly on the jet constituents; the other gives a much better match
    /// to typical experimental methods.
    ///
    /// @note The decision is made by first trying to find a charm-flavoured particle
    /// in the particles list. Most likely this will fail unless charmed hadrons
    /// are set stable. If @a include_decay_products is true (the default), a
    /// fallback is attempted, using the post-hadronization ancestor history of
    /// all constituents.
    DEPRECATED("Prefer the cTags() or cTagged() function")
    bool containsCharm(bool include_decay_products=true) const;

    //@}


    /// @name Effective jet 4-vector properties
    //@{

    /// Get equivalent single momentum four-vector.
    const FourMomentum& momentum() const { return _momentum; }

    /// Apply an active Lorentz transform to this jet
    /// @note The Rivet jet momentum, constituent particles, and tag particles will be modified.
    /// @warning The FastJet cluster sequence and pseudojets will not be modified: don't use them after transformation!
    Jet& transformBy(const LorentzTransform& lt);

    /// Get the total energy of this jet.
    double totalEnergy() const { return momentum().E(); }

    /// Get the energy carried in this jet by neutral particles.
    double neutralEnergy() const;

    /// Get the energy carried in this jet by hadrons.
    double hadronicEnergy() const;

    //@}


    /// @name Interaction with FastJet
    //@{

    /// Access the internal FastJet3 PseudoJet (as a const reference)
    const fastjet::PseudoJet& pseudojet() const { return _pseudojet; }

    /// Cast operator to FastJet3 PseudoJet (as a const reference)
    operator const fastjet::PseudoJet& () const { return pseudojet(); }

    //@}


    /// @name Set the jet constituents and properties
    //@{

    /// @brief Set the jet data from a FastJet PseudoJet, with optional particle constituents and tags lists.
    ///
    /// @note The particles() list will be extracted from PseudoJet constituents
    /// by default, making use of an attached user info if one is found.
    Jet& setState(const fastjet::PseudoJet& pj, const Particles& particles=Particles(), const Particles& tags=Particles());

    /// Set all the jet data, with optional full particle constituent and tag information.
    Jet& setState(const FourMomentum& mom, const Particles& particles, const Particles& tags=Particles());

    /// @brief Set the particles collection with full particle information.
    ///
    /// If set, this overrides particle info extracted from the PseudoJet
    Jet& setParticles(const Particles& particles);
    Jet& setConstituents(const Particles& particles) { return setParticles(particles); }

    /// Reset this jet as empty.
    Jet& clear();

    //@}


  private:

    /// FJ3 PseudoJet member to unify PseudoJet and Jet
    fastjet::PseudoJet _pseudojet;

    /// Full constituent particle information. (Filled from PseudoJet if possible.)
    /// @todo Make these mutable or similar? Add a flag to force a cache rebuild?
    Particles _particles;

    /// Particles used to tag this jet (can be anything, but c and b hadrons are the most common)
    Particles _tags;

    /// Effective jet 4-vector (just for caching)
    mutable FourMomentum _momentum;

  };


  /// @name String representation and streaming support
  //@{

  /// Allow a Jet to be passed to an ostream.
  std::ostream& operator << (std::ostream& os, const Jet& j);

  //@}


}


#include "Rivet/Tools/JetUtils.hh"

#endif
