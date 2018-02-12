// -*- C++ -*-
#ifndef RIVET_ParticleFinder_HH
#define RIVET_ParticleFinder_HH

#include "Rivet/Projection.hh"

namespace Rivet {


  /// @brief Base class for projections which return subsets of an event's particles
  class ParticleFinder : public Projection {
  public:

    /// @name Object lifetime management
    //@{

    /// Construction using Cuts object
    ParticleFinder(const Cut& c=Cuts::OPEN)
      : _cuts(c), _theParticles()
    { }

    /// Virtual destructor for inheritance
    virtual ~ParticleFinder() {}

    /// Clone on the heap.
    virtual unique_ptr<Projection> clone() const = 0;

    //@}


    /// @name Particle accessors
    //@{

    /// Count the final-state particles
    size_t size() const { return particles().size(); }
    /// Count the final-state particles after a Cut is applied
    size_t size(const Cut& c) const { return particles(c).size(); }
    /// Count the final-state particles after a selection functor is applied
    size_t size(const ParticleSelector& s) const { return particles(s).size(); }

    /// Is this final state empty?
    bool empty() const { return size() == 0; }
    /// Is this final state empty after a Cut is applied?
    bool empty(const Cut& c) const { return size(c) == 0; }
    /// Is this final state empty after a selection functor is applied?
    bool empty(const ParticleSelector& s) const { return size(s) == 0; }

    /// Get the particles in no particular order, with no cuts
    virtual const Particles& particles() const { return _theParticles; }

    /// Get the raw particles in no particular order, with no cuts
    ///
    /// @note Raw particles are the final-state constituents, as opposed to
    /// potentially composite particles returned as the finder's particles()
    Particles rawParticles() const {
      Particles rtn;
      for (const Particle& p : particles()) rtn += p.rawConstituents();
      return rtn;
    }

    /// @brief Get the particles with selection cuts
    /// @note Returns a copy rather than a reference, due to the cuts.
    Particles particles(const Cut& c) const {
      return filter_select(particles(), c);
    }

    /// @brief Get the particles with selection cuts via a functor
    /// @note Returns a copy rather than a reference, due to the cuts.
    Particles particles(const ParticleSelector& selector) const {
      return filter_select(particles(), selector);
    }

    /// Get the particles, ordered by supplied sorting function object
    /// @note Returns a copy rather than a reference, due to cuts and sorting.
    Particles particles(const ParticleSorter& sorter, const Cut& c=Cuts::open()) const {
      return sortBy(particles(c), sorter);
    }

    /// Get the particles, ordered by supplied sorting function object
    /// @note Returns a copy rather than a reference, due to cuts and sorting.
    Particles particles(const Cut& c, const ParticleSorter& sorter) const {
      return sortBy(particles(c), sorter);
    }

    /// Get the particles, ordered by a sorting functor and filtered by a selection functor
    /// @note Returns a copy rather than a reference, due to cuts and sorting.
    Particles particles(const ParticleSelector& selector, const ParticleSorter& sorter) const {
      return sortBy(particles(selector), sorter);
    }

    /// Get the particles, ordered by a sorting functor and filtered by a selection functor
    /// @note Returns a copy rather than a reference, due to cuts and sorting.
    Particles particles(const ParticleSorter& sorter, const ParticleSelector& selector) const {
      return sortBy(particles(selector), sorter);
    }

    /// Get the particles, ordered by decreasing \f$ p_T \f$ and with optional cuts
    ///
    /// This is a very common use-case, so is available as syntatic sugar for particles(c, cmpMomByPt).
    Particles particlesByPt(const Cut& c=Cuts::open()) const {
      return particles(c, cmpMomByPt);
    }

    /// Get the particles, ordered by decreasing \f$ p_T \f$ and with optional cuts
    ///
    /// This is a very common use-case, so is available as syntatic sugar for particles(f, cmpMomByPt).
    Particles particlesByPt(const ParticleSelector& selector) const {
      return particles(selector, cmpMomByPt);
    }

    /// Get the particles, ordered by decreasing \f$ p_T \f$ and with a cut on minimum \f$ p_T \f$
    ///
    /// This is a very common use-case, so is available as syntatic sugar for particles(Cuts::pT >= ptmin, cmpMomByPt).
    Particles particlesByPt(double ptmin) const {
      return particles(Cuts::pT >= ptmin, cmpMomByPt);
    }

    //@}


    /// @todo Replace with cuts() accessor
    ///virtual Cut cuts() const { return _cuts; }


    /// @name For JetAlg compatibility
    //@{

    typedef Particle entity_type;
    typedef Particles collection_type;

    /// Template-usable interface common to JetAlg
    const collection_type& entities() const { return particles(); }

    //@}


  protected:

    /// Apply the projection to the event
    virtual void project(const Event& e) = 0;

    /// Compare projections
    virtual int compare(const Projection& p) const;

    /// The kinematic cuts cuts
    Cut _cuts;

    /// The found particles returned by the particles() methods
    Particles _theParticles;

  };


}

#endif
