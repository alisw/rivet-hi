// -*- C++ -*-
#ifndef RIVET_UnstableFinalState_HH
#define RIVET_UnstableFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Project out all physical-but-decayed particles in an event.
  ///
  /// The particles returned by are unique unstable particles, such as hadrons
  /// which are decayed by the generator. If, for example, you set Ks and Lambda
  /// particles stable in the generator, they will not be returned. Also, you
  /// should be aware that all unstable particles in a decay chain are returned:
  /// if you are looking for something like the number of B hadrons in an event
  /// and there is a decay chain from e.g. B** -> B, you will count both B
  /// mesons unless you are careful to check for ancestor/descendent relations
  /// between the particles. Duplicate particles in the event record, i.e. those
  /// which differ only in bookkeeping details or photon emissions, are stripped
  /// from the returned particles collection.
  ///
  /// @todo Rename header, with fallback
  /// @todo Convert to a general ParticleFinder since this is explicitly not a final state... but needs care
  /// @todo Make TauFinder inherit/use
  class UnstableParticles : public FinalState {
  public:

    /// @name Standard constructors and destructors.
    //@{

    /// Cut-based / default constructor
    UnstableParticles(const Cut& c=Cuts::open())
      : FinalState(c)
    {
      setName("UnstableParticles");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(UnstableParticles);

    //@}

  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e);

  };


  // Backward compatibility alias
  using UnstableFinalState = UnstableParticles;


}


#endif
