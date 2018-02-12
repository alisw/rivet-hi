// -*- C++ -*-
#ifndef RIVET_FinalState_HH
#define RIVET_FinalState_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Project out all final-state particles in an event.
  /// Probably the most important projection in Rivet!
  class FinalState : public ParticleFinder {
  private:

    // Hide lossy copy constructors for all classes derived from FinalState
    template<typename T> FinalState(const T& rhs);
    template<typename T> FinalState const& operator=(T const& rhs);


  public:

    /// @name Standard constructors etc.
    //@{

    /// Construction using Cuts object
    FinalState(const Cut& c=Cuts::open());

    // /// Construction using Cuts object and another FinalState
    // FinalState(const Cut& c=Cuts::open(), const FinalState& fsp=FinalState());

    /// Old constructor with numeric cut arguments, retained for compatibility
    /// @deprecated Use the versions with Cut arguments
    FinalState(double mineta, double maxeta, double minpt=0.0*GeV);

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(FinalState);

    //@}


    /// Apply the projection to the event.
    virtual void project(const Event& e);

    /// Compare projections.
    virtual int compare(const Projection& p) const;

    /// Decide if a particle is to be accepted or not.
    /// @todo Rename to _accept or acceptFinal?
    virtual bool accept(const Particle& p) const;

  };


}

#endif
