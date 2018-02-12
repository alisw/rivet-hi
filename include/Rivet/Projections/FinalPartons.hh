// -*- C++ -*-
#ifndef RIVET_FinalPartons_HH
#define RIVET_FinalPartons_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class FinalPartons : public FinalState {
  public:

    /// Constructor
    FinalPartons(const Cut& c=Cuts::open())
      : FinalState(c) { }

    /// Clone method
    DEFAULT_RIVET_PROJ_CLONE(FinalPartons);

    /// Do the calculation
    void project(const Event& e);


  protected:

    /// Cut-applying method overload
    bool accept(const Particle& p) const;

  };


}

#endif
