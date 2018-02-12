// -*- C++ -*-
#ifndef RIVET_ChargedFinalState_HH
#define RIVET_ChargedFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Project only charged final state particles.
  class ChargedFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Construction from another FinalState
    ChargedFinalState(const FinalState& fsp);

    /// Construction using Cuts object
    ChargedFinalState(const Cut& c=Cuts::open());

    /// Single eta-range constructor.
    ChargedFinalState(double mineta, double maxeta, double minpt=0*GeV);

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(ChargedFinalState);

    //@}


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;
  };


}


#endif
