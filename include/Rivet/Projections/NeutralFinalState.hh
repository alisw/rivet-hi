// -*- C++ -*-
#ifndef RIVET_NeutralFinalState_HH
#define RIVET_NeutralFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Project only neutral final state particles.
  class NeutralFinalState : public FinalState {

  public:

    /// @name Constructors
    //@{

    /// Construction from another FinalState
    NeutralFinalState(const FinalState& fsp, double etmin=0*GeV)
      : _Etmin(etmin)
    {
      setName("NeutralFinalState");
      addProjection(fsp, "FS");
    }

    /// Construction using Cuts object
    NeutralFinalState(const Cut& c=Cuts::open()) : _Etmin(0.0*GeV) {
      setName("NeutralFinalState");
      addProjection(FinalState(c), "FS");
    }

    /// Construction from explicit eta range and min ET cut values
    NeutralFinalState(double mineta, double maxeta, double etmin=0*GeV)
      : _Etmin(etmin)
    {
      setName("NeutralFinalState");
      addProjection(FinalState(mineta, maxeta, 0.0*GeV), "FS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(NeutralFinalState);

    //@}


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// The minimum allowed transverse energy.
    double _Etmin;

    /// Compare projections.
    int compare(const Projection& p) const;
  };


}


#endif
