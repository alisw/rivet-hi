// -*- C++ -*-
#ifndef RIVET_PrimaryHadrons_HH
#define RIVET_PrimaryHadrons_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Project out the first hadrons from hadronisation.
  ///
  /// @todo Also be able to return taus? Prefer a separate tau finder.
  /// @todo This assumes that the primary hadrons are unstable... should we also look for stable primary hadrons?
  class PrimaryHadrons : public FinalState {
  public:

    /// @name Constructors and destructors.
    //@{

    /// Constructor with cuts argument
    PrimaryHadrons(const Cut& c=Cuts::open()) {
      setName("PrimaryHadrons");
      addProjection(UnstableFinalState(c), "UFS");
    }

    /// Constructor with specification of the minimum and maximum pseudorapidity
    /// \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
    PrimaryHadrons(double mineta, double maxeta, double minpt=0.0*GeV) {
      setName("PrimaryHadrons");
      addProjection(UnstableFinalState(Cuts::etaIn(mineta, maxeta) && Cuts::pT > minpt), "UFS");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PrimaryHadrons);

    //@}


  protected:

    /// Apply the projection to the event.
    virtual void project(const Event& e);

  };


}


#endif
