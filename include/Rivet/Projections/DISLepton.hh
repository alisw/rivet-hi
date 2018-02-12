// -*- C++ -*-
#ifndef RIVET_DISLepton_HH
#define RIVET_DISLepton_HH

#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  /// @brief Get the incoming and outgoing leptons in a DIS event.
  class DISLepton : public Projection {
  public:

    /// @name Constructors.
    //@{

    DISLepton(){
      setName("DISLepton");
      addProjection(Beam(), "Beam");
      addProjection(PromptFinalState(), "PromptFS");
    }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(DISLepton);

    //@}


  protected:

    /// Perform the projection operation on the supplied event.
    virtual void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection& p) const;


  public:

    /// The incoming lepton
    const Particle& in() const { return _incoming; }

    /// The outgoing lepton
    const Particle& out() const { return _outgoing; }

    /// Sign of the incoming lepton pz component
    int pzSign() const { return sign(_incoming.pz()); }


  private:

    /// The incoming lepton
    Particle _incoming;

    /// The outgoing lepton
    Particle _outgoing;

    // /// The charge sign of the DIS current
    // double _charge;

  };

}


#endif
