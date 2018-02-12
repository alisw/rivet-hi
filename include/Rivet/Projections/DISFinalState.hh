// -*- C++ -*-
#ifndef RIVET_DISFinalState_HH
#define RIVET_DISFinalState_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief Final state particles boosted to the hadronic center of mass system.
  ///
  /// NB. The DIS scattered lepton is not included in the final state particles.
  class DISFinalState: public FinalState {
  public:

    /// Type of DIS boost to apply
    enum BoostType { HCM, BREIT, LAB };


    /// @name Constructors
    //@{

    /// Constructor with explicit FinalState
    /// @note The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(const FinalState& fs, BoostType boosttype, const DISKinematics& kinematicsp=DISKinematics())
      : _boosttype(boosttype)
    {
      setName("DISFinalState");
      declare(fs, "FS");
      declare(kinematicsp, "Kinematics");
    }

    /// Constructor with optional FinalState
    /// @note The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(BoostType boosttype, const FinalState& fs=FinalState(), const DISKinematics& kinematicsp=DISKinematics())
      : DISFinalState(fs, boosttype, kinematicsp)
    {    }

    /// Constructor with explicit cuts to define final-state particles
    /// @note The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(const Cut& c, BoostType boosttype, const DISKinematics& kinematicsp=DISKinematics())
      : DISFinalState(FinalState(c), boosttype, kinematicsp)
    {    }

    /// Constructor with explicit cuts to define final-state particles
    /// @note The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    DISFinalState(BoostType boosttype, const Cut& c, const DISKinematics& kinematicsp=DISKinematics())
      : DISFinalState(FinalState(c), boosttype, kinematicsp)
    {    }

    // /// @brief Constructor with default FinalState
    // /// @note The DISKinematics has no parameters, hence explicitly passing it as an arg shouldn't be necessary.
    // DISFinalState(BoostType boosttype, const DISKinematics& kinematicsp=DISKinematics())
    //   : DISFinalState(FinalState(), boosttype, kinematicsp)
    // {    }

    /// Backward compatible constructor with default FinalState
    /// @deprecated Prefer a version that doesn't need a DISKinematics argument
    DISFinalState(const DISKinematics& kinematicsp, BoostType boosttype)
      : DISFinalState(FinalState(), boosttype, kinematicsp)
    {    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(DISFinalState);

    //@}


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const {
      const DISFinalState& other = dynamic_cast<const DISFinalState&>(p);
      return mkNamedPCmp(p, "Kinematics") || mkNamedPCmp(p, "FS") || cmp(_boosttype, other._boosttype);
    }


  private:

    BoostType _boosttype;

  };


}

#endif
