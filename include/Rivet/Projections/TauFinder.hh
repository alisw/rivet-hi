#ifndef RIVET_TauFinder_HH
#define RIVET_TauFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of unstable taus
  ///
  /// @todo Convert to a general ParticleFinder, since it's not a true final state? Needs some care...
  class TauFinder : public FinalState {
  public:

    enum DecayType { ANY=0, LEPTONIC=1, HADRONIC };

    static bool isHadronic(const Particle& tau) {
      assert(tau.abspid() == PID::TAU);
      return any(tau.stableDescendants(), isHadron);
    }

    static bool isLeptonic(const Particle& tau) {
      return !isHadronic(tau);
    }


    TauFinder(DecayType decaytype, const Cut& cut=Cuts::open()) {
      /// @todo What about directness/promptness?
      setName("TauFinder");
      _dectype = decaytype;
      addProjection(UnstableFinalState(cut), "UFS");
    }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(TauFinder);


    const Particles& taus() const { return _theParticles; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare with other projections.
    virtual int compare(const Projection& p) const;


  private:

    /// The decaytype enum
    DecayType _dectype;

  };


  /// @todo Make this the canonical name in future
  using Taus = TauFinder;


}


#endif
