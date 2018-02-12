// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief OPAL b-fragmentation measurement for weak B-hadron decays
  /// @author Simone Amoroso
  class OPAL_2003_I599181 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_2003_I599181);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableFinalState(), "UFS");

      // Book histograms
      _histXbweak     = bookHisto1D(1, 1, 1);
      _histMeanXbweak = bookProfile1D(2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Get event weight for histo filling
      const double weight = event.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const UnstableFinalState& ufs = apply<UnstableFinalState>(event, "UFS");
      // Get Bottom hadrons
      const Particles bhads = filter_select(ufs.particles(), isBottomHadron);

      for (const Particle& bhad : bhads) {
        // Check for weak decay, i.e. no more bottom present in children
        if (bhad.children(lastParticleWith(hasBottom)).empty()) {
          const double xp = bhad.E()/meanBeamMom;
          _histXbweak->fill(xp, weight);
          _histMeanXbweak->fill(_histMeanXbweak->bin(0).xMid(), xp, weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_histXbweak);
    }

    //@}


  private:

    Histo1DPtr _histXbweak;
    Profile1DPtr _histMeanXbweak;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2003_I599181);


}
