// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"



namespace Rivet {


  class DELPHI_2011_I890503 : public Analysis {
  public:

    /// Constructor
    DELPHI_2011_I890503()
      : Analysis("DELPHI_2011_I890503")
    {
    }


    /// Book projections and histograms
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableFinalState(), "UFS");

      _histXbweak     = bookHisto1D(1, 1, 1);
      _histMeanXbweak = bookProfile1D(2, 1, 1);
    }


    void analyze(const Event& e) {

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (apply<FinalState>(e, "FS").particles().size() < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get event weight for histo filling
      const double weight = e.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
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


    // Finalize
    void finalize() {
      normalize(_histXbweak);
    }


  private:

    Histo1DPtr _histXbweak;
    Profile1DPtr _histMeanXbweak;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(DELPHI_2011_I890503);

}
