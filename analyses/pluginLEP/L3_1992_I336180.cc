// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief L3 inclusive eta production in hadronic Z0 decays
  /// @author Simone Amoroso 
  class L3_1992_I336180 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(L3_1992_I336180);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableFinalState(), "UFS");

      // Book histograms
      _histXpEta = bookHisto1D( 1, 1, 1);
      _histLnXpEta = bookHisto1D( 2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.                                                    
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.particles().size() < 2) {
	MSG_DEBUG("Failed ncharged cut");
	vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get event weight for histo filling                                                                                                    
      const double weight = event.weight();

      // Get beams and average beam momentum                                                                                                
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const Particles& etas = apply<UnstableFinalState>(event, "UFS").particles(Cuts::abspid==PID::ETA);

      foreach (const Particle& p, etas) {
	double xp = p.p3().mod()/meanBeamMom;
        MSG_DEBUG("Eta xp = " << xp);
        _histXpEta->fill(xp, weight);
        _histLnXpEta->fill(log(1./xp), weight);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histXpEta, 1./sumOfWeights());
      scale(_histLnXpEta, 1./sumOfWeights());
    }
    
    //@}
    
    
  private:
    
    Histo1DPtr _histXpEta;
    Histo1DPtr _histLnXpEta;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(L3_1992_I336180);


}
