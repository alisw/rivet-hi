// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief A Measurement of K*+- (892) production in hadronic Z0 decays
  /// @author Simone Amoroso
  class OPAL_1993_I342766 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(OPAL_1993_I342766);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableFinalState(), "UFS");
      // Book histograms
      _histXeKStar892   = bookHisto1D( 1, 1, 1);
      _histMeanKStar892   = bookHisto1D( 2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
	MSG_DEBUG("Failed leptonic event cut");
	vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get event weight for histo filling
      const double weight = event.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
				   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableFinalState& ufs = apply<UnstableFinalState>(event, "UFS");

      foreach (const Particle& p, ufs.particles(Cuts::abspid==323)) {
        double xp = p.p3().mod()/meanBeamMom;
        _histXeKStar892->fill(xp, weight);
        _histMeanKStar892->fill(_histMeanKStar892->bin(0).xMid(), weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histXeKStar892, 1./sumOfWeights());
      scale(_histMeanKStar892, 1./sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    Histo1DPtr _histXeKStar892;
    Histo1DPtr _histMeanKStar892;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1993_I342766);


}
