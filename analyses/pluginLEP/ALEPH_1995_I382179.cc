// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief ALEPH pi+-, K+-, p and pbar differential cross-sections at the Z peak
  class ALEPH_1995_I382179 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALEPH_1995_I382179);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      _histXpPion = bookHisto1D( 1, 1, 1);
      _histXpKaon = bookHisto1D( 2, 1, 1);
      _histXpProton = bookHisto1D( 3, 1, 1);

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

      foreach (const Particle& p, fs.particles()) {
	int id = p.abspid();
	// charged pions
	if (id == PID::PIPLUS || id == PID::PIMINUS) {
	  _histXpPion->fill(p.p3().mod()/meanBeamMom, weight);
	} else if(id == PID::KPLUS || id == PID::KMINUS) {
	  _histXpKaon->fill(p.p3().mod()/meanBeamMom, weight);
	} else if(id == PID::PROTON || id == PID::ANTIPROTON) {
	  _histXpProton->fill(p.p3().mod()/meanBeamMom, weight);
	}
      }


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histXpPion, 1./sumOfWeights());
      scale(_histXpKaon, 1./sumOfWeights());
      scale(_histXpProton, 1./sumOfWeights());
    }

    //@}


  private:


    /// @name Histograms
    Histo1DPtr _histXpPion;
    Histo1DPtr _histXpKaon;
    Histo1DPtr _histXpProton;


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALEPH_1995_I382179);


}
