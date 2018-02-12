// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class PLUTO_1980_I154270 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(PLUTO_1980_I154270);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      if(fuzzyEquals(sqrtS()/GeV,30.75)) {
	_hist=bookProfile1D(1, 2, 1);
      }
      else if (fuzzyEquals(sqrtS()/GeV,9.4 ) ||
	       fuzzyEquals(sqrtS()/GeV,12.0) ||
	       fuzzyEquals(sqrtS()/GeV,13.0) ||
	       fuzzyEquals(sqrtS()/GeV,17.0) ||
	       fuzzyEquals(sqrtS()/GeV,22.0) ||
	       fuzzyEquals(sqrtS()/GeV,27.6) ||
	       fuzzyEquals(sqrtS()/GeV,30.2) ||
	       fuzzyEquals(sqrtS()/GeV,30.7) ||
	       fuzzyEquals(sqrtS()/GeV,31.3)) {
	_hist=bookProfile1D(1, 1, 1);
      }
      else {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      unsigned int nPart(0);
      foreach (const Particle& p, cfs.particles()) {
        // check if prompt or not
        const GenParticle* pmother = p.genParticle();
        const GenVertex* ivertex = pmother->production_vertex();
        bool prompt = true;
        while (ivertex) {
          int n_inparts = ivertex->particles_in_size();
          if (n_inparts < 1) break;
          pmother = particles(ivertex, HepMC::parents)[0]; // first mother particle
          int mother_pid = abs(pmother->pdg_id());
          if (mother_pid==PID::K0S || mother_pid==PID::LAMBDA) {
            prompt = false;
            break;
          }
          else if (mother_pid<6) {
            break;
          }
          ivertex = pmother->production_vertex();
        }
	if(prompt) ++nPart;
      }
      _hist->fill(sqrtS(),nPart,event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() { }

    //@}


  private:

    // Histogram
    Profile1DPtr _hist;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PLUTO_1980_I154270);


}
