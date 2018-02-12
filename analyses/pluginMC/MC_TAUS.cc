// -*- C++ -*-
#include "Rivet/Analyses/MC_ParticleAnalysis.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  /// @brief MC validation analysis for taus
  class MC_TAUS : public MC_ParticleAnalysis {
  public:

    /// Constructor
    MC_TAUS()
      : MC_ParticleAnalysis("MC_TAUS", 2, "tau")
    {    }


    /// Book projections and histograms
    void init() {
      TauFinder taus(TauFinder::ANY);
      declare(taus, "Taus");

      MC_ParticleAnalysis::init();
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const Particles taus = apply<TauFinder>(event, "Taus").particlesByPt(0.5*GeV);
      MC_ParticleAnalysis::_analyze(event, taus);
    }


    /// Normalisations etc.
    void finalize() {
      MC_ParticleAnalysis::finalize();
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_TAUS);

}
