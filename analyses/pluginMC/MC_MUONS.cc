// -*- C++ -*-
#include "Rivet/Analyses/MC_ParticleAnalysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for muons
  class MC_MUONS : public MC_ParticleAnalysis {
  public:

    MC_MUONS()
      : MC_ParticleAnalysis("MC_MUONS", 2, "muon")
    {    }


  public:

    void init() {
      IdentifiedFinalState muons;
      muons.acceptIdPair(PID::MUON);
      declare(muons, "Muons");

      MC_ParticleAnalysis::init();
    }


    void analyze(const Event& event) {
      const Particles mus = apply<FinalState>(event, "Muons").particlesByPt(Cuts::pT > 0.5*GeV);
      MC_ParticleAnalysis::_analyze(event, mus);
    }


    void finalize() {
      MC_ParticleAnalysis::finalize();
    }

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_MUONS);

}
