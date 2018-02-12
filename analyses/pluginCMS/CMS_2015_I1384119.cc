// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2015_I1384119 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1384119);


    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fsa(Cuts::abseta < 20);
      declare(fsa, "FSA");
      const ChargedFinalState cfs(Cuts::abseta < 2);
      declare(cfs, "CFS");

      _hist_dNch_dEta_inel = bookHisto1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Apply inelastic selection (veto pp -> pp elastic events)
      const FinalState& fsa = apply<FinalState>(event, "FSA");
      if (fsa.size() <= 2) vetoEvent;

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      foreach (const Particle& p, cfs.particles()) {
        const int id = p.abspid();
        // continue if particle is a proton, a kaon or a pion
        if (id == 211 || id == 321 || id == 2212) ///< @todo Use PID:: ID constants
          _hist_dNch_dEta_inel->fill(p.eta(), event.weight());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_dNch_dEta_inel,  1/sumOfWeights());
    }


  private:

    /// Histograms
    Histo1DPtr _hist_dNch_dEta_inel;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1384119);

}
