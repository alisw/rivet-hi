// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Inclusive jet cross sections using early 13 TeV data
  /// @author Stefan von Buddenbrock <stef.von.b@cern.ch>
  class ATLAS_2016_CONF_2016_092 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_CONF_2016_092);


    /// Bookings
    void init() {

      // Define the jets
      FastJets antiKT04Jets(FinalState(Cuts::open()),  FastJets::ANTIKT, 0.4);
      antiKT04Jets.useInvisibles();
      declare(antiKT04Jets, "antiKT04Jets");

      // Book histograms
      const double y_bins[] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
      for (size_t i = 0; i < 6; i++)
        _h_pT.addHistogram(y_bins[i], y_bins[i+1], bookHisto1D(i+1, 1, 1));

    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<FastJets>(event, "antiKT04Jets")
        .jetsByPt(Cuts::pT > 100*GeV && Cuts::absrap < 3.0);
      for (const Jet& j : jets)
        _h_pT.fill(j.absrap(), j.pT()/GeV, event.weight());
    }


    /// Post-run scaling
    void finalize() {
      // Divide by 2 to only get positive rapidity values
      _h_pT.scale(0.5*crossSection()/picobarn/sumOfWeights(), this);
    }


    /// Histograms
    BinnedHistogram<double> _h_pT;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_CONF_2016_092);

}
