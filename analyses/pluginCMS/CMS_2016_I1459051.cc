// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Inclusive jet pT at 13 TeV
  class CMS_2016_I1459051 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1459051);


    /// Book histograms and initialize projections:
    void init() {

      // Initialize the projections
      const FinalState fs;
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "JetsAK7");

      // Book sets of histograms, binned in absolute rapidity
      // AK7
      _hist_sigmaAK7.addHistogram(0.0, 0.5, bookHisto1D(1, 1, 1));
      _hist_sigmaAK7.addHistogram(0.5, 1.0, bookHisto1D(2, 1, 1));
      _hist_sigmaAK7.addHistogram(1.0, 1.5, bookHisto1D(3, 1, 1));
      _hist_sigmaAK7.addHistogram(1.5, 2.0, bookHisto1D(4, 1, 1));
      _hist_sigmaAK7.addHistogram(2.0, 2.5, bookHisto1D(5, 1, 1));
      _hist_sigmaAK7.addHistogram(2.5, 3.0, bookHisto1D(6, 1, 1));
      _hist_sigmaAK7Forward = bookHisto1D(7, 1, 1);
      // AK4
      _hist_sigmaAK4.addHistogram(0.0, 0.5, bookHisto1D(8, 1, 1));
      _hist_sigmaAK4.addHistogram(0.5, 1.0, bookHisto1D(9, 1, 1));
      _hist_sigmaAK4.addHistogram(1.0, 1.5, bookHisto1D(10, 1, 1));
      _hist_sigmaAK4.addHistogram(1.5, 2.0, bookHisto1D(11, 1, 1));
      _hist_sigmaAK4.addHistogram(2.0, 2.5, bookHisto1D(12, 1, 1));
      _hist_sigmaAK4.addHistogram(2.5, 3.0, bookHisto1D(13, 1, 1));
      _hist_sigmaAK4Forward = bookHisto1D(14, 1, 1);

    }


    /// Per-event analysis
    void analyze(const Event &event) {

      const double weight = event.weight();

      // AK4 jets
      const FastJets& fjAK4 = applyProjection<FastJets>(event, "JetsAK4");
      const Jets& jetsAK4 = fjAK4.jets(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      for (const Jet& j : jetsAK4) {
        _hist_sigmaAK4.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK4Forward->fill(j.pT(), weight);
      }

      // AK7 jets
      const FastJets& fjAK7 = applyProjection<FastJets>(event, "JetsAK7");
      const Jets& jetsAK7 = fjAK7.jets(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      for (const Jet& j : jetsAK7) {
        _hist_sigmaAK7.fill(j.absrap(), j.pT(), weight);
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK7Forward->fill(j.pT(), weight);
      }

    }


    // Finalize
    void finalize() {
      /// @todo What is the cross-section unit?
      _hist_sigmaAK4.scale(crossSection()/sumOfWeights()/2.0, this);
      _hist_sigmaAK7.scale(crossSection()/sumOfWeights()/2.0, this);
      scale(_hist_sigmaAK4Forward,crossSection()/sumOfWeights()/3.0);
      scale(_hist_sigmaAK7Forward,crossSection()/sumOfWeights()/3.0);
    }


    /// @name Histograms
    //@{
    BinnedHistogram<double> _hist_sigmaAK4;
    BinnedHistogram<double> _hist_sigmaAK7;
    Histo1DPtr _hist_sigmaAK4Forward;
    Histo1DPtr _hist_sigmaAK7Forward;
    //@}

  };


  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2016_I1459051);

}
