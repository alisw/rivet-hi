// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {

  // This analysis is a derived from the class Analysis:
  class CMS_2013_I1208923 : public Analysis {

  public:
    // Constructor
    CMS_2013_I1208923()
      : Analysis("CMS_2013_I1208923") {
      //setNeedsCrossSection(true);
    }

    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;
      declare(fs, "FS");

      // Initialize the projections
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "Jets");

      // Book histograms
      _h_sigma.addHistogram(0.0, 0.5, bookHisto1D(1, 1, 1));
      _h_sigma.addHistogram(0.5, 1.0, bookHisto1D(1, 1, 2));
      _h_sigma.addHistogram(1.0, 1.5, bookHisto1D(1, 1, 3));
      _h_sigma.addHistogram(1.5, 2.0, bookHisto1D(1, 1, 4));
      _h_sigma.addHistogram(2.0, 2.5, bookHisto1D(1, 1, 5));
      
      _h_invMass.addHistogram(0.0, 0.5, bookHisto1D(2, 1, 1));
      _h_invMass.addHistogram(0.5, 1.0, bookHisto1D(2, 1, 2));
      _h_invMass.addHistogram(1.0, 1.5, bookHisto1D(2, 1, 3));
      _h_invMass.addHistogram(1.5, 2.0, bookHisto1D(2, 1, 4));
      _h_invMass.addHistogram(2.0, 2.5, bookHisto1D(2, 1, 5));
    }

    // Analysis
    void analyze(const Event &event) {
      const double weight = event.weight();
      const FastJets &fJets = apply<FastJets>(event, "Jets");
      
      // Fill the jet pT spectra
      const Jets& jets = fJets.jetsByPt(Cuts::pt>100.*GeV && Cuts::absrap <2.5);
      foreach (const Jet &j, jets) {
        _h_sigma.fill(fabs(j.momentum().rapidity()), j.momentum().pT() / GeV, weight);
      }

      // Require two jets
      const Jets& dijets = fJets.jetsByPt(Cuts::pt>30.*GeV && Cuts::absrap < 2.5);
      if (dijets.size() > 1) {
        if (dijets[0].momentum().pT() / GeV > 60.) {
          // Fill the invariant mass histogram
          double ymax = max(dijets[0].momentum().absrapidity(), dijets[1].momentum().absrapidity());
          double invMass = FourMomentum(dijets[0].momentum() + dijets[1].momentum()).mass();
          _h_invMass.fill(fabs(ymax), invMass, weight);
        }
      } 

    }


    // Scale histograms by the production cross section
    void finalize() {
      _h_sigma.scale(  crossSection() / sumOfWeights() / 2.0, this);
      _h_invMass.scale(crossSection() / sumOfWeights() / 2.0, this);
    }

  private:
    BinnedHistogram<double> _h_sigma;
    BinnedHistogram<double> _h_invMass;
  };

  // This global object acts as a hook for the plugin system.
  DECLARE_RIVET_PLUGIN(CMS_2013_I1208923);
}
