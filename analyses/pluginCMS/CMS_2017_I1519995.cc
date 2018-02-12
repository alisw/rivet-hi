// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Search for new physics with dijet angular distributions in proton-proton collisions at sqrt{(s) = 13 TeV
  class CMS_2017_I1519995 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1519995);


    /// Book projections and histograms
    void init() {
      FastJets antikt(FinalState(), FastJets::ANTIKT, 0.4);
      declare(antikt, "ANTIKT");
      _h_chi_dijet.addHistogram(4800., 8000., bookHisto1D(1, 1, 1));
      _h_chi_dijet.addHistogram(4200., 4800., bookHisto1D(2, 1, 1));
      _h_chi_dijet.addHistogram(3600., 4200., bookHisto1D(3, 1, 1));
      _h_chi_dijet.addHistogram(3000., 3600., bookHisto1D(4, 1, 1));
      _h_chi_dijet.addHistogram(2400., 3000., bookHisto1D(5, 1, 1));
      _h_chi_dijet.addHistogram(1900., 2400., bookHisto1D(6, 1, 1));
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      const FourMomentum j0(jets[0].mom()), j1(jets[1].mom());
      if (fabs(j0.rap()+j1.rap())/2 > 1.11) vetoEvent;

      const double mjj = (j0+j1).mass();
      const double chi = exp(fabs(j0.rap()-j1.rap()));
      if (chi < 16) _h_chi_dijet.fill(mjj/GeV, chi, event.weight());
    }


    /// Normalize histograms
    void finalize() {
      for (Histo1DPtr hist : _h_chi_dijet.getHistograms()) normalize(hist);
    }


  private:

    BinnedHistogram<double> _h_chi_dijet;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1519995);

}
