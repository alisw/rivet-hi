// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  class CMS_2012_I1090423 : public Analysis {
  public:

    CMS_2012_I1090423()
      : Analysis("CMS_2012_I1090423")
    { }


    void init() {
      FinalState fs;
      FastJets antikt(fs, FastJets::ANTIKT, 0.5);
      declare(antikt, "ANTIKT");
      _h_chi_dijet.addHistogram(3000, 7000, bookHisto1D(1, 1, 1));
      _h_chi_dijet.addHistogram(2400, 3000, bookHisto1D(2, 1, 1));
      _h_chi_dijet.addHistogram(1900, 2400, bookHisto1D(3, 1, 1));
      _h_chi_dijet.addHistogram(1500, 1900, bookHisto1D(4, 1, 1));
      _h_chi_dijet.addHistogram(1200, 1500, bookHisto1D(5, 1, 1));
      _h_chi_dijet.addHistogram(1000, 1200, bookHisto1D(6, 1, 1));
      _h_chi_dijet.addHistogram( 800, 1000, bookHisto1D(7, 1, 1));
      _h_chi_dijet.addHistogram( 600,  800, bookHisto1D(8, 1, 1));
      _h_chi_dijet.addHistogram( 400,  600, bookHisto1D(9, 1, 1));
    }


    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "ANTIKT").jetsByPt();
      if (jets.size() < 2) vetoEvent;

      const double y0 = jets[0].rapidity();
      const double y1 = jets[1].rapidity();
      if (fabs(y0+y1)/2 > 1.11) vetoEvent;

      const double chi = exp(fabs(y0-y1));
      if (chi > 16) vetoEvent;

      const FourMomentum jj = jets[0].momentum() + jets[1].momentum();
       _h_chi_dijet.fill(jj.mass(), chi, event.weight());
    }


    void finalize() {
      foreach (Histo1DPtr hist, _h_chi_dijet.getHistograms()) {
        normalize(hist);
      }
    }


  private:

    BinnedHistogram<double> _h_chi_dijet;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1090423);

}
