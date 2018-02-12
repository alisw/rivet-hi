// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// Rivet analysis class for ATLAS_2015_I1387176 dataset
  class ATLAS_2015_I1387176 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1387176);


    /// Initialization, called once before running
    void init() {
      // Projections
      FastJets jets(FinalState(), FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "Jets");

      // Book histograms
      _hist_EEC  = bookHisto1D(1, 1, 1);
      _hist_AEEC = bookScatter2D(2, 1, 1);

      // add dummy histogram for heterogenous merging
      string hname = "d01-x01-y01";
      const Scatter2D& ref = refData(hname);
      hname = "d01-x01-y02";
      _hist_dummy = bookHisto1D(hname, ref);
    }

    void analyze(const Event& event) {

      const double evtWeight = event.weight();
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 50.0*GeV && Cuts::abseta < 2.5);

      if (jets.size() < 2)  vetoEvent;
      if (jets[0].pT() + jets[1].pT() < 500*GeV)  vetoEvent;

      double sumEt = 0.0;
      foreach (Jet j, jets)  sumEt += j.E() / cosh(j.eta());

      foreach (Jet j1, jets) {
        double et1 = j1.E() / cosh(j1.eta());

        foreach (Jet j2, jets) {
          double et2 = j2.E() / cosh(j2.eta());
          double etWeight = et1 * et2 / ( sumEt * sumEt );
          double dPhi = deltaPhi(j1, j2);
          double cosPhi = cos(dPhi);
          if (cosPhi == 1.0)  cosPhi = 0.9999;

          _hist_EEC->fill(cosPhi, etWeight * evtWeight);
          _hist_dummy->fill(cosPhi, etWeight * evtWeight);
	      }
      }
    }

    void finalize() {

      scale(_hist_dummy, crossSectionPerEvent());
      normalize(_hist_EEC);

      vector<Point2D> points;
      size_t nBins = _hist_EEC->numBins();
      for (size_t k = 0; k < nBins/2; ++k) {
        double x = _hist_EEC->bin(k).midpoint();
        double y = _hist_EEC->bin(k).height() - _hist_EEC->bin(nBins-(k+1)).height();
        double ex = _hist_EEC->bin(k).xWidth()/2;
        double e1 = _hist_EEC->bin(k).heightErr();
        double e2 = _hist_EEC->bin(nBins-(k+1)).heightErr();
        double ey = sqrt( e1 * e1 + e2 * e2 );
        points.push_back(Point2D(x, y, ex, ey));
      }

      _hist_AEEC->addPoints(points);
    }

  private:
    Histo1DPtr _hist_EEC;
    Histo1DPtr _hist_dummy;
    Scatter2DPtr _hist_AEEC;
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1387176);

}
