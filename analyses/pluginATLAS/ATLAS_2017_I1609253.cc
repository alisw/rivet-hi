// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Multijet transverse energy-energy correlations (TEEC) at 8 TeV
  class ATLAS_2017_I1609253 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1609253);


    /// Initialization, called once before running
    void init() {

      // Projections
      const FastJets jets(FinalState(), FastJets::ANTIKT, 0.4, JetAlg::ALL_MUONS, JetAlg::ALL_INVISIBLES);
      addProjection(jets, "Jets");

      // Book histograms
      _hist_EEC1  = bookHisto1D(   1, 1, 1);
      _hist_AEEC1 = bookScatter2D( 2, 1, 1);
      _hist_EEC2  = bookHisto1D(   3, 1, 1);
      _hist_AEEC2 = bookScatter2D( 4, 1, 1);
      _hist_EEC3  = bookHisto1D(   5, 1, 1);
      _hist_AEEC3 = bookScatter2D( 6, 1, 1);
      _hist_EEC4  = bookHisto1D(   7, 1, 1);
      _hist_AEEC4 = bookScatter2D( 8, 1, 1);
      _hist_EEC5  = bookHisto1D(   9, 1, 1);
      _hist_AEEC5 = bookScatter2D(10, 1, 1);
      _hist_EEC6  = bookHisto1D(  11, 1, 1);
      _hist_AEEC6 = bookScatter2D(12, 1, 1);
    }


    void analyze(const Event& event) {

      const double evtWeight = event.weight();
      const Jets& jets = applyProjection<FastJets>(event, "Jets").jetsByPt(Cuts::abseta < 2.5 && Cuts::pT > 100*GeV);
      if (jets.size() < 2)  vetoEvent;

      double sumPt12 = jets[0].pt() + jets[1].pt();
      if (sumPt12 < 800*GeV)  vetoEvent;

      double sumEt = 0.;
      for (const Jet& j : jets) sumEt += j.Et();

      for (const Jet& j1 : jets) {
        double et1 = j1.Et();

        for (const Jet& j2 : jets) {
          double et2 = j2.Et();

          double etWeight = et1*et2/(sumEt*sumEt);
          double dPhi = deltaPhi(j1, j2);
          double cosPhi = cos(dPhi);
          if (cos(dPhi) == 1.0)  cosPhi = 0.9999;

          if (sumPt12 >  800*GeV && sumPt12 <=  850*GeV)  _hist_EEC1->fill(cosPhi, etWeight*evtWeight);
          if (sumPt12 >  850*GeV && sumPt12 <=  900*GeV)  _hist_EEC2->fill(cosPhi, etWeight*evtWeight);
          if (sumPt12 >  900*GeV && sumPt12 <= 1000*GeV)  _hist_EEC3->fill(cosPhi, etWeight*evtWeight);
          if (sumPt12 > 1000*GeV && sumPt12 <= 1100*GeV)  _hist_EEC4->fill(cosPhi, etWeight*evtWeight);
          if (sumPt12 > 1100*GeV && sumPt12 <= 1400*GeV)  _hist_EEC5->fill(cosPhi, etWeight*evtWeight);
          if (sumPt12 > 1400*GeV)  _hist_EEC6->fill(cosPhi, etWeight*evtWeight);
        }
      }
    }


    void finalize() {

      normalize(_hist_EEC1);
      normalize(_hist_EEC2);
      normalize(_hist_EEC3);
      normalize(_hist_EEC4);
      normalize(_hist_EEC5);
      normalize(_hist_EEC6);

      vector<Point2D> points1, points2, points3, points4, points5, points6;
      size_t nBins = _hist_EEC1->numBins();
      for (size_t k = 0; k < nBins/2; ++k) {
        double x = _hist_EEC1->bin(k).midpoint(); double ex = _hist_EEC1->bin(k).xWidth()/2;

        double y1 = _hist_EEC1->bin(k).height() - _hist_EEC1->bin(nBins-(k+1)).height();
        double ey1 = sqrt( pow(_hist_EEC1->bin(k).heightErr(),2) + pow(_hist_EEC1->bin(nBins-(k+1)).heightErr(),2) );
        points1.push_back(Point2D(x,y1,ex,ey1));

        double y2 = _hist_EEC2->bin(k).height() - _hist_EEC2->bin(nBins-(k+1)).height();
        double ey2 = sqrt( pow(_hist_EEC2->bin(k).heightErr(),2) + pow(_hist_EEC2->bin(nBins-(k+1)).heightErr(),2) );
        points2.push_back(Point2D(x,y2,ex,ey2));

        double y3 = _hist_EEC3->bin(k).height() - _hist_EEC3->bin(nBins-(k+1)).height();
        double ey3 = sqrt( pow(_hist_EEC3->bin(k).heightErr(),2) + pow(_hist_EEC3->bin(nBins-(k+1)).heightErr(),2) );
        points3.push_back(Point2D(x,y3,ex,ey3));

        double y4 = _hist_EEC4->bin(k).height() - _hist_EEC4->bin(nBins-(k+1)).height();
        double ey4 = sqrt( pow(_hist_EEC4->bin(k).heightErr(),2) + pow(_hist_EEC4->bin(nBins-(k+1)).heightErr(),2) );
        points4.push_back(Point2D(x,y4,ex,ey4));

        double y5 = _hist_EEC5->bin(k).height() - _hist_EEC5->bin(nBins-(k+1)).height();
        double ey5 = sqrt( pow(_hist_EEC5->bin(k).heightErr(),2) + pow(_hist_EEC5->bin(nBins-(k+1)).heightErr(),2) );
        points5.push_back(Point2D(x,y5,ex,ey5));

        double y6 = _hist_EEC6->bin(k).height() - _hist_EEC6->bin(nBins-(k+1)).height();
        double ey6 = sqrt( pow(_hist_EEC6->bin(k).heightErr(),2) + pow(_hist_EEC6->bin(nBins-(k+1)).heightErr(),2) );
        points6.push_back(Point2D(x,y6,ex,ey6));
      }

      _hist_AEEC1->addPoints(points1);
      _hist_AEEC2->addPoints(points2);
      _hist_AEEC3->addPoints(points3);
      _hist_AEEC4->addPoints(points4);
      _hist_AEEC5->addPoints(points5);
      _hist_AEEC6->addPoints(points6);
    }


  private:

    Histo1DPtr _hist_EEC1, _hist_EEC2, _hist_EEC3, _hist_EEC4, _hist_EEC5, _hist_EEC6;
    Scatter2DPtr _hist_AEEC1, _hist_AEEC2, _hist_AEEC3, _hist_AEEC4, _hist_AEEC5, _hist_AEEC6;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1609253);


}
