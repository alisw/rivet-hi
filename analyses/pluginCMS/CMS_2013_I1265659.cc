// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class CMS_2013_I1265659 : public Analysis {
  public:

    /// Constructor
    CMS_2013_I1265659()
      : Analysis("CMS_2013_I1265659")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {
      const FastJets jets(FinalState(-10, 10, 0.0*GeV), FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      _h_hTotD = bookHisto1D(1, 1, 1);
      _h_hTotDF = bookHisto1D(1, 1, 2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      if (jets.size() < 3) vetoEvent;

      const FourMomentum jet1 = jets[0].momentum();
      const FourMomentum jet2 = jets[1].momentum();
      const FourMomentum jet3 = jets[2].momentum();

      // Cut on lead jet pT and lead/sublead jet centrality
      if (jet1.pT() < 100*GeV) vetoEvent;
      if (jet1.abseta() > 2.5 || jet2.abseta() > 2.5) vetoEvent;

      // Construct eta & phi distances between 2nd and 3rd jets
      double dEta23 = jet3.eta() - jet2.eta(); ///< Note not abs
      double dPhi23 = jet3.phi() - jet2.phi(); ///< Note not abs
      if (dPhi23 > M_PI)  dPhi23 -= 2*M_PI; ///< @todo Use mapTo... functions?
      if (dPhi23 < -M_PI) dPhi23 += 2*M_PI; ///< @todo Use mapTo... functions?

      // Cut on distance between 2nd and 3rd jets
      const double R23 = add_quad(dPhi23, dEta23);
      if (!inRange(R23, 0.5, 1.5)) vetoEvent;

      // Cut on dijet mass
      const FourMomentum diJet = jet1 + jet2;
      if (diJet.mass() < 220*GeV) vetoEvent;

      // Calc beta and fill histogram (choose central or fwd histo inline)
      double beta = fabs(atan2(dPhi23, sign(jet2.eta())*dEta23));
      ((jet2.abseta() < 0.8) ? _h_hTotD : _h_hTotDF)->fill(beta, event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double width = _h_hTotD->bin(0).xWidth();
      normalize(_h_hTotD, width);
      normalize(_h_hTotDF, width);
    }


  private:

    /// @name Histograms
    Histo1DPtr _h_hTotD;
    Histo1DPtr _h_hTotDF;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1265659);

}
