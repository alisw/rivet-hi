// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// Jet mass as a function of ystar
  class ATLAS_2014_I1268975 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1268975()
      : Analysis("ATLAS_2014_I1268975")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs,"FinalState");

      FastJets fj04(fs,  FastJets::ANTIKT, 0.4);
      fj04.useInvisibles();
      declare(fj04, "AntiKT04");

      FastJets fj06(fs,  FastJets::ANTIKT, 0.6);
      fj06.useInvisibles();
      declare(fj06, "AntiKT06");

      double ystarbins[] = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};

      size_t massDsOffset(0);
      for (size_t alg = 0; alg < 2; ++alg) {
        for (size_t i = 0; i < 6; ++i) {
          _mass[alg].addHistogram(ystarbins[i], ystarbins[i+1], bookHisto1D(1 + massDsOffset, 1, i+1));
        }
        massDsOffset += 1;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Jets jetAr[2];
      jetAr[AKT4] = apply<FastJets>(event, "AntiKT04").jetsByPt(50*GeV);
      jetAr[AKT6] = apply<FastJets>(event, "AntiKT06").jetsByPt(50*GeV);

      // Loop over jet "radii" used in analysis
      for (size_t alg = 0; alg < 2; ++alg) {

        // Identify dijets
        vector<FourMomentum> leadjets;
        foreach (const Jet& jet, jetAr[alg]) {
          if (jet.absrap() < 3.0 && leadjets.size() < 2) {
            if (leadjets.empty() && jet.pT() < 100*GeV) continue;
            leadjets.push_back(jet.momentum());
          }
        }

        // Make sure we have a leading jet with pT > 100 GeV and a second to leading jet with pT > 50 GeV
        if (leadjets.size() < 2) {
          MSG_DEBUG("Could not find two suitable leading jets");
          continue;
        }

        const double y1    = leadjets[0].rapidity();
        const double y2    = leadjets[1].rapidity();
        const double ystar = fabs(y1-y2) / 2.;
        const double m     = (leadjets[0] + leadjets[1]).mass();

        // Fill mass histogram
        _mass[alg].fill(ystar, m/TeV, event.weight());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t alg = 0; alg < 2; ++alg) {
        _mass[alg].scale(crossSectionPerEvent()/picobarn, this);
      }
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here
    enum Alg { AKT4=0, AKT6=1 };

    /// The di-jet mass spectrum binned in rapidity for akt6 and akt4 jets (array index is jet type from enum above)
    BinnedHistogram<double> _mass[2];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1268975);

}
