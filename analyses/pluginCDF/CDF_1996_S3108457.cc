// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/SmearedJets.hh"

namespace Rivet {


  /// @brief CDF properties of high-mass multi-jet events
  class CDF_1996_S3108457 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    CDF_1996_S3108457()
      : Analysis("CDF_1996_S3108457")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      const FinalState fs(-4.2, 4.2);

      FastJets fj(fs, FastJets::CDFJETCLU, 0.7);
      declare(fj, "Jets");

      // Smear Energy and mass with the 10% uncertainty quoted in the paper
      SmearedJets sj_E(fj, [](const Jet& jet){ return P4_SMEAR_MASS_GAUSS(P4_SMEAR_E_GAUSS(jet, 0.1*jet.E()), 0.1*jet.mass()); });
      declare(sj_E, "SmearedJets_E");


      /// Book histograms here, e.g.:
      for (size_t i=0; i<5; ++i) {
        _h_m[i] = bookHisto1D(1+i, 1, 1);
        _h_costheta[i] = bookHisto1D(10+i, 1, 1);
        _h_pT[i] = bookHisto1D(15+i, 1, 1);
      }
      /// @todo Ratios of mass histograms left out: Binning doesn't work out
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get the smeared jets
      Jets SJets = apply<JetAlg>(event, "SmearedJets_E").jets(Cuts::Et > 20.0*GeV, cmpMomByEt);
      if (SJets.size() < 2 || SJets.size() > 6) vetoEvent;

      // Calculate Et, total jet 4 Momentum
      double sumEt(0), sumE(0);
      FourMomentum JS(0,0,0,0);

      foreach(const Jet& jet, SJets) {
        sumEt += jet.Et()*GeV;
        sumE  += jet.E()*GeV;
        JS+=jet.momentum();
      }

      if (sumEt < 420*GeV || sumE > 2000*GeV) vetoEvent;

      double mass = JS.mass();

      LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(JS.betaVec());
      FourMomentum jet0boosted(cms_boost.transform(SJets[0].momentum()));
      double costheta0 = fabs(cos(jet0boosted.theta()));

      if (costheta0 < 2.0/3.0) {
        _h_m[SJets.size()-2]->fill(mass, weight);
      }
      if (mass > 600.0*GeV) _h_costheta[JS.size()-2]->fill(costheta0, weight);
      if (costheta0 < 2.0/3.0 && mass > 600.0*GeV) {
        foreach (const Jet& jet, SJets) {
          _h_pT[SJets.size()-2]->fill(jet.pT(), weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      /// Normalise, scale and otherwise manipulate histograms here
      for (size_t i=0; i<5; ++i) {
        normalize(_h_m[i], 40.0);
        normalize(_h_costheta[i], 2.0);
        normalize(_h_pT[i], 20.0);
      }

    }

    //@}


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_m[5];
    Histo1DPtr _h_costheta[5];
    Histo1DPtr _h_pT[5];

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_1996_S3108457);

}
