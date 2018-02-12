// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Ratios of jet pT spectra, related to ratios of differential jet cross sections
  class CMS_2014_I1298810 : public Analysis {
  public:

    /// Constructor
    CMS_2014_I1298810()
      : Analysis("CMS_2014_I1298810")
    {    }


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      FastJets jetsak5(FinalState(), FastJets::ANTIKT, 0.5);
      declare(jetsak5, "JetsAK5");
      FastJets jetsak7(FinalState(), FastJets::ANTIKT, 0.7);
      declare(jetsak7, "JetsAK7");

      // Histograms
      _h_pt_05_ak5    = bookHisto1D(1, 1, 1);
      _h_pt_05_10_ak5 = bookHisto1D(2, 1, 1);
      _h_pt_10_15_ak5 = bookHisto1D(3, 1, 1);
      _h_pt_15_20_ak5 = bookHisto1D(4, 1, 1);
      _h_pt_20_25_ak5 = bookHisto1D(5, 1, 1);
      _h_pt_25_30_ak5 = bookHisto1D(6, 1, 1);

      _h_pt_05_ak7    = bookHisto1D(7, 1, 1);
      _h_pt_05_10_ak7 = bookHisto1D(8, 1, 1);
      _h_pt_10_15_ak7 = bookHisto1D(9, 1, 1);
      _h_pt_15_20_ak7 = bookHisto1D(10, 1, 1);
      _h_pt_20_25_ak7 = bookHisto1D(11, 1, 1);
      _h_pt_25_30_ak7 = bookHisto1D(12, 1, 1);

      _h_pt_05_ratio    = bookScatter2D(13, 1, 1);
      _h_pt_05_10_ratio = bookScatter2D(14, 1, 1);
      _h_pt_10_15_ratio = bookScatter2D(15, 1, 1);
      _h_pt_15_20_ratio = bookScatter2D(16, 1, 1);
      _h_pt_20_25_ratio = bookScatter2D(17, 1, 1);
      _h_pt_25_30_ratio = bookScatter2D(18, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& jetsak5 = apply<FastJets>(event, "JetsAK5").jetsByPt(56*GeV);
      const Jets& jetsak7 = apply<FastJets>(event, "JetsAK7").jetsByPt(56*GeV);
      if (jetsak5.size() < 1 && jetsak7.size() < 1) vetoEvent;

      const double weight = event.weight();

      // Filling R = 0.5 jets
      foreach(const Jet& jet, jetsak5) {
        if (jet.absrapidity() < 0.5) {
          _h_pt_05_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 1.0) {
          _h_pt_05_10_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 1.5) {
          _h_pt_10_15_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 2.0) {
          _h_pt_15_20_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 2.5) {
          _h_pt_20_25_ak5->fill(jet.pT()/GeV, weight);
        } else if (jet.absrapidity() < 3.0) {
          _h_pt_25_30_ak5->fill(jet.pT()/GeV, weight);
        }
      }


      // Filling R = 0.7 jets
      foreach(const Jet& jet, jetsak7) {
        if (jet.absrapidity() < 0.5) {
          _h_pt_05_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 1.0) {
          _h_pt_05_10_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 1.5) {
          _h_pt_10_15_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 2.0) {
          _h_pt_15_20_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 2.5) {
          _h_pt_20_25_ak7->fill(jet.pT() * GeV, weight);
        } else if (jet.absrapidity() < 3.0) {
          _h_pt_25_30_ak7->fill(jet.pT() * GeV, weight);
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pt_05_ak5,    crossSection()/sumOfWeights());
      scale(_h_pt_05_10_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_10_15_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_15_20_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_20_25_ak5, crossSection()/sumOfWeights());
      scale(_h_pt_25_30_ak5, crossSection()/sumOfWeights());

      scale(_h_pt_05_ak7,    crossSection()/sumOfWeights());
      scale(_h_pt_05_10_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_10_15_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_15_20_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_20_25_ak7, crossSection()/sumOfWeights());
      scale(_h_pt_25_30_ak7, crossSection()/sumOfWeights());

      divide(_h_pt_05_ak5,    _h_pt_05_ak7,    _h_pt_05_ratio);
      divide(_h_pt_05_10_ak5, _h_pt_05_10_ak7, _h_pt_05_10_ratio);
      divide(_h_pt_10_15_ak5, _h_pt_10_15_ak7, _h_pt_10_15_ratio);
      divide(_h_pt_15_20_ak5, _h_pt_15_20_ak7, _h_pt_15_20_ratio);
      divide(_h_pt_20_25_ak5, _h_pt_20_25_ak7, _h_pt_20_25_ratio);
      divide(_h_pt_25_30_ak5, _h_pt_25_30_ak7, _h_pt_25_30_ratio);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_pt_05_ak5, _h_pt_05_10_ak5, _h_pt_10_15_ak5, _h_pt_15_20_ak5, _h_pt_20_25_ak5, _h_pt_25_30_ak5;
    Histo1DPtr _h_pt_05_ak7, _h_pt_05_10_ak7, _h_pt_10_15_ak7, _h_pt_15_20_ak7, _h_pt_20_25_ak7, _h_pt_25_30_ak7;
    Scatter2DPtr _h_pt_05_ratio, _h_pt_05_10_ratio, _h_pt_10_15_ratio, _h_pt_15_20_ratio, _h_pt_20_25_ratio, _h_pt_25_30_ratio;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2014_I1298810);

}
