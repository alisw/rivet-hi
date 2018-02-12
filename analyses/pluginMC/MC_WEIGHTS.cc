// -*- C++ -*-
#include "Rivet/Analysis.hh"

namespace Rivet {

  /// @brief Analysis for the generated cross section
  class MC_WEIGHTS : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_WEIGHTS()
      : Analysis("MC_WEIGHTS")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      /// @todo Convert to Scatter1D or Counter
      _h_weight_100 = bookHisto1D("weight_100", 200, -100.0, 100.0);
      _h_weight_10  = bookHisto1D("weight_10", 200, -10.0, 10.0);
      _h_logweight_pos  = bookHisto1D("logweight_pos", logspace(100, 0.1, 10000.0));
      _h_logweight_neg  = bookHisto1D("logweight_neg", logspace(100, 0.1, 10000.0));

      _h_xsfraction_neg   = bookScatter2D("xsfraction_neg");

      _sow_pos = _sow_neg = _nevts = 0.;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double w = event.weight();

      _nevts += 1.0;
      _h_weight_100->fill(w, 1.0);
      _h_weight_10->fill(w, 1.0);
      if (w<0.0) {
        _h_logweight_neg->fill(fabs(w), 1.0);
        _sow_neg += fabs(w);
      }
      else {
        _h_logweight_pos->fill(w, 1.0);
        _sow_pos += w;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_weight_100, 1.0/_nevts);
      scale(_h_weight_10, 1.0/_nevts);
      scale(_h_logweight_pos, 1.0/_nevts);
      scale(_h_logweight_neg, 1.0/_nevts);
      /// @todo correct unc estimate:
      _h_xsfraction_neg->addPoint(0, _sow_neg/(_sow_neg+_sow_pos), 0.5, 0.0);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Scatter2DPtr _h_xsfraction_neg;
    Histo1DPtr _h_weight_100, _h_weight_10, _h_logweight_pos, _h_logweight_neg;
    double _sow_pos, _sow_neg, _nevts;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WEIGHTS);

}
