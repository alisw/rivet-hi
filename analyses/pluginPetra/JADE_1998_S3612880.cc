// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class JADE_1998_S3612880 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1998_S3612880);


    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::pT > 0.1*GeV);
      declare(cfs, "CFS");
      declare(FastJets(cfs, FastJets::DURHAM, 0.7), "DurhamJets");

      // Thrust
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Histos
      int offset = 0;
      switch (int(sqrtS()/GeV)) {

        case 44:
          offset = 0;
          _h_thrust  = bookHisto1D( 2+offset, 1, 1);
          _h_MH = bookHisto1D( 3 + offset, 1, 1);
          _h_BT = bookHisto1D( 4 + offset, 1, 1);
          _h_BW = bookHisto1D( 5 + offset, 1, 1);
          _h_y23 = bookHisto1D(10, 1, 1);
          break;
        case 35:
          offset = 4;
          _h_thrust  = bookHisto1D( 2+offset, 1, 1);
          _h_MH = bookHisto1D( 3 + offset, 1, 1);
          _h_BT = bookHisto1D( 4 + offset, 1, 1);
          _h_BW = bookHisto1D( 5 + offset, 1, 1);
          _h_y23 = bookHisto1D(11, 1, 1);
          break;
        case 22:
          _h_y23 = bookHisto1D(12, 1, 1);
          break;
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // JADE hadronic event selection
      if (cfs.particles().size() < 3 ) vetoEvent;

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const Vector3 & thrustAxis = thrust.thrustAxis ();
      double theta = thrustAxis.theta();
      if ( fabs(cos(theta)) >= 0.8 ) {
        MSG_DEBUG("Failed thrust angle cut: " << fabs(cos(theta)));
        vetoEvent;
      }
      /// @todo Evis, pmiss, pbal

      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");

      const double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);

      // Make sure we don't run into a segfault by trying to fill non-existing histos
      int s = int(sqrtS()/GeV);
      if (s == 44 || s == 35) {
        _h_thrust->fill(1. - thrust.thrust(), weight);
        _h_MH->fill(sqrt(hemi.scaledM2high()), weight);
        _h_BT->fill(hemi.Bsum(), weight);
        _h_BW->fill(hemi.Bmax(), weight);
      }
      _h_y23->fill(y23, weight);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Make sure we don't try to normalise non-existing histos
      const int s = int(sqrtS()/GeV);
      if (s == 44 || s == 35) {
        normalize(_h_thrust);
        normalize(_h_MH);
        normalize(_h_BT);
        normalize(_h_BW);
      }
      normalize(_h_y23);
    }

    //@}


  private:

    Histo1DPtr _h_thrust, _h_MH, _h_BT, _h_BW, _h_y23;


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1998_S3612880);

}
