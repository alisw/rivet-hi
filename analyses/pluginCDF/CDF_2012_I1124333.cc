// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  /// @ CDF Run II Z \f$ p_\perp \f$ in Drell-Yan events
  /// @author Simone Amoroso
  class CDF_2012_I1124333 : public Analysis {
  public:

    /// Constructor
    CDF_2012_I1124333()
      : Analysis("CDF_2012_I1124333")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      ///  Initialise and register projections here
      ZFinder zfinder(FinalState(), Cuts::open(), PID::ELECTRON, 66*GeV, 116*GeV, 0.0, ZFinder::NOCLUSTER);
      declare(zfinder, "ZFinder");


      ///  Book histograms here, e.g.:
      //      _hist_z_xs = bookHisto1D(1, 1, 1);
      _hist_zpt = bookHisto1D(2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      /// @todo Do the event by event analysis here
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1) {
        MSG_DEBUG("Num e+ e- pairs found = " << zfinder.bosons().size());
	vetoEvent;
      }
      const FourMomentum& pZ = zfinder.bosons()[0].momentum();
      if (pZ.mass2() < 0) {
	MSG_DEBUG("Negative Z mass**2 = " << pZ.mass2()/GeV2 << "!");
	vetoEvent;
      }

      MSG_DEBUG("Dilepton mass = " << pZ.mass()/GeV << " GeV");
      _hist_zpt->fill(pZ.pT(), weight);
      //      _hist_z_xs->fill(1, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_zpt, crossSection()/picobarn/sumOfWeights());
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here


    /// @name Histograms
    //@{
    Histo1DPtr _hist_zpt;
    //    Histo1DPtr _hist_z_xs;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2012_I1124333);

}
