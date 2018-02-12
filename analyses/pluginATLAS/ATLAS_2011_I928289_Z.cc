// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class ATLAS_2011_I928289_Z : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I928289_Z()
      : Analysis("ATLAS_2011_I928289_Z")
    {
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;

      Cut cut = (Cuts::pT >= 20.0*GeV);

      ZFinder zfinder_ee_bare(   fs, cut, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      ZFinder zfinder_ee_dressed(fs, cut, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      ZFinder zfinder_mm_bare(   fs, cut, PID::MUON    , 66.0*GeV, 116.0*GeV, 0.0, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      ZFinder zfinder_mm_dressed(fs, cut, PID::MUON    , 66.0*GeV, 116.0*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);

      declare(zfinder_ee_bare   , "ZFinder_ee_bare"   );
      declare(zfinder_ee_dressed, "ZFinder_ee_dressed");
      declare(zfinder_mm_bare   , "ZFinder_mm_bare"   );
      declare(zfinder_mm_dressed, "ZFinder_mm_dressed");

      // y(Z) cross-section dependence
      _h_Z_y_ee_bare     = bookHisto1D(1, 1, 1);
      _h_Z_y_ee_dressed  = bookHisto1D(1, 1, 2);
      _h_Z_y_mm_bare     = bookHisto1D(1, 1, 3);
      _h_Z_y_mm_dressed  = bookHisto1D(1, 1, 4);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder_ee_bare     = apply<ZFinder>(event, "ZFinder_ee_bare"   );
      const ZFinder& zfinder_ee_dressed  = apply<ZFinder>(event, "ZFinder_ee_dressed");
      const ZFinder& zfinder_mm_bare     = apply<ZFinder>(event, "ZFinder_mm_bare"   );
      const ZFinder& zfinder_mm_dressed  = apply<ZFinder>(event, "ZFinder_mm_dressed");

      const double weight = event.weight();
      fillPlots1D(zfinder_ee_bare   , _h_Z_y_ee_bare   , weight);
      fillPlots1D(zfinder_ee_dressed, _h_Z_y_ee_dressed, weight);
      fillPlots1D(zfinder_mm_bare   , _h_Z_y_mm_bare   , weight);
      fillPlots1D(zfinder_mm_dressed, _h_Z_y_mm_dressed, weight);

    }


    void fillPlots1D(const ZFinder& zfinder, Histo1DPtr hist, double weight) {
      if (zfinder.bosons().size() != 1) return;
      const FourMomentum zmom = zfinder.bosons()[0].momentum();
      hist->fill(zmom.absrap(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double sumw(sumOfWeights());
      MSG_DEBUG("Cross-Section/pb: " << xs_pb      );
      MSG_DEBUG("Sum of weights  : " << sumw       );
      MSG_DEBUG("nEvents         : " << numEvents());

      // Normalise, scale and otherwise manipulate histograms here
      const double sf(0.5 * xs_pb / sumw); // 0.5 accounts for rapidity bin width
      scale(_h_Z_y_ee_bare   , sf);
      scale(_h_Z_y_ee_dressed, sf);
      scale(_h_Z_y_mm_bare   , sf);
      scale(_h_Z_y_mm_dressed, sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Z_y_ee_bare;
    Histo1DPtr _h_Z_y_ee_dressed;
    Histo1DPtr _h_Z_y_mm_bare;
    Histo1DPtr _h_Z_y_mm_dressed;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I928289_Z);

}
