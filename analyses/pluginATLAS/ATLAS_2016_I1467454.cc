// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief High-mass Drell-Yan at 8 TeV
  class ATLAS_2016_I1467454 : public Analysis {
  public:

    /// Constructor
    ATLAS_2016_I1467454(const string& name="ATLAS_2016_I1467454")
      : Analysis(name)
    {
      _mode = 0; // use electron channel by default
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 30*GeV;
      ZFinder zfinder(fs, cuts, _mode? PID::MUON : PID::ELECTRON, 116*GeV, 1500*GeV, 0.1);
      declare(zfinder, "ZFinder");

      size_t ch = _mode? 11 : 0; // offset
      _hist_mll = bookHisto1D(18 + ch, 1, 1);

      vector<double> mll_bins = { 116., 150., 200., 300., 500., 1500. };
      for (size_t i = 0; i < (mll_bins.size() - 1); ++i) {
        _hist_rap.addHistogram( mll_bins[i], mll_bins[i+1], bookHisto1D(19 + ch + i, 1, 1));
        _hist_deta.addHistogram(mll_bins[i], mll_bins[i+1], bookHisto1D(24 + ch + i, 1, 1));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1)  vetoEvent;

      const Particle z0  = zfinder.bosons()[0];
      /// @todo Could use z0.constituents()
      const Particle el1 = zfinder.constituentLeptons()[0];
      const Particle el2 = zfinder.constituentLeptons()[1];

      if (el1.pT() > 40*GeV || el2.pT() > 40*GeV) {
        const double mass = z0.mass();
        const double weight = event.weight();
        _hist_mll->fill(mass/GeV, weight);
        _hist_rap. fill(mass/GeV, z0.absrap(), weight);
        _hist_deta.fill(mass/GeV, deltaEta(el1,el2), weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection()/sumOfWeights();
      scale(_hist_mll, sf);
      _hist_rap.scale(sf*0.5,  this);
      _hist_deta.scale(sf*0.5, this);

    }

    //@}


    /// Choose to work in electron or muon mode
    size_t _mode;


    /// @name Histograms
    //@{
    Histo1DPtr _hist_mll;
    BinnedHistogram<double> _hist_rap, _hist_deta;
    //@}

  };


  /// High-mass Drell-Yan at 8 TeV, electron channel
  struct ATLAS_2016_I1467454_EL : public ATLAS_2016_I1467454 {
    ATLAS_2016_I1467454_EL() : ATLAS_2016_I1467454("ATLAS_2016_I1467454_EL") { _mode = 0; }
  };


  /// High-mass Drell-Yan at 8 TeV, muon channel
  struct ATLAS_2016_I1467454_MU : public ATLAS_2016_I1467454 {
    ATLAS_2016_I1467454_MU() : ATLAS_2016_I1467454("ATLAS_2016_I1467454_MU") { _mode = 1; }
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1467454);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1467454_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1467454_MU);

}
