// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {



  /// @brief MC validation analysis for Z events
  class MC_ZINC : public Analysis {
  public:

    /// Default constructor
    MC_ZINC(string name="MC_ZINC")
		 : Analysis(name) {
		 _dR=0.2;
		 _lepton=PID::ELECTRON;
	 }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      ZFinder zfinder(fs, cut, _lepton, 65.0*GeV, 115.0*GeV, _dR, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
      declare(zfinder, "ZFinder");

      _h_Z_mass = bookHisto1D("Z_mass", 50, 66.0, 116.0);
      _h_Z_pT = bookHisto1D("Z_pT", logspace(100, 1.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));
      _h_Z_pT_peak = bookHisto1D("Z_pT_peak", 25, 0.0, 25.0);
      _h_Z_y = bookHisto1D("Z_y", 40, -4.0, 4.0);
      _h_Z_phi = bookHisto1D("Z_phi", 25, 0.0, TWOPI);
      _h_lepton_pT = bookHisto1D("lepton_pT", logspace(100, 10.0, 0.25*(sqrtS()>0.?sqrtS():14000.)/GeV));
      _h_lepton_eta = bookHisto1D("lepton_eta", 40, -4.0, 4.0);

    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      const double weight = e.weight();

      FourMomentum zmom(zfinder.bosons()[0].momentum());
      _h_Z_mass->fill(zmom.mass()/GeV, weight);
      _h_Z_pT->fill(zmom.pT()/GeV, weight);
      _h_Z_pT_peak->fill(zmom.pT()/GeV, weight);
      _h_Z_y->fill(zmom.rapidity(), weight);
      _h_Z_phi->fill(zmom.phi(), weight);
      for (const Particle& l : zfinder.constituents()) {
        _h_lepton_pT->fill(l.pT()/GeV, weight);
        _h_lepton_eta->fill(l.eta(), weight);
      }
    }


    /// Finalize
    void finalize() {
      const double s = crossSection()/picobarn/sumOfWeights();
      scale(_h_Z_mass, s);
      scale(_h_Z_pT, s);
      scale(_h_Z_pT_peak, s);
      scale(_h_Z_y, s);
      scale(_h_Z_phi, s);
      scale(_h_lepton_pT, s);
      scale(_h_lepton_eta, s);
    }

    //@}


  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    //@{
    double _dR;
    PdgId _lepton;
    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Z_mass;
    Histo1DPtr _h_Z_pT;
    Histo1DPtr _h_Z_pT_peak;
    Histo1DPtr _h_Z_y;
    Histo1DPtr _h_Z_phi;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;
    //@}

  };



  struct MC_ZINC_EL : public MC_ZINC {
    MC_ZINC_EL() : MC_ZINC("MC_ZINC_EL") {
      _dR = 0.2;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_ZINC_EL_BARE : public MC_ZINC {
    MC_ZINC_EL_BARE() : MC_ZINC("MC_ZINC_EL_BARE") {
      _dR = 0;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_ZINC_MU : public MC_ZINC {
    MC_ZINC_MU() : MC_ZINC("MC_ZINC_MU") {
      _dR = 0.2;
      _lepton = PID::MUON;
    }
  };

  struct MC_ZINC_MU_BARE : public MC_ZINC {
    MC_ZINC_MU_BARE() : MC_ZINC("MC_ZINC_MU_BARE") {
      _dR = 0;
      _lepton = PID::MUON;
    }
  };



  // The hooks for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ZINC);
  DECLARE_RIVET_PLUGIN(MC_ZINC_EL);
  DECLARE_RIVET_PLUGIN(MC_ZINC_EL_BARE);
  DECLARE_RIVET_PLUGIN(MC_ZINC_MU);
  DECLARE_RIVET_PLUGIN(MC_ZINC_MU_BARE);


}
