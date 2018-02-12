// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {



  /// @brief MC validation analysis for inclusive W events
  class MC_WINC : public Analysis {
  public:

    /// Default constructor
    MC_WINC(string name="MC_WINC")
      : Analysis(name)
    {
		 _dR=0.2;
		 _lepton=PID::ELECTRON;
	 }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      WFinder wfinder(fs, Cuts::abseta < 3.5 && Cuts::pT > 25*GeV, _lepton, 60.0*GeV, 100.0*GeV, 25.0*GeV, _dR);
      declare(wfinder, "WFinder");

      double sqrts = sqrtS()>0. ? sqrtS() : 14000.;
      _h_W_mass = bookHisto1D("W_mass", 50, 55.0, 105.0);
      _h_W_pT = bookHisto1D("W_pT", logspace(100, 1.0, 0.5*sqrts));
      _h_W_pT_peak = bookHisto1D("W_pT_peak", 25, 0.0, 125.0);
      _h_W_y = bookHisto1D("W_y", 40, -4.0, 4.0);
      _h_W_phi = bookHisto1D("W_phi", 25, 0.0, TWOPI);
      _h_Wplus_pT = bookHisto1D("Wplus_pT", logspace(25, 10.0, 0.5*sqrts));
      _h_Wminus_pT = bookHisto1D("Wminus_pT", logspace(25, 10.0, 0.5*sqrts));
      _h_lepton_pT = bookHisto1D("lepton_pT", logspace(100, 10.0, 0.25*sqrts));
      _h_lepton_eta = bookHisto1D("lepton_eta", 40, -4.0, 4.0);
      _htmp_dsigminus_deta = bookHisto1D("lepton_dsigminus_deta", 20, 0.0, 4.0);
      _htmp_dsigplus_deta  = bookHisto1D("lepton_dsigplus_deta", 20, 0.0, 4.0);

      _h_asym = bookScatter2D("W_chargeasymm_eta");
      _h_asym_pT = bookScatter2D("W_chargeasymm_pT");
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const WFinder& wfinder = apply<WFinder>(e, "WFinder");
      if (wfinder.bosons().size() != 1) {
        vetoEvent;
      }

      double charge3_x_eta = 0;
      int charge3 = 0;
      FourMomentum emom;
      FourMomentum wmom(wfinder.bosons().front().momentum());
      _h_W_mass->fill(wmom.mass(), weight);
      _h_W_pT->fill(wmom.pT(), weight);
      _h_W_pT_peak->fill(wmom.pT(), weight);
      _h_W_y->fill(wmom.rapidity(), weight);
      _h_W_phi->fill(wmom.phi(), weight);
      Particle l=wfinder.constituentLeptons()[0];
      _h_lepton_pT->fill(l.pT(), weight);
      _h_lepton_eta->fill(l.eta(), weight);
      if (PID::threeCharge(l.pid()) != 0) {
        emom = l.momentum();
        charge3_x_eta = PID::threeCharge(l.pid()) * emom.eta();
        charge3 = PID::threeCharge(l.pid());
      }
      assert(charge3_x_eta != 0);
      assert(charge3!=0);
      if (emom.Et() > 30/GeV) {
        if (charge3_x_eta < 0) {
          _htmp_dsigminus_deta->fill(emom.eta(), weight);
        } else {
          _htmp_dsigplus_deta->fill(emom.eta(), weight);
        }
      }
      if (charge3 < 0) {
        _h_Wminus_pT->fill(wmom.pT(), weight);
      } else {
        _h_Wplus_pT->fill(wmom.pT(), weight);
      }
    }


    /// Finalize
    void finalize() {
      scale(_h_W_mass, crossSection()/sumOfWeights());
      scale(_h_W_pT, crossSection()/sumOfWeights());
      scale(_h_W_pT_peak, crossSection()/sumOfWeights());
      scale(_h_W_y, crossSection()/sumOfWeights());
      scale(_h_W_phi, crossSection()/sumOfWeights());
      scale(_h_lepton_pT, crossSection()/sumOfWeights());
      scale(_h_lepton_eta, crossSection()/sumOfWeights());

      // Construct asymmetry: (dsig+/deta - dsig-/deta) / (dsig+/deta + dsig-/deta) for each Et region
      divide(*_htmp_dsigplus_deta - *_htmp_dsigminus_deta,
             *_htmp_dsigplus_deta + *_htmp_dsigminus_deta,
             _h_asym);

      // // W charge asymmetry vs. pTW: dsig+/dpT / dsig-/dpT
      divide(_h_Wplus_pT, _h_Wminus_pT,
             _h_asym_pT);

      scale(_h_Wplus_pT, crossSection()/picobarn/sumOfWeights());
      scale(_h_Wminus_pT, crossSection()/picobarn/sumOfWeights());
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
    Histo1DPtr _h_W_mass;
    Histo1DPtr _h_W_pT;
    Histo1DPtr _h_W_pT_peak;
    Histo1DPtr _h_W_y;
    Histo1DPtr _h_W_phi;
    Histo1DPtr _h_Wplus_pT;
    Histo1DPtr _h_Wminus_pT;
    Histo1DPtr _h_lepton_pT;
    Histo1DPtr _h_lepton_eta;

    Histo1DPtr _htmp_dsigminus_deta;
    Histo1DPtr _htmp_dsigplus_deta;

    Scatter2DPtr _h_asym;
    Scatter2DPtr _h_asym_pT;
    //@}

  };



  struct MC_WINC_EL : public MC_WINC {
    MC_WINC_EL() : MC_WINC("MC_WINC_EL") {
      _dR = 0.2;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_WINC_EL_BARE : public MC_WINC {
    MC_WINC_EL_BARE() : MC_WINC("MC_WINC_EL_BARE") {
      _dR = 0;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_WINC_MU : public MC_WINC {
    MC_WINC_MU() : MC_WINC("MC_WINC_MU") {
      _dR = 0.2;
      _lepton = PID::MUON;
    }
  };

  struct MC_WINC_MU_BARE : public MC_WINC {
    MC_WINC_MU_BARE() : MC_WINC("MC_WINC_MU_BARE") {
      _dR = 0;
      _lepton = PID::MUON;
    }
  };



  // The hooks for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WINC);
  DECLARE_RIVET_PLUGIN(MC_WINC_EL);
  DECLARE_RIVET_PLUGIN(MC_WINC_EL_BARE);
  DECLARE_RIVET_PLUGIN(MC_WINC_MU);
  DECLARE_RIVET_PLUGIN(MC_WINC_MU_BARE);


}
