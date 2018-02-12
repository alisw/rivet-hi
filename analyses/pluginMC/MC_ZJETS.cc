// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  class MC_ZJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_ZJETS(string name = "MC_ZJETS")
      : MC_JetAnalysis(name, 4, "Jets")
	  {
		  _dR=0.2;
		  _lepton=PID::ELECTRON;
	  }


    /// @name Analysis methods
    //@{

    /// Initialize
    void init() {
      FinalState fs;
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      ZFinder zfinder(fs, cut, _lepton, 65*GeV, 115*GeV, _dR, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
      declare(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      _h_Z_jet1_deta = bookHisto1D("Z_jet1_deta", 50, -5, 5);
      _h_Z_jet1_dR = bookHisto1D("Z_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      const FourMomentum& zmom = zfinder.bosons()[0].momentum();

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(_jetptcut);
      if (jets.size() > 0) {
        const double weight = e.weight();
        _h_Z_jet1_deta->fill(zmom.eta()-jets[0].eta(), weight);
        _h_Z_jet1_dR->fill(deltaR(zmom, jets[0].momentum()), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_Z_jet1_deta, crossSection()/picobarn/sumOfWeights());
      scale(_h_Z_jet1_dR, crossSection()/picobarn/sumOfWeights());
      MC_JetAnalysis::finalize();
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
    Histo1DPtr _h_Z_jet1_deta;
    Histo1DPtr _h_Z_jet1_dR;
    //@}

  };



  struct MC_ZJETS_EL : public MC_ZJETS {
    MC_ZJETS_EL() : MC_ZJETS("MC_ZJETS_EL") {
      _dR = 0.2;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_ZJETS_EL_BARE : public MC_ZJETS {
    MC_ZJETS_EL_BARE() : MC_ZJETS("MC_ZJETS_EL_BARE") {
      _dR = 0;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_ZJETS_MU : public MC_ZJETS {
    MC_ZJETS_MU() : MC_ZJETS("MC_ZJETS_MU") {
      _dR = 0.2;
      _lepton = PID::MUON;
    }
  };

  struct MC_ZJETS_MU_BARE : public MC_ZJETS {
    MC_ZJETS_MU_BARE() : MC_ZJETS("MC_ZJETS_MU_BARE") {
      _dR = 0;
      _lepton = PID::MUON;
    }
  };



  // The hooks for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ZJETS);
  DECLARE_RIVET_PLUGIN(MC_ZJETS_EL);
  DECLARE_RIVET_PLUGIN(MC_ZJETS_EL_BARE);
  DECLARE_RIVET_PLUGIN(MC_ZJETS_MU);
  DECLARE_RIVET_PLUGIN(MC_ZJETS_MU_BARE);

}
