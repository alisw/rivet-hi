// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Analysis.hh"

namespace Rivet {


  /// @brief MC validation analysis for W + jets events
  class MC_WJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_WJETS(string name="MC_WJETS")
      : MC_JetAnalysis(name, 4, "Jets")
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
      FastJets jetpro(wfinder.remainingFinalState(), FastJets::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      _h_W_jet1_deta = bookHisto1D("W_jet1_deta", 50, -5.0, 5.0);
      _h_W_jet1_dR = bookHisto1D("W_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const double weight = e.weight();

      const WFinder& wfinder = apply<WFinder>(e, "WFinder");
      if (wfinder.bosons().size() != 1) {
        vetoEvent;
      }
      FourMomentum wmom(wfinder.bosons().front().momentum());

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(_jetptcut);
      if (jets.size() > 0) {
        _h_W_jet1_deta->fill(wmom.eta()-jets[0].eta(), weight);
        _h_W_jet1_dR->fill(deltaR(wmom, jets[0].momentum()), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_W_jet1_deta, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_jet1_dR, crossSection()/picobarn/sumOfWeights());
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
    Histo1DPtr _h_W_jet1_deta;
    Histo1DPtr _h_W_jet1_dR;
    //@}

  };



  struct MC_WJETS_EL : public MC_WJETS {
    MC_WJETS_EL() : MC_WJETS("MC_WJETS_EL") {
      _dR = 0.2;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_WJETS_EL_BARE : public MC_WJETS {
    MC_WJETS_EL_BARE() : MC_WJETS("MC_WJETS_EL_BARE") {
      _dR = 0;
      _lepton = PID::ELECTRON;
    }
  };

  struct MC_WJETS_MU : public MC_WJETS {
    MC_WJETS_MU() : MC_WJETS("MC_WJETS_MU") {
      _dR = 0.2;
      _lepton = PID::MUON;
    }
  };

  struct MC_WJETS_MU_BARE : public MC_WJETS {
    MC_WJETS_MU_BARE() : MC_WJETS("MC_WJETS_MU_BARE") {
      _dR = 0;
      _lepton = PID::MUON;
    }
  };



  // The hooks for the plugin system
  DECLARE_RIVET_PLUGIN(MC_WJETS);
  DECLARE_RIVET_PLUGIN(MC_WJETS_EL);
  DECLARE_RIVET_PLUGIN(MC_WJETS_EL_BARE);
  DECLARE_RIVET_PLUGIN(MC_WJETS_MU);
  DECLARE_RIVET_PLUGIN(MC_WJETS_MU_BARE);

}
