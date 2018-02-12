// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"

namespace Rivet {


  struct Plots {
    string label;

    Histo1DPtr h_dy;
    Histo1DPtr h_mjj;
    Histo1DPtr h_njets;
    Histo1DPtr h_dphijj;
    Histo1DPtr h_ptbal;

    Histo1DPtr h_jetveto_mjj_veto;
    Histo1DPtr h_jetveto_mjj_inc;
    Histo1DPtr h_jetveto_dy_veto;
    Histo1DPtr h_jetveto_dy_inc;

    Histo1DPtr h_ptbaleff_mjj_veto;
    Histo1DPtr h_ptbaleff_mjj_inc;
    Histo1DPtr h_ptbaleff_dy_veto;
    Histo1DPtr h_ptbaleff_dy_inc;

    Profile1DPtr p_avgnjets_dy;
    Profile1DPtr p_avgnjets_mjj;
  };


  struct Variables {

    Variables(const vector<const Jet*>& jets, const Particle* lep1, const Particle* lep2) {
      FourMomentum j1 = jets.at(0)->momentum();
      FourMomentum j2 = jets.at(1)->momentum();
      jet1pt = j1.pT();
      jet2pt = j2.pT();
      assert(jet1pt > jet2pt);

      zpt = (lep1->mom() + lep2->mom()).pT();

      deltay = fabs(j1.rapidity() - j2.rapidity());
      mjj = (j1 + j2).mass();
      deltaphijj = deltaPhi(j1, j2) / PI;

      FourMomentum gapjet(0., 0., 0., 0.);
      ngapjets = _getNumGapJets(jets, gapjet);

      double ptbal_vec = (j1 + j2 + lep1->mom() + lep2->mom()).pT();
      double ptbal_sc = j1.pT() + j2.pT() + lep1->pT() + lep2->pT();
      ptbalance2 = ptbal_vec / ptbal_sc;

      double ptbal3_vec = (j1 + j2 + gapjet + lep1->mom() + lep2->mom()).pT();
      double ptbal3_sc = j1.pT() + j2.pT() + gapjet.pT() + lep1->pT() + lep2->pT();
      ptbalance3 = ptbal3_vec / ptbal3_sc;

      pass_jetveto = gapjet.pT() < 25.0*GeV;
      pass_ptbaleff = ptbalance2 < 0.15;
    }


    double jet1pt;
    double jet2pt;
    double zpt;

    double deltay;
    double mjj;
    double deltaphijj;
    double ptbalance2;
    double ptbalance3;
    int ngapjets;

    double dilepton_dr;

    bool pass_jetveto;
    bool pass_ptbaleff;


  private:

    bool _isBetween(const Jet* probe, const Jet* boundary1, const Jet* boundary2) {
      double y_p = probe->rapidity();
      double y_b1 = boundary1->rapidity();
      double y_b2 = boundary2->rapidity();

      double y_min = std::min(y_b1, y_b2);
      double y_max = std::max(y_b1, y_b2);

      if (y_p > y_min && y_p < y_max) return true;
      else return false;
    }

    int _getNumGapJets(const vector<const Jet*>& jets, FourMomentum& thirdJet) {
      if (jets.size() < 2) return 0;
      // The vector of jets is already sorted by pT. So the boundary jets will be the first two.
      const Jet* bj1 = jets.at(0);
      const Jet* bj2 = jets.at(1);

      int n_between = 0;
      // Start loop at the 3rd hardest pT jet
      for (size_t i = 2; i < jets.size(); ++i) {
        const Jet* j = jets.at(i);
        // If this jet is between the boundary jets and is hard enough, increment counter
        if (_isBetween(j, bj1, bj2)) {
          if (n_between == 0) thirdJet = j->momentum();
          ++n_between;
        }
      }
      return n_between;
    }

  };



  class ATLAS_2014_I1279489 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1279489()
      : Analysis("ATLAS_2014_I1279489")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs(-5.0, 5.0);

      IdentifiedFinalState photon_fs(fs);
      photon_fs.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState electron_fs(fs);
      electron_fs.acceptIdPair(PID::ELECTRON);

      IdentifiedFinalState muon_fs(fs);
      muon_fs.acceptIdPair(PID::MUON);

      DressedLeptons dressed_electrons(photon_fs, electron_fs, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(dressed_electrons, "DressedElectrons");

      DressedLeptons dressed_muons(photon_fs, muon_fs, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 25*GeV);
      declare(dressed_muons, "DressedMuons");

      FastJets jets(fs, FastJets::ANTIKT, 0.4);
      declare(jets, "Jets");

      initialisePlots(baseline_plots, "baseline");
      initialisePlots(highpt_plots, "highpt");
      initialisePlots(search_plots, "search");
      initialisePlots(control_plots, "control");
      initialisePlots(highmass_plots, "highmass");
    }


    void initialisePlots(Plots& plots, const string& phase_space){
      /****************************************
       * Plot labeling:                       *
       * format = d0_-x0_-y0_                 *
       * d01 = baseline fiducial region       *
       * d02 = high-pt fiducial region        *
       * d03 = search fiducial region         *
       * d04 = control fiducial region        *
       * d05 = high-mass fiducial region      *
       *                                      *
       * x01 = mjj on x-axis                  *
       * x02 = delta-y on x-axis              *
       * x03 = njets on x-axis                *
       * x04 = dphijj on x-axis               *
       * x05 = ptbalance on x-axis            *
       *                                      *
       * y01 = differential cross-section     *
       * y02 = jet veto efficiency            *
       * y03 = ptbalance efficiency           *
       * y04 = average njets                  *
       ****************************************/
      plots.label = phase_space;

      if (phase_space=="baseline") {
        plots.h_mjj = bookHisto1D(1, 1, 1);
        plots.h_dy = bookHisto1D(1, 2, 1);

        plots.h_jetveto_mjj_veto = bookHisto1D("jetveto_mjj_baseline_veto", refData(1,1,2));
        plots.h_jetveto_mjj_inc = bookHisto1D("jetveto_mjj_baseline_inc", refData(1,1,2));
        plots.h_jetveto_dy_veto = bookHisto1D("jetveto_dy_baseline_veto", refData(1,2,2));
        plots.h_jetveto_dy_inc = bookHisto1D("jetveto_dy_baseline_inc", refData(1,2,2));

        plots.h_ptbaleff_mjj_veto = bookHisto1D("ptbaleff_mjj_baseline_veto", refData(1,1,3));
        plots.h_ptbaleff_mjj_inc = bookHisto1D("ptbaleff_mjj_baseline_inc", refData(1,1,3));
        plots.h_ptbaleff_dy_veto = bookHisto1D("ptbaleff_dy_baseline_veto", refData(1,2,3));
        plots.h_ptbaleff_dy_inc = bookHisto1D("ptbaleff_dy_baseline_inc", refData(1,2,3));

        plots.p_avgnjets_mjj = bookProfile1D(1,1,4);
        plots.p_avgnjets_dy = bookProfile1D(1,2,4);
      }

      if (phase_space=="highpt") {
        plots.h_mjj = bookHisto1D(2, 1, 1);
        plots.h_dy = bookHisto1D(2, 2, 1);

        plots.h_jetveto_mjj_veto = bookHisto1D("jetveto_mjj_highpt_veto", refData(2,1,2));
        plots.h_jetveto_mjj_inc = bookHisto1D("jetveto_mjj_highpt_inc", refData(2,1,2));
        plots.h_jetveto_dy_veto = bookHisto1D("jetveto_dy_highpt_veto", refData(2,2,2));
        plots.h_jetveto_dy_inc = bookHisto1D("jetveto_dy_highpt_inc", refData(2,2,2));

        plots.h_ptbaleff_mjj_veto = bookHisto1D("ptbaleff_mjj_highpt_veto", refData(2,1,3));
        plots.h_ptbaleff_mjj_inc = bookHisto1D("ptbaleff_mjj_highpt_inc", refData(2,1,3));
        plots.h_ptbaleff_dy_veto = bookHisto1D("ptbaleff_dy_highpt_veto", refData(2,2,3));
        plots.h_ptbaleff_dy_inc = bookHisto1D("ptbaleff_dy_highpt_inc", refData(2,2,3));

        plots.p_avgnjets_mjj = bookProfile1D(2,1,4);
        plots.p_avgnjets_dy = bookProfile1D(2,2,4);
      }

      if (phase_space=="search") {
        plots.h_mjj = bookHisto1D(3,1,1);
        plots.h_dy = bookHisto1D(3,2,1);
      }

      if (phase_space=="control") {
        plots.h_mjj = bookHisto1D(4,1,1);
        plots.h_dy = bookHisto1D(4,2,1);
      }

      if (phase_space=="highmass") {
        plots.h_njets = bookHisto1D(5, 3, 1);
        plots.h_dphijj = bookHisto1D(5, 4, 1);
        plots.h_ptbal = bookHisto1D(5, 5, 1);
      }
    }



    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Make sure that we have a Z-candidate:
      const Particle *lep1 = NULL, *lep2 = NULL;
      //
      const vector<DressedLepton>& muons = apply<DressedLeptons>(event, "DressedMuons").dressedLeptons();
      if (muons.size() == 2) {
        const FourMomentum dimuon = muons[0].mom() + muons[1].mom();
        if ( inRange(dimuon.mass()/GeV, 81.0, 101.0) && muons[0].threeCharge() != muons[1].threeCharge() ) {
          lep1 = &muons[0];
          lep2 = &muons[1];
        }
      }
      //
      const vector<DressedLepton>& electrons = apply<DressedLeptons>(event, "DressedElectrons").dressedLeptons();
      if (electrons.size() == 2) {
        const FourMomentum dielectron = electrons[0].mom() + electrons[1].mom();
        if ( inRange(dielectron.mass()/GeV, 81.0, 101.0) && electrons[0].threeCharge() != electrons[1].threeCharge() ) {
          if (lep1 && lep2) {
            MSG_INFO("Found Z candidates using both electrons and muons! Continuing with the muon-channel candidate");
          } else {
            lep1 = &electrons[0];
            lep2 = &electrons[1];
          }
        }
      }
      // If there's no Z-candidate, we won't use this event:
      if (!lep1 || !lep2) vetoEvent;


      // Do lepton-jet overlap removal:
      vector<const Jet*> good_jets;
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::absrap < 4.4);
      foreach(const Jet& j, jets) {
        bool nearby_lepton = false;
        foreach (const Particle& m, muons)
          if (deltaR(j, m) < 0.3) nearby_lepton = true;
        foreach (const Particle& e, electrons)
          if (deltaR(j, e) < 0.3) nearby_lepton = true;
        if (!nearby_lepton)
          good_jets.push_back(&j);
      }
      // If we don't have at least 2 good jets, we won't use this event.
      if (good_jets.size() < 2) vetoEvent;


      // Plotting, using variables and histo classes calculated by the Variables object constructor
      Variables vars(good_jets, lep1, lep2);
      bool pass_baseline = (vars.jet1pt > 55.0*GeV && vars.jet2pt > 45.0*GeV);
      bool pass_highpt = (vars.jet1pt > 85.0*GeV && vars.jet2pt > 75.0*GeV);
      bool pass_highmass = (pass_baseline && vars.mjj > 1000.0*GeV);
      bool pass_search = (pass_baseline && vars.zpt > 20.0*GeV && vars.ngapjets == 0 && vars.ptbalance2 < 0.15 && vars.mjj > 250.0*GeV);
      bool pass_control = (pass_baseline && vars.zpt > 20.0*GeV && vars.ngapjets > 0 && vars.ptbalance3 < 0.15 && vars.mjj > 250.0*GeV);
      //
      const double weight = event.weight();
      if (pass_baseline) fillPlots(vars, baseline_plots, "baseline", weight);
      if (pass_highpt) fillPlots(vars, highpt_plots, "highpt", weight);
      if (pass_highmass) fillPlots(vars, highmass_plots, "highmass", weight);
      if (pass_search) fillPlots(vars, search_plots, "search", weight);
      if (pass_control) fillPlots(vars, control_plots, "control", weight);
    }


    void fillPlots(const Variables& vars, Plots& plots, string phase_space, double weight) {
      if (phase_space == "baseline" || phase_space == "highpt" || phase_space == "search" || phase_space == "control") {
        plots.h_dy->fill(vars.deltay, weight);
        plots.h_mjj->fill(vars.mjj, weight);
      }

      if (phase_space == "baseline" || phase_space == "highpt") {
        if (vars.pass_jetveto) {
          plots.h_jetveto_dy_veto->fill(vars.deltay, weight);
          plots.h_jetveto_mjj_veto->fill(vars.mjj, weight);
        }
        plots.h_jetveto_dy_inc->fill(vars.deltay, weight);
        plots.h_jetveto_mjj_inc->fill(vars.mjj, weight);

        if (vars.pass_ptbaleff) {
          plots.h_ptbaleff_mjj_veto->fill(vars.mjj, weight);
          plots.h_ptbaleff_dy_veto->fill(vars.deltay, weight);
        }
        plots.h_ptbaleff_mjj_inc->fill(vars.mjj, weight);
        plots.h_ptbaleff_dy_inc->fill(vars.deltay, weight);

        plots.p_avgnjets_dy->fill(vars.deltay, vars.ngapjets, weight);
        plots.p_avgnjets_mjj->fill(vars.mjj, vars.ngapjets, weight);
      }

      if (phase_space == "highmass") {
        plots.h_njets->fill(vars.ngapjets, weight);
        plots.h_dphijj->fill(vars.deltaphijj, weight);
        plots.h_ptbal->fill(vars.ptbalance2, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      finalizePlots(baseline_plots);
      finalizePlots(highpt_plots);
      finalizePlots(search_plots);
      finalizePlots(control_plots);
      finalizePlots(highmass_plots);
      finalizeEfficiencies(baseline_plots);
      finalizeEfficiencies(highpt_plots);
    }

    void finalizePlots(Plots& plots) {
      if (plots.h_dy) normalize(plots.h_dy);
      if (plots.h_mjj) normalize(plots.h_mjj);
      if (plots.h_dphijj) normalize(plots.h_dphijj);
      if (plots.h_njets) normalize(plots.h_njets);
      if (plots.h_ptbal) normalize(plots.h_ptbal);
    }

    void finalizeEfficiencies(Plots& plots) {
      int region_index = 0;
      if (plots.label=="baseline") region_index = 1;
      else if (plots.label=="highpt") region_index = 2;
      else return;

      if (plots.h_jetveto_mjj_veto && plots.h_jetveto_mjj_inc) divide(plots.h_jetveto_mjj_veto, plots.h_jetveto_mjj_inc, bookScatter2D(region_index, 1, 2));
      getScatter2D(region_index, 1, 2)->addAnnotation("InclusiveSumWeights", plots.h_jetveto_mjj_inc->integral());
      removeAnalysisObject(plots.h_jetveto_mjj_veto); removeAnalysisObject(plots.h_jetveto_mjj_inc);

      if (plots.h_jetveto_dy_veto && plots.h_jetveto_dy_inc) divide(plots.h_jetveto_dy_veto, plots.h_jetveto_dy_inc, bookScatter2D(region_index, 2, 2));
      getScatter2D(region_index, 2, 2)->addAnnotation("InclusiveSumWeights", plots.h_jetveto_dy_inc->integral());
      removeAnalysisObject(plots.h_jetveto_dy_veto); removeAnalysisObject(plots.h_jetveto_dy_inc);

      if (plots.h_ptbaleff_mjj_veto && plots.h_ptbaleff_mjj_inc) divide(plots.h_ptbaleff_mjj_veto, plots.h_ptbaleff_mjj_inc, bookScatter2D(region_index, 1, 3));
      getScatter2D(region_index, 1, 3)->addAnnotation("InclusiveSumWeights", plots.h_ptbaleff_mjj_inc->integral());
      removeAnalysisObject(plots.h_ptbaleff_mjj_veto); removeAnalysisObject(plots.h_ptbaleff_mjj_inc);

      if (plots.h_ptbaleff_dy_veto && plots.h_ptbaleff_dy_inc) divide(plots.h_ptbaleff_dy_veto, plots.h_ptbaleff_dy_inc, bookScatter2D(region_index, 2, 3));
      getScatter2D(region_index, 2, 3)->addAnnotation("InclusiveSumWeights", plots.h_ptbaleff_dy_inc->integral());
      removeAnalysisObject(plots.h_ptbaleff_dy_veto); removeAnalysisObject(plots.h_ptbaleff_dy_inc);
    }

    //@}


  private:

    //Variables* vars;

    Plots baseline_plots;
    Plots highpt_plots;
    Plots search_plots;
    Plots control_plots;
    Plots highmass_plots;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1279489);

}
