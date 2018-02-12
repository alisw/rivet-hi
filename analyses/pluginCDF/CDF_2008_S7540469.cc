// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement differential Z/\f$ \gamma^* \f$ + jet + \f$ X \f$ cross sections
  /// @author Frank Siegert
  class CDF_2008_S7540469 : public Analysis {

  public:

    /// Constructor
    CDF_2008_S7540469()
      : Analysis("CDF_2008_S7540469")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      // Full final state
      FinalState fs(-5.0, 5.0);
      declare(fs, "FS");

      // Leading electrons in tracking acceptance
      IdentifiedFinalState elfs(Cuts::abseta < 5 && Cuts::pT > 25*GeV);
      elfs.acceptIdPair(PID::ELECTRON);
      declare(elfs, "LeadingElectrons");

      _h_jet_multiplicity = bookHisto1D(1, 1, 1);
      _h_jet_pT_cross_section_incl_1jet = bookHisto1D(2, 1, 1);
      _h_jet_pT_cross_section_incl_2jet = bookHisto1D(3, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event & event) {
      const double weight = event.weight();

      // Skip if the event is empty
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.empty()) {
        MSG_DEBUG("Skipping event " << numEvents() << " because no final state pair found");
        vetoEvent;
      }

      // Find the Z candidates
      const FinalState & electronfs = apply<FinalState>(event, "LeadingElectrons");
      std::vector<std::pair<Particle, Particle> > Z_candidates;
      Particles all_els=electronfs.particles();
      for (size_t i=0; i<all_els.size(); ++i) {
        for (size_t j=i+1; j<all_els.size(); ++j) {
          bool candidate=true;
          double mZ = FourMomentum(all_els[i].momentum()+all_els[j].momentum()).mass()/GeV;
          if (mZ < 66.0 || mZ > 116.0) {
            candidate = false;
          }
          double abs_eta_0 = fabs(all_els[i].eta());
          double abs_eta_1 = fabs(all_els[j].eta());
          if (abs_eta_1 < abs_eta_0) {
            double tmp = abs_eta_0;
            abs_eta_0 = abs_eta_1;
            abs_eta_1 = tmp;
          }
          if (abs_eta_0 > 1.0) {
            candidate = false;
          }
          if (!(abs_eta_1 < 1.0 || (inRange(abs_eta_1, 1.2, 2.8)))) {
            candidate = false;
          }
          if (candidate) {
            Z_candidates.push_back(make_pair(all_els[i], all_els[j]));
          }
        }
      }
      if (Z_candidates.size() != 1) {
        MSG_DEBUG("Skipping event " << numEvents() << " because no unique electron pair found ");
        vetoEvent;
      }

      // Now build the jets on a FS without the electrons from the Z (including QED radiation)
      Particles jetparts;
      for (const Particle& p : fs.particles()) {
        bool copy = true;
        if (p.pid() == PID::PHOTON) {
          FourMomentum p_e0 = Z_candidates[0].first.momentum();
          FourMomentum p_e1 = Z_candidates[0].second.momentum();
          FourMomentum p_P = p.momentum();
          if (deltaR(p_e0, p_P) < 0.2) copy = false;
          if (deltaR(p_e1, p_P) < 0.2) copy = false;
        } else {
          if (p.genParticle()->barcode() == Z_candidates[0].first.genParticle()->barcode()) copy = false;
          if (p.genParticle()->barcode() == Z_candidates[0].second.genParticle()->barcode()) copy = false;
        }
        if (copy) jetparts.push_back(p);
      }

      // Proceed to lepton dressing
      const PseudoJets pjs = mkPseudoJets(jetparts);
      const auto jplugin = make_shared<fastjet::CDFMidPointPlugin>(0.7, 0.5, 1.0);
      const Jets jets_all = mkJets(fastjet::ClusterSequence(pjs, jplugin.get()).inclusive_jets());
      const Jets jets_cut = sortByPt(filterBy(jets_all, Cuts::pT > 30*GeV && Cuts::abseta < 2.1));
      // FastJets jetpro(FastJets::CDFMIDPOINT, 0.7);
      // jetpro.calc(jetparts);
      // // Take jets with pt > 30, |eta| < 2.1:
      // const Jets& jets = jetpro.jets();
      // Jets jets_cut;
      // foreach (const Jet& j, jets) {
      //   if (j.pT()/GeV > 30.0 && j.abseta() < 2.1) {
      //     jets_cut.push_back(j);
      //   }
      // }
      // // Sort by pT:
      // sort(jets_cut.begin(), jets_cut.end(), cmpMomByPt);

      // Return if there are no jets:
      MSG_DEBUG("Num jets above 30 GeV = " << jets_cut.size());
      if (jets_cut.empty()) {
        MSG_DEBUG("No jets pass cuts ");
        vetoEvent;
      }

      // Cut on Delta R between Z electrons and *all* jets
      for (const Jet& j : jets_cut) {
        if (deltaR(Z_candidates[0].first, j) < 0.7) vetoEvent;
        if (deltaR(Z_candidates[0].second, j) < 0.7) vetoEvent;
      }

      // Fill histograms
      for (size_t njet=1; njet<=jets_cut.size(); ++njet) {
        _h_jet_multiplicity->fill(njet, weight);
      }
      for (const Jet& j : jets_cut) {
        if (jets_cut.size() > 0) {
          _h_jet_pT_cross_section_incl_1jet->fill(j.pT(), weight);
        }
        if (jets_cut.size() > 1) {
          _h_jet_pT_cross_section_incl_2jet->fill(j.pT(), weight);
        }
      }
    }


    /// Rescale histos
    void finalize() {
      const double invlumi = crossSection()/femtobarn/sumOfWeights();
      scale(_h_jet_multiplicity, invlumi);
      scale(_h_jet_pT_cross_section_incl_1jet, invlumi);
      scale(_h_jet_pT_cross_section_incl_2jet, invlumi);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_jet_multiplicity;
    Histo1DPtr _h_jet_pT_cross_section_incl_1jet;
    Histo1DPtr _h_jet_pT_cross_section_incl_2jet;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_2008_S7540469);

}
