// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement of isolated gamma + jet + X differential cross-sections
  class ATLAS_2013_I1244522 : public Analysis {
  public:

    // Constructor
    ATLAS_2013_I1244522()
      : Analysis("ATLAS_2013_I1244522")
    {     }


    // Book histograms and initialise projections before the run
    void init() {
      FinalState fs;

      // Voronoi eta-phi tassellation with KT jets, for ambient energy density calculation
      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      declare(fj, "KtJetsD05");

      // Leading photon
      LeadingParticlesFinalState photonfs(PromptFinalState(FinalState(-2.37, 2.37, 45.0*GeV)));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "JetFS");

      // Jets
      FastJets jetpro(vfs, FastJets::ANTIKT, 0.6);
      jetpro.useInvisibles();
      declare(jetpro, "Jets");

      // Histograms
      _h_ph_pt      = bookHisto1D(1, 1, 1);
      _h_jet_pt     = bookHisto1D(2, 1, 1);
      _h_jet_rap    = bookHisto1D(3, 1, 1);
      _h_dphi_phjet = bookHisto1D(4, 1, 1);
      _h_costheta_biased_phjet = bookHisto1D(5, 1, 1);
      _h_mass_phjet            = bookHisto1D(6, 1, 1);
      _h_costheta_phjet        = bookHisto1D(7, 1, 1);

    }


    size_t getEtaBin(double eta) const {
      const double aeta = fabs(eta);
      return binIndex(aeta, _eta_bins_areaoffset);
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the photon
      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1 )  vetoEvent;
      const Particle& photon = photons[0];

      if (inRange(photon.abseta(), 1.37, 1.52))  vetoEvent;

      //Compute isolation energy in cone of radius .4 around photon (all particles)
      FourMomentum mom_in_EtCone;
      const Particles& fs = apply<VetoedFinalState>(event, "JetFS").particles();
      for (const Particle& p : fs) {
        // Check if it's outside the cone of 0.4
        if (deltaR(photon, p) >= 0.4) continue;
        // Increment isolation energy
        mom_in_EtCone += p.momentum();
      }

      // Get the jets
      Jets alljets = apply<FastJets>(event, "Jets").jetsByPt(40.0*GeV);
      Jets jets;
      for (const Jet& jet : alljets)
        if (deltaR(photon, jet) > 1.0) jets += jet;
      if (jets.empty())  vetoEvent;
      Jet leadingJet = jets[0];
      if (leadingJet.absrap() > 2.37) vetoEvent;

      // Get the area-filtered jet inputs for computing median energy density, etc.
      vector<double> ptDensity;
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      FastJets fast_jets = apply<FastJets>(event, "KtJetsD05");
      const auto clust_seq_area = fast_jets.clusterSeqArea();
      foreach (const Jet& jet, fast_jets.jets()) {
        const double area = clust_seq_area->area(jet);
        if (area > 1e-4 && jet.abseta() < _eta_bins_areaoffset.back())
          ptDensities.at( getEtaBin(jet.abseta()) ).push_back(jet.pT()/area);
      }

      // Compute the median energy density, etc.
      for (size_t b = 0; b < _eta_bins_areaoffset.size() - 1; ++b) {
        const int njets = ptDensities[b].size();
        ptDensity += (njets > 0) ? median(ptDensities[b]) : 0;
      }

      // Compute the isolation energy correction (cone area*energy density)
      const double etCone_area = PI*sqr(0.4) - (5.0*.025)*(7.0*PI/128.);
      const double correction = ptDensity[getEtaBin(photon.abseta())] * etCone_area;

      // Apply isolation cut on area-corrected value
      if (mom_in_EtCone.Et() - correction >= 4*GeV)  vetoEvent;

      // Fill histos
      const double weight = event.weight();
      const double dy = deltaRap(photon, leadingJet);
      const double costheta_yj = tanh(dy/2);
      _h_ph_pt->fill(photon.pT()/GeV, weight);
      _h_jet_pt->fill(leadingJet.pT()/GeV, weight);
      _h_jet_rap->fill(leadingJet.absrap(), weight);
      _h_dphi_phjet->fill(deltaPhi(photon, leadingJet), weight);
      _h_costheta_biased_phjet->fill(costheta_yj, weight);
      if (costheta_yj < 0.829022) {
        const FourMomentum yj = photon.momentum() + leadingJet.momentum();
        if (yj.mass() > 160.939*GeV) {
          if (fabs(photon.eta() + leadingJet.rap()) < 2.37) {
            _h_mass_phjet->fill(yj.mass()/GeV, weight);
            _h_costheta_phjet->fill(costheta_yj, weight);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / picobarn / sumOfWeights();
      scale(_h_ph_pt,                 sf);
      scale(_h_jet_pt,                sf);
      scale(_h_jet_rap,               sf);
      scale(_h_dphi_phjet,            sf);
      scale(_h_costheta_biased_phjet, sf);
      scale(_h_mass_phjet,            sf);
      scale(_h_costheta_phjet,        sf);
    }


  private:

    Histo1DPtr _h_ph_pt, _h_jet_pt, _h_jet_rap, _h_dphi_phjet, _h_costheta_biased_phjet, _h_mass_phjet, _h_costheta_phjet;

    const vector<double> _eta_bins_areaoffset = {0.0, 1.5, 3.0};

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1244522);


}
