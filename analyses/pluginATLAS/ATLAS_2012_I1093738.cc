// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement of isolated gamma + jet + X differential cross-sections
  ///
  /// Inclusive isolated gamma + jet cross-sections, differential in pT(gamma), for
  /// various photon and jet rapidity configurations.
  ///
  /// @author Giovanni Marchiori
  class ATLAS_2012_I1093738 : public Analysis {
  public:

    // Constructor
    ATLAS_2012_I1093738()
      : Analysis("ATLAS_2012_I1093738")
    {    }


    // Book histograms and initialise projections before the run
    void init() {
      // Final state
      FinalState fs;
      declare(fs, "FS");

      // Voronoi eta-phi tessellation with KT jets, for ambient energy density calculation
      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      declare(fj, "KtJetsD05");

      // Leading photon
      LeadingParticlesFinalState photonfs(FinalState(-1.37, 1.37, 25.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "JetFS");

      // Jets
      FastJets jetpro(vfs, FastJets::ANTIKT, 0.4);
      jetpro.useInvisibles();
      declare(jetpro, "Jets");

      _h_phbarrel_jetcentral_SS = bookHisto1D(1, 1, 1);
      _h_phbarrel_jetmedium_SS  = bookHisto1D(2, 1, 1);
      _h_phbarrel_jetforward_SS = bookHisto1D(3, 1, 1);

      _h_phbarrel_jetcentral_OS = bookHisto1D(4, 1, 1);
      _h_phbarrel_jetmedium_OS  = bookHisto1D(5, 1, 1);
      _h_phbarrel_jetforward_OS = bookHisto1D(6, 1, 1);
    }


    int getEtaBin(double eta, int what) const {
      const double aeta = fabs(eta);
      if (what == 0) return binIndex(aeta, _eta_bins_ph);
      if (what == 1) return binIndex(aeta, _eta_bins_jet);
      return binIndex(aeta, _eta_bins_areaoffset);
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the photon
      const FinalState& photonfs = apply<FinalState>(event, "LeadingPhoton");
      if (photonfs.particles().size() < 1) vetoEvent;
      const FourMomentum photon = photonfs.particles().front().momentum();

      // Get the jet
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(20.0*GeV);
      if (jets.empty()) vetoEvent;
      FourMomentum leadingJet = jets[0].momentum();

      // Require jet separated from photon
      if (deltaR(photon, leadingJet) < 1.0) vetoEvent;

      // Veto if leading jet is outside plotted rapidity regions
      if (leadingJet.absrap() > 4.4) vetoEvent;

      // Compute the jet pT densities
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      FastJets fastjets = apply<FastJets>(event, "KtJetsD05");
      const shared_ptr<fastjet::ClusterSequenceArea> clust_seq_area = fastjets.clusterSeqArea();
      for (const Jet& jet : fastjets.jets()) {
        const double area = clust_seq_area->area(jet); //< Implicit call to pseudojet()
        if (area > 1e-4 && jet.abseta() < _eta_bins_areaoffset.back()) {
          ptDensities.at(getEtaBin(jet.abseta(), 2)) += jet.pT()/area;
        }
      }

      // Compute the median event energy density
      /// @todo This looks equivalent to median(ptDensities[b]) -- isn't SKIPNHARDJETS meant to be used as an offset?
      const unsigned int SKIPNHARDJETS = 0;
      vector<double> ptDensity;
      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; b++) {
        double median = 0.0;
        if (ptDensities[b].size() > SKIPNHARDJETS) {
          std::sort(ptDensities[b].begin(), ptDensities[b].end());
          const int nDens = ptDensities[b].size() - SKIPNHARDJETS;
          if (nDens % 2 == 0) {
            median = (ptDensities[b][nDens/2]+ptDensities[b][(nDens-2)/2])/2;
          } else {
            median = ptDensities[b][(nDens-1)/2];
          }
        }
        ptDensity.push_back(median);
      }

      // Compute photon isolation with a standard ET cone
      const Particles fs = apply<FinalState>(event, "JetFS").particles();
      FourMomentum mom_in_EtCone;
      const double ISO_DR = 0.4;
      const double CLUSTER_ETA_WIDTH = 0.25*5.0;
      const double CLUSTER_PHI_WIDTH = (PI/128.)*7.0;
      for (const Particle& p : fs) {
        // Check if it's in the cone of .4
        if (deltaR(photon, p) >= ISO_DR) continue;
        // Check if it's in the 5x7 central core
        if (fabs(deltaEta(photon, p)) < CLUSTER_ETA_WIDTH*0.5 &&
            fabs(deltaPhi(photon, p)) < CLUSTER_PHI_WIDTH*0.5) continue;
        // Increment sum
        mom_in_EtCone += p.momentum();
      }

      // Figure out the correction (area*density)
      const double ETCONE_AREA = PI*ISO_DR*ISO_DR - CLUSTER_ETA_WIDTH*CLUSTER_PHI_WIDTH;
      const double correction = ptDensity[getEtaBin(photon.abseta(),2)] * ETCONE_AREA;

      // Require photon to be isolated
      if (mom_in_EtCone.Et()-correction > 4.0*GeV) vetoEvent;

      const int photon_jet_sign = sign( leadingJet.rapidity() * photon.rapidity() );

      // Fill histos
      const double abs_jet_rapidity = fabs(leadingJet.rapidity());
      const double photon_pt = photon.pT()/GeV;
      const double abs_photon_eta = fabs(photon.eta());
      const double weight = event.weight();
      if (abs_photon_eta < 1.37) {
        if (abs_jet_rapidity < 1.2) {
          if (photon_jet_sign >= 1) {
            _h_phbarrel_jetcentral_SS->fill(photon_pt, weight);
          } else {
            _h_phbarrel_jetcentral_OS->fill(photon_pt, weight);
          }
        } else if (abs_jet_rapidity < 2.8) {
          if (photon_jet_sign >= 1) {
            _h_phbarrel_jetmedium_SS->fill(photon_pt, weight);
          } else {
            _h_phbarrel_jetmedium_OS->fill(photon_pt, weight);
          }
        } else if (abs_jet_rapidity < 4.4) {
          if (photon_jet_sign >= 1) {
            _h_phbarrel_jetforward_SS->fill(photon_pt, weight);
          } else {
            _h_phbarrel_jetforward_OS->fill(photon_pt, weight);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_phbarrel_jetcentral_SS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetcentral_OS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetmedium_SS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetmedium_OS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetforward_SS, crossSection()/sumOfWeights());
      scale(_h_phbarrel_jetforward_OS, crossSection()/sumOfWeights());
    }


  private:

    Histo1DPtr _h_phbarrel_jetcentral_SS, _h_phbarrel_jetmedium_SS, _h_phbarrel_jetforward_SS;
    Histo1DPtr _h_phbarrel_jetcentral_OS, _h_phbarrel_jetmedium_OS, _h_phbarrel_jetforward_OS;

    const vector<double> _eta_bins_ph = {0.0, 1.37, 1.52, 2.37};
    const vector<double> _eta_bins_jet = {0.0, 1.2, 2.8, 4.4};
    const vector<double> _eta_bins_areaoffset = {0.0, 1.5, 3.0};

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1093738);


}
