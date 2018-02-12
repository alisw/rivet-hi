// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

namespace Rivet {


  /// @brief Measurement of prompt isolated photon + b/c-jet + X differential cross-sections
  class ATLAS_2017_I1632756 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1632756);


    /// Book histograms and initialise projections before the run
    void init() {
      // particles for photon isolation: no muons, no neutrinos
      declare(VisibleFinalState(Cuts::abspid != PID::MUON), "caloParticles");

      // Voronoi eta-phi tessellation with KT jets, for ambient energy density calculation
      FastJets fj(FinalState(), FastJets::KT, 0.5, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      declare(fj, "KtJetsD05");

      // Leading photon
      LeadingParticlesFinalState photonfs(PromptFinalState(Cuts::abseta < 2.37 && Cuts::pT > 25*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // Jets
      FastJets jetpro(FinalState(), FastJets::ANTIKT, 0.4, JetAlg::DECAY_MUONS, JetAlg::DECAY_INVISIBLES);
      declare(jetpro, "Jets");

      // Heavy hadrons
      declare(HeavyHadrons(), "HeavyHadrons");

      // Book the dsigma/dEt (in eta bins) histograms
      // d02 and d03 are for photon+b; d04 and d05 are for photon+c
      for (size_t i = 0; i < _eta_bins.size() - 1; ++i) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue; // This region is not used
        int offset = i > 1? 1 : 2;
        _h_Et_photonb[i] = bookHisto1D(i + offset, 1, 1);
        _h_Et_photonc[i] = bookHisto1D(i + offset + 2, 1, 1);
      }

    }


    /// Return eta bin for either dsigma/dET histogram (area_eta=false) or energy density correction (area_eta=true)
    size_t _getEtaBin(double eta_w, bool area_eta) const {
      const double eta = fabs(eta_w);
      if (!area_eta) {
        return binIndex(eta, _eta_bins);
      } else {
        return binIndex(eta, _eta_bins_areaoffset);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the leading photon
      const Particles& photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particlesByPt();
      if (photons.empty())  vetoEvent;
      const Particle& leadingPhoton = photons[0];

      // Veto events with leading photon in ECAL crack
      if (inRange(leadingPhoton.abseta(), 1.37, 1.56))  vetoEvent;

      // Compute isolation energy in cone of radius .4 around photon (all particles except muons, neutrinos and leading photon)
      FourMomentum mom_in_EtCone;
      const Particles& fs = apply<FinalState>(event, "caloParticles").particles();
      for (const Particle& p : fs) {
        // increment if inside cone of 0.4
        if (deltaR(leadingPhoton, p) < 0.4)  mom_in_EtCone += p.momentum();
      }
      // Remove the photon energy from the isolation
      mom_in_EtCone -= leadingPhoton.momentum();

      // Get the area-filtered jet inputs for computing median energy density, etc.
      vector<double> ptDensity;
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      const FastJets& fast_jets = apply<FastJets>(event, "KtJetsD05");
      const auto clust_seq_area = fast_jets.clusterSeqArea();
      for (const Jet& jet : fast_jets.jets()) {
        const double area = clust_seq_area->area(jet);
        if (area > 10e-4 && jet.abseta() < _eta_bins_areaoffset.back())
          ptDensities.at( _getEtaBin(jet.abseta(), true) ).push_back(jet.pT()/area);
      }

      // Compute the median energy density, etc.
      for (size_t b = 0; b < _eta_bins_areaoffset.size() - 1; ++b) {
        const double ptmedian = (ptDensities[b].size() > 0) ? median(ptDensities[b]) : 0;
        ptDensity.push_back(ptmedian);
      }

      // Compute the isolation energy correction (cone area*energy density)
      const double etCone_area = PI * sqr(0.4);
      const double correction = ptDensity[_getEtaBin(leadingPhoton.abseta(), true)] * etCone_area;

      // Apply isolation cut on area-corrected value
      // cut is Etiso < 4.8GeV + 4.2E-03 * Et_gamma.
      if (mom_in_EtCone.Et() - correction > 4.8*GeV + 0.0042*leadingPhoton.Et())  vetoEvent;


      // Get the leading jet
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV);
      ifilter_discard(jets, deltaRLess(leadingPhoton, 0.4));
      if (jets.empty())  vetoEvent;
      const Jet& leadingJet = jets[0];

      // Veto events with leading jet outside |y|<2.5
      if (leadingJet.absrap() > 2.5)  vetoEvent;

      // Veto events with leading photon and leading jet with deltaR < 1.0
      if (deltaR(leadingPhoton, leadingJet) < 1.0)  vetoEvent;

      // Veto events with leading jet not b-tagged (deltaR match with a B-hadron) nor c-tagged (deltaR match with a C-hadron)
      const Particles& allBs = apply<HeavyHadrons>(event, "HeavyHadrons").bHadrons(5*GeV);
      bool bTagged = false;
      for (const Particle& thisB : allBs) {
        if(deltaR(thisB, leadingJet) < 0.3) {
          bTagged = true;
          break;
        }
      }

      bool cTagged = false;
      if (!bTagged) {
        const Particles& allCs = apply<HeavyHadrons>(event, "HeavyHadrons").cHadrons(5*GeV);
        for (const Particle& thisC : allCs) {
          if (deltaR(thisC, leadingJet) < 0.3) {
            cTagged = true;
            break;
          }
        }
        if (!cTagged) vetoEvent;
      }

      // Fill histograms
      const size_t eta_bin = _getEtaBin(leadingPhoton.abseta(), false);
      if (bTagged) _h_Et_photonb[eta_bin]->fill(leadingPhoton.Et()/GeV, event.weight());
      if (cTagged) _h_Et_photonc[eta_bin]->fill(leadingPhoton.Et()/GeV, event.weight());

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / (picobarn * sumOfWeights());
      for (size_t i = 0; i < _eta_bins.size() - 1; ++i) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue; // This region is not used
        scale(_h_Et_photonb[i], sf);
        scale(_h_Et_photonc[i], sf);
      }
    }


  private:

      Histo1DPtr _h_Et_photonb[3];
      Histo1DPtr _h_Et_photonc[3];

      const vector<double> _eta_bins = { 0.00, 1.37, 1.56, 2.37 };
      const vector<double> _eta_bins_areaoffset = { 0.0, 1.5, 3.0 };

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1632756);


}
