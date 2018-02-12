// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2010_S8914702 : public Analysis {
  public:

    /// Constructor
    ATLAS_2010_S8914702()
      : Analysis("ATLAS_2010_S8914702")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      declare(fs, "FS");

      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      declare(fj, "KtJetsD05");

      LeadingParticlesFinalState photonfs(FinalState(Cuts::abseta < 1.81 && Cuts::pT > 15*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      size_t hist_bin = 0;
      for (size_t i = 0; i < _eta_bins.size()-1; ++i) {
        if (fabs(_eta_bins[i] - 1.37) < .0001) continue;
        _h_Et_photon[i] = bookHisto1D(1, 1, hist_bin+1);
        hist_bin += 1;
      }
    }


    size_t getEtaBin(double eta, bool area_eta) const {
      return (!area_eta) ? binIndex(fabs(eta), _eta_bins) : binIndex(fabs(eta), _eta_bins_areaoffset);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1) vetoEvent;

      const Particle& leadingPhoton = photons[0];
      if (inRange(leadingPhoton.abseta(), 1.37, 1.52)) vetoEvent;
      const int eta_bin = getEtaBin(leadingPhoton.abseta(), false);

      const Particles& fs = apply<FinalState>(event, "FS").particles();
      FourMomentum mom_in_EtCone;
      for (const Particle& p : fs) {
        // Check if it's in the cone of .4
        if (deltaR(leadingPhoton, p) >= 0.4) continue;
        // Check if it's in the 5x7 central core
        if (fabs(deltaEta(leadingPhoton, p)) < .025*5.0*0.5 &&
            fabs(deltaPhi(leadingPhoton, p)) < (PI/128.)*7.0*0.5) continue;
        // Increment
        mom_in_EtCone += p.momentum();
      }
      MSG_DEBUG("Done with initial Et cone");

      // Get the jet pT densities
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      FastJets fastjets = apply<FastJets>(event, "KtJetsD05");
      const shared_ptr<fastjet::ClusterSequenceArea> clust_seq_area = fastjets.clusterSeqArea();
      for (const Jet& jet : fastjets.jets()) {
        const double area = clust_seq_area->area(jet); //< Implicit call to pseudojet()
        if (area > 1e-4 && jet.abseta() < _eta_bins_areaoffset.back()) {
          ptDensities.at(getEtaBin(jet.abseta(), true)) += jet.pT()/area;
        }
      }

      // Now compute the median energy densities
      vector<double> ptDensity;
      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; ++b) {
        ptDensity += ptDensities[b].empty() ? 0 : median(ptDensities[b]);
      }

      // Now figure out the correction
      const double ETCONE_AREA = PI*.4*.4 - (7.0*.025)*(5.0*PI/128.);
      const double correction = ptDensity[getEtaBin(leadingPhoton.abseta(), true)]*ETCONE_AREA;
      MSG_DEBUG("Jet area correction done");

      // Shouldn't need to subtract photon
      // NB. Using expected cut at hadron/particle level, not cut at reco level
      if (mom_in_EtCone.Et() - correction/*-leadingPhoton.Et()*/ > 4.0*GeV) vetoEvent;
      MSG_DEBUG("Passed isolation cut");

      // Fill histogram
      _h_Et_photon[eta_bin]->fill(leadingPhoton.Et()/GeV, event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i = 0; i < _eta_bins.size()-1; ++i) {
        if (fabs(_eta_bins[i] - 1.37) < .0001) continue;
        scale(_h_Et_photon[i], crossSection()/sumOfWeights());
      }
    }


  private:

    Histo1DPtr _h_Et_photon[6];

    const vector<double> _eta_bins = {0.00, 0.60, 1.37, 1.52, 1.81};
    const vector<double> _eta_bins_areaoffset = {0.0, 1.5, 3.0};

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2010_S8914702);

}
