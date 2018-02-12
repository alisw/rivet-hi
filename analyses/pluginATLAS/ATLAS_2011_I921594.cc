// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Inclusive isolated prompt photon analysis with full 2010 LHC data
  class ATLAS_2011_I921594 : public Analysis {
  public:

    /// Constructor
    ATLAS_2011_I921594()
      : Analysis("ATLAS_2011_I921594")
    {    }


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      declare(fs, "FS");

      // Consider the final state jets for the energy density calculation
      FastJets fj(fs, FastJets::KT, 0.5);
      fj.useJetArea(new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec()));
      declare(fj, "KtJetsD05");

      // Consider the leading pt photon with |eta|<2.37 and pT>45 GeV
      LeadingParticlesFinalState photonfs(FinalState(-2.37, 2.37, 45*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // Book the dsigma/dEt (in eta bins) histograms
      for (size_t i = 0; i < _eta_bins.size()-1; i++) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue; // skip this bin
        _h_Et_photon[i] = bookHisto1D(1, 1, i+1);
      }
    }


    /// Return eta bin for either dsigma/dET histogram (area_eta=false) or energy density correction (area_eta=true)
    size_t _getEtaBin(double eta, bool area_eta) const {
      const double aeta = fabs(eta);
      return (!area_eta) ? binIndex(aeta, _eta_bins) : binIndex(aeta, _eta_bins_areaoffset);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve leading photon
      const Particles& photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1) vetoEvent;
      const Particle& leadingPhoton = photons[0];

      // Veto events with photon in ECAL crack
      if (inRange(leadingPhoton.abseta(), 1.37, 1.52)) vetoEvent;

      // Compute isolation energy in cone of radius .4 around photon (all particles)
      FourMomentum mom_in_EtCone;
      Particles fs = apply<FinalState>(event, "FS").particles();
      for (const Particle& p : fs) {
        // Check if it's outside the cone of 0.4
        if (deltaR(leadingPhoton, p) >= 0.4) continue;
        // Don't count particles in the 5x7 central core
        if (deltaEta(leadingPhoton, p) < .025*5.0*0.5 &&
            deltaPhi(leadingPhoton, p) < (PI/128.)*7.0*0.5) continue;
        // Increment isolation energy
        mom_in_EtCone += p.momentum();
      }

      // Get the area-filtered jet inputs for computing median energy density, etc.
      vector< vector<double> > ptDensities(_eta_bins_areaoffset.size()-1);
      FastJets fast_jets = apply<FastJets>(event, "KtJetsD05");
      const shared_ptr<fastjet::ClusterSequenceArea> clust_seq_area = fast_jets.clusterSeqArea();
      for (const Jet& jet : fast_jets.jets()) {
        const double area = clust_seq_area->area(jet); //< Implicit call to .pseudojet()
        if (area > 1e-4 && jet.abseta() < _eta_bins_areaoffset.back())
          ptDensities.at( _getEtaBin(jet.abseta(), true) ).push_back(jet.pT()/area);
      }

      // Compute the median energy density, etc.
      vector<double> ptDensity;
      for (size_t b = 0; b < _eta_bins_areaoffset.size()-1; b++) {
        ptDensity += ptDensities[b].empty() ? 0 : median(ptDensities[b]);
      }

      // Compute the isolation energy correction (cone area*energy density)
      const double ETCONE_AREA = M_PI*sqr(0.4) - (7.0*.025)*(5.0*PI/128.);
      const double correction = ptDensity[_getEtaBin(leadingPhoton.abseta(), true)] * ETCONE_AREA;

      // Apply isolation cut on area-corrected value
      if (mom_in_EtCone.Et() - correction > 4*GeV) vetoEvent;

      // Fill histograms
      const size_t eta_bin = _getEtaBin(leadingPhoton.abseta(), false);
      _h_Et_photon[eta_bin]->fill(leadingPhoton.Et(), event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (size_t i = 0; i < _eta_bins.size()-1; i++) {
        if (fuzzyEquals(_eta_bins[i], 1.37)) continue;
        scale(_h_Et_photon[i], crossSection()/picobarn/sumOfWeights());
      }
    }


  private:

    Histo1DPtr _h_Et_photon[5];

    const vector<double> _eta_bins = {0.00, 0.60, 1.37, 1.52, 1.81, 2.37};
    const vector<double> _eta_bins_areaoffset = {0.0, 1.5, 3.0};

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2011_I921594);

}
