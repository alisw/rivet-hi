// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Isolated diphoton + X differential cross-sections
  class ATLAS_2017_I1591327 : public Analysis {
  public:

    // Constructor
    ATLAS_2017_I1591327() : Analysis("ATLAS_2017_I1591327") {
      setNeedsCrossSection(true);
    }


    // Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs, "FS");

      FastJets fj(fs, FastJets::KT, 0.5);
      _area_def = new fastjet::AreaDefinition(fastjet::VoronoiAreaSpec());
      fj.useJetArea(_area_def);
      declare(fj, "KtJetsD05");

      IdentifiedFinalState photonfs(Cuts::abseta < 2.37 && Cuts::pT > 30*GeV);
      photonfs.acceptId(PID::PHOTON);
      declare(photonfs, "Photon");

      // Histograms
      _h_M       = bookHisto1D(2, 1, 1);
      _h_pT      = bookHisto1D(3, 1, 1);
      _h_at      = bookHisto1D(4, 1, 1);
      _h_phistar = bookHisto1D(5, 1, 1);
      _h_costh   = bookHisto1D(6, 1, 1);
      _h_dPhi    = bookHisto1D(7, 1, 1);
    }


    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Require at least 2 photons in final state
      const Particles photons = apply<IdentifiedFinalState>(event, "Photon").particlesByPt();
      if (photons.size() < 2) vetoEvent;

      // Compute the median energy density
      _ptDensity.clear();
      _sigma.clear();
      _Njets.clear();
      vector<vector<double> > ptDensities;
      vector<double> emptyVec;
      ptDensities.assign(ETA_BINS.size()-1, emptyVec);

      // Get jets, and corresponding jet areas
      const shared_ptr<fastjet::ClusterSequenceArea> clust_seq_area = applyProjection<FastJets>(event, "KtJetsD05").clusterSeqArea();
      for (const fastjet::PseudoJet& jet : apply<FastJets>(event, "KtJetsD05").pseudoJets(0.0*GeV)) {
        const double aeta = fabs(jet.eta());
        const double pt = jet.perp();
        const double area = clust_seq_area->area(jet);
        if (area < 1e-3) continue;
        const int ieta = binIndex(aeta, ETA_BINS);
        if (ieta != -1) ptDensities[ieta].push_back(pt/area);
      }

      // Compute median jet properties over the jets in the event
      for (size_t b = 0; b < ETA_BINS.size()-1; ++b) {
        double median = 0.0, sigma = 0.0;
        int Njets = 0;
        if (ptDensities[b].size() > 0) {
          std::sort(ptDensities[b].begin(), ptDensities[b].end());
          int nDens = ptDensities[b].size();
          median = (nDens % 2 == 0) ? (ptDensities[b][nDens/2]+ptDensities[b][(nDens-2)/2])/2 : ptDensities[b][(nDens-1)/2];
          sigma = ptDensities[b][(int)(.15865*nDens)];
          Njets = nDens;
        }
        _ptDensity.push_back(median);
        _sigma.push_back(sigma);
        _Njets.push_back(Njets);
      }

      // Loop over photons and fill vector of isolated ones
      Particles isolated_photons;
      for (const Particle& photon : photons) {
        // Check if it's a prompt photon (needed for SHERPA 2->5 sample, otherwise I also get photons from hadron decays in jets)
        if (photon.fromDecay()) continue;

        // Remove photons in ECAL crack region
        if (inRange(photon.abseta(), 1.37, 1.56))  continue;
        const double eta_P = photon.eta();
        const double phi_P = photon.phi();

        // Compute isolation via particles within an R=0.4 cone of the photon
        const Particles fs = apply<FinalState>(event, "FS").particles();
        FourMomentum mom_in_EtCone;
        for (const Particle& p : fs) {
          // Reject if not in cone
          if (deltaR(photon.momentum(), p.momentum()) > 0.4)  continue;
          // Reject if in the 5x7 cell central core
          if (fabs(eta_P - p.eta()) < 0.025 * 5 * 0.5 &&
              fabs(phi_P - p.phi()) < PI/128. * 7 * 0.5)  continue;
          // Sum momentum
          mom_in_EtCone += p.momentum();
        }
        // Now figure out the correction (area*density)
        const double EtCone_area = M_PI*sqr(0.4) - (7*.025)*(5*M_PI/128.); // cone area - central core rectangle
        const double correction = _ptDensity[binIndex(fabs(eta_P), ETA_BINS)] * EtCone_area;

        // Discard the photon if there is more than 11 GeV of cone activity
        // NOTE: Shouldn't need to subtract photon itself (it's in the central core)
        if (mom_in_EtCone.Et() - correction > 11*GeV)  continue;
        // Add isolated photon to list
        isolated_photons.push_back(photon);
      }

      // Require at least two isolated photons
      if (isolated_photons.size() < 2) vetoEvent;

      // Select leading pT pair
      sortByPt(isolated_photons);
      const FourMomentum y1 = isolated_photons[0];
      const FourMomentum y2 = isolated_photons[1];

      // Leading photon should have pT > 40 GeV, subleading > 30 GeV
      if (y1.pT() < 40.*GeV) vetoEvent;
      if (y2.pT() < 30.*GeV) vetoEvent;

      // Require the two photons to be separated (dR>0.4)
      if (deltaR(y1,y2) < 0.4) vetoEvent;

      const FourMomentum yy = y1 + y2;
      const double Myy = yy.mass();
      const double pTyy = yy.pT();
      const double dPhiyy = mapAngle0ToPi(y1.phi() - y2.phi());

      // phi*
      const double costhetastar_ = fabs(tanh(( y1.eta() - y2.eta() ) / 2.));
      const double sinthetastar_ = sqrt(1. - pow(costhetastar_, 2));
      const double phistar = tan(0.5 * (PI - dPhiyy)) * sinthetastar_;

      // a_t
      const Vector3 t_hat(y1.x()-y2.x(), y1.y()-y2.y(), 0.);
      const double factor = t_hat.mod();
      const Vector3 t_hatx(t_hat.x()/factor, t_hat.y()/factor, t_hat.z()/factor);
      const Vector3 At(y1.x()+y2.x(), y1.y()+y2.y(), 0.);
      // Compute a_t transverse component with respect to t_hat
      const double at = At.cross(t_hatx).mod();

      // Fill histograms
      const double weight = event.weight();
      _h_M->fill(Myy, weight);
      _h_pT->fill(pTyy, weight);
      _h_dPhi->fill(dPhiyy, weight);
      _h_costh->fill(costhetastar_, weight);
      _h_phistar->fill(phistar, weight);
      _h_at->fill(at, weight);
    }


    // Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() / (femtobarn * sumOfWeights());
      scale({_h_M, _h_pT, _h_dPhi, _h_costh, _h_phistar, _h_at}, sf);
    }


  private:

    Histo1DPtr _h_M, _h_pT, _h_dPhi, _h_costh, _h_phistar, _h_at;

    fastjet::AreaDefinition* _area_def;

    const vector<double> ETA_BINS = {0.0, 1.5, 3.0};
    vector<double> _ptDensity, _sigma, _Njets;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1591327);

}
