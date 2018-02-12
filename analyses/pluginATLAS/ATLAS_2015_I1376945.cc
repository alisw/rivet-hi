// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Colour flow in hadronic top decay at 8 TeV
  class ATLAS_2015_I1376945 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1376945);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;

      PromptFinalState promptFs(fs);
      promptFs.acceptTauDecays(true);
      promptFs.acceptMuonDecays(false);

      IdentifiedFinalState neutrino_fs(promptFs);
      neutrino_fs.acceptNeutrinos();
      declare(neutrino_fs, "NEUTRINO_FS");

      IdentifiedFinalState Photon(fs);
      Photon.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState bare_muons_fs(promptFs);
      bare_muons_fs.acceptIdPair(PID::MUON);

      IdentifiedFinalState bare_elecs_fs(promptFs);
      bare_elecs_fs.acceptIdPair(PID::ELECTRON);

      Cut lep_cuts = (Cuts::abseta < 2.5) & (Cuts::pT > 1*MeV);
      DressedLeptons muons(Photon, bare_muons_fs, 0.1, lep_cuts);
      declare(muons, "MUONS");

      DressedLeptons elecs(Photon, bare_elecs_fs, 0.1, lep_cuts);
      declare(elecs, "ELECS");

      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(muons);
      vfs.addVetoOnThisFinalState(elecs);
      vfs.addVetoOnThisFinalState(neutrino_fs);

      FastJets fjets(vfs, FastJets::ANTIKT, 0.4);
      fjets.useInvisibles();
      declare(fjets, "jets");

      h_pull_all     = bookHisto1D(4,1,1);
      h_pull_charged = bookHisto1D(5,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      /**************
       *    JETS    *
       **************/
      const Jets& allJets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25.0*GeV && Cuts::absrap < 2.5);
      const vector<DressedLepton>& all_elecs = apply<DressedLeptons>(event, "ELECS").dressedLeptons();
      const vector<DressedLepton>& all_muons = apply<DressedLeptons>(event, "MUONS").dressedLeptons();
      Jets goodJets;
      foreach (const Jet j, allJets) {
        bool keep = true;
        foreach (const DressedLepton el, all_elecs)  keep &= deltaR(j, el) >= 0.2;
        if (keep)  goodJets += j;
      }
      if ( goodJets.size() < 4 )  vetoEvent;

      /****************
       *    LEPTONS   *
       ****************/
      vector<DressedLepton> muons, vetoMuons;
      foreach (const DressedLepton mu, all_muons) {
        bool keep = true;
        foreach (const Jet j, goodJets)  keep &= deltaR(j, mu) >= 0.4;
        if (keep && mu.pt() > 15*GeV) {
          vetoMuons.push_back(mu);
          if (mu.pt() > 25*GeV)  muons.push_back(mu);
        }
      }

      vector<DressedLepton> elecs, vetoElecs;
      foreach (const DressedLepton el, all_elecs) {
        bool keep = true;
        foreach (const Jet j, goodJets)  keep &= deltaR(j, el) >= 0.4;
        if (keep && el.pt() > 15*GeV) {
          vetoElecs.push_back(el);
          if (el.pt() > 25*GeV)  elecs.push_back(el);
        }
      }

      if (muons.empty() && elecs.empty())  vetoEvent;

      bool muCandidate = !( muons.size() < 1 || vetoMuons.size() > 1 || vetoElecs.size() > 0 );
      bool elCandidate = !( elecs.size() < 1 || vetoElecs.size() > 1 || vetoMuons.size() > 0 );

      if (!elCandidate && !muCandidate)  vetoEvent;

      /******************************
       *    ELECTRON-MUON OVERLAP   *
       ******************************/
      foreach (const DressedLepton electron, elecs) {
        foreach (const DressedLepton muon, muons) {
          double d_theta = fabs(muon.theta() - electron.theta());
          double d_phi = deltaPhi(muon.phi(), electron.phi());
          if (d_theta < 0.005 && d_phi < 0.005)  vetoEvent;
        }
      }

      /****************
       *  NEUTRINOS   *
       ****************/
      const Particles& neutrinos = apply<IdentifiedFinalState>(event, "NEUTRINO_FS").particlesByPt();
      FourMomentum metVector = FourMomentum(0.,0.,0.,0.);
      foreach (const Particle& n, neutrinos) {
        metVector += n.momentum();
      }
      double met = metVector.pt();
      if (met <= 20*GeV)  vetoEvent;

      if ( (_mT(muCandidate? muons[0] : elecs[0], metVector) + met) <= 60. )  vetoEvent;

      /****************
       *    B-JETS    *
       ****************/
      Jets bJets, wJets;
      foreach(Jet j, goodJets) {
        bool b_tagged = false;
        Particles bTags = j.bTags();
        foreach ( Particle b, bTags ) {
          b_tagged |= b.pT() > 5*GeV;
        }
        if (b_tagged)  bJets += j;
        if (!b_tagged && j.abseta() < 2.1)  wJets += j;
      }

      if ( bJets.size() < 2 || wJets.size() < 2 )  vetoEvent;

      double pull_angle = fabs(CalculatePullAngle(wJets[0], wJets[1], 0));
      h_pull_all->fill(pull_angle / Rivet::PI, weight);

      double pull_angle_charged = fabs(CalculatePullAngle(wJets[0], wJets[1], 1));
      h_pull_charged->fill(pull_angle_charged / Rivet::PI, weight);

    }

    Vector3 CalculatePull(Jet& jet, bool &isCharged) {
      Vector3 pull(0.0, 0.0, 0.0);
      double PT = jet.pT();
      Particles& constituents = jet.particles();
      Particles charged_constituents;
      if (isCharged) {
        foreach (Particle p, constituents) {
          if (p.threeCharge() != 0)  charged_constituents += p;
        }
        constituents = charged_constituents;
      }
      // calculate axis
      FourMomentum axis;
      foreach (Particle p, constituents)  axis += p.momentum();
      Vector3 J(axis.rap(), axis.phi(MINUSPI_PLUSPI), 0.0);
      // calculate pull
      foreach (Particle p, constituents) {
        Vector3 ri = Vector3(p.rap(), p.phi(MINUSPI_PLUSPI), 0.0) - J;
        while (ri.y() >  Rivet::PI) ri.setY(ri.y() - Rivet::TWOPI);
        while (ri.y() < -Rivet::PI) ri.setY(ri.y() + Rivet::TWOPI);
        pull.setX(pull.x() + (ri.mod() * ri.x() * p.pT()) / PT);
        pull.setY(pull.y() + (ri.mod() * ri.y() * p.pT()) / PT);
      }
      return pull;
    }

    double CalculatePullAngle(Jet& jet1, Jet& axisjet, bool isCharged) {
      Vector3 pull_vector = CalculatePull(jet1, isCharged);
      pull_vector = Vector3(1000.*pull_vector.x(), 1000.*pull_vector.y(), 0.);
      double drap = axisjet.rap() - jet1.rap();
      double dphi = axisjet.phi(MINUSPI_PLUSPI) - jet1.phi(MINUSPI_PLUSPI);
      Vector3 j2_vector(drap, dphi, 0.0);
      return mapAngleMPiToPi(deltaPhi(pull_vector, j2_vector));
    }

    double _mT(const FourMomentum &l, FourMomentum &nu) const {
      return sqrt( 2 * l.pT() * nu.pT() * (1 - cos(deltaPhi(l, nu))) );
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(h_pull_all);
      normalize(h_pull_charged);
    }

    //@}


  private:

    Histo1DPtr h_pull_all;
    Histo1DPtr h_pull_charged;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1376945);

}
