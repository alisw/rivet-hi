// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/RivetMT2.hh"

namespace Rivet {


  /// @author Peter Richardson
  class ATLAS_2012_I1112263 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1112263()
      : Analysis("ATLAS_2012_I1112263")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // projection to find the electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 10*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // projection to find the muons
      IdentifiedFinalState muons(Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");

      // for pTmiss
      declare(VisibleFinalState(Cuts::abseta < 4.9),"vfs");

      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);

      /// Jet finder
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // all tracks (to do deltaR with leptons)
      declare(ChargedFinalState(Cuts::abseta < 3.0),"cfs");

      // Book histograms
      _hist_leptonpT_SR1.push_back(bookHisto1D("hist_lepton_pT_1_SR1",11,0.,220.));
      _hist_leptonpT_SR1.push_back(bookHisto1D("hist_lepton_pT_2_SR1", 7,0.,140.));
      _hist_leptonpT_SR1.push_back(bookHisto1D("hist_lepton_pT_3_SR1", 8,0.,160.));
      _hist_leptonpT_SR2.push_back(bookHisto1D("hist_lepton_pT_1_SR2",11,0.,220.));
      _hist_leptonpT_SR2.push_back(bookHisto1D("hist_lepton_pT_2_SR2", 7,0.,140.));
      _hist_leptonpT_SR2.push_back(bookHisto1D("hist_lepton_pT_3_SR2", 8,0.,160.));
      _hist_etmiss_SR1_A = bookHisto1D("hist_etmiss_SR1_A",15,10.,310.);
      _hist_etmiss_SR1_B = bookHisto1D("hist_etmiss_SR1_B", 9,10.,190.);
      _hist_etmiss_SR2_A = bookHisto1D("hist_etmiss_SR2_A",15,10.,310.);
      _hist_etmiss_SR2_B = bookHisto1D("hist_etmiss_SR2_B", 9,10.,190.);
      _hist_mSFOS= bookHisto1D("hist_mSFOF",9,0.,180.);

      _count_SR1 = bookHisto1D("count_SR1", 1, 0., 1.);
      _count_SR2 = bookHisto1D("count_SR2", 1, 0., 1.);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the jet candidates
      Jets cand_jets;
      foreach (const Jet& jet, apply<FastJets>(event, "AntiKtJets04").jetsByPt(20.0*GeV) ) {
        if ( fabs( jet.eta() ) < 2.8 ) {
          cand_jets.push_back(jet);
        }
      }

      // Candidate muons
      Particles cand_mu;
      Particles chg_tracks =
        apply<ChargedFinalState>(event, "cfs").particles();
      foreach ( const Particle & mu, apply<IdentifiedFinalState>(event, "muons").particlesByPt() ) {
        double pTinCone = -mu.pT();
        foreach ( const Particle & track, chg_tracks ) {
          if ( deltaR(mu.momentum(),track.momentum()) <= 0.2 )
            pTinCone += track.pT();
        }
        if ( pTinCone < 1.8*GeV )
          cand_mu.push_back(mu);
      }

      // Candidate electrons
      Particles cand_e;
      foreach ( const Particle & e, apply<IdentifiedFinalState>(event, "elecs").particlesByPt() ) {
        double eta = e.eta();
        // Remove electrons with pT<15 in old veto region
        // (NOT EXPLICIT IN THIS PAPER BUT IN SIMILAR 4 LEPTON PAPER and THIS DESCRPITION
        //  IS MUCH WORSE SO ASSUME THIS IS DONE)
        if ( fabs(eta)>1.37 && fabs(eta) < 1.52 && e.perp()< 15.*GeV)
          continue;
        double pTinCone = -e.perp();
        foreach ( const Particle & track, chg_tracks ) {
          if ( deltaR(e.momentum(),track.momentum()) <= 0.2 )
            pTinCone += track.pT();
        }
        if (pTinCone/e.perp()<0.1) {
          cand_e.push_back(e);
        }
      }

      // Resolve jet/lepton ambiguity
      // (NOT EXPLICIT IN THIS PAPER BUT IN SIMILAR 4 LEPTON PAPER and THIS DESCRPITION
      //  IS MUCH WORSE SO ASSUME THIS IS DONE)
      Jets recon_jets;
      foreach ( const Jet& jet, cand_jets ) {
        bool away_from_e = true;
        foreach ( const Particle & e, cand_e ) {
          if ( deltaR(e.momentum(),jet.momentum()) <= 0.2 ) {
            away_from_e = false;
            break;
          }
        }
        if ( away_from_e )
          recon_jets.push_back( jet );
      }

      // Only keep electrons more than R=0.4 from jets
      Particles recon_e;
      foreach ( const Particle & e, cand_e ) {
        bool away = true;
        foreach ( const Jet& jet, recon_jets ) {
          if ( deltaR(e.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        // ... and 0.1 from any muons
        if ( ! away ) {
          foreach ( const Particle & mu, cand_e ) {
            if ( deltaR(mu.momentum(),e.momentum()) < 0.1 ) {
              away = false;
              break;
            }
          }
        }
        if ( away )
          recon_e.push_back( e );
      }
      // Only keep muons more than R=0.4 from jets
      Particles recon_mu;
      foreach ( const Particle & mu, cand_mu ) {
        bool away = true;
        foreach ( const Jet& jet, recon_jets ) {
          if ( deltaR(mu.momentum(),jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        // ... and 0.1 from any electrona
        if ( ! away ) {
          foreach ( const Particle & e, cand_e ) {
            if ( deltaR(mu.momentum(),e.momentum()) < 0.1 ) {
              away = false;
              break;
            }
          }
        }
        if ( away )
          recon_mu.push_back( mu );
      }

      // pTmiss
      Particles vfs_particles =
        apply<VisibleFinalState>(event, "vfs").particles();
      FourMomentum pTmiss;
      foreach ( const Particle & p, vfs_particles ) {
        pTmiss -= p.momentum();
      }
      double eTmiss = pTmiss.pT();

      // Now only use recon_jets, recon_mu, recon_e

      // Reject events with wrong number of leptons
      if ( recon_mu.size() + recon_e.size() != 3 ) {
        MSG_DEBUG("To few charged leptons left after selection");
        vetoEvent;
      }

      // ATLAS calo problem
      if (rand()/static_cast<double>(RAND_MAX) <= 0.42) {
        foreach ( const Particle & e, recon_e ) {
          double eta = e.eta();
          double phi = e.azimuthalAngle(MINUSPI_PLUSPI);
          if (inRange(eta, -0.1, 1.5) && inRange(phi, -0.9, -0.5)) vetoEvent;
        }
        foreach ( const Jet & jet, recon_jets ) {
          const double eta = jet.rapidity();
          const double phi = jet.azimuthalAngle(MINUSPI_PLUSPI);
          if (jet.perp() > 40*GeV && inRange(eta, -0.1, 1.5) && inRange(phi, -0.9, -0.5)) vetoEvent;
        }
      }

      if ( !( !recon_e .empty() && recon_e[0] .perp() > 25*GeV) &&
           !( !recon_mu.empty() && recon_mu[0].perp() > 20*GeV) ) {
        MSG_DEBUG("Hardest lepton fails trigger");
        vetoEvent;
      }

      // eTmiss cut
      if (eTmiss < 50*GeV) vetoEvent;

      // Check at least 1 SFOS pair
      double mSFOS=1e30, mdiff=1e30*GeV;
      size_t nSFOS=0;
      for (size_t ix = 0; ix < recon_e.size(); ++ix) {
        for (size_t iy = ix+1; iy < recon_e.size(); ++iy) {
          if (recon_e[ix].pid()*recon_e[iy].pid() > 0) continue;
          ++nSFOS;
          double mtest = (recon_e[ix].momentum() + recon_e[iy].momentum()).mass();
          // Veto is mass<20
          if (mtest < 20*GeV) vetoEvent;
          if (fabs(mtest - 90*GeV) < mdiff) {
            mSFOS = mtest;
            mdiff = fabs(mtest - 90*GeV);
          }
        }
      }
      for (size_t ix = 0; ix < recon_mu.size(); ++ix) {
        for (size_t iy = ix+1; iy < recon_mu.size(); ++iy) {
          if (recon_mu[ix].pid()*recon_mu[iy].pid() > 0) continue;
          ++nSFOS;
          double mtest = (recon_mu[ix].momentum() + recon_mu[iy].momentum()).mass();
          // Veto is mass < 20*GeV
          if (mtest < 20*GeV) vetoEvent;
          if (fabs(mtest - 90*GeV) < mdiff) {
            mSFOS = mtest;
            mdiff = fabs(mtest - 90*GeV);
          }
        }
      }
      // Require at least 1 SFOS pair
      if (nSFOS == 0) vetoEvent;
      // b-jet veto in SR!
      if (mdiff > 10*GeV) {
        foreach (const Jet & jet, recon_jets ) {
          if (jet.bTagged() && rand()/static_cast<double>(RAND_MAX) <= 0.60) vetoEvent;
        }
      }

	  // Histogram filling
	  const double weight = event.weight();

      // Region SR1, Z depleted
      if (mdiff > 10*GeV) {
        _count_SR1->fill(0.5, weight);
        _hist_etmiss_SR1_A->fill(eTmiss, weight);
        _hist_etmiss_SR1_B->fill(eTmiss, weight);
        _hist_mSFOS->fill(mSFOS, weight);
      }
      // Region SR2, Z enriched
      else {
        _count_SR2->fill(0.5, weight);
        _hist_etmiss_SR2_A->fill(eTmiss, weight);
        _hist_etmiss_SR2_B->fill(eTmiss, weight);
      }
      // Make the control plots
      // lepton pT
      size_t ie=0, imu=0;
      for (size_t ix = 0; ix < 3; ++ix) {
        Histo1DPtr hist = (mdiff > 10*GeV) ? _hist_leptonpT_SR1[ix] :  _hist_leptonpT_SR2[ix];
        double pTe  = (ie  < recon_e .size()) ? recon_e [ie ].perp() : -1*GeV;
        double pTmu = (imu < recon_mu.size()) ? recon_mu[imu].perp() : -1*GeV;
        if (pTe > pTmu) {
          hist->fill(pTe, weight);
          ++ie;
        } else {
          hist->fill(pTmu, weight);
          ++imu;
        }
      }

    }


    void finalize() {
      const double norm = crossSection()/femtobarn*2.06/sumOfWeights();
      // These are number of events at 2.06fb^-1 per 20 GeV
      for (size_t ix = 0; ix < 3; ++ix) {
        scale(_hist_leptonpT_SR1[ix], norm*20.);
        scale(_hist_leptonpT_SR2[ix], norm*20.);
      }
      scale(_hist_etmiss_SR1_A, norm*20.);
      scale(_hist_etmiss_SR1_B, norm*20.);
      scale(_hist_etmiss_SR2_A, norm*20.);
      scale(_hist_etmiss_SR2_B, norm*20.);
      scale(_hist_mSFOS, norm*20.);
      // These are number of events at 2.06fb^-1
      scale(_count_SR1, norm);
      scale(_count_SR2, norm);
    }

    //@}


  private:

  /// @name Histograms
  //@{
  vector<Histo1DPtr> _hist_leptonpT_SR1, _hist_leptonpT_SR2;
  Histo1DPtr _hist_etmiss_SR1_A, _hist_etmiss_SR1_B, _hist_etmiss_SR2_A, _hist_etmiss_SR2_B;
  Histo1DPtr _hist_mSFOS, _count_SR1, _count_SR2;
  //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1112263);

}
