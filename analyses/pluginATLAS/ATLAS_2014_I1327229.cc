// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"


namespace Rivet {


  class ATLAS_2014_I1327229 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1327229()
      : Analysis("ATLAS_2014_I1327229") {    }


    /// Book histograms and initialise projections before the run
    void init() {

      // To calculate the acceptance without having the fiducial lepton efficiencies included, this part can be turned off
      _use_fiducial_lepton_efficiency = true;

      // Random numbers for simulation of ATLAS detector reconstruction efficiency
      /// @todo Replace with SmearedParticles etc.
      srand(160385);

      // Read in all signal regions
      _signal_regions = getSignalRegions();

      // Set number of events per signal region to 0
      for (size_t i = 0; i < _signal_regions.size(); i++)
        _eventCountsPerSR[_signal_regions[i]] = 0.0;

      // Final state including all charged and neutral particles
      const FinalState fs(-5.0, 5.0, 1*GeV);
      declare(fs, "FS");

      // Final state including all charged particles
      declare(ChargedFinalState(-2.5, 2.5, 1*GeV), "CFS");

      // Final state including all visible particles (to calculate MET, Jets etc.)
      declare(VisibleFinalState(-5.0,5.0),"VFS");

      // Final state including all AntiKt 04 Jets
      VetoedFinalState vfs;
      vfs.addVetoPairId(PID::MUON);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "AntiKtJets04");

      // Final state including all unstable particles (including taus)
      declare(UnstableFinalState(Cuts::abseta < 5.0 && Cuts::pT > 5*GeV),"UFS");

      // Final state including all electrons
      IdentifiedFinalState elecs(Cuts::abseta < 2.47 && Cuts::pT > 10*GeV);
      elecs.acceptIdPair(PID::ELECTRON);
      declare(elecs, "elecs");

      // Final state including all muons
      IdentifiedFinalState muons(Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "muons");



      /// Book histograms:
      _h_HTlep_all = bookHisto1D("HTlep_all", 30,0,3000);
      _h_HTjets_all = bookHisto1D("HTjets_all", 30,0,3000);
      _h_MET_all = bookHisto1D("MET_all", 30,0,1500);
      _h_Meff_all = bookHisto1D("Meff_all", 50,0,5000);
      _h_min_pT_all = bookHisto1D("min_pT_all", 50, 0, 2000);
      _h_mT_all = bookHisto1D("mT_all", 50, 0, 2000);

      _h_e_n = bookHisto1D("e_n", 10, -0.5, 9.5);
      _h_mu_n = bookHisto1D("mu_n", 10, -0.5, 9.5);
      _h_tau_n = bookHisto1D("tau_n", 10, -0.5, 9.5);

      _h_pt_1_3l = bookHisto1D("pt_1_3l", 100, 0, 2000);
      _h_pt_2_3l = bookHisto1D("pt_2_3l", 100, 0, 2000);
      _h_pt_3_3l = bookHisto1D("pt_3_3l", 100, 0, 2000);
      _h_pt_1_2ltau = bookHisto1D("pt_1_2ltau", 100, 0, 2000);
      _h_pt_2_2ltau = bookHisto1D("pt_2_2ltau", 100, 0, 2000);
      _h_pt_3_2ltau = bookHisto1D("pt_3_2ltau", 100, 0, 2000);

      _h_excluded = bookHisto1D("excluded", 2, -0.5, 1.5);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Muons
      Particles muon_candidates;
      const Particles charged_tracks = apply<ChargedFinalState>(event, "CFS").particles();
      const Particles visible_particles = apply<VisibleFinalState>(event, "VFS").particles();
      for (const Particle& mu : apply<IdentifiedFinalState>(event, "muons").particlesByPt() ) {

        // Calculate pTCone30 variable (pT of all tracks within dR<0.3 - pT of muon itself)
        double pTinCone = -mu.pT();
        for (const Particle& track : charged_tracks ) {
          if (deltaR(mu.momentum(),track.momentum()) < 0.3 )
            pTinCone += track.pT();
        }

        // Calculate eTCone30 variable (pT of all visible particles within dR<0.3)
        double eTinCone = 0.;
        for (const Particle& visible_particle : visible_particles) {
          if (visible_particle.abspid() != PID::MUON && inRange(deltaR(mu.momentum(),visible_particle.momentum()), 0.1, 0.3))
            eTinCone += visible_particle.pT();
        }

        // Apply reconstruction efficiency and simulate reconstruction
        int muon_id = 13;
        if (mu.hasAncestor(PID::TAU) || mu.hasAncestor(-PID::TAU)) muon_id = 14;
        const double eff = (_use_fiducial_lepton_efficiency) ? apply_reco_eff(muon_id,mu) : 1.0;
        const bool keep_muon = rand()/static_cast<double>(RAND_MAX)<=eff;

        // Keep muon if pTCone30/pT < 0.15 and eTCone30/pT < 0.2 and reconstructed
        if (keep_muon && pTinCone/mu.pT() <= 0.1 && eTinCone/mu.pT() < 0.1)
          muon_candidates.push_back(mu);
      }

      // Electrons
      Particles electron_candidates;
      for (const Particle& e : apply<IdentifiedFinalState>(event, "elecs").particlesByPt() ) {
        // Neglect electrons in crack regions
        if (inRange(e.abseta(), 1.37, 1.52)) continue;

        // Calculate pTCone30 variable (pT of all tracks within dR<0.3 - pT of electron itself)
        double pTinCone = -e.pT();
        for (const Particle& track : charged_tracks) {
          if (deltaR(e.momentum(), track.momentum()) < 0.3 ) pTinCone += track.pT();
        }

        // Calculate eTCone30 variable (pT of all visible particles (except muons) within dR<0.3)
        double eTinCone = 0.;
        for (const Particle& visible_particle : visible_particles) {
          if (visible_particle.abspid() != PID::MUON && inRange(deltaR(e.momentum(),visible_particle.momentum()), 0.1, 0.3))
            eTinCone += visible_particle.pT();
        }

        // Apply reconstruction efficiency and simulate reconstruction
        int elec_id = 11;
        if (e.hasAncestor(15) || e.hasAncestor(-15)) elec_id = 12;
        const double eff = (_use_fiducial_lepton_efficiency) ? apply_reco_eff(elec_id,e) : 1.0;
        const bool keep_elec = rand()/static_cast<double>(RAND_MAX)<=eff;

        // Keep electron if pTCone30/pT < 0.13 and eTCone30/pT < 0.2 and reconstructed
        if (keep_elec && pTinCone/e.pT() <= 0.1  && eTinCone/e.pT() < 0.1)
          electron_candidates.push_back(e);
      }


      // Taus
      Particles tau_candidates;
      for (const Particle& tau : apply<UnstableFinalState>(event, "UFS").particles() ) {
        // Only pick taus out of all unstable particles
        if ( tau.abspid() != PID::TAU) continue;
        // Check that tau has decayed into daughter particles
        if (tau.genParticle()->end_vertex() == 0) continue;
       // Calculate visible tau momentum using the tau neutrino momentum in the tau decay
        FourMomentum daughter_tau_neutrino_momentum = get_tau_neutrino_momentum(tau);
        Particle tau_vis = tau;
        tau_vis.setMomentum(tau.momentum()-daughter_tau_neutrino_momentum);
        // keep only taus in certain eta region and above 15 GeV of visible tau pT
        if ( tau_vis.pT()/GeV <= 15.0 || tau_vis.abseta() > 2.5) continue;

        // Get prong number (number of tracks) in tau decay and check if tau decays leptonically
        unsigned int nprong = 0;
        bool lep_decaying_tau = false;
        get_prong_number(tau.genParticle(),nprong,lep_decaying_tau);

        // Apply reconstruction efficiency and simulate reconstruction
        int tau_id = 15;
        if (nprong == 1) tau_id = 15;
        else if (nprong == 3) tau_id = 16;


        const double eff = (_use_fiducial_lepton_efficiency) ? apply_reco_eff(tau_id,tau_vis) : 1.0;
        const bool keep_tau = rand()/static_cast<double>(RAND_MAX)<=eff;

        // Keep tau if nprong = 1, it decays hadronically and it is reconstructed
        if ( !lep_decaying_tau && nprong == 1 && keep_tau) tau_candidates.push_back(tau_vis);
      }

      // Jets (all anti-kt R=0.4 jets with pT > 30 GeV and eta < 4.9
      Jets jet_candidates;
      for (const Jet& jet : apply<FastJets>(event, "AntiKtJets04").jetsByPt(30.0*GeV) ) {
        if (jet.abseta() < 4.9 ) jet_candidates.push_back(jet);
      }

      // ETmiss
      Particles vfs_particles = apply<VisibleFinalState>(event, "VFS").particles();
      FourMomentum pTmiss;
      for (const Particle& p : vfs_particles)
        pTmiss -= p.momentum();
      double eTmiss = pTmiss.pT()/GeV;


      // -------------------------
      // Overlap removal

      // electron - electron
      Particles electron_candidates_2;
      for(size_t ie = 0; ie < electron_candidates.size(); ++ie) {
        const Particle& e = electron_candidates[ie];
        bool away = true;
        // If electron pair within dR < 0.1: remove electron with lower pT
        for(size_t ie2 = 0; ie2 < electron_candidates_2.size(); ++ie2) {
          if (deltaR(e.momentum(),electron_candidates_2[ie2].momentum()) < 0.1 ) {
            away = false;
            break;
          }
        }
        // If isolated keep it
        if ( away )
          electron_candidates_2.push_back( e );
      }

      // jet - electron
      Jets recon_jets;
      for (const Jet& jet : jet_candidates) {
        bool away = true;
        // If jet within dR < 0.2 of electron: remove jet
        for (const Particle& e : electron_candidates_2) {
          if (deltaR(e.momentum(), jet.momentum()) < 0.2 ) {
            away = false;
            break;
          }
        }

        // jet - tau
        if ( away )  {
          // If jet within dR < 0.2 of tau: remove jet
          for (const Particle& tau : tau_candidates) {
            if (deltaR(tau.momentum(), jet.momentum()) < 0.2 ) {
              away = false;
              break;
            }
          }
        }
        // If isolated keep it
        if ( away )
          recon_jets.push_back( jet );
      }

      // electron - jet
      Particles recon_leptons, recon_e;
      for (size_t ie = 0; ie < electron_candidates_2.size(); ++ie) {
        const Particle& e = electron_candidates_2[ie];
        // If electron within 0.2 < dR < 0.4 from any jets: remove electron
        bool away = true;
        for (const Jet& jet : recon_jets) {
          if (deltaR(e.momentum(), jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        // electron - muon
        // If electron within dR < 0.1 of a muon: remove electron
        if (away) {
          for (const Particle& mu : muon_candidates) {
            if (deltaR(mu.momentum(),e.momentum()) < 0.1) {
              away = false;
              break;
            }
          }
        }
        // If isolated keep it
        if ( away )  {
          recon_e.push_back( e );
          recon_leptons.push_back( e );
          }
      }

      // tau - electron
      Particles recon_tau;
      for (const Particle& tau : tau_candidates) {
        bool away = true;
        // If tau within dR < 0.2 of an electron: remove tau
        for (const Particle & e : recon_e) {
          if (deltaR(tau.momentum(),e.momentum()) < 0.2 ) {
            away = false;
            break;
          }
        }
        // tau - muon
        // If tau within dR < 0.2 of a muon: remove tau
        if (away)  {
          for (const Particle& mu : muon_candidates) {
            if (deltaR(tau.momentum(), mu.momentum()) < 0.2 ) {
              away = false;
              break;
            }
          }
        }
        // If isolated keep it
        if (away) recon_tau.push_back( tau );
      }

      // muon - jet
      Particles recon_mu, trigger_mu;
      // If muon within dR < 0.4 of a jet: remove muon
      for (const Particle& mu : muon_candidates ) {
        bool away = true;
        for (const Jet& jet : recon_jets) {
          if (deltaR(mu.momentum(), jet.momentum()) < 0.4 ) {
            away = false;
            break;
          }
        }
        if (away)  {
          recon_mu.push_back( mu );
          recon_leptons.push_back( mu );
          if (mu.abseta() < 2.4) trigger_mu.push_back( mu );
        }
      }

      // End overlap removal
      // ---------------------

      // Jet cleaning
      if (rand()/static_cast<double>(RAND_MAX) <= 0.42) {
        for (const Jet& jet : recon_jets ) {
          const double eta = jet.rapidity();
          const double phi = jet.azimuthalAngle(MINUSPI_PLUSPI);
          if(jet.pT() > 25*GeV && inRange(eta,-0.1,1.5) && inRange(phi,-0.9,-0.5)) vetoEvent;
        }
      }

      // Event selection
      // Require at least 3 charged tracks in event
      if (charged_tracks.size() < 3) vetoEvent;

      // And at least one e/mu passing trigger
      if( !( !recon_e.empty() && recon_e[0].pT()>26.*GeV)  &&
          !( !trigger_mu.empty() && trigger_mu[0].pT()>26.*GeV) ) {
        MSG_DEBUG("Hardest lepton fails trigger");
        vetoEvent;
      }

      // And only accept events with at least 2 electrons and muons and at least 3 leptons in total
      if (recon_mu.size() + recon_e.size() + recon_tau.size() < 3 || recon_leptons.size() < 2) vetoEvent;

      // Getting the event weight
      const double weight = event.weight();

      // Sort leptons by decreasing pT
      sortByPt(recon_leptons);
      sortByPt(recon_tau);

      // Calculate HTlep, fill lepton pT histograms & store chosen combination of 3 leptons
      double HTlep = 0.;
      Particles chosen_leptons;
      if (recon_leptons.size() > 2) {
        _h_pt_1_3l->fill(recon_leptons[0].pT()/GeV, weight);
        _h_pt_2_3l->fill(recon_leptons[1].pT()/GeV, weight);
        _h_pt_3_3l->fill(recon_leptons[2].pT()/GeV, weight);
        HTlep = (recon_leptons[0].pT() + recon_leptons[1].pT() + recon_leptons[2].pT())/GeV;
        chosen_leptons.push_back( recon_leptons[0] );
        chosen_leptons.push_back( recon_leptons[1] );
        chosen_leptons.push_back( recon_leptons[2] );
      }
      else {
        _h_pt_1_2ltau->fill(recon_leptons[0].pT()/GeV, weight);
        _h_pt_2_2ltau->fill(recon_leptons[1].pT()/GeV, weight);
        _h_pt_3_2ltau->fill(recon_tau[0].pT()/GeV, weight);
        HTlep = recon_leptons[0].pT()/GeV + recon_leptons[1].pT()/GeV + recon_tau[0].pT()/GeV;
        chosen_leptons.push_back( recon_leptons[0] );
        chosen_leptons.push_back( recon_leptons[1] );
        chosen_leptons.push_back( recon_tau[0] );
      }

      // Calculate mT and mTW variable
      Particles mT_leptons;
      Particles mTW_leptons;
      for (size_t i1 = 0; i1 < 3; i1 ++)  {
        for (size_t i2 = i1+1; i2 < 3; i2 ++)  {
          double OSSF_inv_mass = isOSSF_mass(chosen_leptons[i1],chosen_leptons[i2]);
          if (OSSF_inv_mass != 0.)  {
            for (size_t i3 = 0; i3 < 3 ; i3 ++)  {
              if (i3 != i2 && i3 != i1)  {
                mT_leptons.push_back(chosen_leptons[i3]);
                if ( fabs(91.0 - OSSF_inv_mass) < 20. )
                  mTW_leptons.push_back(chosen_leptons[i3]);
              }
            }
          }
          else  {
            mT_leptons.push_back(chosen_leptons[0]);
            mTW_leptons.push_back(chosen_leptons[0]);
          }
        }
      }

      sortByPt(mT_leptons);
      sortByPt(mTW_leptons);

      double mT = sqrt(2*pTmiss.pT()/GeV*mT_leptons[0].pT()/GeV*(1-cos(pTmiss.phi()-mT_leptons[0].phi())));
      double mTW = sqrt(2*pTmiss.pT()/GeV*mTW_leptons[0].pT()/GeV*(1-cos(pTmiss.phi()-mTW_leptons[0].phi())));

      // Calculate Min pT variable
      double min_pT = chosen_leptons[2].pT()/GeV;

      // Number of prompt e/mu and had taus
      _h_e_n->fill(recon_e.size(),weight);
      _h_mu_n->fill(recon_mu.size(),weight);
      _h_tau_n->fill(recon_tau.size(),weight);

      // Calculate HTjets variable
      double HTjets = 0.;
      for (const Jet& jet : recon_jets)
        HTjets += jet.pT()/GeV;

      // Calculate meff variable
      double meff = eTmiss + HTjets;
      Particles all_leptons;
      for (const Particle& e : recon_e ) {
        meff += e.pT()/GeV;
        all_leptons.push_back( e );
      }
      for (const Particle& mu : recon_mu) {
        meff += mu.pT()/GeV;
        all_leptons.push_back( mu );
      }
      for (const Particle& tau : recon_tau) {
        meff += tau.pT()/GeV;
        all_leptons.push_back( tau );
      }

      // Fill histograms of kinematic variables
      _h_HTlep_all->fill(HTlep,weight);
      _h_HTjets_all->fill(HTjets,weight);
      _h_MET_all->fill(eTmiss,weight);
      _h_Meff_all->fill(meff,weight);
      _h_min_pT_all->fill(min_pT,weight);
      _h_mT_all->fill(mT,weight);

      // Determine signal region (3l / 2ltau , onZ / offZ OSSF / offZ no-OSSF)
      // 3l vs. 2ltau
      string basic_signal_region;
      if (recon_mu.size() + recon_e.size() > 2)
        basic_signal_region += "3l_";
      else if ( (recon_mu.size() + recon_e.size() == 2) && (recon_tau.size() > 0))
        basic_signal_region += "2ltau_";

      // Is there an OSSF pair or a three lepton combination with an invariant mass close to the Z mass
      int onZ = isonZ(chosen_leptons);
      if (onZ == 1) basic_signal_region += "onZ";
      else if (onZ == 0)  {
        bool OSSF = isOSSF(chosen_leptons);
        if (OSSF) basic_signal_region += "offZ_OSSF";
        else basic_signal_region += "offZ_noOSSF";
        }

      // Check in which signal regions this event falls and adjust event counters
      // INFO: The b-jet signal regions of the paper are not included in this Rivet implementation
      fillEventCountsPerSR(basic_signal_region,onZ,HTlep,eTmiss,HTjets,meff,min_pT,mTW,weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Normalize to an integrated luminosity of 1 fb-1
      double norm = crossSection()/femtobarn/sumOfWeights();

      string best_signal_region = "";
      double ratio_best_SR = 0.;

      // Loop over all signal regions and find signal region with best sensitivity (ratio signal events/visible cross-section)
      for (size_t i = 0; i < _signal_regions.size(); i++) {
        double signal_events = _eventCountsPerSR[_signal_regions[i]] * norm;
        // Use expected upper limits to find best signal region:
        double UL95 = getUpperLimit(_signal_regions[i],false);
        double ratio = signal_events / UL95;
        if (ratio > ratio_best_SR)  {
          best_signal_region = _signal_regions.at(i);
          ratio_best_SR = ratio;
        }
      }

      double signal_events_best_SR = _eventCountsPerSR[best_signal_region] * norm;
      double exp_UL_best_SR = getUpperLimit(best_signal_region, false);
      double obs_UL_best_SR = getUpperLimit(best_signal_region, true);


      // Print out result
      cout << "----------------------------------------------------------------------------------------" << endl;
      cout << "Number of total events: " << sumOfWeights() << endl;
      cout << "Best signal region: " << best_signal_region << endl;
      cout << "Normalized number of signal events in this best signal region (per fb-1): " << signal_events_best_SR << endl;
      cout << "Efficiency*Acceptance: " << _eventCountsPerSR[best_signal_region]/sumOfWeights() << endl;
      cout << "Cross-section [fb]: " << crossSection()/femtobarn << endl;
      cout << "Expected visible cross-section (per fb-1): " << exp_UL_best_SR << endl;
      cout << "Ratio (signal events / expected visible cross-section): " << ratio_best_SR << endl;
      cout << "Observed visible cross-section (per fb-1): " << obs_UL_best_SR << endl;
      cout << "Ratio (signal events / observed visible cross-section): " <<  signal_events_best_SR/obs_UL_best_SR << endl;
      cout << "----------------------------------------------------------------------------------------" << endl;

      cout << "Using the EXPECTED limits (visible cross-section) of the analysis: " << endl;
      if (signal_events_best_SR > exp_UL_best_SR)  {
        cout << "Since the number of signal events > the visible cross-section, this model/grid point is EXCLUDED with 95% C.L." << endl;
        _h_excluded->fill(1);
      }
      else  {
        cout << "Since the number of signal events < the visible cross-section, this model/grid point is NOT EXCLUDED." << endl;
        _h_excluded->fill(0);
      }
      cout << "----------------------------------------------------------------------------------------" << endl;

      cout << "Using the OBSERVED limits (visible cross-section) of the analysis: " << endl;
      if (signal_events_best_SR > obs_UL_best_SR)  {
        cout << "Since the number of signal events > the visible cross-section, this model/grid point is EXCLUDED with 95% C.L." << endl;
        _h_excluded->fill(1);
      }
      else  {
        cout << "Since the number of signal events < the visible cross-section, this model/grid point is NOT EXCLUDED." << endl;
        _h_excluded->fill(0);
      }
      cout << "----------------------------------------------------------------------------------------" << endl;
      cout << "INFO: The b-jet signal regions of the paper are not included in this Rivet implementation." << endl;
      cout << "----------------------------------------------------------------------------------------" << endl;


      /// Normalize to cross section

      if (norm != 0)  {
        scale(_h_HTlep_all, norm);
        scale(_h_HTjets_all, norm);
        scale(_h_MET_all, norm);
        scale(_h_Meff_all, norm);
        scale(_h_min_pT_all, norm);
        scale(_h_mT_all, norm);

        scale(_h_pt_1_3l, norm);
        scale(_h_pt_2_3l, norm);
        scale(_h_pt_3_3l, norm);
        scale(_h_pt_1_2ltau, norm);
        scale(_h_pt_2_2ltau, norm);
        scale(_h_pt_3_2ltau, norm);

        scale(_h_e_n, norm);
        scale(_h_mu_n, norm);
        scale(_h_tau_n, norm);

        scale(_h_excluded, norm);
      }

    }


    /// Helper functions
    //@{
    /// Function giving a list of all signal regions
    vector<string> getSignalRegions()  {

      // List of basic signal regions
      vector<string> basic_signal_regions;
      basic_signal_regions.push_back("3l_offZ_OSSF");
      basic_signal_regions.push_back("3l_offZ_noOSSF");
      basic_signal_regions.push_back("3l_onZ");
      basic_signal_regions.push_back("2ltau_offZ_OSSF");
      basic_signal_regions.push_back("2ltau_offZ_noOSSF");
      basic_signal_regions.push_back("2ltau_onZ");

      // List of kinematic variables
      vector<string> kinematic_variables;
      kinematic_variables.push_back("HTlep");
      kinematic_variables.push_back("METStrong");
      kinematic_variables.push_back("METWeak");
      kinematic_variables.push_back("Meff");
      kinematic_variables.push_back("MeffStrong");
      kinematic_variables.push_back("MeffMt");
      kinematic_variables.push_back("MinPt");

      vector<string> signal_regions;
      // Loop over all kinematic variables and basic signal regions
      for (size_t i0 = 0; i0 < kinematic_variables.size(); i0++)  {
        for (size_t i1 = 0; i1 < basic_signal_regions.size(); i1++)  {
          // Is signal region onZ?
          int onZ = (basic_signal_regions[i1].find("onZ") != string::npos) ? 1 : 0;
          // Get cut values for this kinematic variable
          vector<int> cut_values = getCutsPerSignalRegion(kinematic_variables[i0], onZ);
          // Loop over all cut values
          for (size_t i2 = 0; i2 < cut_values.size(); i2++)  {
            // Push signal region into vector
            signal_regions.push_back( kinematic_variables[i0] + "_" + basic_signal_regions[i1] + "_cut_" + toString(cut_values[i2]) );
          }
        }
      }
      return signal_regions;
    }



    /// Function giving all cut values per kinematic variable
    vector<int> getCutsPerSignalRegion(const string& signal_region, int onZ = 0)  {
      vector<int> cutValues;

      // Cut values for HTlep
      if (signal_region.compare("HTlep") == 0)  {
        cutValues.push_back(0);
        cutValues.push_back(200);
        cutValues.push_back(500);
        cutValues.push_back(800);
        }
      // Cut values for MinPt
      else if (signal_region.compare("MinPt") == 0)  {
        cutValues.push_back(0);
        cutValues.push_back(50);
        cutValues.push_back(100);
        cutValues.push_back(150);
        }
      // Cut values for METStrong (HTjets > 150 GeV) and METWeak (HTjets < 150 GeV)
      else if (signal_region.compare("METStrong") == 0 || signal_region.compare("METWeak") == 0)  {
        cutValues.push_back(0);
        cutValues.push_back(100);
        cutValues.push_back(200);
        cutValues.push_back(300);
        }
      // Cut values for Meff
      if (signal_region.compare("Meff") == 0)  {
        cutValues.push_back(0);
        cutValues.push_back(600);
        cutValues.push_back(1000);
        cutValues.push_back(1500);
        }
      // Cut values for MeffStrong (MET > 100 GeV)
      if ((signal_region.compare("MeffStrong") == 0 || signal_region.compare("MeffMt") == 0) && onZ ==1)  {
        cutValues.push_back(0);
        cutValues.push_back(600);
        cutValues.push_back(1200);
        }

      return cutValues;
    }

    /// function fills map _eventCountsPerSR by looping over all signal regions
    /// and looking if the event falls into this signal region
    void fillEventCountsPerSR(const string& basic_signal_region, int onZ,
                              double HTlep, double eTmiss, double HTjets,
                              double meff, double min_pT, double mTW,
                              double weight)  {

      // Get cut values for HTlep, loop over them and add event if cut is passed
      vector<int> cut_values = getCutsPerSignalRegion("HTlep", onZ);
      for (size_t i = 0; i < cut_values.size(); i++)  {
        if (HTlep > cut_values[i])
          _eventCountsPerSR[("HTlep_" + basic_signal_region + "_cut_" + toString(cut_values[i]))] += weight;
      }

      // Get cut values for MinPt, loop over them and add event if cut is passed
      cut_values = getCutsPerSignalRegion("MinPt", onZ);
      for (size_t i = 0; i < cut_values.size(); i++)  {
        if (min_pT > cut_values[i])
          _eventCountsPerSR[("MinPt_" + basic_signal_region + "_cut_" + toString(cut_values[i]))] += weight;
      }

      // Get cut values for METStrong, loop over them and add event if cut is passed
      cut_values = getCutsPerSignalRegion("METStrong", onZ);
      for (size_t i = 0; i < cut_values.size(); i++)  {
        if (eTmiss > cut_values[i] && HTjets > 150.)
          _eventCountsPerSR[("METStrong_" + basic_signal_region + "_cut_" + toString(cut_values[i]))] += weight;
      }

      // Get cut values for METWeak, loop over them and add event if cut is passed
      cut_values = getCutsPerSignalRegion("METWeak", onZ);
      for (size_t i = 0; i < cut_values.size(); i++)  {
        if (eTmiss > cut_values[i] && HTjets <= 150.)
          _eventCountsPerSR[("METWeak_" + basic_signal_region + "_cut_" + toString(cut_values[i]))] += weight;
      }

      // Get cut values for Meff, loop over them and add event if cut is passed
      cut_values = getCutsPerSignalRegion("Meff", onZ);
      for (size_t i = 0; i < cut_values.size(); i++)  {
        if (meff > cut_values[i])
          _eventCountsPerSR[("Meff_" + basic_signal_region + "_cut_" + toString(cut_values[i]))] += weight;
      }

      // Get cut values for MeffStrong, loop over them and add event if cut is passed
      cut_values = getCutsPerSignalRegion("MeffStrong", onZ);
      for (size_t i = 0; i < cut_values.size(); i++)  {
        if (meff > cut_values[i] && eTmiss > 100.)
          _eventCountsPerSR[("MeffStrong_" + basic_signal_region + "_cut_" + toString(cut_values[i]))] += weight;
      }

      // Get cut values for MeffMt, loop over them and add event if cut is passed
      cut_values = getCutsPerSignalRegion("MeffMt", onZ);
      for (size_t i = 0; i < cut_values.size(); i++)  {
        if (meff > cut_values[i] && mTW > 100. && onZ == 1)
          _eventCountsPerSR[("MeffMt_" + basic_signal_region + "_cut_" + toString(cut_values[i]))] += weight;
      }

    }

    /// Function returning 4-momentum of daughter-particle if it is a tau neutrino
    FourMomentum get_tau_neutrino_momentum(const Particle& p)  {
      assert(p.abspid() == PID::TAU);
      const GenVertex* dv = p.genParticle()->end_vertex();
      assert(dv != NULL);
      // Loop over all daughter particles
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        if (abs((*pp)->pdg_id()) == PID::NU_TAU) return FourMomentum((*pp)->momentum());
      }
      return FourMomentum();
    }

    /// Function calculating the prong number of taus
    void get_prong_number(const GenParticle* p, unsigned int& nprong, bool& lep_decaying_tau)  {
      assert(p != NULL);
      const GenVertex* dv = p->end_vertex();
      assert(dv != NULL);
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        // If they have status 1 and are charged they will produce a track and the prong number is +1
        if ((*pp)->status() == 1 )  {
          const int id = (*pp)->pdg_id();
          if (Rivet::PID::charge(id) != 0 ) ++nprong;
          // Check if tau decays leptonically
          if (( abs(id) == PID::ELECTRON || abs(id) == PID::MUON || abs(id) == PID::TAU ) && abs(p->pdg_id()) == PID::TAU) lep_decaying_tau = true;
        }
        // If the status of the daughter particle is 2 it is unstable and the further decays are checked
        else if ((*pp)->status() == 2 )  {
          get_prong_number((*pp),nprong,lep_decaying_tau);
        }
      }
    }

    /// Function giving fiducial lepton efficiency
    double apply_reco_eff(int flavor, const Particle& p) {
      double pt = p.pT()/GeV;
      double eta = p.eta();

      double eff = 0.;

      if (flavor == 11) { // weight prompt electron -- now including data/MC ID SF in eff.
        double avgrate = 0.685;
        const static double wz_ele[] =  {0.0256,0.522,0.607,0.654,0.708,0.737,0.761,0.784,0.815,0.835,0.851,0.841,0.898};
        // double ewz_ele[] = {0.000257,0.00492,0.00524,0.00519,0.00396,0.00449,0.00538,0.00513,0.00773,0.00753,0.0209,0.0964,0.259};
        int ibin = 0;
        if(pt > 10  && pt < 15) ibin = 0;
        if(pt > 15  && pt < 20) ibin = 1;
        if(pt > 20  && pt < 25) ibin = 2;
        if(pt > 25  && pt < 30) ibin = 3;
        if(pt > 30  && pt < 40) ibin = 4;
        if(pt > 40  && pt < 50) ibin = 5;
        if(pt > 50  && pt < 60) ibin = 6;
        if(pt > 60  && pt < 80) ibin = 7;
        if(pt > 80  && pt < 100) ibin = 8;
        if(pt > 100 && pt < 200) ibin = 9;
        if(pt > 200 && pt < 400) ibin = 10;
        if(pt > 400 && pt < 600) ibin = 11;
        if(pt > 600) ibin = 12;
        double eff_pt = 0.;
        eff_pt = wz_ele[ibin];

        eta = fabs(eta);

        const static double wz_ele_eta[] =  {0.65,0.714,0.722,0.689,0.635,0.615};
        // double ewz_ele_eta[] = {0.00642,0.00355,0.00335,0.004,0.00368,0.00422};
        ibin = 0;
        if(eta > 0 && eta < 0.1) ibin = 0;
        if(eta > 0.1 && eta < 0.5) ibin = 1;
        if(eta > 0.5 && eta < 1.0) ibin = 2;
        if(eta > 1.0 && eta < 1.5) ibin = 3;
        if(eta > 1.5 && eta < 2.0) ibin = 4;
        if(eta > 2.0 && eta < 2.5) ibin = 5;
        double eff_eta = 0.;
        eff_eta = wz_ele_eta[ibin];

        eff = (eff_pt * eff_eta) / avgrate;
      }

      if (flavor == 12) { // weight electron from tau
        double avgrate = 0.476;
        const static double wz_ele[] =  {0.00855,0.409,0.442,0.55,0.632,0.616,0.615,0.642,0.72,0.617};
        // double ewz_ele[] = {0.000573,0.0291,0.0366,0.0352,0.0363,0.0474,0.0628,0.0709,0.125,0.109};
        int ibin = 0;
        if(pt > 10  && pt < 15) ibin = 0;
        if(pt > 15  && pt < 20) ibin = 1;
        if(pt > 20  && pt < 25) ibin = 2;
        if(pt > 25  && pt < 30) ibin = 3;
        if(pt > 30  && pt < 40) ibin = 4;
        if(pt > 40  && pt < 50) ibin = 5;
        if(pt > 50  && pt < 60) ibin = 6;
        if(pt > 60  && pt < 80) ibin = 7;
        if(pt > 80  && pt < 100) ibin = 8;
        if(pt > 100)           ibin = 9;
        double eff_pt = 0.;
        eff_pt = wz_ele[ibin];

        eta = fabs(eta);

        const static double wz_ele_eta[] =  {0.546,0.5,0.513,0.421,0.47,0.433};
        //double ewz_ele_eta[] = {0.0566,0.0257,0.0263,0.0263,0.0303,0.0321};
        ibin = 0;
        if(eta > 0 && eta < 0.1) ibin = 0;
        if(eta > 0.1 && eta < 0.5) ibin = 1;
        if(eta > 0.5 && eta < 1.0) ibin = 2;
        if(eta > 1.0 && eta < 1.5) ibin = 3;
        if(eta > 1.5 && eta < 2.0) ibin = 4;
        if(eta > 2.0 && eta < 2.5) ibin = 5;
        double eff_eta = 0.;
        eff_eta = wz_ele_eta[ibin];

        eff = (eff_pt * eff_eta) / avgrate;
      }

      if (flavor == 13) { // weight prompt muon
        int ibin = 0;
        if(pt > 10  && pt < 15) ibin = 0;
        if(pt > 15  && pt < 20) ibin = 1;
        if(pt > 20  && pt < 25) ibin = 2;
        if(pt > 25  && pt < 30) ibin = 3;
        if(pt > 30  && pt < 40) ibin = 4;
        if(pt > 40  && pt < 50) ibin = 5;
        if(pt > 50  && pt < 60) ibin = 6;
        if(pt > 60  && pt < 80) ibin = 7;
        if(pt > 80  && pt < 100) ibin = 8;
        if(pt > 100 && pt < 200) ibin = 9;
        if(pt > 200 && pt < 400) ibin = 10;
        if(pt > 400) ibin = 11;
        if(fabs(eta) < 0.1) {
          const static double wz_mu[] =  {0.00705,0.402,0.478,0.49,0.492,0.499,0.527,0.512,0.53,0.528,0.465,0.465};
          //double ewz_mu[] = {0.000298,0.0154,0.017,0.0158,0.0114,0.0123,0.0155,0.0133,0.0196,0.0182,0.0414,0.0414};
          double eff_pt = 0.;
          eff_pt = wz_mu[ibin];
          eff = eff_pt;
        }
        if(fabs(eta) > 0.1) {
          const static double wz_mu[] =  {0.0224,0.839,0.887,0.91,0.919,0.923,0.925,0.925,0.922,0.918,0.884,0.834};
          //double ewz_mu[] = {0.000213,0.00753,0.0074,0.007,0.00496,0.00534,0.00632,0.00583,0.00849,0.00804,0.0224,0.0963};
          double eff_pt = 0.;
          eff_pt = wz_mu[ibin];
          eff = eff_pt;
        }
      }

      if (flavor == 14) { // weight muon from tau
        int ibin = 0;
        if(pt > 10  && pt < 15) ibin = 0;
        if(pt > 15  && pt < 20) ibin = 1;
        if(pt > 20  && pt < 25) ibin = 2;
        if(pt > 25  && pt < 30) ibin = 3;
        if(pt > 30  && pt < 40) ibin = 4;
        if(pt > 40  && pt < 50) ibin = 5;
        if(pt > 50  && pt < 60) ibin = 6;
        if(pt > 60  && pt < 80) ibin = 7;
        if(pt > 80  && pt < 100) ibin = 8;
        if(pt > 100) ibin = 9;

        if(fabs(eta) < 0.1) {
          const static double wz_mu[] =  {0.0,0.664,0.124,0.133,0.527,0.283,0.495,0.25,0.5,0.331};
          //double ewz_mu[] = {0.0,0.192,0.0437,0.0343,0.128,0.107,0.202,0.125,0.25,0.191};
          double eff_pt = 0.;
          eff_pt = wz_mu[ibin];
          eff = eff_pt;
        }
        if(fabs(eta) > 0.1) {
          const static double wz_mu[] =  {0.0,0.617,0.655,0.676,0.705,0.738,0.712,0.783,0.646,0.745};
          //double ewz_mu[] = {0.0,0.043,0.0564,0.0448,0.0405,0.0576,0.065,0.0825,0.102,0.132};
          double eff_pt = 0.;
          eff_pt = wz_mu[ibin];
          eff = eff_pt;
        }
      }

      if (flavor == 15) { // weight hadronic tau 1p
        double avgrate = 0.16;
        const static double wz_tau1p[] =  {0.0,0.0311,0.148,0.229,0.217,0.292,0.245,0.307,0.227,0.277};
        //double ewz_tau1p[] = {0.0,0.00211,0.0117,0.0179,0.0134,0.0248,0.0264,0.0322,0.0331,0.0427};
        int ibin = 0;
        if(pt > 10  && pt < 15) ibin = 0;
        if(pt > 15  && pt < 20) ibin = 1;
        if(pt > 20  && pt < 25) ibin = 2;
        if(pt > 25  && pt < 30) ibin = 3;
        if(pt > 30  && pt < 40) ibin = 4;
        if(pt > 40  && pt < 50) ibin = 5;
        if(pt > 50  && pt < 60) ibin = 6;
        if(pt > 60  && pt < 80) ibin = 7;
        if(pt > 80  && pt < 100) ibin = 8;
        if(pt > 100) ibin = 9;
        double eff_pt = 0.;
        eff_pt = wz_tau1p[ibin];

        const static double wz_tau1p_eta[] = {0.166,0.15,0.188,0.175,0.142,0.109};
        //double ewz_tau1p_eta[] ={0.0166,0.00853,0.0097,0.00985,0.00949,0.00842};
        ibin = 0;
        if(eta > 0.0 && eta < 0.1) ibin = 0;
        if(eta > 0.1 && eta < 0.5) ibin = 1;
        if(eta > 0.5 && eta < 1.0) ibin = 2;
        if(eta > 1.0 && eta < 1.5) ibin = 3;
        if(eta > 1.5 && eta < 2.0) ibin = 4;
        if(eta > 2.0 && eta < 2.5) ibin = 5;
        double eff_eta = 0.;
        eff_eta = wz_tau1p_eta[ibin];

        eff = (eff_pt * eff_eta) / avgrate;
      }

      return eff;
    }


    /// Function giving observed and expected upper limits (on the visible cross-section)
    double getUpperLimit(const string& signal_region, bool observed)  {

      map<string,double> upperLimitsObserved;
      map<string,double> upperLimitsExpected;

      upperLimitsObserved["HTlep_3l_offZ_OSSF_cut_0"] = 2.435;
      upperLimitsObserved["HTlep_3l_offZ_OSSF_cut_200"] = 0.704;
      upperLimitsObserved["HTlep_3l_offZ_OSSF_cut_500"] = 0.182;
      upperLimitsObserved["HTlep_3l_offZ_OSSF_cut_800"] = 0.147;
      upperLimitsObserved["HTlep_2ltau_offZ_OSSF_cut_0"] = 13.901;
      upperLimitsObserved["HTlep_2ltau_offZ_OSSF_cut_200"] = 1.677;
      upperLimitsObserved["HTlep_2ltau_offZ_OSSF_cut_500"] = 0.141;
      upperLimitsObserved["HTlep_2ltau_offZ_OSSF_cut_800"] = 0.155;
      upperLimitsObserved["HTlep_3l_offZ_noOSSF_cut_0"] = 1.054;
      upperLimitsObserved["HTlep_3l_offZ_noOSSF_cut_200"] = 0.341;
      upperLimitsObserved["HTlep_3l_offZ_noOSSF_cut_500"] = 0.221;
      upperLimitsObserved["HTlep_3l_offZ_noOSSF_cut_800"] = 0.140;
      upperLimitsObserved["HTlep_2ltau_offZ_noOSSF_cut_0"] = 4.276;
      upperLimitsObserved["HTlep_2ltau_offZ_noOSSF_cut_200"] = 0.413;
      upperLimitsObserved["HTlep_2ltau_offZ_noOSSF_cut_500"] = 0.138;
      upperLimitsObserved["HTlep_2ltau_offZ_noOSSF_cut_800"] = 0.150;
      upperLimitsObserved["HTlep_3l_onZ_cut_0"] = 29.804;
      upperLimitsObserved["HTlep_3l_onZ_cut_200"] = 3.579;
      upperLimitsObserved["HTlep_3l_onZ_cut_500"] = 0.466;
      upperLimitsObserved["HTlep_3l_onZ_cut_800"] = 0.298;
      upperLimitsObserved["HTlep_2ltau_onZ_cut_0"] = 205.091;
      upperLimitsObserved["HTlep_2ltau_onZ_cut_200"] = 3.141;
      upperLimitsObserved["HTlep_2ltau_onZ_cut_500"] = 0.290;
      upperLimitsObserved["HTlep_2ltau_onZ_cut_800"] = 0.157;
      upperLimitsObserved["METStrong_3l_offZ_OSSF_cut_0"] = 1.111;
      upperLimitsObserved["METStrong_3l_offZ_OSSF_cut_100"] = 0.354;
      upperLimitsObserved["METStrong_3l_offZ_OSSF_cut_200"] = 0.236;
      upperLimitsObserved["METStrong_3l_offZ_OSSF_cut_300"] = 0.150;
      upperLimitsObserved["METStrong_2ltau_offZ_OSSF_cut_0"] = 1.881;
      upperLimitsObserved["METStrong_2ltau_offZ_OSSF_cut_100"] = 0.406;
      upperLimitsObserved["METStrong_2ltau_offZ_OSSF_cut_200"] = 0.194;
      upperLimitsObserved["METStrong_2ltau_offZ_OSSF_cut_300"] = 0.134;
      upperLimitsObserved["METStrong_3l_offZ_noOSSF_cut_0"] = 0.770;
      upperLimitsObserved["METStrong_3l_offZ_noOSSF_cut_100"] = 0.295;
      upperLimitsObserved["METStrong_3l_offZ_noOSSF_cut_200"] = 0.149;
      upperLimitsObserved["METStrong_3l_offZ_noOSSF_cut_300"] = 0.140;
      upperLimitsObserved["METStrong_2ltau_offZ_noOSSF_cut_0"] = 2.003;
      upperLimitsObserved["METStrong_2ltau_offZ_noOSSF_cut_100"] = 0.806;
      upperLimitsObserved["METStrong_2ltau_offZ_noOSSF_cut_200"] = 0.227;
      upperLimitsObserved["METStrong_2ltau_offZ_noOSSF_cut_300"] = 0.138;
      upperLimitsObserved["METStrong_3l_onZ_cut_0"] = 6.383;
      upperLimitsObserved["METStrong_3l_onZ_cut_100"] = 0.959;
      upperLimitsObserved["METStrong_3l_onZ_cut_200"] = 0.549;
      upperLimitsObserved["METStrong_3l_onZ_cut_300"] = 0.182;
      upperLimitsObserved["METStrong_2ltau_onZ_cut_0"] = 10.658;
      upperLimitsObserved["METStrong_2ltau_onZ_cut_100"] = 0.637;
      upperLimitsObserved["METStrong_2ltau_onZ_cut_200"] = 0.291;
      upperLimitsObserved["METStrong_2ltau_onZ_cut_300"] = 0.227;
      upperLimitsObserved["METWeak_3l_offZ_OSSF_cut_0"] = 1.802;
      upperLimitsObserved["METWeak_3l_offZ_OSSF_cut_100"] = 0.344;
      upperLimitsObserved["METWeak_3l_offZ_OSSF_cut_200"] = 0.189;
      upperLimitsObserved["METWeak_3l_offZ_OSSF_cut_300"] = 0.148;
      upperLimitsObserved["METWeak_2ltau_offZ_OSSF_cut_0"] = 12.321;
      upperLimitsObserved["METWeak_2ltau_offZ_OSSF_cut_100"] = 0.430;
      upperLimitsObserved["METWeak_2ltau_offZ_OSSF_cut_200"] = 0.137;
      upperLimitsObserved["METWeak_2ltau_offZ_OSSF_cut_300"] = 0.134;
      upperLimitsObserved["METWeak_3l_offZ_noOSSF_cut_0"] = 0.562;
      upperLimitsObserved["METWeak_3l_offZ_noOSSF_cut_100"] = 0.153;
      upperLimitsObserved["METWeak_3l_offZ_noOSSF_cut_200"] = 0.154;
      upperLimitsObserved["METWeak_3l_offZ_noOSSF_cut_300"] = 0.141;
      upperLimitsObserved["METWeak_2ltau_offZ_noOSSF_cut_0"] = 2.475;
      upperLimitsObserved["METWeak_2ltau_offZ_noOSSF_cut_100"] = 0.244;
      upperLimitsObserved["METWeak_2ltau_offZ_noOSSF_cut_200"] = 0.141;
      upperLimitsObserved["METWeak_2ltau_offZ_noOSSF_cut_300"] = 0.142;
      upperLimitsObserved["METWeak_3l_onZ_cut_0"] = 24.769;
      upperLimitsObserved["METWeak_3l_onZ_cut_100"] = 0.690;
      upperLimitsObserved["METWeak_3l_onZ_cut_200"] = 0.198;
      upperLimitsObserved["METWeak_3l_onZ_cut_300"] = 0.138;
      upperLimitsObserved["METWeak_2ltau_onZ_cut_0"] = 194.360;
      upperLimitsObserved["METWeak_2ltau_onZ_cut_100"] = 0.287;
      upperLimitsObserved["METWeak_2ltau_onZ_cut_200"] = 0.144;
      upperLimitsObserved["METWeak_2ltau_onZ_cut_300"] = 0.130;
      upperLimitsObserved["Meff_3l_offZ_OSSF_cut_0"] = 2.435;
      upperLimitsObserved["Meff_3l_offZ_OSSF_cut_600"] = 0.487;
      upperLimitsObserved["Meff_3l_offZ_OSSF_cut_1000"] = 0.156;
      upperLimitsObserved["Meff_3l_offZ_OSSF_cut_1500"] = 0.140;
      upperLimitsObserved["Meff_2ltau_offZ_OSSF_cut_0"] = 13.901;
      upperLimitsObserved["Meff_2ltau_offZ_OSSF_cut_600"] = 0.687;
      upperLimitsObserved["Meff_2ltau_offZ_OSSF_cut_1000"] = 0.224;
      upperLimitsObserved["Meff_2ltau_offZ_OSSF_cut_1500"] = 0.155;
      upperLimitsObserved["Meff_3l_offZ_noOSSF_cut_0"] = 1.054;
      upperLimitsObserved["Meff_3l_offZ_noOSSF_cut_600"] = 0.249;
      upperLimitsObserved["Meff_3l_offZ_noOSSF_cut_1000"] = 0.194;
      upperLimitsObserved["Meff_3l_offZ_noOSSF_cut_1500"] = 0.145;
      upperLimitsObserved["Meff_2ltau_offZ_noOSSF_cut_0"] = 4.276;
      upperLimitsObserved["Meff_2ltau_offZ_noOSSF_cut_600"] = 0.772;
      upperLimitsObserved["Meff_2ltau_offZ_noOSSF_cut_1000"] = 0.218;
      upperLimitsObserved["Meff_2ltau_offZ_noOSSF_cut_1500"] = 0.204;
      upperLimitsObserved["Meff_3l_onZ_cut_0"] = 29.804;
      upperLimitsObserved["Meff_3l_onZ_cut_600"] = 2.933;
      upperLimitsObserved["Meff_3l_onZ_cut_1000"] = 0.912;
      upperLimitsObserved["Meff_3l_onZ_cut_1500"] = 0.225;
      upperLimitsObserved["Meff_2ltau_onZ_cut_0"] = 205.091;
      upperLimitsObserved["Meff_2ltau_onZ_cut_600"] = 1.486;
      upperLimitsObserved["Meff_2ltau_onZ_cut_1000"] = 0.641;
      upperLimitsObserved["Meff_2ltau_onZ_cut_1500"] = 0.204;
      upperLimitsObserved["MeffStrong_3l_offZ_OSSF_cut_0"] = 0.479;
      upperLimitsObserved["MeffStrong_3l_offZ_OSSF_cut_600"] = 0.353;
      upperLimitsObserved["MeffStrong_3l_offZ_OSSF_cut_1200"] = 0.187;
      upperLimitsObserved["MeffStrong_2ltau_offZ_OSSF_cut_0"] = 0.617;
      upperLimitsObserved["MeffStrong_2ltau_offZ_OSSF_cut_600"] = 0.320;
      upperLimitsObserved["MeffStrong_2ltau_offZ_OSSF_cut_1200"] = 0.281;
      upperLimitsObserved["MeffStrong_3l_offZ_noOSSF_cut_0"] = 0.408;
      upperLimitsObserved["MeffStrong_3l_offZ_noOSSF_cut_600"] = 0.240;
      upperLimitsObserved["MeffStrong_3l_offZ_noOSSF_cut_1200"] = 0.150;
      upperLimitsObserved["MeffStrong_2ltau_offZ_noOSSF_cut_0"] = 0.774;
      upperLimitsObserved["MeffStrong_2ltau_offZ_noOSSF_cut_600"] = 0.417;
      upperLimitsObserved["MeffStrong_2ltau_offZ_noOSSF_cut_1200"] = 0.266;
      upperLimitsObserved["MeffStrong_3l_onZ_cut_0"] = 1.208;
      upperLimitsObserved["MeffStrong_3l_onZ_cut_600"] = 0.837;
      upperLimitsObserved["MeffStrong_3l_onZ_cut_1200"] = 0.269;
      upperLimitsObserved["MeffStrong_2ltau_onZ_cut_0"] = 0.605;
      upperLimitsObserved["MeffStrong_2ltau_onZ_cut_600"] = 0.420;
      upperLimitsObserved["MeffStrong_2ltau_onZ_cut_1200"] = 0.141;
      upperLimitsObserved["MeffMt_3l_onZ_cut_0"] = 1.832;
      upperLimitsObserved["MeffMt_3l_onZ_cut_600"] = 0.862;
      upperLimitsObserved["MeffMt_3l_onZ_cut_1200"] = 0.222;
      upperLimitsObserved["MeffMt_2ltau_onZ_cut_0"] = 1.309;
      upperLimitsObserved["MeffMt_2ltau_onZ_cut_600"] = 0.481;
      upperLimitsObserved["MeffMt_2ltau_onZ_cut_1200"] = 0.146;
      upperLimitsObserved["MinPt_3l_offZ_OSSF_cut_0"] = 2.435;
      upperLimitsObserved["MinPt_3l_offZ_OSSF_cut_50"] = 0.500;
      upperLimitsObserved["MinPt_3l_offZ_OSSF_cut_100"] = 0.203;
      upperLimitsObserved["MinPt_3l_offZ_OSSF_cut_150"] = 0.128;
      upperLimitsObserved["MinPt_2ltau_offZ_OSSF_cut_0"] = 13.901;
      upperLimitsObserved["MinPt_2ltau_offZ_OSSF_cut_50"] = 0.859;
      upperLimitsObserved["MinPt_2ltau_offZ_OSSF_cut_100"] = 0.158;
      upperLimitsObserved["MinPt_2ltau_offZ_OSSF_cut_150"] = 0.155;
      upperLimitsObserved["MinPt_3l_offZ_noOSSF_cut_0"] = 1.054;
      upperLimitsObserved["MinPt_3l_offZ_noOSSF_cut_50"] = 0.295;
      upperLimitsObserved["MinPt_3l_offZ_noOSSF_cut_100"] = 0.148;
      upperLimitsObserved["MinPt_3l_offZ_noOSSF_cut_150"] = 0.137;
      upperLimitsObserved["MinPt_2ltau_offZ_noOSSF_cut_0"] = 4.276;
      upperLimitsObserved["MinPt_2ltau_offZ_noOSSF_cut_50"] = 0.314;
      upperLimitsObserved["MinPt_2ltau_offZ_noOSSF_cut_100"] = 0.134;
      upperLimitsObserved["MinPt_2ltau_offZ_noOSSF_cut_150"] = 0.140;
      upperLimitsObserved["MinPt_3l_onZ_cut_0"] = 29.804;
      upperLimitsObserved["MinPt_3l_onZ_cut_50"] = 1.767;
      upperLimitsObserved["MinPt_3l_onZ_cut_100"] = 0.690;
      upperLimitsObserved["MinPt_3l_onZ_cut_150"] = 0.301;
      upperLimitsObserved["MinPt_2ltau_onZ_cut_0"] = 205.091;
      upperLimitsObserved["MinPt_2ltau_onZ_cut_50"] = 1.050;
      upperLimitsObserved["MinPt_2ltau_onZ_cut_100"] = 0.155;
      upperLimitsObserved["MinPt_2ltau_onZ_cut_150"] = 0.146;
      upperLimitsObserved["nbtag_3l_offZ_OSSF_cut_0"] = 2.435;
      upperLimitsObserved["nbtag_3l_offZ_OSSF_cut_1"] = 0.865;
      upperLimitsObserved["nbtag_3l_offZ_OSSF_cut_2"] = 0.474;
      upperLimitsObserved["nbtag_2ltau_offZ_OSSF_cut_0"] = 13.901;
      upperLimitsObserved["nbtag_2ltau_offZ_OSSF_cut_1"] = 1.566;
      upperLimitsObserved["nbtag_2ltau_offZ_OSSF_cut_2"] = 0.426;
      upperLimitsObserved["nbtag_3l_offZ_noOSSF_cut_0"] = 1.054;
      upperLimitsObserved["nbtag_3l_offZ_noOSSF_cut_1"] = 0.643;
      upperLimitsObserved["nbtag_3l_offZ_noOSSF_cut_2"] = 0.321;
      upperLimitsObserved["nbtag_2ltau_offZ_noOSSF_cut_0"] = 4.276;
      upperLimitsObserved["nbtag_2ltau_offZ_noOSSF_cut_1"] = 2.435;
      upperLimitsObserved["nbtag_2ltau_offZ_noOSSF_cut_2"] = 1.073;
      upperLimitsObserved["nbtag_3l_onZ_cut_0"] = 29.804;
      upperLimitsObserved["nbtag_3l_onZ_cut_1"] = 3.908;
      upperLimitsObserved["nbtag_3l_onZ_cut_2"] = 0.704;
      upperLimitsObserved["nbtag_2ltau_onZ_cut_0"] = 205.091;
      upperLimitsObserved["nbtag_2ltau_onZ_cut_1"] = 9.377;
      upperLimitsObserved["nbtag_2ltau_onZ_cut_2"] = 0.657;
      upperLimitsExpected["HTlep_3l_offZ_OSSF_cut_0"] = 2.893;
      upperLimitsExpected["HTlep_3l_offZ_OSSF_cut_200"] = 1.175;
      upperLimitsExpected["HTlep_3l_offZ_OSSF_cut_500"] = 0.265;
      upperLimitsExpected["HTlep_3l_offZ_OSSF_cut_800"] = 0.155;
      upperLimitsExpected["HTlep_2ltau_offZ_OSSF_cut_0"] = 14.293;
      upperLimitsExpected["HTlep_2ltau_offZ_OSSF_cut_200"] = 1.803;
      upperLimitsExpected["HTlep_2ltau_offZ_OSSF_cut_500"] = 0.159;
      upperLimitsExpected["HTlep_2ltau_offZ_OSSF_cut_800"] = 0.155;
      upperLimitsExpected["HTlep_3l_offZ_noOSSF_cut_0"] = 0.836;
      upperLimitsExpected["HTlep_3l_offZ_noOSSF_cut_200"] = 0.340;
      upperLimitsExpected["HTlep_3l_offZ_noOSSF_cut_500"] = 0.218;
      upperLimitsExpected["HTlep_3l_offZ_noOSSF_cut_800"] = 0.140;
      upperLimitsExpected["HTlep_2ltau_offZ_noOSSF_cut_0"] = 4.132;
      upperLimitsExpected["HTlep_2ltau_offZ_noOSSF_cut_200"] = 0.599;
      upperLimitsExpected["HTlep_2ltau_offZ_noOSSF_cut_500"] = 0.146;
      upperLimitsExpected["HTlep_2ltau_offZ_noOSSF_cut_800"] = 0.148;
      upperLimitsExpected["HTlep_3l_onZ_cut_0"] = 32.181;
      upperLimitsExpected["HTlep_3l_onZ_cut_200"] = 4.879;
      upperLimitsExpected["HTlep_3l_onZ_cut_500"] = 0.473;
      upperLimitsExpected["HTlep_3l_onZ_cut_800"] = 0.266;
      upperLimitsExpected["HTlep_2ltau_onZ_cut_0"] = 217.801;
      upperLimitsExpected["HTlep_2ltau_onZ_cut_200"] = 3.676;
      upperLimitsExpected["HTlep_2ltau_onZ_cut_500"] = 0.235;
      upperLimitsExpected["HTlep_2ltau_onZ_cut_800"] = 0.150;
      upperLimitsExpected["METStrong_3l_offZ_OSSF_cut_0"] = 1.196;
      upperLimitsExpected["METStrong_3l_offZ_OSSF_cut_100"] = 0.423;
      upperLimitsExpected["METStrong_3l_offZ_OSSF_cut_200"] = 0.208;
      upperLimitsExpected["METStrong_3l_offZ_OSSF_cut_300"] = 0.158;
      upperLimitsExpected["METStrong_2ltau_offZ_OSSF_cut_0"] = 2.158;
      upperLimitsExpected["METStrong_2ltau_offZ_OSSF_cut_100"] = 0.461;
      upperLimitsExpected["METStrong_2ltau_offZ_OSSF_cut_200"] = 0.186;
      upperLimitsExpected["METStrong_2ltau_offZ_OSSF_cut_300"] = 0.138;
      upperLimitsExpected["METStrong_3l_offZ_noOSSF_cut_0"] = 0.495;
      upperLimitsExpected["METStrong_3l_offZ_noOSSF_cut_100"] = 0.284;
      upperLimitsExpected["METStrong_3l_offZ_noOSSF_cut_200"] = 0.150;
      upperLimitsExpected["METStrong_3l_offZ_noOSSF_cut_300"] = 0.146;
      upperLimitsExpected["METStrong_2ltau_offZ_noOSSF_cut_0"] = 1.967;
      upperLimitsExpected["METStrong_2ltau_offZ_noOSSF_cut_100"] = 0.732;
      upperLimitsExpected["METStrong_2ltau_offZ_noOSSF_cut_200"] = 0.225;
      upperLimitsExpected["METStrong_2ltau_offZ_noOSSF_cut_300"] = 0.147;
      upperLimitsExpected["METStrong_3l_onZ_cut_0"] = 7.157;
      upperLimitsExpected["METStrong_3l_onZ_cut_100"] = 1.342;
      upperLimitsExpected["METStrong_3l_onZ_cut_200"] = 0.508;
      upperLimitsExpected["METStrong_3l_onZ_cut_300"] = 0.228;
      upperLimitsExpected["METStrong_2ltau_onZ_cut_0"] = 12.441;
      upperLimitsExpected["METStrong_2ltau_onZ_cut_100"] = 0.534;
      upperLimitsExpected["METStrong_2ltau_onZ_cut_200"] = 0.243;
      upperLimitsExpected["METStrong_2ltau_onZ_cut_300"] = 0.218;
      upperLimitsExpected["METWeak_3l_offZ_OSSF_cut_0"] = 2.199;
      upperLimitsExpected["METWeak_3l_offZ_OSSF_cut_100"] = 0.391;
      upperLimitsExpected["METWeak_3l_offZ_OSSF_cut_200"] = 0.177;
      upperLimitsExpected["METWeak_3l_offZ_OSSF_cut_300"] = 0.144;
      upperLimitsExpected["METWeak_2ltau_offZ_OSSF_cut_0"] = 12.431;
      upperLimitsExpected["METWeak_2ltau_offZ_OSSF_cut_100"] = 0.358;
      upperLimitsExpected["METWeak_2ltau_offZ_OSSF_cut_200"] = 0.150;
      upperLimitsExpected["METWeak_2ltau_offZ_OSSF_cut_300"] = 0.135;
      upperLimitsExpected["METWeak_3l_offZ_noOSSF_cut_0"] = 0.577;
      upperLimitsExpected["METWeak_3l_offZ_noOSSF_cut_100"] = 0.214;
      upperLimitsExpected["METWeak_3l_offZ_noOSSF_cut_200"] = 0.155;
      upperLimitsExpected["METWeak_3l_offZ_noOSSF_cut_300"] = 0.140;
      upperLimitsExpected["METWeak_2ltau_offZ_noOSSF_cut_0"] = 2.474;
      upperLimitsExpected["METWeak_2ltau_offZ_noOSSF_cut_100"] = 0.382;
      upperLimitsExpected["METWeak_2ltau_offZ_noOSSF_cut_200"] = 0.144;
      upperLimitsExpected["METWeak_2ltau_offZ_noOSSF_cut_300"] = 0.146;
      upperLimitsExpected["METWeak_3l_onZ_cut_0"] = 26.305;
      upperLimitsExpected["METWeak_3l_onZ_cut_100"] = 1.227;
      upperLimitsExpected["METWeak_3l_onZ_cut_200"] = 0.311;
      upperLimitsExpected["METWeak_3l_onZ_cut_300"] = 0.188;
      upperLimitsExpected["METWeak_2ltau_onZ_cut_0"] = 205.198;
      upperLimitsExpected["METWeak_2ltau_onZ_cut_100"] = 0.399;
      upperLimitsExpected["METWeak_2ltau_onZ_cut_200"] = 0.166;
      upperLimitsExpected["METWeak_2ltau_onZ_cut_300"] = 0.140;
      upperLimitsExpected["Meff_3l_offZ_OSSF_cut_0"] = 2.893;
      upperLimitsExpected["Meff_3l_offZ_OSSF_cut_600"] = 0.649;
      upperLimitsExpected["Meff_3l_offZ_OSSF_cut_1000"] = 0.252;
      upperLimitsExpected["Meff_3l_offZ_OSSF_cut_1500"] = 0.150;
      upperLimitsExpected["Meff_2ltau_offZ_OSSF_cut_0"] = 14.293;
      upperLimitsExpected["Meff_2ltau_offZ_OSSF_cut_600"] = 0.657;
      upperLimitsExpected["Meff_2ltau_offZ_OSSF_cut_1000"] = 0.226;
      upperLimitsExpected["Meff_2ltau_offZ_OSSF_cut_1500"] = 0.154;
      upperLimitsExpected["Meff_3l_offZ_noOSSF_cut_0"] = 0.836;
      upperLimitsExpected["Meff_3l_offZ_noOSSF_cut_600"] = 0.265;
      upperLimitsExpected["Meff_3l_offZ_noOSSF_cut_1000"] = 0.176;
      upperLimitsExpected["Meff_3l_offZ_noOSSF_cut_1500"] = 0.146;
      upperLimitsExpected["Meff_2ltau_offZ_noOSSF_cut_0"] = 4.132;
      upperLimitsExpected["Meff_2ltau_offZ_noOSSF_cut_600"] = 0.678;
      upperLimitsExpected["Meff_2ltau_offZ_noOSSF_cut_1000"] = 0.243;
      upperLimitsExpected["Meff_2ltau_offZ_noOSSF_cut_1500"] = 0.184;
      upperLimitsExpected["Meff_3l_onZ_cut_0"] = 32.181;
      upperLimitsExpected["Meff_3l_onZ_cut_600"] = 3.219;
      upperLimitsExpected["Meff_3l_onZ_cut_1000"] = 0.905;
      upperLimitsExpected["Meff_3l_onZ_cut_1500"] = 0.261;
      upperLimitsExpected["Meff_2ltau_onZ_cut_0"] = 217.801;
      upperLimitsExpected["Meff_2ltau_onZ_cut_600"] = 1.680;
      upperLimitsExpected["Meff_2ltau_onZ_cut_1000"] = 0.375;
      upperLimitsExpected["Meff_2ltau_onZ_cut_1500"] = 0.178;
      upperLimitsExpected["MeffStrong_3l_offZ_OSSF_cut_0"] = 0.571;
      upperLimitsExpected["MeffStrong_3l_offZ_OSSF_cut_600"] = 0.386;
      upperLimitsExpected["MeffStrong_3l_offZ_OSSF_cut_1200"] = 0.177;
      upperLimitsExpected["MeffStrong_2ltau_offZ_OSSF_cut_0"] = 0.605;
      upperLimitsExpected["MeffStrong_2ltau_offZ_OSSF_cut_600"] = 0.335;
      upperLimitsExpected["MeffStrong_2ltau_offZ_OSSF_cut_1200"] = 0.249;
      upperLimitsExpected["MeffStrong_3l_offZ_noOSSF_cut_0"] = 0.373;
      upperLimitsExpected["MeffStrong_3l_offZ_noOSSF_cut_600"] = 0.223;
      upperLimitsExpected["MeffStrong_3l_offZ_noOSSF_cut_1200"] = 0.150;
      upperLimitsExpected["MeffStrong_2ltau_offZ_noOSSF_cut_0"] = 0.873;
      upperLimitsExpected["MeffStrong_2ltau_offZ_noOSSF_cut_600"] = 0.428;
      upperLimitsExpected["MeffStrong_2ltau_offZ_noOSSF_cut_1200"] = 0.210;
      upperLimitsExpected["MeffStrong_3l_onZ_cut_0"] = 2.034;
      upperLimitsExpected["MeffStrong_3l_onZ_cut_600"] = 1.093;
      upperLimitsExpected["MeffStrong_3l_onZ_cut_1200"] = 0.293;
      upperLimitsExpected["MeffStrong_2ltau_onZ_cut_0"] = 0.690;
      upperLimitsExpected["MeffStrong_2ltau_onZ_cut_600"] = 0.392;
      upperLimitsExpected["MeffStrong_2ltau_onZ_cut_1200"] = 0.156;
      upperLimitsExpected["MeffMt_3l_onZ_cut_0"] = 2.483;
      upperLimitsExpected["MeffMt_3l_onZ_cut_600"] = 0.845;
      upperLimitsExpected["MeffMt_3l_onZ_cut_1200"] = 0.255;
      upperLimitsExpected["MeffMt_2ltau_onZ_cut_0"] = 1.448;
      upperLimitsExpected["MeffMt_2ltau_onZ_cut_600"] = 0.391;
      upperLimitsExpected["MeffMt_2ltau_onZ_cut_1200"] = 0.146;
      upperLimitsExpected["MinPt_3l_offZ_OSSF_cut_0"] = 2.893;
      upperLimitsExpected["MinPt_3l_offZ_OSSF_cut_50"] = 0.703;
      upperLimitsExpected["MinPt_3l_offZ_OSSF_cut_100"] = 0.207;
      upperLimitsExpected["MinPt_3l_offZ_OSSF_cut_150"] = 0.143;
      upperLimitsExpected["MinPt_2ltau_offZ_OSSF_cut_0"] = 14.293;
      upperLimitsExpected["MinPt_2ltau_offZ_OSSF_cut_50"] = 0.705;
      upperLimitsExpected["MinPt_2ltau_offZ_OSSF_cut_100"] = 0.149;
      upperLimitsExpected["MinPt_2ltau_offZ_OSSF_cut_150"] = 0.155;
      upperLimitsExpected["MinPt_3l_offZ_noOSSF_cut_0"] = 0.836;
      upperLimitsExpected["MinPt_3l_offZ_noOSSF_cut_50"] = 0.249;
      upperLimitsExpected["MinPt_3l_offZ_noOSSF_cut_100"] = 0.135;
      upperLimitsExpected["MinPt_3l_offZ_noOSSF_cut_150"] = 0.136;
      upperLimitsExpected["MinPt_2ltau_offZ_noOSSF_cut_0"] = 4.132;
      upperLimitsExpected["MinPt_2ltau_offZ_noOSSF_cut_50"] = 0.339;
      upperLimitsExpected["MinPt_2ltau_offZ_noOSSF_cut_100"] = 0.149;
      upperLimitsExpected["MinPt_2ltau_offZ_noOSSF_cut_150"] = 0.145;
      upperLimitsExpected["MinPt_3l_onZ_cut_0"] = 32.181;
      upperLimitsExpected["MinPt_3l_onZ_cut_50"] = 2.260;
      upperLimitsExpected["MinPt_3l_onZ_cut_100"] = 0.438;
      upperLimitsExpected["MinPt_3l_onZ_cut_150"] = 0.305;
      upperLimitsExpected["MinPt_2ltau_onZ_cut_0"] = 217.801;
      upperLimitsExpected["MinPt_2ltau_onZ_cut_50"] = 1.335;
      upperLimitsExpected["MinPt_2ltau_onZ_cut_100"] = 0.162;
      upperLimitsExpected["MinPt_2ltau_onZ_cut_150"] = 0.149;
      upperLimitsExpected["nbtag_3l_offZ_OSSF_cut_0"] = 2.893;
      upperLimitsExpected["nbtag_3l_offZ_OSSF_cut_1"] = 0.923;
      upperLimitsExpected["nbtag_3l_offZ_OSSF_cut_2"] = 0.452;
      upperLimitsExpected["nbtag_2ltau_offZ_OSSF_cut_0"] = 14.293;
      upperLimitsExpected["nbtag_2ltau_offZ_OSSF_cut_1"] = 1.774;
      upperLimitsExpected["nbtag_2ltau_offZ_OSSF_cut_2"] = 0.549;
      upperLimitsExpected["nbtag_3l_offZ_noOSSF_cut_0"] = 0.836;
      upperLimitsExpected["nbtag_3l_offZ_noOSSF_cut_1"] = 0.594;
      upperLimitsExpected["nbtag_3l_offZ_noOSSF_cut_2"] = 0.298;
      upperLimitsExpected["nbtag_2ltau_offZ_noOSSF_cut_0"] = 4.132;
      upperLimitsExpected["nbtag_2ltau_offZ_noOSSF_cut_1"] = 2.358;
      upperLimitsExpected["nbtag_2ltau_offZ_noOSSF_cut_2"] = 0.958;
      upperLimitsExpected["nbtag_3l_onZ_cut_0"] = 32.181;
      upperLimitsExpected["nbtag_3l_onZ_cut_1"] = 3.868;
      upperLimitsExpected["nbtag_3l_onZ_cut_2"] = 0.887;
      upperLimitsExpected["nbtag_2ltau_onZ_cut_0"] = 217.801;
      upperLimitsExpected["nbtag_2ltau_onZ_cut_1"] = 9.397;
      upperLimitsExpected["nbtag_2ltau_onZ_cut_2"] = 0.787;



      if (observed) return upperLimitsObserved[signal_region];
      else          return upperLimitsExpected[signal_region];
    }


    /// Function checking if there is an OSSF lepton pair or a combination of 3 leptons with an invariant mass close to the Z mass
    int isonZ (const Particles& particles) {
      int onZ = 0;
      double best_mass_2 = 999.;
      double best_mass_3 = 999.;

      // Loop over all 2 particle combinations to find invariant mass of OSSF pair closest to Z mass
      for (const Particle& p1 : particles)  {
        for (const Particle& p2 : particles)  {
          double mass_difference_2_old = fabs(91.0 - best_mass_2);
          double mass_difference_2_new = fabs(91.0 - (p1.momentum() + p2.momentum()).mass()/GeV);

          // If particle combination is OSSF pair calculate mass difference to Z mass
          if ((p1.pid()*p2.pid() == -121 || p1.pid()*p2.pid() == -169))  {

            // Get invariant mass closest to Z mass
            if (mass_difference_2_new < mass_difference_2_old)
              best_mass_2 = (p1.momentum() + p2.momentum()).mass()/GeV;
          // In case there is an OSSF pair take also 3rd lepton into account (e.g. from FSR and photon to electron conversion)
            for (const Particle& p3 : particles  )  {
            double mass_difference_3_old = fabs(91.0 - best_mass_3);
            double mass_difference_3_new = fabs(91.0 - (p1.momentum() + p2.momentum() + p3.momentum()).mass()/GeV);
            if (mass_difference_3_new < mass_difference_3_old)
              best_mass_3 = (p1.momentum() + p2.momentum() + p3.momentum()).mass()/GeV;
            }
          }
        }
      }

      // Pick the minimum invariant mass of the best OSSF pair combination and the best 3 lepton combination
      double best_mass = min(best_mass_2,best_mass_3);
      // if this mass is in a 20 GeV window around the Z mass, the event is classified as onZ
      if ( fabs(91.0 - best_mass) < 20. ) onZ = 1;
      return onZ;
    }

    /// function checking if two leptons are an OSSF lepton pair and giving out the invariant mass (0 if no OSSF pair)
    double isOSSF_mass (const Particle& p1, const Particle& p2) {
      double inv_mass = 0.;
      // Is particle combination OSSF pair?
      if ((p1.pid()*p2.pid() == -121 || p1.pid()*p2.pid() == -169))  {
        // Get invariant mass
        inv_mass = (p1.momentum() + p2.momentum()).mass()/GeV;
      }
      return inv_mass;
    }

    /// Function checking if there is an OSSF lepton pair
    bool isOSSF (const Particles& particles)  {
      for (size_t i1=0 ; i1 < 3 ; i1 ++)  {
        for (size_t i2 = i1+1 ; i2 < 3 ; i2 ++)  {
          if ((particles[i1].pid()*particles[i2].pid() == -121 || particles[i1].pid()*particles[i2].pid() == -169))  {
            return true;
          }
        }
      }
      return false;
    }

    //@}

  private:

    /// Histograms
    //@{
    Histo1DPtr _h_HTlep_all, _h_HTjets_all, _h_MET_all, _h_Meff_all, _h_min_pT_all, _h_mT_all;
    Histo1DPtr _h_pt_1_3l, _h_pt_2_3l, _h_pt_3_3l, _h_pt_1_2ltau, _h_pt_2_2ltau, _h_pt_3_2ltau;
    Histo1DPtr _h_e_n, _h_mu_n, _h_tau_n;
    Histo1DPtr _h_excluded;
    //@}

    /// Fiducial efficiencies to model the effects of the ATLAS detector
    bool _use_fiducial_lepton_efficiency;

    /// List of signal regions and event counts per signal region
    vector<string> _signal_regions;
    map<string, double> _eventCountsPerSR;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1327229);

}
