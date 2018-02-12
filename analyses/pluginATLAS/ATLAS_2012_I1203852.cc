// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/InvMassFinalState.hh"

namespace Rivet {

  /// Generic Z candidate
  struct Zstate : public ParticlePair {
    Zstate() { }
    Zstate(ParticlePair _particlepair) : ParticlePair(_particlepair) { }
    FourMomentum mom() const { return first.momentum() + second.momentum(); }
    operator FourMomentum() const { return mom(); }
    static bool cmppT(const Zstate& lx, const Zstate& rx) { return lx.mom().pT() < rx.mom().pT(); }
  };



  /// @name ZZ analysis
  class ATLAS_2012_I1203852 : public Analysis {
  public:

    /// Default constructor
    ATLAS_2012_I1203852()
      : Analysis("ATLAS_2012_I1203852")
    {    }

    void init() {

      // NB Missing ET is not required to be neutrinos
      FinalState fs(-5.0, 5.0, 0.0*GeV);

      // Final states to form Z bosons
      vids.push_back(make_pair(PID::ELECTRON, PID::POSITRON));
      vids.push_back(make_pair(PID::MUON, PID::ANTIMUON));

      IdentifiedFinalState Photon(fs);
      Photon.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState bare_EL(fs);
      bare_EL.acceptIdPair(PID::ELECTRON);

      IdentifiedFinalState bare_MU(fs);
      bare_MU.acceptIdPair(PID::MUON);


      // Selection 1: ZZ-> llll selection
      Cut etaranges_lep = Cuts::abseta < 3.16 && Cuts::pT > 7*GeV;

      DressedLeptons electron_sel4l(Photon, bare_EL, 0.1, etaranges_lep);
      declare(electron_sel4l, "ELECTRON_sel4l");
      DressedLeptons muon_sel4l(Photon, bare_MU, 0.1, etaranges_lep);
      declare(muon_sel4l, "MUON_sel4l");


      // Selection 2: ZZ-> llnunu selection
      Cut etaranges_lep2 = Cuts::abseta < 2.5 && Cuts::pT > 10*GeV;

      DressedLeptons electron_sel2l2nu(Photon, bare_EL, 0.1, etaranges_lep2);
      declare(electron_sel2l2nu, "ELECTRON_sel2l2nu");
      DressedLeptons muon_sel2l2nu(Photon, bare_MU, 0.1, etaranges_lep2);
      declare(muon_sel2l2nu, "MUON_sel2l2nu");


      /// Get all neutrinos. These will not be used to form jets.
      IdentifiedFinalState neutrino_fs(Cuts::abseta < 4.5);
      neutrino_fs.acceptNeutrinos();
      declare(neutrino_fs, "NEUTRINO_FS");

      // Calculate missing ET from the visible final state, not by requiring neutrinos
      addProjection(MissingMomentum(Cuts::abseta < 4.5), "MISSING");

      VetoedFinalState jetinput;
      jetinput.addVetoOnThisFinalState(bare_MU);
      jetinput.addVetoOnThisFinalState(neutrino_fs);

      FastJets jetpro(fs, FastJets::ANTIKT, 0.4);
      declare(jetpro, "jet");

      // Both ZZ on-shell histos
      _h_ZZ_xsect = bookHisto1D(1, 1, 1);
      _h_ZZ_ZpT   = bookHisto1D(3, 1, 1);
      _h_ZZ_phill = bookHisto1D(5, 1, 1);
      _h_ZZ_mZZ   = bookHisto1D(7, 1, 1);

      // One Z off-shell (ZZstar) histos
      _h_ZZs_xsect = bookHisto1D(1, 1, 2);

      // ZZ -> llnunu histos
      _h_ZZnunu_xsect = bookHisto1D(1, 1, 3);
      _h_ZZnunu_ZpT   = bookHisto1D(4, 1, 1);
      _h_ZZnunu_phill = bookHisto1D(6, 1, 1);
      _h_ZZnunu_mZZ   = bookHisto1D(8, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      const double weight = e.weight();

      ////////////////////////////////////////////////////////////////////
      // preselection of leptons for ZZ-> llll final state
      ////////////////////////////////////////////////////////////////////

      Particles leptons_sel4l;

      const vector<DressedLepton>& mu_sel4l = apply<DressedLeptons>(e, "MUON_sel4l").dressedLeptons();
      const vector<DressedLepton>& el_sel4l = apply<DressedLeptons>(e, "ELECTRON_sel4l").dressedLeptons();

      vector<DressedLepton> leptonsFS_sel4l;
      leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), mu_sel4l.begin(), mu_sel4l.end() );
      leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), el_sel4l.begin(), el_sel4l.end() );

      ////////////////////////////////////////////////////////////////////
      // OVERLAP removal dR(l,l)>0.2
      ////////////////////////////////////////////////////////////////////
      foreach ( const DressedLepton& l1, leptonsFS_sel4l) {
        bool isolated = true;
        foreach (DressedLepton& l2, leptonsFS_sel4l) {
          const double dR = deltaR(l1, l2);
          if (dR < 0.2 && l1 != l2) { isolated = false; break; }
        }
        if (isolated) leptons_sel4l.push_back(l1);
      }

      //////////////////////////////////////////////////////////////////
      // Exactly two opposite charged leptons
      //////////////////////////////////////////////////////////////////

      // calculate total 'flavour' charge
      double totalcharge = 0;
      foreach (Particle& l, leptons_sel4l) totalcharge += l.pid();

      // Analyze 4 lepton events
      if (leptons_sel4l.size() == 4 && totalcharge == 0  ) {
        Zstate Z1, Z2;

        // Identifies Z states from 4 lepton pairs
        identifyZstates(Z1, Z2,leptons_sel4l);

        ////////////////////////////////////////////////////////////////////////////
        // Z MASS WINDOW
        //  -ZZ: for both Z: 66<mZ<116 GeV
        //  -ZZ*: one Z on-shell: 66<mZ<116 GeV, one Z off-shell: mZ>20 GeV
        ///////////////////////////////////////////////////////////////////////////

        Zstate leadPtZ = std::max(Z1, Z2, Zstate::cmppT);

        double mZ1   = Z1.mom().mass();
        double mZ2   = Z2.mom().mass();
        double ZpT   = leadPtZ.mom().pT();
        double phill = fabs(deltaPhi(leadPtZ.first, leadPtZ.second));
        if (phill > M_PI) phill = 2*M_PI-phill;
        double mZZ   = (Z1.mom() + Z2.mom()).mass();

        if (mZ1 > 20*GeV && mZ2 > 20*GeV) {
          // ZZ* selection
          if (inRange(mZ1, 66*GeV, 116*GeV) || inRange(mZ2, 66*GeV, 116*GeV)) {
            _h_ZZs_xsect  -> fill(sqrtS()*GeV,  weight);
          }

          // ZZ selection
          if (inRange(mZ1, 66*GeV, 116*GeV) && inRange(mZ2, 66*GeV, 116*GeV)) {
            _h_ZZ_xsect  -> fill(sqrtS()*GeV,  weight);
            _h_ZZ_ZpT    -> fill(ZpT   , weight);
            _h_ZZ_phill  -> fill(phill , weight);
            _h_ZZ_mZZ    -> fill(mZZ   , weight);
          }
        }
      }

      ////////////////////////////////////////////////////////////////////
      /// preselection of leptons for ZZ-> llnunu final state
      ////////////////////////////////////////////////////////////////////

      Particles leptons_sel2l2nu; // output
      const vector<DressedLepton>& mu_sel2l2nu = apply<DressedLeptons>(e, "MUON_sel2l2nu").dressedLeptons();
      const vector<DressedLepton>& el_sel2l2nu = apply<DressedLeptons>(e, "ELECTRON_sel2l2nu").dressedLeptons();

      vector<DressedLepton> leptonsFS_sel2l2nu;
      leptonsFS_sel2l2nu.insert( leptonsFS_sel2l2nu.end(), mu_sel2l2nu.begin(), mu_sel2l2nu.end() );
      leptonsFS_sel2l2nu.insert( leptonsFS_sel2l2nu.end(), el_sel2l2nu.begin(), el_sel2l2nu.end() );

      // Lepton preselection for ZZ-> llnunu
      if ((mu_sel2l2nu.empty() || el_sel2l2nu.empty()) // cannot have opposite flavour
           && (leptonsFS_sel2l2nu.size() == 2) // exactly two leptons
           && (leptonsFS_sel2l2nu[0].charge() * leptonsFS_sel2l2nu[1].charge() < 1 ) // opposite charge
           && (deltaR(leptonsFS_sel2l2nu[0], leptonsFS_sel2l2nu[1]) > 0.3) // overlap removal
           && (leptonsFS_sel2l2nu[0].pT() > 20*GeV && leptonsFS_sel2l2nu[1].pT() > 20*GeV)) { // trigger requirement
        leptons_sel2l2nu.insert(leptons_sel2l2nu.end(), leptonsFS_sel2l2nu.begin(), leptonsFS_sel2l2nu.end());
      }
      if (leptons_sel2l2nu.empty()) vetoEvent; // no further analysis, fine to veto

      Particles leptons_sel2l2nu_jetveto;
      foreach (const DressedLepton& l, mu_sel2l2nu) leptons_sel2l2nu_jetveto.push_back(l.constituentLepton());
      foreach (const DressedLepton& l, el_sel2l2nu) leptons_sel2l2nu_jetveto.push_back(l.constituentLepton());
      double ptll = (leptons_sel2l2nu[0].momentum() + leptons_sel2l2nu[1].momentum()).pT();

      // Find Z1-> ll
      FinalState fs2(-3.2, 3.2);
      InvMassFinalState imfs(fs2, vids, 20*GeV, sqrtS());
      imfs.calc(leptons_sel2l2nu);
      if (imfs.particlePairs().size() != 1) vetoEvent;
      const ParticlePair& Z1constituents = imfs.particlePairs()[0];
      FourMomentum Z1 = Z1constituents.first.momentum() + Z1constituents.second.momentum();

      // Z to neutrinos candidate from missing ET
      const MissingMomentum & missmom = applyProjection<MissingMomentum>(e, "MISSING");
      const FourMomentum Z2 = missmom.missingMomentum(ZMASS);
      double met_Znunu = missmom.missingEt(); //Z2.pT();

      // mTZZ
      const double mT2_1st_term = add_quad(ZMASS, ptll) + add_quad(ZMASS, met_Znunu);
      const double mT2_2nd_term = Z1.px() + Z2.px();
      const double mT2_3rd_term = Z1.py() + Z2.py();
      const double mTZZ = sqrt(sqr(mT2_1st_term) - sqr(mT2_2nd_term) - sqr(mT2_3rd_term));

      if (!inRange(Z2.mass(), 66*GeV, 116*GeV)) vetoEvent;
      if (!inRange(Z1.mass(), 76*GeV, 106*GeV)) vetoEvent;

      /////////////////////////////////////////////////////////////
      // AXIAL MET < 75 GeV
      ////////////////////////////////////////////////////////////

      double dPhiZ1Z2 = fabs(deltaPhi(Z1, Z2));
      if (dPhiZ1Z2 > M_PI) dPhiZ1Z2 = 2*M_PI - dPhiZ1Z2;
      const double axialEtmiss = -Z2.pT()*cos(dPhiZ1Z2);
      if (axialEtmiss < 75*GeV) vetoEvent;

      const double ZpT   = Z1.pT();
      double phill = fabs(deltaPhi(Z1constituents.first, Z1constituents.second));
      if (phill > M_PI) phill = 2*M_PI - phill;


      ////////////////////////////////////////////////////////////////////////////
      // JETS
      //    -"j": found by "jetpro" projection && pT() > 25 GeV && |eta| < 4.5
      //    -"goodjets": "j"  && dR(electron/muon,jet) > 0.3
      //
      // JETVETO: veto all events with at least one good jet
      ///////////////////////////////////////////////////////////////////////////
      vector<Jet> good_jets;
      foreach (const Jet& j, apply<FastJets>(e, "jet").jetsByPt(25)) {
        if (j.abseta() > 4.5) continue;
        bool isLepton = 0;
        foreach (const Particle& l, leptons_sel2l2nu_jetveto) {
          const double dR = deltaR(l.momentum(), j.momentum());
          if (dR < 0.3) { isLepton = true; break; }
        }
        if (!isLepton) good_jets.push_back(j);
      }
      size_t n_sel_jets = good_jets.size();
      if (n_sel_jets != 0) vetoEvent;


      /////////////////////////////////////////////////////////////
      // Fractional MET and lepton pair difference: "RatioMet"< 0.4
      ////////////////////////////////////////////////////////////
      double ratioMet = fabs(Z2.pT() - Z1.pT()) / Z1.pT();
      if (ratioMet  > 0.4 ) vetoEvent;


      // End of ZZllnunu selection: now fill histograms
      _h_ZZnunu_xsect->fill(sqrtS()/GeV, weight);
      _h_ZZnunu_ZpT  ->fill(ZpT, weight);
      _h_ZZnunu_phill->fill(phill, weight);
      _h_ZZnunu_mZZ  ->fill(mTZZ, weight);

    }


    /// Finalize
    void finalize() {
      const double norm = crossSection()/sumOfWeights()/femtobarn;

      scale(_h_ZZ_xsect, norm);
      normalize(_h_ZZ_ZpT);
      normalize(_h_ZZ_phill);
      normalize(_h_ZZ_mZZ);

      scale(_h_ZZs_xsect, norm);
      scale(_h_ZZnunu_xsect, norm);
      normalize(_h_ZZnunu_ZpT);
      normalize(_h_ZZnunu_phill);
      normalize(_h_ZZnunu_mZZ);
    }


  private:

    void identifyZstates(Zstate& Z1, Zstate& Z2, const Particles& leptons_sel4l);
    Histo1DPtr _h_ZZ_xsect, _h_ZZ_ZpT, _h_ZZ_phill, _h_ZZ_mZZ;
    Histo1DPtr _h_ZZs_xsect;
    Histo1DPtr _h_ZZnunu_xsect, _h_ZZnunu_ZpT, _h_ZZnunu_phill, _h_ZZnunu_mZZ;
    vector< pair<PdgId,PdgId> > vids;
    const double ZMASS = 91.1876; // GeV


  };


  /// 4l to ZZ assignment -- algorithm
  void ATLAS_2012_I1203852::identifyZstates(Zstate& Z1, Zstate& Z2, const Particles& leptons_sel4l) {

    /////////////////////////////////////////////////////////////////////////////
    /// ZZ->4l pairing
    /// - Exactly two same flavour opposite charged leptons
    /// - Ambiguities in pairing are resolved by choosing the combination
    ///     that results in the smaller value of the sum |mll - mZ| for the two pairs
    /////////////////////////////////////////////////////////////////////////////

    Particles part_pos_el, part_neg_el, part_pos_mu, part_neg_mu;
    foreach (const Particle& l , leptons_sel4l) {
      if (l.abspid() == PID::ELECTRON) {
        if (l.pid() < 0) part_neg_el.push_back(l);
        if (l.pid() > 0) part_pos_el.push_back(l);
      }
      else if (l.abspid() == PID::MUON) {
        if (l.pid() < 0) part_neg_mu.push_back(l);
        if (l.pid() > 0) part_pos_mu.push_back(l);
      }
    }

    // ee/mm channel
    if ( part_neg_el.size() == 2 || part_neg_mu.size() == 2) {

      Zstate Zcand_1, Zcand_2, Zcand_3, Zcand_4;
      if (part_neg_el.size() == 2) { // ee
        Zcand_1 = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[0] ) );
        Zcand_2 = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[1] ) );
        Zcand_3 = Zstate( ParticlePair( part_neg_el[1],  part_pos_el[0] ) );
        Zcand_4 = Zstate( ParticlePair( part_neg_el[1],  part_pos_el[1] ) );
      } else { // mumu
        Zcand_1 = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[0] ) );
        Zcand_2 = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[1] ) );
        Zcand_3 = Zstate( ParticlePair( part_neg_mu[1],  part_pos_mu[0] ) );
        Zcand_4 = Zstate( ParticlePair( part_neg_mu[1],  part_pos_mu[1] ) );
      }

      // We can have the following pairs: (Z1 + Z4) || (Z2 + Z3)
      double minValue_1, minValue_2;
      minValue_1 = fabs( Zcand_1.mom().mass() - ZMASS ) + fabs( Zcand_4.mom().mass() - ZMASS);
      minValue_2 = fabs( Zcand_2.mom().mass() - ZMASS ) + fabs( Zcand_3.mom().mass() - ZMASS);
      if (minValue_1 < minValue_2 ) {
        Z1 = Zcand_1;
        Z2 = Zcand_4;
      } else {
        Z1 = Zcand_2;
        Z2 = Zcand_3;
      }

    // emu channel
    } else if (part_neg_mu.size() == 1 && part_neg_el.size() == 1) {
      Z1 = Zstate ( ParticlePair (part_neg_mu[0],  part_pos_mu[0] ) );
      Z2 = Zstate ( ParticlePair (part_neg_el[0],  part_pos_el[0] ) );
    }

  }


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1203852);

}
