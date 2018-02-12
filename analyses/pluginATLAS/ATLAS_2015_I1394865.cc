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


  /// Inclusive 4-lepton lineshape
  class ATLAS_2015_I1394865 : public Analysis {
  public:

    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1394865);


    void init() {
      FinalState fs(Cuts::abseta < 5.0);

      IdentifiedFinalState photon(fs, PID::PHOTON);
      IdentifiedFinalState bare_EL(fs, {PID::ELECTRON, -PID::ELECTRON});
      IdentifiedFinalState bare_MU(fs, {PID::MUON, -PID::MUON});

      // Selection 1: ZZ-> llll selection
      Cut etaranges_el = Cuts::abseta < 2.5 && Cuts::pT > 7*GeV;
      Cut etaranges_mu = Cuts::abseta < 2.7 && Cuts::pT > 6*GeV;

      DressedLeptons electron_sel4l(photon, bare_EL, 0.1, etaranges_el);
      declare(electron_sel4l, "ELECTRON_sel4l");
      DressedLeptons muon_sel4l(photon, bare_MU, 0.1, etaranges_mu);
      declare(muon_sel4l, "MUON_sel4l");


      // Both ZZ on-shell histos
      _h_ZZ_mZZ  = bookHisto1D(1, 1, 1);
      _h_ZZ_pTZZ = bookHisto1D(2, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      const double weight = e.weight();

      ////////////////////////////////////////////////////////////////////
      // Preselection of leptons for ZZ-> llll final state
      ////////////////////////////////////////////////////////////////////

      Particles leptons_sel4l;
      const vector<DressedLepton>& mu_sel4l = apply<DressedLeptons>(e, "MUON_sel4l").dressedLeptons();
      const vector<DressedLepton>& el_sel4l = apply<DressedLeptons>(e, "ELECTRON_sel4l").dressedLeptons();
      const vector<DressedLepton> leptonsFS_sel4l = mu_sel4l + el_sel4l;
      // leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), mu_sel4l.begin(), mu_sel4l.end() );
      // leptonsFS_sel4l.insert( leptonsFS_sel4l.end(), el_sel4l.begin(), el_sel4l.end() );

      // mu: pT > 6 GeV, eta < 2.7; ele: pT > 7 GeV, eta < 2.5
      for (const DressedLepton& l : leptonsFS_sel4l) {
        if (l.abspid() == PID::ELECTRON) leptons_sel4l.push_back(l);  // REDUNDANT: if (l.pT() > 7*GeV && l.abseta() < 2.5)
        else if (l.abspid() == PID::MUON) leptons_sel4l.push_back(l); // REDUNDANT: if (l.pT() > 6*GeV && l.abseta() < 2.7)
      }

      //////////////////////////////////////////////////////////////////
      // Exactly two opposite charged leptons
      //////////////////////////////////////////////////////////////////

      // Calculate total 'flavour' charge
      double totalcharge = 0;
      for (const Particle& l : leptons_sel4l)  totalcharge += l.pid();

      // Analyze 4 lepton events
      if (leptons_sel4l.size() != 4 || totalcharge != 0) vetoEvent;

      // Identify Z states from 4 lepton pairs
      Zstate Z1, Z2, Z1_alt, Z2_alt;
      if ( !identifyZstates(Z1, Z2, Z1_alt, Z2_alt, leptons_sel4l) )  vetoEvent;

      const double mZ1 = Z1.mom().mass();
      const double mZ2 = Z2.mom().mass();
      const double mZ1_alt = Z1_alt.mom().mass();
      const double mZ2_alt = Z2_alt.mom().mass();
      const double pTZ1 = Z1.mom().pT();
      const double pTZ2 = Z2.mom().pT();
      const double mZZ = (Z1.mom() + Z2.mom()).mass();
      const double pTZZ = (Z1.mom() + Z2.mom()).pT();

      // Event selections
      // pT(Z) > 2 GeV
      bool pass = pTZ1 > 2*GeV && pTZ2 > 2*GeV;
      if (!pass) vetoEvent;

      // Lepton kinematics: pT > 20, 15, 10 (8 if muon) GeV
      int n1 = 0, n2 = 0, n3 = 0;
      for (Particle& l : leptons_sel4l) {
        if (l.pT() > 20*GeV) ++n1;
        if (l.pT() > 15*GeV) ++n2;
        if (l.pT() > 10*GeV && l.abspid() == PID::ELECTRON) ++n3;
        if (l.pT() >  8*GeV && l.abspid() == PID::MUON) ++n3;
      }
      pass = pass && n1>=1 && n2>=2 && n3>=3;
      if (!pass) vetoEvent;

      // Dilepton mass: 50 < mZ1 < 120 GeV, 12 < mZ2 < 120 GeV
      pass = pass && mZ1 > 50*GeV && mZ1 < 120*GeV;
      pass = pass && mZ2 > 12*GeV && mZ2 < 120*GeV;
      if (!pass) vetoEvent;

      // Lepton separation: deltaR(l, l') > 0.1 (0.2) for same- (different-) flavor leptons
      for (size_t i = 0; i < leptons_sel4l.size(); ++i) {
        for (size_t j = i + 1; j < leptons_sel4l.size(); ++j) {
          const Particle& l1 = leptons_sel4l[i];
          const Particle& l2 = leptons_sel4l[j];
          pass = pass && deltaR(l1, l2) > (l1.abspid() == l2.abspid() ? 0.1 : 0.2);
          if (!pass) vetoEvent;
        }
      }

      // J/Psi veto: m(l+l-) > 5 GeV
      pass = pass && mZ1 > 5*GeV && mZ2 > 5*GeV && mZ1_alt > 5*GeV && mZ2_alt > 5*GeV;
      if (!pass) vetoEvent;

      // 80 < m4l < 1000 GeV
      pass = pass && mZZ > 80*GeV && mZZ < 1000*GeV;
      if (!pass) vetoEvent;

      // Fill histograms
      _h_ZZ_mZZ->fill(mZZ, weight);
      _h_ZZ_pTZZ->fill(pTZZ, weight);
    }


    /// Finalize
    void finalize() {
      const double norm = crossSection()/sumOfWeights()/femtobarn;
      scale(_h_ZZ_mZZ,  norm);
      scale(_h_ZZ_pTZZ, norm);
    }


    /// Generic Z candidate
    struct Zstate : public ParticlePair {
      Zstate() { }
      Zstate(ParticlePair _particlepair) : ParticlePair(_particlepair) { }
      FourMomentum mom() const { return first.momentum() + second.momentum(); }
      operator FourMomentum() const { return mom(); }
      static bool cmppT(const Zstate& lx, const Zstate& rx) { return lx.mom().pT() < rx.mom().pT(); }
    };


    /// @brief 4l to ZZ assignment algorithm
    ///
    /// ZZ->4l pairing
    /// - At least two same flavour opposite sign (SFOS) lepton pairs
    /// - Ambiguities in pairing are resolved following the procedure
    ///   1. the leading Z (Z1) is choosen as the SFOS with dilepton mass closet to Z mass
    ///   2. the subleading Z (Z2) is choosen as the remaining SFOS dilepton pair
    ///
    /// Z1, Z2: the selected pairing
    /// Z1_alt, Z2_alt: the alternative pairing (the same as Z1, Z2 in 2e2m case)
    bool identifyZstates(Zstate& Z1, Zstate& Z2, Zstate& Z1_alt, Zstate& Z2_alt, const Particles& leptons_sel4l) {
      const double ZMASS = 91.1876*GeV;
      bool findZZ = false;

      Particles part_pos_el, part_neg_el, part_pos_mu, part_neg_mu;
      for (const Particle& l : leptons_sel4l) {
        if (l.abspid() == PID::ELECTRON) {
          if (l.pid() < 0) part_neg_el.push_back(l);
          if (l.pid() > 0) part_pos_el.push_back(l);
        }
        else if (l.abspid() == PID::MUON) {
          if (l.pid() < 0) part_neg_mu.push_back(l);
          if (l.pid() > 0) part_pos_mu.push_back(l);
        }
      }

      // eeee/mmmm channel
      if ((part_neg_el.size() == 2 && part_pos_el.size() == 2) || (part_neg_mu.size() == 2 && part_pos_mu.size() == 2)) {
        findZZ = true;

        Zstate Zcand_1, Zcand_2, Zcand_3, Zcand_4;
        Zstate Zcand_1_tmp, Zcand_2_tmp, Zcand_3_tmp, Zcand_4_tmp;
        if (part_neg_el.size() == 2) { // eeee
          Zcand_1_tmp = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[0] ) );
          Zcand_2_tmp = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[1] ) );
          Zcand_3_tmp = Zstate( ParticlePair( part_neg_el[1],  part_pos_el[0] ) );
          Zcand_4_tmp = Zstate( ParticlePair( part_neg_el[1],  part_pos_el[1] ) );
        }
        else { // mmmm
          Zcand_1_tmp = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[0] ) );
          Zcand_2_tmp = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[1] ) );
          Zcand_3_tmp = Zstate( ParticlePair( part_neg_mu[1],  part_pos_mu[0] ) );
          Zcand_4_tmp = Zstate( ParticlePair( part_neg_mu[1],  part_pos_mu[1] ) );
        }

        // We can have the following pairs: (Z1 + Z4) || (Z2 + Z3)
        // Firstly, reorder withing each quadruplet to have
        //  - fabs(mZ1 - ZMASS) < fabs(mZ4 - ZMASS)
        //  - fabs(mZ2 - ZMASS) < fabs(mZ3 - ZMASS)
        if (fabs(Zcand_1_tmp.mom().mass() - ZMASS) < fabs(Zcand_4_tmp.mom().mass() - ZMASS)) {
          Zcand_1 = Zcand_1_tmp;
          Zcand_4 = Zcand_4_tmp;
        } else {
          Zcand_1 = Zcand_4_tmp;
          Zcand_4 = Zcand_1_tmp;
        }
        if (fabs(Zcand_2_tmp.mom().mass() - ZMASS) < fabs(Zcand_3_tmp.mom().mass() - ZMASS)) {
          Zcand_2 = Zcand_2_tmp;
          Zcand_3 = Zcand_3_tmp;
        } else {
          Zcand_2 = Zcand_3_tmp;
          Zcand_3 = Zcand_2_tmp;
        }

        // We can have the following pairs: (Z1 + Z4) || (Z2 + Z3)
        // Secondly, select the leading and subleading Z following
        //   1. the leading Z (Z1) is choosen as the SFOS with dilepton mass closet to Z mass
        //   2. the subleading Z (Z2) is choosen as the remaining SFOS dilepton pair
        if (fabs(Zcand_1.mom().mass() - ZMASS) < fabs(Zcand_2.mom().mass() - ZMASS)) {
          Z1 = Zcand_1;
          Z2 = Zcand_4;
          Z1_alt = Zcand_2;
          Z2_alt = Zcand_3;
        } else {
          Z1 = Zcand_2;
          Z2 = Zcand_3;
          Z1_alt = Zcand_1;
          Z2_alt = Zcand_4;
        }
      } // end of eeee/mmmm channel
      else if (part_neg_el.size() == 1 && part_pos_el.size() == 1 && part_neg_mu.size() == 1 && part_pos_mu.size() == 1) { // 2e2m channel
        findZZ = true;

        Zstate Zcand_1, Zcand_2;

        Zcand_1 = Zstate( ParticlePair( part_neg_mu[0],  part_pos_mu[0] ) );
        Zcand_2 = Zstate( ParticlePair( part_neg_el[0],  part_pos_el[0] ) );

        if (fabs(Zcand_1.mom().mass() - ZMASS) < fabs(Zcand_2.mom().mass() - ZMASS)) {
          Z1 = Zcand_1;
          Z2 = Zcand_2;
        } else {
          Z1 = Zcand_2;
          Z2 = Zcand_1;
        }
        Z1_alt = Z1;
        Z2_alt = Z2;
      }

      return findZZ;
    }


  private:

    Histo1DPtr _h_ZZ_pTZZ, _h_ZZ_mZZ;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1394865);

}
