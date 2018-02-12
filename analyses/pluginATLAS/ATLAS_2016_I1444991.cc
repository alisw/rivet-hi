// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VisibleFinalState.hh"

namespace Rivet {


  class ATLAS_2016_I1444991 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1444991);


  public:

    /// Book histograms and initialise projections before the run
    void init() {

      // All particles within |eta| < 5.0
      const FinalState FS(Cuts::abseta < 5.0);

      // Project photons for dressing
      IdentifiedFinalState photon_id(FS);
      photon_id.acceptIdPair(PID::PHOTON);

      // Project dressed electrons with pT > 15 GeV and |eta| < 2.47
      IdentifiedFinalState el_id(FS);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState el_bare(el_id);
      Cut cuts = (Cuts::abseta < 2.47) && ( (Cuts::abseta <= 1.37) || (Cuts::abseta >= 1.52) ) && (Cuts::pT > 15*GeV);
      DressedLeptons el_dressed_FS(photon_id, el_bare, 0.1, cuts, true);
      declare(el_dressed_FS,"EL_DRESSED_FS");

      // Project dressed muons with pT > 15 GeV and |eta| < 2.5
      IdentifiedFinalState mu_id(FS);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState mu_bare(mu_id);
      DressedLeptons mu_dressed_FS(photon_id, mu_bare, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 15*GeV, true);
      declare(mu_dressed_FS,"MU_DRESSED_FS");

      // get MET from generic invisibles
      VetoedFinalState inv_fs(FS);
      inv_fs.addVetoOnThisFinalState(VisibleFinalState(FS));
      declare(inv_fs, "InvisibleFS");

      // Project jets
      FastJets jets(FS, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(JetAlg::NO_INVISIBLES);
      jets.useMuons(JetAlg::NO_MUONS);
      declare(jets, "jets");

      // Book histograms
      _h_Njets        = bookHisto1D( 2,1,1);
      _h_PtllMET      = bookHisto1D( 3,1,1);
      _h_Yll          = bookHisto1D( 4,1,1);
      _h_PtLead       = bookHisto1D( 5,1,1);
      _h_Njets_norm   = bookHisto1D( 6,1,1);
      _h_PtllMET_norm = bookHisto1D( 7,1,1);
      _h_Yll_norm     = bookHisto1D( 8,1,1);
      _h_PtLead_norm  = bookHisto1D( 9,1,1);
      _h_JetVeto      = bookScatter2D(10, 1, 1, true);

      //histos for jetveto
      std::vector<double> ptlead25_bins = { 0., 25., 300. };
      std::vector<double> ptlead40_bins = { 0., 40., 300. };
      _h_pTj1_sel25 = bookHisto1D( "pTj1_sel25", ptlead25_bins, "", "", "" );
      _h_pTj1_sel40 = bookHisto1D( "pTj1_sel40", ptlead40_bins, "", "", "" );
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Get final state particles
      const FinalState& ifs = applyProjection<FinalState>(event, "InvisibleFS");
      const vector<DressedLepton>& good_mu = applyProjection<DressedLeptons>(event, "MU_DRESSED_FS").dressedLeptons();
      const vector<DressedLepton>& el_dressed = applyProjection<DressedLeptons>(event, "EL_DRESSED_FS").dressedLeptons();
      const Jets& jets = applyProjection<FastJets>(event, "jets").jetsByPt(Cuts::pT>25*GeV && Cuts::abseta < 4.5);

      //find good electrons
      vector<DressedLepton> good_el;
      for (const DressedLepton& el : el_dressed){
        bool keep = true;
        for (const DressedLepton& mu : good_mu) {
          keep &= deltaR(el, mu) >= 0.1;
        }
        if (keep)  good_el += el;
      }

      // select only emu events
      if ((good_el.size() != 1) || good_mu.size() != 1)  vetoEvent;

      //built dilepton
      FourMomentum dilep = good_el[0].momentum() + good_mu[0].momentum();
      double Mll = dilep.mass();
      double Yll = dilep.rapidity();
      double DPhill = fabs(deltaPhi(good_el[0], good_mu[0]));
      double pTl1 = (good_el[0].pT() > good_mu[0].pT())? good_el[0].pT() : good_mu[0].pT();

      //get MET
      FourMomentum met;
      foreach (const Particle& p, ifs.particles())  met += p.momentum();

      // do a few cuts before looking at jets
      if (pTl1 <= 22. || DPhill >= 1.8 || met.pT() <= 20.)  vetoEvent;
      if (Mll <= 10. || Mll >= 55.)  vetoEvent;

      Jets jets_selected;
      foreach (const Jet &j, jets) {
        if( j.abseta() > 2.4 && j.pT()<=30*GeV ) continue;
        bool keep = true;
        foreach(DressedLepton el, good_el) {
          keep &= deltaR(j, el) >= 0.3;
        }
        if (keep)  jets_selected += j;
      }

      double PtllMET = (met + good_el[0].momentum() + good_mu[0].momentum()).pT();

      double Njets = jets_selected.size() > 2 ? 2 : jets_selected.size();
      double pTj1 = jets_selected.size()? jets_selected[0].pT() : 0.1;

      // Fill histograms
      _h_Njets->fill(Njets, weight);
      _h_PtllMET->fill(PtllMET, weight);
      _h_Yll->fill(fabs(Yll), weight);
      _h_PtLead->fill(pTj1, weight);
      _h_Njets_norm->fill(Njets, weight);
      _h_PtllMET_norm->fill(PtllMET, weight);
      _h_Yll_norm->fill(fabs(Yll), weight);
      _h_PtLead_norm->fill(pTj1, weight);
      _h_pTj1_sel25->fill(pTj1, weight);
      _h_pTj1_sel40->fill(pTj1, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xs = crossSectionPerEvent()/femtobarn;

      /// @todo Normalise, scale and otherwise manipulate histograms here
      scale(_h_Njets, xs);
      scale(_h_PtllMET, xs);
      scale(_h_Yll, xs);
      scale(_h_PtLead, xs);
      normalize(_h_Njets_norm);
      normalize(_h_PtllMET_norm);
      normalize(_h_Yll_norm);
      normalize(_h_PtLead_norm);
      scale(_h_pTj1_sel25, xs);
      scale(_h_pTj1_sel40, xs);
      normalize(_h_pTj1_sel25);
      normalize(_h_pTj1_sel40);
      // fill jet veto efficiency histogram
      _h_JetVeto->point(0).setY(_h_pTj1_sel25->bin(0).sumW(), sqrt(_h_pTj1_sel25->bin(0).sumW2()));
      _h_JetVeto->point(1).setY(_h_PtLead_norm->bin(0).sumW(), sqrt(_h_PtLead_norm->bin(0).sumW2()));
      _h_JetVeto->point(2).setY(_h_pTj1_sel40->bin(0).sumW(), sqrt(_h_pTj1_sel25->bin(0).sumW2()));

    }

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_Njets;
    Histo1DPtr _h_PtllMET;
    Histo1DPtr _h_Yll;
    Histo1DPtr _h_PtLead;
    Histo1DPtr _h_Njets_norm;
    Histo1DPtr _h_PtllMET_norm;
    Histo1DPtr _h_Yll_norm;
    Histo1DPtr _h_PtLead_norm;

    Scatter2DPtr _h_JetVeto;

    Histo1DPtr _h_pTj1_sel25;
    Histo1DPtr _h_pTj1_sel40;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1444991);

}
