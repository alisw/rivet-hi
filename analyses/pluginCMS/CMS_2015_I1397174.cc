#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Fully leptonic partonic ttbar analysis
  class CMS_2015_I1397174 : public Analysis {
  public:

    /// Minimal constructor
    CMS_2015_I1397174()
      : Analysis("CMS_2015_I1397174") { }


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {

      // Parton level top quarks
      addProjection(PartonicTops(PartonicTops::E_MU, false), "PartonTops");

      // Find jets not related to the top/W decays
      VetoedFinalState vfs;
      vfs.addDecayProductsVeto(PID::WPLUSBOSON);
      vfs.addDecayProductsVeto(PID::WMINUSBOSON);
      FastJets fj(vfs, FastJets::ANTIKT, 0.5, JetAlg::ALL_MUONS, JetAlg::ALL_INVISIBLES);
      addProjection(fj, "Jets");

      // Book histograms
      _hVis_nJet30_abs       = bookHisto1D( 1, 1, 1);
      _hVis_nJet30           = bookHisto1D( 2, 1, 1);
      _hVis_nJet60_abs       = bookHisto1D( 3, 1, 1);
      _hVis_nJet60           = bookHisto1D( 4, 1, 1);
      _hVis_nJet100_abs      = bookHisto1D( 5, 1, 1);
      _hVis_nJet100          = bookHisto1D( 6, 1, 1);

      _hVis_addJet1Pt_abs    = bookHisto1D( 7, 1, 1);
      _hVis_addJet1Pt        = bookHisto1D( 8, 1, 1);
      _hVis_addJet1Eta_abs   = bookHisto1D( 9, 1, 1);
      _hVis_addJet1Eta       = bookHisto1D(10, 1, 1);
      _hVis_addJet2Pt_abs    = bookHisto1D(11, 1, 1);
      _hVis_addJet2Pt        = bookHisto1D(12, 1, 1);
      _hVis_addJet2Eta_abs   = bookHisto1D(13, 1, 1);
      _hVis_addJet2Eta       = bookHisto1D(14, 1, 1);
      _hVis_addJJMass_abs    = bookHisto1D(15, 1, 1);
      _hVis_addJJMass        = bookHisto1D(16, 1, 1);
      _hVis_addJJDR_abs      = bookHisto1D(17, 1, 1);
      _hVis_addJJDR          = bookHisto1D(18, 1, 1);
      _hVis_addJJHT_abs      = bookHisto1D(19, 1, 1);
      _hVis_addJJHT          = bookHisto1D(20, 1, 1);

      _hFull_addJet1Pt_abs   = bookHisto1D(21, 1, 1);
      _hFull_addJet1Pt       = bookHisto1D(22, 1, 1);
      _hFull_addJet1Eta_abs  = bookHisto1D(23, 1, 1);
      _hFull_addJet1Eta      = bookHisto1D(24, 1, 1);
      _hFull_addJet2Pt_abs   = bookHisto1D(25, 1, 1);
      _hFull_addJet2Pt       = bookHisto1D(26, 1, 1);
      _hFull_addJet2Eta_abs  = bookHisto1D(27, 1, 1);
      _hFull_addJet2Eta      = bookHisto1D(28, 1, 1);
      _hFull_addJJMass_abs   = bookHisto1D(29, 1, 1);
      _hFull_addJJMass       = bookHisto1D(30, 1, 1);
      _hFull_addJJDR_abs     = bookHisto1D(31, 1, 1);
      _hFull_addJJDR         = bookHisto1D(32, 1, 1);
      _hFull_addJJHT_abs     = bookHisto1D(33, 1, 1);
      _hFull_addJJHT         = bookHisto1D(34, 1, 1);

      _hVis_addBJet1Pt_abs   = bookHisto1D(35, 1, 1);
      _hVis_addBJet1Pt       = bookHisto1D(36, 1, 1);
      _hVis_addBJet1Eta_abs  = bookHisto1D(37, 1, 1);
      _hVis_addBJet1Eta      = bookHisto1D(38, 1, 1);
      _hVis_addBJet2Pt_abs   = bookHisto1D(39, 1, 1);
      _hVis_addBJet2Pt       = bookHisto1D(40, 1, 1);
      _hVis_addBJet2Eta_abs  = bookHisto1D(41, 1, 1);
      _hVis_addBJet2Eta      = bookHisto1D(42, 1, 1);
      _hVis_addBBMass_abs    = bookHisto1D(43, 1, 1);
      _hVis_addBBMass        = bookHisto1D(44, 1, 1);
      _hVis_addBBDR_abs      = bookHisto1D(45, 1, 1);
      _hVis_addBBDR          = bookHisto1D(46, 1, 1);

      _hFull_addBJet1Pt_abs  = bookHisto1D(47, 1, 1);
      _hFull_addBJet1Pt      = bookHisto1D(48, 1, 1);
      _hFull_addBJet1Eta_abs = bookHisto1D(49, 1, 1);
      _hFull_addBJet1Eta     = bookHisto1D(50, 1, 1);
      _hFull_addBJet2Pt_abs  = bookHisto1D(51, 1, 1);
      _hFull_addBJet2Pt      = bookHisto1D(52, 1, 1);
      _hFull_addBJet2Eta_abs = bookHisto1D(53, 1, 1);
      _hFull_addBJet2Eta     = bookHisto1D(54, 1, 1);
      _hFull_addBBMass_abs   = bookHisto1D(55, 1, 1);
      _hFull_addBBMass       = bookHisto1D(56, 1, 1);
      _hFull_addBBDR_abs     = bookHisto1D(57, 1, 1);
      _hFull_addBBDR         = bookHisto1D(58, 1, 1);

      _h_gap_addJet1Pt       = bookProfile1D(59, 1, 1);
      _h_gap_addJet1Pt_eta0  = bookProfile1D(60, 1, 1);
      _h_gap_addJet1Pt_eta1  = bookProfile1D(61, 1, 1);
      _h_gap_addJet1Pt_eta2  = bookProfile1D(62, 1, 1);
      _h_gap_addJet2Pt       = bookProfile1D(63, 1, 1);
      _h_gap_addJet2Pt_eta0  = bookProfile1D(64, 1, 1);
      _h_gap_addJet2Pt_eta1  = bookProfile1D(65, 1, 1);
      _h_gap_addJet2Pt_eta2  = bookProfile1D(66, 1, 1);
      _h_gap_addJetHT        = bookProfile1D(67, 1, 1);
      _h_gap_addJetHT_eta0   = bookProfile1D(68, 1, 1);
      _h_gap_addJetHT_eta1   = bookProfile1D(69, 1, 1);
      _h_gap_addJetHT_eta2   = bookProfile1D(70, 1, 1);
    }


    void analyze(const Event& event) {

      // The objects used in the PAPER 12-041 are defined as follows (see p.16 for details):
      //
      //   * Leptons    : from the W boson decays after FSR
      //   * Jets       : anti-kT R=0.5 to all stable particles
      //                               exclude W->enu, munu, taunu
      //   * B jet      : B-Ghost matched
      //   * B from top : B hadron from top->b decay
      //
      // Visible phase space definition:
      //
      //   * Leptons         : pT > 20, |eta| < 2.4
      //   * B jets from top : pT > 30, |eta| < 2.4
      //     Additional jets : pT > 20, |eta| < 2.4
      //   *
      // Full phase space definition:
      //
      //   * Correction to dilepton BR from W boson BR
      //   * No cut on top decay products
      //   * Additional jets : pT > 20, |eta| < 2.4

      // Do the analysis only for the ttbar full leptonic channel, removing tau decays
      const Particles partontops = apply<ParticleFinder>(event, "PartonTops").particlesByPt();
      if (partontops.size() != 2) vetoEvent;
      const Particle& t1 = partontops[0];
      const Particle& t2 = partontops[1];

      // Apply acceptance cuts on top-decay leptons (existence should be guaranteed)
      const auto isPromptChLepton = [](const Particle& p){return isChargedLepton(p) && !fromDecay(p);};
      const Particle lep1 = t1.allDescendants(lastParticleWith(isPromptChLepton)).front();
      const Particle lep2 = t2.allDescendants(lastParticleWith(isPromptChLepton)).front();
      if (lep1.pT() < 1e-9*GeV || lep2.pT() < 1e-9*GeV) vetoEvent; // sanity check?

      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);
      int nJet30 = 0, nJet60 = 0, nJet100 = 0;
      Jets topBJets, addJets, addBJets, addJets_eta0, addJets_eta1, addJets_eta2;
      for (const Jet& jet : jets) {
        if (jet.pT() >  30*GeV) nJet30 += 1;
        if (jet.pT() >  60*GeV) nJet60 += 1;
        if (jet.pT() > 100*GeV) nJet100 += 1;

        const bool isBtagged = jet.bTagged();
        const bool isBFromTop = any(jet.bTags(), hasParticleAncestorWith(Cuts::abspid == PID::TQUARK));

        if (isBFromTop) {
          if (jet.pT() > 30*GeV) topBJets.push_back(jet);
        } else {
          addJets.push_back(jet);
          if (isBtagged) addBJets.push_back(jet);
          if      (jet.abseta() < 0.8 ) addJets_eta0.push_back(jet);
          else if (jet.abseta() < 1.5 ) addJets_eta1.push_back(jet);
          else if (jet.abseta() < 2.4 ) addJets_eta2.push_back(jet);
        }
      }


      const bool isVisiblePS = topBJets.size() >= 2
        && lep1.pT() > 20*GeV && lep1.abseta() < 2.4 && lep2.pT() > 20*GeV && lep2.abseta() < 2.4;
      MSG_DEBUG(isVisiblePS << ": #b(top) = " << topBJets.size()
                << "; l1 = " << lep1.pT() << ", " << lep1.abseta()
                << "; l2 = " << lep2.pT() << ", " << lep2.abseta());

      const double weight = event.weight();


      if (isVisiblePS) {
        fillWithOF(_hVis_nJet30_abs,  nJet30, weight);
        fillWithOF(_hVis_nJet30,      nJet30, weight);
        fillWithOF(_hVis_nJet60_abs,  nJet60, weight);
        fillWithOF(_hVis_nJet60,      nJet60, weight);
        fillWithOF(_hVis_nJet100_abs, nJet100, weight);
        fillWithOF(_hVis_nJet100,     nJet100, weight);

        fillGapFractions(addJets, _h_gap_addJet1Pt, _h_gap_addJet2Pt, _h_gap_addJetHT, weight);
        fillGapFractions(addJets_eta0, _h_gap_addJet1Pt_eta0, _h_gap_addJet2Pt_eta0, _h_gap_addJetHT_eta0, weight);
        fillGapFractions(addJets_eta1, _h_gap_addJet1Pt_eta1, _h_gap_addJet2Pt_eta1, _h_gap_addJetHT_eta1, weight);
        fillGapFractions(addJets_eta2, _h_gap_addJet1Pt_eta2, _h_gap_addJet2Pt_eta2, _h_gap_addJetHT_eta2, weight);
      }

      // Plots with two additional jets
      if (addJets.size() >= 1) {
        const double ht = sum(addJets, pT, 0.0);
        _hFull_addJJHT_abs->fill(ht/GeV, weight);
        _hFull_addJJHT    ->fill(ht/GeV, weight);
        if (isVisiblePS) {
          _hVis_addJJHT_abs->fill(ht/GeV, weight);
          _hVis_addJJHT    ->fill(ht/GeV, weight);
        }

        const Jet& j1 = addJets[0];
        _hFull_addJet1Pt_abs ->fill(j1.pT()/GeV, weight);
        _hFull_addJet1Pt     ->fill(j1.pT()/GeV, weight);
        _hFull_addJet1Eta_abs->fill(j1.abseta(), weight);
        _hFull_addJet1Eta    ->fill(j1.abseta(), weight);
        if (isVisiblePS) {
          _hVis_addJet1Pt_abs ->fill(j1.pT()/GeV, weight);
          _hVis_addJet1Pt     ->fill(j1.pT()/GeV, weight);
          _hVis_addJet1Eta_abs->fill(j1.abseta(), weight);
          _hVis_addJet1Eta    ->fill(j1.abseta(), weight);
        }

        if (addJets.size() >= 2) {
          const Jet& j2 = addJets[1];

          _hFull_addJet2Pt_abs ->fill(j2.pT()/GeV, weight);
          _hFull_addJet2Pt     ->fill(j2.pT()/GeV, weight);
          _hFull_addJet2Eta_abs->fill(j2.abseta(), weight);
          _hFull_addJet2Eta    ->fill(j2.abseta(), weight);
          if (isVisiblePS) {
            _hVis_addJet2Pt_abs ->fill(j2.pT()/GeV, weight);
            _hVis_addJet2Pt     ->fill(j2.pT()/GeV, weight);
            _hVis_addJet2Eta_abs->fill(j2.abseta(), weight);
            _hVis_addJet2Eta    ->fill(j2.abseta(), weight);
          }

          const double jjmass = (j1.mom() + j2.mom()).mass();
          const double jjdR = deltaR(j1, j2);
          _hFull_addJJMass_abs->fill(jjmass/GeV, weight);
          _hFull_addJJMass    ->fill(jjmass/GeV, weight);
          _hFull_addJJDR_abs  ->fill(jjdR, weight);
          _hFull_addJJDR      ->fill(jjdR, weight);
          if (isVisiblePS) {
            _hVis_addJJMass_abs->fill(jjmass/GeV, weight);
            _hVis_addJJMass    ->fill(jjmass/GeV, weight);
            _hVis_addJJDR_abs  ->fill(jjdR, weight);
            _hVis_addJJDR      ->fill(jjdR, weight);
          }
        }
      }


      // Same set of plots if there are additional b-jets
      if (addBJets.size() >= 1) {
        const Jet& b1 = addBJets[0];
        _hFull_addBJet1Pt_abs ->fill(b1.pT()/GeV, weight);
        _hFull_addBJet1Pt     ->fill(b1.pT()/GeV, weight);
        _hFull_addBJet1Eta_abs->fill(b1.abseta(), weight);
        _hFull_addBJet1Eta    ->fill(b1.abseta(), weight);
        if (isVisiblePS) {
          _hVis_addBJet1Pt_abs ->fill(b1.pT()/GeV, weight);
          _hVis_addBJet1Pt     ->fill(b1.pT()/GeV, weight);
          _hVis_addBJet1Eta_abs->fill(b1.abseta(), weight);
          _hVis_addBJet1Eta    ->fill(b1.abseta(), weight);
        }

        if (addBJets.size() >= 2) {
          const Jet& b2 = addBJets[1];

          _hFull_addBJet2Pt_abs ->fill(b2.pT()/GeV, weight);
          _hFull_addBJet2Pt     ->fill(b2.pT()/GeV, weight);
          _hFull_addBJet2Eta_abs->fill(b2.abseta(), weight);
          _hFull_addBJet2Eta    ->fill(b2.abseta(), weight);
          if (isVisiblePS) {
            _hVis_addBJet2Pt_abs ->fill(b2.pT()/GeV, weight);
            _hVis_addBJet2Pt     ->fill(b2.pT()/GeV, weight);
            _hVis_addBJet2Eta_abs->fill(b2.abseta(), weight);
            _hVis_addBJet2Eta    ->fill(b2.abseta(), weight);
          }

          const double bbmass = (b1.mom() + b2.mom()).mass();
          const double bbdR = deltaR(b1, b2);
          _hFull_addBBMass_abs->fill(bbmass/GeV, weight);
          _hFull_addBBMass    ->fill(bbmass/GeV, weight);
          _hFull_addBBDR_abs  ->fill(bbdR, weight);
          _hFull_addBBDR      ->fill(bbdR, weight);
          if (isVisiblePS) {
            _hVis_addBBMass_abs->fill(bbmass/GeV, weight);
            _hVis_addBBMass    ->fill(bbmass/GeV, weight);
            _hVis_addBBDR_abs  ->fill(bbdR, weight);
            _hVis_addBBDR      ->fill(bbdR, weight);
          }
        }
      }

    }


    void finalize() {
      const double ttbarXS = !std::isnan(crossSectionPerEvent()) ? crossSection() : 252.89*picobarn;
      if (std::isnan(crossSectionPerEvent()))
        MSG_INFO("No valid cross-section given, using NNLO (arXiv:1303.6254; sqrt(s)=8 TeV, m_t=172.5 GeV): " << ttbarXS/picobarn << " pb");

      normalize({_hVis_nJet30,_hVis_nJet60, _hVis_nJet100,
            _hVis_addJet1Pt, _hVis_addJet1Eta, _hVis_addJet2Pt, _hVis_addJet2Eta,
            _hVis_addJJMass, _hVis_addJJDR, _hVis_addJJHT,
            _hFull_addJet1Pt, _hFull_addJet1Eta, _hFull_addJet2Pt, _hFull_addJet2Eta,
            _hFull_addJJMass, _hFull_addJJDR, _hFull_addJJHT,
            _hVis_addBJet1Pt, _hVis_addBJet1Eta, _hVis_addBJet2Pt, _hVis_addBJet2Eta,
            _hVis_addBBMass, _hVis_addBBDR,
            _hFull_addBJet1Pt, _hFull_addBJet1Eta, _hFull_addBJet2Pt, _hFull_addBJet2Eta,
            _hFull_addBBMass, _hFull_addBBDR});

      const double xsPerWeight = ttbarXS/picobarn / sumOfWeights();
      scale({_hVis_nJet30_abs, _hVis_nJet60_abs, _hVis_nJet100_abs,
            _hVis_addJet1Pt_abs, _hVis_addJet1Eta_abs, _hVis_addJet2Pt_abs, _hVis_addJet2Eta_abs,
            _hVis_addJJMass_abs, _hVis_addJJDR_abs, _hVis_addJJHT_abs,
            _hVis_addBJet1Pt_abs, _hVis_addBJet1Eta_abs, _hVis_addBJet2Pt_abs, _hVis_addBJet2Eta_abs,
            _hVis_addBBMass_abs, _hVis_addBBDR_abs}, xsPerWeight);

      const double sfull = xsPerWeight / 0.0454; //< correct for dilepton branching fraction
      scale({_hFull_addJet1Pt_abs, _hFull_addJet1Eta_abs, _hFull_addJet2Pt_abs, _hFull_addJet2Eta_abs,
            _hFull_addJJMass_abs, _hFull_addJJDR_abs, _hFull_addJJHT_abs,
            _hFull_addBJet1Pt_abs, _hFull_addBJet1Eta_abs, _hFull_addBJet2Pt_abs, _hFull_addBJet2Eta_abs,
            _hFull_addBBMass_abs, _hFull_addBBDR_abs}, sfull);
    }

    //@}


    void fillWithOF(Histo1DPtr h, double x, double w) {
      h->fill(std::min(x, h->xMax()-1e-9), w);
    }


    void fillGapFractions(const Jets& addJets, Profile1DPtr h_gap_addJet1Pt, Profile1DPtr h_gap_addJet2Pt, Profile1DPtr h_gap_addJetHT, double weight) {
      const double j1pt = (addJets.size() > 0) ? addJets[0].pT() : 0;
      for (size_t i = 0; i < h_gap_addJet1Pt->numBins(); ++i) {
        const double binCenter = h_gap_addJet1Pt->bin(i).xMid();
        h_gap_addJet1Pt->fillBin(i, int(j1pt/GeV < binCenter), weight);
      }

      const double j2pt = (addJets.size() > 1) ? addJets[1].pT() : 0;
      for (size_t i = 0; i < h_gap_addJet2Pt->numBins(); ++i) {
        const double binCenter = h_gap_addJet2Pt->bin(i).xMid();
        h_gap_addJet2Pt->fillBin(i, int(j2pt/GeV < binCenter), weight);
      }

      const double ht = sum(addJets, pT, 0.);
      for (size_t i = 0; i < h_gap_addJetHT->numBins(); ++i) {
        const double binCenter = h_gap_addJetHT->bin(i).xMid();
        h_gap_addJetHT->fillBin(i, int(ht/GeV < binCenter) , weight);
      }
    }


    // @name Histogram data members
    //@{

    Histo1DPtr _hVis_nJet30_abs, _hVis_nJet60_abs, _hVis_nJet100_abs;
    Histo1DPtr _hVis_addJet1Pt_abs, _hVis_addJet1Eta_abs, _hVis_addJet2Pt_abs, _hVis_addJet2Eta_abs;
    Histo1DPtr _hVis_addJJMass_abs, _hVis_addJJDR_abs, _hVis_addJJHT_abs;
    Histo1DPtr _hFull_addJet1Pt_abs, _hFull_addJet1Eta_abs, _hFull_addJet2Pt_abs, _hFull_addJet2Eta_abs;
    Histo1DPtr _hFull_addJJMass_abs, _hFull_addJJDR_abs, _hFull_addJJHT_abs;
    Histo1DPtr _hVis_addBJet1Pt_abs, _hVis_addBJet1Eta_abs, _hVis_addBJet2Pt_abs, _hVis_addBJet2Eta_abs;
    Histo1DPtr _hVis_addBBMass_abs, _hVis_addBBDR_abs;
    Histo1DPtr _hFull_addBJet1Pt_abs, _hFull_addBJet1Eta_abs, _hFull_addBJet2Pt_abs, _hFull_addBJet2Eta_abs;
    Histo1DPtr _hFull_addBBMass_abs, _hFull_addBBDR_abs;

    Histo1DPtr _hVis_nJet30, _hVis_nJet60, _hVis_nJet100;
    Histo1DPtr _hVis_addJet1Pt, _hVis_addJet1Eta, _hVis_addJet2Pt, _hVis_addJet2Eta;
    Histo1DPtr _hVis_addJJMass, _hVis_addJJDR, _hVis_addJJHT;
    Histo1DPtr _hFull_addJet1Pt, _hFull_addJet1Eta, _hFull_addJet2Pt, _hFull_addJet2Eta;
    Histo1DPtr _hFull_addJJMass, _hFull_addJJDR, _hFull_addJJHT;
    Histo1DPtr _hVis_addBJet1Pt, _hVis_addBJet1Eta, _hVis_addBJet2Pt, _hVis_addBJet2Eta;
    Histo1DPtr _hVis_addBBMass, _hVis_addBBDR;
    Histo1DPtr _hFull_addBJet1Pt, _hFull_addBJet1Eta, _hFull_addBJet2Pt, _hFull_addBJet2Eta;
    Histo1DPtr _hFull_addBBMass, _hFull_addBBDR;

    Profile1DPtr _h_gap_addJet1Pt, _h_gap_addJet1Pt_eta0, _h_gap_addJet1Pt_eta1, _h_gap_addJet1Pt_eta2;
    Profile1DPtr _h_gap_addJet2Pt, _h_gap_addJet2Pt_eta0, _h_gap_addJet2Pt_eta1, _h_gap_addJet2Pt_eta2;
    Profile1DPtr _h_gap_addJetHT, _h_gap_addJetHT_eta0, _h_gap_addJetHT_eta1, _h_gap_addJetHT_eta2;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1397174);


}
