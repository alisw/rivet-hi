// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Charge asymmetry in top quark pair production in dilepton channel
  class ATLAS_2016_I1449082 : public Analysis {
  public:

    const double MW = 80.300*GeV;
    const double MTOP = 172.5*GeV;

    enum MeasureType { kInclMeas, kmttMeas, kbetaMeas, kptMeas, kNmeas };
    const size_t kNbins = 2;
    const double bins[kNmeas][3];
    const string measStr[kNmeas] = {"incl","mtt","beta","ptt"};
    const string rangeStr[kNmeas][2];


    /// Constructor
    //DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1449082);
    ATLAS_2016_I1449082() : Analysis("ATLAS_2016_I1449082"),
                            // inclusive (dummy), mtt [GeV], beta, pTtt
                            bins{ { 0., 1., 2. }, { 0., 500., 2000.}, { 0., 0.6 , 1.0}, { 0., 30. , 1000.} },
                            rangeStr{ { "0_1", "1_2"}, { "0_500", "500_2000"}, { "0_0.6", "0.6_1.0"}, { "0_30" , "30_1000"} }
    {  }

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Cuts
      const Cut eta_full = Cuts::abseta < 5.0;
      const Cut lep_cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
      // All final state particles
      FinalState fs(eta_full);
      // Get photons to dress leptons
      IdentifiedFinalState photons(fs, PID::PHOTON);



      // Electron projections
      // ---------------------
      // Electron/muons are defined from electron/muon and photons within a cone of DR = 0.1.
      // No isolation condition is imposed. The parent of the electron/muon is required not to be a hadron or quark.
      IdentifiedFinalState el_id(fs, {PID::ELECTRON,-PID::ELECTRON});
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      // Electron dressing
      DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts, true);
      declare(dressedelectrons, "dressedelectrons");
      DressedLeptons dressedelectrons_full(photons, electrons, 0.1, eta_full, true);

      // Muon projections
      // ---------------------
      IdentifiedFinalState mu_id(fs, {PID::MUON,-PID::MUON});
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      // Muon dressing
      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true);
      declare(dressedmuons, "dressedmuons");
      DressedLeptons dressedmuons_full(photons, muons, 0.1, eta_full, true);

      // Neutrino projections
      // ---------------------
      // Missing ET is calculated as the 4–vector sum of neutrinos from W/Z-boson decays. Tau decays are
      // included. A neutrino is treated as a detectable particle and is selected for consideration in the same
      // way as electrons or muons, i.e. the parent is required not to be a hadron or quark (u − b).
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      declare(neutrinos, "neutrinos");

      // Jets projections
      // ---------------------
      // Jets are defined with the anti-kt algorithm, clustering all stable particles excluding the electrons,
      // muons, neutrinos, and photons used in the definition of the selected leptons.
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(dressedelectrons_full);
      vfs.addVetoOnThisFinalState(dressedmuons_full);
      vfs.addVetoOnThisFinalState(neutrinos);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "Jets");

      // Book histograms
      _h_dEta      = bookHisto1D(1, 1, 1);
      _h_dY        = bookHisto1D(2, 1, 1);
      for (size_t iM = 0; iM < kNmeas; ++iM) {
        _h_Acll[iM] = bookScatter2D(3+iM, 1, 1);
        _h_Actt[iM] = bookScatter2D(7+iM, 1, 1);
      }
      for (size_t iM = 0; iM < kNmeas; ++iM) {
        for (size_t iB = 0; iB < kNbins; ++iB) {
          _h_dEta_asym[iM][iB] = bookHisto1D( "dEta_asym_" + measStr[iM] + "_bin" + rangeStr[iM][iB],  2,  -10., 10.);
          _h_dY_asym  [iM][iB] = bookHisto1D( "dY_asym_"   + measStr[iM] + "_bin" + rangeStr[iM][iB],  2,  -10., 10.);
        }
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the electrons and muons
      const vector<DressedLepton> dressedelectrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      const vector<DressedLepton> dressedmuons     = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      const vector<DressedLepton> leptons = dressedelectrons + dressedmuons;
      // Require at least 2 leptons in the event
      if (leptons.size() < 2) vetoEvent;

      // Get the neutrinos
      const Particles neutrinos = apply<PromptFinalState>(event, "neutrinos").particlesByPt();
      // Require at least 2 neutrinos in the event (ick)
      if (neutrinos.size() < 2)  vetoEvent;

      // Get jets and apply selection
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      // Require at least 2 jets in the event
      if (jets.size() < 2)  vetoEvent;

      // Remaining selections
      // Events where leptons and jets overlap, within dR = 0.4, are rejected.
      for (const DressedLepton& lepton : leptons) {
        if (any(jets, deltaRLess(lepton, 0.4))) vetoEvent;
      }

      // Construct pseudo-tops
      // Exactly 2 opposite-sign leptons are required (e/mu)
      if (leptons.size() != 2) vetoEvent;
      if ( (leptons[0].charge() * leptons[1].charge()) > 0.) vetoEvent;
      const FourMomentum lep_p = (leptons[0].charge3() > 0) ? leptons[0] : leptons[1];
      const FourMomentum lep_n = (leptons[0].charge3() > 0) ? leptons[1] : leptons[0];

      // Only the 2 leading pT selected neutrinos are considered
      const FourMomentum nu1 = neutrinos[0].momentum();
      const FourMomentum nu2 = neutrinos[1].momentum();

      // Two jets correspond to the two leading jets in the event.
      // If there is any b-tagged jet in the event, then the b-tagged jets
      // are preferentially selected over the non-tagged jets without taking into account its pT.
      // A jet is a b–jet if any B–hadron is included in the jet.
      // Only B-hadrons with an initial pT > 5 GeV are considered.
      Jets bjets, lightjets;
      for (const Jet& jet : jets) {
        (jet.bTagged(Cuts::pT > 5*GeV) ? bjets : lightjets) += jet;
      }
      // Already sorted by construction, since jets is sorted by decreasing pT
      // std::sort(bjets.begin()    , bjets.end()    , cmpMomByPt);
      // std::sort(lightjets.begin(), lightjets.end(), cmpMomByPt);

      // Initially take 2 highest pT jets
      FourMomentum bjet1 = jets[0];
      FourMomentum bjet2 = jets[1];
      if (!bjets.empty()) {
        bjet1 = bjets[0];
        bjet2 = (bjets.size() > 1) ? bjets[1] : lightjets[0]; //< We should have a light jet because >=2 jets requirement
      } else {
        // No btagged jets --> should have >= 2 light jets
        bjet1 = lightjets[0];
        bjet2 = lightjets[1];
      }

      // Construct pseudo-W bosons from lepton-neutrino combinations
      // Minimize the difference between the mass computed from each lepton-neutrino combination and the W boson mass
      const double massDiffW1 = fabs( (nu1 + lep_p).mass() - MW ) + fabs( (nu2 + lep_n).mass() - MW );
      const double massDiffW2 = fabs( (nu1 + lep_n).mass() - MW ) + fabs( (nu2 + lep_p).mass() - MW );
      const FourMomentum Wp = (massDiffW1 < massDiffW2) ? nu1+lep_p : nu2+lep_p;
      const FourMomentum Wn = (massDiffW1 < massDiffW2) ? nu2+lep_n : nu1+lep_n;

      // Construct pseudo-tops from jets and pseudo-W bosons
      // Minimize the difference between the mass computed from each W-boson and b-jet combination and the top mass
      const double massDiffT1 = fabs( (Wp+bjet1).mass()*GeV - MTOP ) + fabs( (Wn+bjet2).mass()*GeV - MTOP );
      const double massDiffT2 = fabs( (Wp+bjet2).mass()*GeV - MTOP ) + fabs( (Wn+bjet1).mass()*GeV - MTOP );
      const FourMomentum top_p = (massDiffT1 < massDiffT2) ? Wp+bjet1 : Wp+bjet2;
      const FourMomentum top_n = (massDiffT1 < massDiffT2) ? Wn+bjet2 : Wn+bjet1;

      // Calculate d|eta|, d|y|, etc.
      double dEta = lep_p.abseta() - lep_n.abseta();
      double dY   = top_p.absrapidity() - top_n.absrapidity();
      double mtt  = (top_p + top_n).mass()*GeV;
      double beta = fabs( (top_p + top_n).pz() ) / (top_p + top_n).E();
      double pttt = (top_p + top_n).pt()*GeV;

      // Fill histos, counters
      const double weight = event.weight();
      _h_dEta->fill(dEta, weight);
      _h_dY  ->fill(dY  , weight);
      // Histos for inclusive and differential asymmetries
      int mttBinID  = getBinID(kmttMeas , mtt);
      int betaBinID = getBinID(kbetaMeas, beta);
      int ptttBinID = getBinID(kptMeas  , pttt);
      for (int iM = 0; iM < kNmeas; ++iM) {
        int binID = -1;
        switch (iM) {
          case kInclMeas : binID = 0;         break;
          case kmttMeas  : binID = mttBinID ; break;
          case kbetaMeas : binID = betaBinID; break;
          case kptMeas   : binID = ptttBinID; break;
          default: binID = -1; break;
        }
        if (binID >= 0) {
          _h_dY_asym  [iM][binID] ->fill(dY  , weight);
          _h_dEta_asym[iM][binID] ->fill(dEta, weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Calculate charge asymmetries and fill them to hists
      // Just for cross-check, calculate asymmetries from original dEta/dY histos
      double asym = 0, err = 0;
      calcAsymAndError(_h_dEta, asym, err);
      MSG_INFO("Lepton inclusive asymmetry from histo:  = " << asym << " +- " << err );
      calcAsymAndError(_h_dY, asym, err);
      MSG_INFO("ttbar inclusive asymmetry from histo:  = "  << asym << " +- " << err );

      // dEta/dY distributions: normalize to unity
      normalize(_h_dEta);
      normalize(_h_dY);

      // Build asymm scatters
      for (size_t iM = 0; iM < kNmeas; ++iM) {
        for (size_t iB = 0; iB < kNbins; ++iB) {
          removeAnalysisObject(_h_dEta_asym[iM][iB]);
          removeAnalysisObject(_h_dY_asym[iM][iB]);

          // Only one bin for inclusive measurement
          if ( (iM == kInclMeas) && (iB != 0)) continue;

          calcAsymAndError(_h_dEta_asym[iM][iB], asym, err);
          _h_Acll[iM]->addPoint( (bins[iM][iB+1] + bins[iM][iB])/2., asym, (bins[iM][iB+1] - bins[iM][iB])/2., err);

          calcAsymAndError(_h_dY_asym[iM][iB], asym, err);
          _h_Actt[iM]->addPoint( (bins[iM][iB+1] + bins[iM][iB])/2., asym, (bins[iM][iB+1] - bins[iM][iB])/2., err);
        }
      }
    }

    //@}


  private:

    void calcAsymAndError(Histo1DPtr hist, double& asym, double& err)  {

      int nBins = hist->numBins();
      if (nBins % 2 != 0) {
      	asym = -999; err = -999.;
        return;
      }

      double Nneg  = 0.;
      double Npos  = 0.;
      double dNneg = 0.;
      double dNpos = 0.;
      for (int iB = 0; iB < nBins; ++iB) {
        if (iB < nBins/2) {
          Nneg  += hist->bin(iB).sumW();
          dNneg += hist->bin(iB).sumW2();
        } else {
          Npos  += hist->bin(iB).sumW();
          dNpos += hist->bin(iB).sumW2();
        }
      }
      dNneg = sqrt(dNneg);
      dNpos = sqrt(dNpos);
      asym = (Npos + Nneg) != 0.0 ? (Npos - Nneg) / (Npos + Nneg) : -999.;

      const double Ntot = Npos + Nneg;
      const double Ntot2 = Ntot * Ntot;
      const double dNpos2 = dNpos * dNpos;
      const double dNneg2 = dNneg * dNneg;
      err = Ntot2 != 0. ? 2. * sqrt( (dNneg2 * Npos * Npos + dNpos2 * Nneg * Nneg) / (Ntot2 * Ntot2)) : -999.;
    }


    int getBinID(MeasureType type, double value)  {
      /// @todo Use Rivet binIndex() function
      for (size_t iBin = 0; iBin < kNbins; ++iBin) {
        if (value <= bins[type][iBin+1]) return iBin;
      }
      return -1;
    }


    /// @name Histograms
    //@{
    Histo1DPtr _h_dEta;
    Histo1DPtr _h_dY;
    Scatter2DPtr _h_Actt[kNmeas];
    Scatter2DPtr _h_Acll[kNmeas];
    // Histograms to calculate the asymmetries from
    /// @todo Use /TMP histos?
    Histo1DPtr _h_dEta_asym[kNmeas][2];
    Histo1DPtr _h_dY_asym  [kNmeas][2];
    //@}

    // Not-scaled histos
    Histo1DPtr _h_dEta_notscaled, _h_dY_notscaled;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1449082);
}
