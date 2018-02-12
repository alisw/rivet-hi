// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include <bitset>

namespace Rivet {


  /// @brief ATLAS 7 TeV jets in ttbar events analysis
  ///
  /// @author W. H. Bell <W.Bell@cern.ch>
  /// @author A. Grohsjean <alexander.grohsjean@desy.de>
  class ATLAS_2014_I1304688 : public Analysis {
  public:

    ATLAS_2014_I1304688():
      Analysis("ATLAS_2014_I1304688"),
      _jet_ntag(0),
      _met_et(0.),
      _met_phi(0.),
      _hMap(),
      //_chanLimit(3),
      _histLimit(6)
    {   }


    void init() {
      // Eta ranges
      /// @todo 1 MeV? Really?
      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;
      Cut eta_lep = Cuts::abseta < 2.5;

      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      declare(electrons, "electrons");
      DressedLeptons dressedelectrons(photons, electrons, 0.1, eta_lep && Cuts::pT > 25*GeV, true);
      declare(dressedelectrons, "dressedelectrons");
      DressedLeptons vetodressedelectrons(photons, electrons, 0.1, eta_lep && Cuts::pT >= 15*GeV, true);
      declare(vetodressedelectrons, "vetodressedelectrons");
      DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full, true);
      declare(ewdressedelectrons, "ewdressedelectrons");

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      declare(muons, "muons");
      vector<pair<double, double> > eta_muon;
      DressedLeptons dressedmuons(photons, muons, 0.1, eta_lep && Cuts::pT >= 25*GeV, true);
      declare(dressedmuons, "dressedmuons");
      DressedLeptons vetodressedmuons(photons, muons, 0.1, eta_lep && Cuts::pT >= 15*GeV, true);
      declare(vetodressedmuons, "vetodressedmuons");
      DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full, true);
      declare(ewdressedmuons, "ewdressedmuons");

      // Projection to find neutrinos and produce MET
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      declare(neutrinos, "neutrinos");

      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(ewdressedelectrons);
      vfs.addVetoOnThisFinalState(ewdressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");

      // Book histograms
      for (unsigned int ihist = 0; ihist < _histLimit ; ihist++) {
        const unsigned int threshLimit = _thresholdLimit(ihist);
        for (unsigned int ithres = 0; ithres < threshLimit; ithres++) {
          _histogram(ihist, ithres); // Create all histograms
        }
      }
    }


    void analyze(const Event& event) {

      // Get the selected objects, using the projections.
      _dressedelectrons = sortByPt(apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons());
      _vetodressedelectrons = apply<DressedLeptons>(event, "vetodressedelectrons").dressedLeptons();

      _dressedmuons = sortByPt(apply<DressedLeptons>(event, "dressedmuons").dressedLeptons());
      _vetodressedmuons = apply<DressedLeptons>(event, "vetodressedmuons").dressedLeptons();

      _neutrinos = apply<PromptFinalState>(event, "neutrinos").particlesByPt();

      _jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);


      // Calculate the missing ET, using the prompt neutrinos only (really?)
      /// @todo Why not use MissingMomentum?
      FourMomentum pmet;
      foreach (const Particle& p, _neutrinos) pmet += p.momentum();
      _met_et = pmet.pT();
      _met_phi = pmet.phi();

      // Check overlap of jets/leptons.
      unsigned int i,j;
      _jet_ntag = 0;
      _overlap = false;
      for (i = 0; i < _jets.size(); i++) {
        const Jet& jet = _jets[i];
        // If dR(el,jet) < 0.4 skip the event
        foreach (const DressedLepton& el, _dressedelectrons) {
          if (deltaR(jet, el) < 0.4) _overlap = true;
        }
        // If dR(mu,jet) < 0.4 skip the event
        foreach (const DressedLepton& mu, _dressedmuons) {
          if (deltaR(jet, mu) < 0.4) _overlap = true;
        }
        // If dR(jet,jet) < 0.5 skip the event
        for (j = 0; j < _jets.size(); j++) {
          const Jet& jet2 = _jets[j];
          if (i == j) continue; // skip the diagonal
          if (deltaR(jet, jet2) < 0.5) _overlap = true;
        }
        // Count the number of b-tags
        if (!jet.bTags().empty()) _jet_ntag += 1;
      }

      // Evaluate basic event selection
      unsigned int ejets_bits = 0, mujets_bits = 0;
      bool pass_ejets = _ejets(ejets_bits);
      bool pass_mujets = _mujets(mujets_bits);

      // Remove events with object overlap
      if (_overlap) vetoEvent;
      // basic event selection requirements
      if (!pass_ejets && !pass_mujets) vetoEvent;

      // Check if the additional pT threshold requirements are passed
      bool pass_jetPt = _additionalJetCuts();

      // Count the jet multiplicity for 25, 40, 60 and 80GeV
      unsigned int ithres, jet_n[4];
      for (ithres = 0; ithres < 4; ithres++) {
        jet_n[ithres] = _countJets(ithres);
      }

      // Fill histograms
      const double weight = event.weight();
      for (unsigned int ihist = 0; ihist < 6; ihist++) {
        if (ihist > 0 && !pass_jetPt) continue; // additional pT threshold cuts for pT plots
        unsigned int threshLimit = _thresholdLimit(ihist);
        for (ithres = 0; ithres < threshLimit; ithres++) {
          if (jet_n[ithres] < 3) continue; // 3 or more jets for ljets
          // Fill
          if (ihist == 0) _histogram(ihist, ithres)->fill(jet_n[ithres], weight); // njets
          else if (ihist == 1) _histogram(ihist, ithres)->fill(_jets[0].pT(), weight); // leading jet pT
          else if (ihist == 2) _histogram(ihist, ithres)->fill(_jets[1].pT(), weight); // 2nd jet pT
          else if (ihist == 3 && jet_n[ithres] >= 3) _histogram(ihist, ithres)->fill(_jets[2].pT(), weight); // 3rd jet pT
          else if (ihist == 4 && jet_n[ithres] >= 4) _histogram(ihist, ithres)->fill(_jets[3].pT(), weight); // 4th jet pT
          else if (ihist == 5 && jet_n[ithres] >= 5) _histogram(ihist, ithres)->fill(_jets[4].pT(), weight); // 5th jet pT
        }
      }
    }


    void finalize() {
      // Normalize to cross-section
      const double norm = crossSection()/sumOfWeights();
      typedef map<unsigned int, Histo1DPtr>::value_type IDtoHisto1DPtr; ///< @todo Remove when C++11 allowed
      foreach (IDtoHisto1DPtr ihpair, _hMap) scale(ihpair.second, norm); ///< @todo Use normalize(ihpair.second, crossSection())
      // Calc averages
      for (unsigned int ihist = 0; ihist < _histLimit ; ihist++) {
        unsigned int threshLimit = _thresholdLimit(ihist);
        for (unsigned int ithres = 0; ithres < threshLimit; ithres++) {
          scale(_histogram(ihist, ithres), 0.5);
        }
      }
    }



  private:


    /// @name Cut helper functions
    //@{

    // Event selection functions
    bool _ejets(unsigned int& cutBits) {
      // 1. More than zero good electrons
      cutBits += 1; if (_dressedelectrons.size() == 0) return false;
      // 2. No additional electrons passing the veto selection
      cutBits += 1 << 1; if (_vetodressedelectrons.size() > 1) return false;
      // 3. No muons passing the veto selection
      cutBits += 1 << 2; if (_vetodressedmuons.size() > 0) return false;
      // 4. total neutrino pT > 30 GeV
      cutBits += 1 << 3; if (_met_et <= 30.0*GeV) return false;
      // 5. MTW > 35 GeV
      cutBits += 1 << 4;
      if (_transMass(_dressedelectrons[0].pT(), _dressedelectrons[0].phi(), _met_et, _met_phi) <= 35*GeV) return false;
      // 6. At least one b-tagged jet
      cutBits += 1 << 5; if (_jet_ntag < 1) return false;
      // 7. At least three good jets
      cutBits += 1 << 6; if (_jets.size() < 3) return false;
      cutBits += 1 << 7;
      return true;
    }

    bool _mujets(unsigned int& cutBits) {
      // 1. More than zero good muons
      cutBits += 1; if (_dressedmuons.size() == 0) return false;
      // 2. No additional muons passing the veto selection
      cutBits += 1 << 1; if (_vetodressedmuons.size() > 1) return false;
      // 3. No electrons passing the veto selection
      cutBits += 1 << 2; if (_vetodressedelectrons.size() > 0) return false;
      // 4. total neutrino pT > 30 GeV
      cutBits += 1 << 3; if (_met_et <= 30*GeV) return false;
      // 5. MTW > 35 GeV
      cutBits += 1 << 4;
      if (_transMass(_dressedmuons[0].pT(), _dressedmuons[0].phi(), _met_et, _met_phi) <= 35*GeV) return false;
      // 6. At least one b-tagged jet
      cutBits += 1 << 5; if (_jet_ntag < 1) return false;
      // 7. At least three good jets
      cutBits += 1 << 6; if (_jets.size() < 3) return false;
      cutBits += 1 << 7;
      return true;
    }

    bool _additionalJetCuts() {
      if (_jets.size() < 2) return false;
      if (_jets[0].pT() <= 50*GeV || _jets[1].pT() <= 35*GeV) return false;
      return true;
    }

    //@}


    /// @name Histogram helper functions
    //@{

    unsigned int _thresholdLimit(unsigned int histId) {
      if (histId == 0) return 4;
      return 1;
    }

    Histo1DPtr _histogram(unsigned int histId, unsigned int thresholdId) {
      assert(histId < _histLimit);
      assert(thresholdId < _thresholdLimit(histId));

      const unsigned int hInd = (histId == 0) ? thresholdId : (_thresholdLimit(0) + (histId-1) + thresholdId);
      if (_hMap.find(hInd) != _hMap.end()) return _hMap[hInd];

      if (histId == 0) _hMap.insert(make_pair(hInd,bookHisto1D(1,thresholdId+1,1)));
      else _hMap.insert(make_pair(hInd,bookHisto1D(2,histId,1)));
      return _hMap[hInd];
    }

    //@}


    /// @name Physics object helper functions
    //@{

    double _transMass(double ptLep, double phiLep, double met, double phiMet) {
      return sqrt(2.0*ptLep*met*(1 - cos(phiLep-phiMet)));
    }

    unsigned int _countJets(unsigned int ithres) {
      if (ithres > 4) assert(0);
      double pTcut[4] = {25.,40.,60.,80.};
      unsigned int i, jet_n = 0;
      for (i = 0; i < _jets.size(); i++) {
        if (_jets[i].pT() > pTcut[ithres]) jet_n++;
      }
      unsigned int ncutoff[4] = {8,7,6,5};
      if (jet_n > ncutoff[ithres]) jet_n = ncutoff[ithres];
      return jet_n;
    }

    //@}


  private:

    /// @name Objects that are used by the event selection decisions
    //@{
    vector<DressedLepton> _dressedelectrons;
    vector<DressedLepton> _vetodressedelectrons;
    vector<DressedLepton> _dressedmuons;
    vector<DressedLepton> _vetodressedmuons;
    Particles _neutrinos;
    Jets _jets;
    unsigned int _jet_ntag;
    /// @todo Why not store the whole MET FourMomentum?
    double _met_et, _met_phi;
    bool _overlap;

    map<unsigned int, Histo1DPtr> _hMap;
    //unsigned int _chanLimit;
    unsigned int _histLimit;
    //@}

  };



  // Declare the class as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1304688);

}
