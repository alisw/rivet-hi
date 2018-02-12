// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief $t\bar{t}$ + jets at 13 TeV
  class ATLAS_2017_I1495243 : public Analysis {
  public:


    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1495243);


    void init() {

      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT > 1.0*MeV;
      Cut eta_lep = Cuts::abseta < 2.5;

      // Collect final state particles
      FinalState FS(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState photons(FS);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      IdentifiedFinalState el_id(FS);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(false);
      DressedLeptons dressedelectrons(photons, electrons, 0.1, Cuts::abseta< 2.5 && Cuts::pT > 25.0*GeV, true);
      addProjection(dressedelectrons, "electrons");
      DressedLeptons fulldressedelectrons(photons, electrons, 0.1, eta_full, true);

      // Projection to find the muons
      IdentifiedFinalState mu_id(FS);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(false);
      DressedLeptons dressedmuons(photons, muons, 0.1, Cuts::abseta < 2.5 && Cuts::pT > 25.0*GeV, true);
      addProjection(dressedmuons, "muons");
      DressedLeptons fulldressedmuons(photons, muons, 0.1, eta_full, true);

      // Projection to find neutrinos to exclude from jets
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(false);

      // Jet clustering
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(fulldressedelectrons);
      vfs.addVetoOnThisFinalState(fulldressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      addProjection(jets, "jets");

      // Book Histograms
      _h["bjet_pt"]  = bookHisto1D(5,1,1);
      _h["2bjet_pt"] = bookHisto1D(6,1,1);
      _h["ljet_pt"]  = bookHisto1D(7,1,1);

      for (size_t i = 0; i < 4; ++i) {
        _h["njet"  + to_str(i)] = bookHisto1D(i+1, 1, 1);
        _h["Q0"    + to_str(i)] = bookHisto1D("_Q0"    + to_str(i+ 7), refData((i>1?"d":"d0") + to_str(i+ 8) + "-x01-y01"));
        _h["MQ0"   + to_str(i)] = bookHisto1D("_MQ0"   + to_str(i+12), refData("d" + to_str(i+12) + "-x01-y01"));
        _h["Qsum"  + to_str(i)] = bookHisto1D("_Qsum"  + to_str(i+16), refData("d" + to_str(i+16) + "-x01-y01"));
        _h["MQsum" + to_str(i)] = bookHisto1D("_MQsum" + to_str(i+20), refData("d" + to_str(i+20) + "-x01-y01"));
        _s["gapFracQ0"    + to_str(i)] = bookScatter2D( 8+i, 1 ,1, true);
        _s["gapFracMQ0"   + to_str(i)] = bookScatter2D(12+i, 1, 1, true);
        _s["gapFracQsum"  + to_str(i)] = bookScatter2D(16+i, 1, 1, true);
        _s["gapFracMQsum" + to_str(i)] = bookScatter2D(20+i, 1, 1, true);
      }
    }


    void analyze(const Event& event) {

      const double weight = event.weight();

      // Get the selected objects, using the projections.
      Jets all_jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      const vector<DressedLepton> electrons = filter_discard(apply<DressedLeptons>(event, "electrons").dressedLeptons(),
        [&](const DressedLepton &e) {
          return any(all_jets, deltaRLess(e, 0.4));
        });

      const vector<DressedLepton> muons = filter_discard(apply<DressedLeptons>(event, "muons").dressedLeptons(),
        [&](const DressedLepton &m) {
          return any(all_jets, deltaRLess(m, 0.4));
        });

      if (electrons.size() != 1 || muons.size() != 1)  vetoEvent;
      if (electrons[0].charge() == muons[0].charge())  vetoEvent;

      Jets bjets, extrajets;
      for (Jet j : all_jets) {
        size_t b_tagged = j.bTags(Cuts::pT > 5*GeV).size();
        if (bjets.size() < 2 && b_tagged)  bjets += j;
        else  extrajets += j;
      }
      if (bjets.size() < 2)  vetoEvent;

      double bjetpt = bjets[0].pt();
      if (bjetpt > 250*GeV)  bjetpt = 275*GeV;
      _h["bjet_pt"]->fill(bjetpt, weight);

      double b2jetpt = bjets[1].pt();
      if (b2jetpt > 150*GeV)  b2jetpt = 175*GeV;
      _h["2bjet_pt"]->fill(b2jetpt, weight);

      if (extrajets.size()) {
        double ljetpt = extrajets[0].pt();
        if (ljetpt > 250*GeV)  ljetpt = 275*GeV;
        _h["ljet_pt"]->fill(ljetpt, weight);
      }

      double Memubb = (electrons[0].momentum() + muons[0].momentum() + bjets[0].momentum() + bjets[1].momentum()).mass();
      vector<double> leadpt = { 0., 0., 0., 0. }, ptsum = { 0., 0., 0., 0. };
      vector<size_t> njetcount = { 0, 0, 0, 0 };
      for (size_t i = 0; i < extrajets.size(); ++i) {
        double absrap = extrajets[i].absrap(), pt = extrajets[i].pT();
        if (pt > 25*GeV)  ++njetcount[0];
        if (pt > 40*GeV)  ++njetcount[1];
        if (pt > 60*GeV)  ++njetcount[2];
        if (pt > 80*GeV)  ++njetcount[3];

        if (absrap < 0.8 && pt > leadpt[0])  leadpt[0] = pt;
        else if (absrap > 0.8 && absrap < 1.5 && pt > leadpt[1])  leadpt[1] = pt;
        else if (absrap > 1.5 && absrap < 2.1 && pt > leadpt[2])  leadpt[2] = pt;
        if (absrap < 2.1 && pt > leadpt[3])  leadpt[3] = pt;

        if (absrap < 0.8)  ptsum[0] += pt;
        else if (absrap > 0.8 && absrap < 1.5)  ptsum[1] += pt;
        else if (absrap > 1.5 && absrap < 2.1)  ptsum[2] += pt;
        if (absrap < 2.1)  ptsum[3] += pt;
      }


      for (size_t i = 0; i < 4; ++i) {
        size_t cutoff = i? 3 : 4;
        if (njetcount[i] > cutoff)  njetcount[i] = cutoff;
        _h["njet" + to_str(i)]->fill(njetcount[i], weight);

        if (leadpt[i] > 305*GeV)  leadpt[i] = 305*GeV;
        _h["Q0" + to_str(i)]->fill(leadpt[i], weight);

        if (ptsum[i] > 505*GeV)  ptsum[i] = 505*GeV;
        _h["Qsum" + to_str(i)]->fill(ptsum[i], weight);
      }


      for (size_t i = 0; i < 4; ++i) {
        if (i == 0 && !(Memubb < 300*GeV))  continue;
        if (i == 1 && !(Memubb > 300*GeV && Memubb < 425*GeV))  continue;
        if (i == 2 && !(Memubb > 425*GeV && Memubb < 600*GeV))  continue;
        if (i == 3 && !(Memubb > 600*GeV))  continue;
        _h["MQ0"   + to_str(i)]->fill(leadpt[3], weight);
        _h["MQsum" + to_str(i)]->fill(ptsum[3],  weight);
      }
    }


    void constructGapFraction(Scatter2DPtr out, Histo1DPtr in) {
      bool hasWeights = in->effNumEntries() != in->numEntries();
      double denW  = in->sumW();
      double denW2 = in->sumW2();
      size_t nEnd = out->numPoints();

      for (size_t i = 0; i < nEnd; ++i) {
          double numW = in->sumW(), numW2 = in->sumW2();
          for (size_t j = i; j < nEnd; ++j) {
            numW  -= in->bin(j).sumW();
            numW2 -= in->bin(j).sumW2();
          }
          double yval = safediv(numW, denW);
          double yerr = sqrt(safediv(yval * (1 - yval), denW));
          if (hasWeights) { // use F. James's approximation for weighted events
            yerr = sqrt( safediv((1 - 2 * yval) * numW2 + yval * yval * denW2, denW * denW) );
          }
          out->point(i).setY(yval, yerr);
      }
    }


    void finalize() {

      // Build gap fraction plots
      for (size_t i = 0; i < 4; ++i) {
        constructGapFraction(_s["gapFracQ0"    + to_str(i)], _h["Q0"    + to_str(i)]);
        constructGapFraction(_s["gapFracMQ0"   + to_str(i)], _h["MQ0"   + to_str(i)]);
        constructGapFraction(_s["gapFracQsum"  + to_str(i)], _h["Qsum"  + to_str(i)]);
        constructGapFraction(_s["gapFracMQsum" + to_str(i)], _h["MQsum" + to_str(i)]);
      }

      // Normalize to cross-section
      for (map<string, Histo1DPtr>::iterator hit = _h.begin(); hit != _h.end(); ++hit) {
        if (hit->first.find("jet") != string::npos)  normalize(hit->second);
      }
    }


  private:

    /// @name Histogram helper functions
    map<string, Histo1DPtr> _h;
    map<string, Scatter2DPtr> _s;
  };


  // Declare the class as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1495243);


}
