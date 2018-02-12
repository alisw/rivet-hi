// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Tools/Cutflow.hh"
#include <tuple>

namespace Rivet {


  /// CMS search for SUSY with multijet + MET signatures in 36/fb of 13 TeV pp data
  class CMS_2017_I1594909 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1594909);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      VisibleFinalState pfall;
      declare(pfall, "PFAll");
      ChargedFinalState pfchg(Cuts::abseta < 2.5);
      declare(pfchg, "PFChg");

      FastJets jets(FinalState(Cuts::abseta < 4.9), FastJets::ANTIKT, 0.4);
      SmearedJets recojets(jets, JET_SMEAR_CMS_RUN2, [](const Jet& j){ return j.bTagged() ? 0.55 : j.cTagged() ? 0.12 : 0.016; });
      declare(recojets, "Jets");

      FinalState electrons(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5);
      SmearedParticles recoelectrons(electrons, ELECTRON_EFF_CMS_RUN2);
      declare(recoelectrons, "Electrons");

      FinalState muons(Cuts::abspid == PID::MUON && Cuts::abseta < 2.4);
      SmearedParticles recomuons(muons, MUON_EFF_CMS_RUN2);
      declare(recomuons, "Muons");

      VisibleFinalState calofs(Cuts::abseta < 4.9 && Cuts::abspid != PID::MUON);
      MissingMomentum met(calofs);
      SmearedMET recomet(met, MET_SMEAR_CMS_RUN2);
      declare(recomet, "MET");


      // Book counters, into a map of 3 indices since the global index is not obvious to calculate
      size_t i = 0;
      for (int j = 1; j <= 5; ++j) {
        for (int b = 1; b <= 4; ++b) {
          if (j == 1 && b == 4) continue;
          for (int k = 1; k <= 10; ++k) {
            if (j > 3 && (k == 1 || k == 4)) continue;
            stringstream s; s << "count_" << (i+1); // << "_" << j << b << k;
            _counts[make_tuple(j,b,k)] = bookCounter(s.str());
            i += 1;
          }
        }
      }
      MSG_DEBUG("Booked " << i << " signal regions (should be 174)");
      // Aggregate SR counters
      for (size_t i = 0; i < 12; ++i)
        _counts_agg[i] = bookCounter("count_agg_" + toString(i+1));


      // Book cut-flow
      _flow = Cutflow("Presel", {"Njet>=2", "HT>300", "HTmiss>300",
            "Nmuon=0", "Nmuisotrk=0", "Nelec=0", "Nelisotrk=0", "Nhadisotrk=0",
            "dPhi_miss,j1>0.5", "dPhi_miss,j2>0.5", "dPhi_miss,j3>0.3", "dPhi_miss,j4>0.3"
            });

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _flow.fillinit();

      // Find leptons and isolation particles
      const Particles elecs = apply<ParticleFinder>(event, "Electrons").particlesByPt();
      const Particles mus = apply<ParticleFinder>(event, "Muons").particlesByPt();
      const Particles pfall = apply<ParticleFinder>(event, "PFAll").particlesByPt();
      const Particles pfiso = filter_select(pfall, [](const Particle& p){ return p.isHadron() || p.pid() == PID::PHOTON; });

      // Find isolated leptons
      const Particles isoleps = filter_select(elecs+mus, [&](const Particle& l){
          const double dR = l.pT() < 50*GeV ? 0.2 : l.pT() < 200*GeV ? 10*GeV/l.pT() : 0.05;
          const double sumpt = sum(filter_select(pfiso, deltaRLess(l, dR)), pT, 0.0);
          return sumpt/l.pT() < (l.abspid() == PID::ELECTRON ? 0.1 : 0.2); //< different I criteria for e and mu
        });

      // Find other isolated tracks
      const Particles pfchg = apply<ParticleFinder>(event, "PFChg").particlesByPt();
      const Particles isochgs = filter_select(pfchg, [&](const Particle& t){
          if (t.abseta() > 2.4) return false;
          if (any(isoleps, deltaRLess(t, 0.01))) return false; //< don't count isolated leptons here
          const double sumpt = sum(filter_select(pfchg, deltaRLess(t, 0.3)), pT, -t.pT());
          return sumpt/t.pT() < ((t.abspid() == PID::ELECTRON || t.abspid() == PID::MUON) ? 0.2 : 0.1);
        });

      // Find and isolate jets
      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV);
      const Jets isojets = filter_select(jets, Cuts::abseta < 2.4); //< @todo Isolation from leptons?!?
      const int njets = isojets.size();
      const Jets isobjets = filter_select(isojets, hasBTag());
      const int nbjets = isobjets.size();
      MSG_DEBUG("Njets = " << jets.size() << ", Nisojets = " << njets << ", Nbjets = " << nbjets);

      // Calculate HT, HTmiss, and pTmiss quantities
      const double ht = sum(jets, pT, 0.0);
      const Vector3 vhtmiss = -sum(jets, pTvec, Vector3());
      const double htmiss = vhtmiss.perp();
      const Vector3& vptmiss = -apply<SmearedMET>(event, "MET").vectorEt();
      const double ptmiss = vptmiss.perp();
      MSG_DEBUG("HT = " << ht/GeV << " GeV, HTmiss = " << htmiss/GeV << " GeV");


      /////////////////////////////////////
      // Event selection

      // Njet cut
      if (njets < 2) vetoEvent;
      _flow.fill(1);
      // HT cut
      if (ht < 300*GeV) vetoEvent;
      _flow.fill(2);
      // HTmiss cut
      if (htmiss < 300*GeV) vetoEvent;
      _flow.fill(3);

      // Isolated leptons cut
      if (!filter_select(isoleps, Cuts::pT > 10*GeV).empty()) vetoEvent;
      // Isolated tracks cut
      for (const Particle& t : isochgs) {
        const double mT = sqrt(2*t.pT()*ptmiss * (1 - cos(deltaPhi(t, vptmiss))) );
        if (mT < 100*GeV) continue;
        const double pTmax = (t.abspid() == PID::ELECTRON || t.abspid() == PID::MUON) ? 5*GeV : 10*GeV;
        if (t.pT() > pTmax) vetoEvent;
      }
      //
      // // Inefficiently separated version of isolation cuts for detailed cutflow debugging
      // // Muon cut
      // if (!filter_select(isoleps, Cuts::pT > 10*GeV && Cuts::abspid == PID::MUON).empty()) vetoEvent;
      // _flow.fill(4);
      // // Muon isotrk cut
      // for (const Particle& t : filter_select(isochgs, Cuts::abspid == PID::MUON)) {
      //   const double mT = sqrt(2*t.pT()*ptmiss * (1 - cos(deltaPhi(t, vptmiss))) );
      //   if (mT > 100*GeV && t.pT() > 5*GeV) vetoEvent;
      // }
      // _flow.fill(5);
      // // Electron cut
      // if (!filter_select(isoleps, Cuts::pT > 10*GeV && Cuts::abspid == PID::ELECTRON).empty()) vetoEvent;
      // _flow.fill(6);
      // // Electron isotrk cut
      // for (const Particle& t : filter_select(isochgs, Cuts::abspid == PID::ELECTRON)) {
      //   const double mT = sqrt(2*t.pT()*ptmiss * (1 - cos(deltaPhi(t, vptmiss))) );
      //   if (mT > 100*GeV && t.pT() > 5*GeV) vetoEvent;
      // }
      // _flow.fill(7);
      // // Hadron isotrk cut
      // for (const Particle& t : filter_select(isochgs, Cuts::abspid != PID::ELECTRON && Cuts::abspid != PID::MUON)) {
      //   const double mT = sqrt(2*t.pT()*ptmiss * (1 - cos(deltaPhi(t, vptmiss))) );
      //   if (mT > 100*GeV && t.pT() > 10*GeV) vetoEvent;
      // }
      _flow.fill(8);


      // dPhi(jet,HTmiss) cuts
      if (deltaPhi(vhtmiss, isojets[0]) < 0.5) vetoEvent;
      _flow.fill(9);
      if (deltaPhi(vhtmiss, isojets[1]) < 0.5) vetoEvent;
      _flow.fill(10);
      if (njets >= 3 && deltaPhi(vhtmiss, isojets[2]) < 0.3) vetoEvent;
      _flow.fill(11);
      if (njets >= 4 && deltaPhi(vhtmiss, isojets[3]) < 0.3) vetoEvent;
      _flow.fill(12);


      /////////////////////////////////////
      // Find SR index and fill counter

      const double w = event.weight();

      const int idx_j = binIndex(njets, vector<int>{2,3,5,7,9}, true);
      const int idx_b = binIndex(nbjets, vector<int>{0,1,2,3}, true);
      int idx_k = -1;
      if (inRange(htmiss/GeV, 300, 350)) {
        idx_k = ht < 500*GeV ? 1 : ht < 1000*GeV ? 2 : 3;
      } else if (inRange(htmiss/GeV, 350, 500) && ht > 350*GeV) {
        idx_k = ht < 500*GeV ? 4 : ht < 1000*GeV ? 5 : 6;
      } else if (inRange(htmiss/GeV, 500, 750) && ht > 500*GeV) {
        idx_k = ht < 1000*GeV ? 7 : 8;
      } else if (htmiss/GeV > 750 && ht > 750*GeV) {
        idx_k = ht < 1500*GeV ? 9 : 10;
      }

      // Fill via 3-tuple index
      if (idx_j >= 0 && idx_b >= 0 && idx_k >= 0) {
        const auto idx = make_tuple(idx_j+1,idx_b+1,idx_k);
        if (has_key(_counts, idx)) _counts[idx]->fill(w);
      }


      /////////////////////////////////////
      // Aggregate SRs

      // Region Njet Nb-jet HT [GeV] HTmiss [GeV] Parton multiplicity Heavy flavor ? ∆m
      if (njets >= 2 && nbjets == 0 && ht >=  500*GeV && htmiss >= 500*GeV) _counts_agg[0]->fill(w);
      if (njets >= 3 && nbjets == 0 && ht >= 1500*GeV && htmiss >= 750*GeV) _counts_agg[1]->fill(w);
      if (njets >= 5 && nbjets == 0 && ht >=  500*GeV && htmiss >= 500*GeV) _counts_agg[2]->fill(w);
      if (njets >= 5 && nbjets == 0 && ht >= 1500*GeV && htmiss >= 750*GeV) _counts_agg[3]->fill(w);
      if (njets >= 9 && nbjets == 0 && ht >= 1500*GeV && htmiss >= 750*GeV) _counts_agg[4]->fill(w);
      if (njets >= 2 && nbjets >= 2 && ht >=  500*GeV && htmiss >= 500*GeV) _counts_agg[5]->fill(w);
      if (njets >= 3 && nbjets >= 1 && ht >=  750*GeV && htmiss >= 750*GeV) _counts_agg[6]->fill(w);
      if (njets >= 5 && nbjets >= 3 && ht >=  500*GeV && htmiss >= 500*GeV) _counts_agg[7]->fill(w);
      if (njets >= 5 && nbjets >= 2 && ht >= 1500*GeV && htmiss >= 750*GeV) _counts_agg[8]->fill(w);
      if (njets >= 9 && nbjets >= 3 && ht >=  750*GeV && htmiss >= 750*GeV) _counts_agg[9]->fill(w);
      if (njets >= 7 && nbjets >= 1 && ht >=  300*GeV && htmiss >= 300*GeV) _counts_agg[10]->fill(w);
      if (njets >= 5 && nbjets >= 1 && ht >=  750*GeV && htmiss >= 750*GeV) _counts_agg[11]->fill(w);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = 35.9*crossSection()/femtobarn;
      const double sf = norm/sumOfWeights();

      for (auto& idx_cptr : _counts)
        scale(idx_cptr.second, sf);
      for (CounterPtr& cptr : _counts_agg)
        scale(cptr, sf);

      _flow.scale(sf);
      MSG_INFO("CUTFLOWS:\n\n" << _flow);

    }

    //@}


  private:

    Cutflow _flow;

    map<tuple<int,int,int>, CounterPtr> _counts;
    CounterPtr _counts_agg[12];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1594909);


}
