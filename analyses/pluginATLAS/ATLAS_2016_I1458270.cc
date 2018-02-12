// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief ATLAS 0-lepton SUSY search with 3.2/fb of 13 TeV pp data
  class ATLAS_2016_I1458270 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1458270);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 4.8);
      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_ATLAS_RUN2, JET_BTAG_ATLAS_RUN2_MV2C20), "RecoJets");

      MissingMomentum mm(calofs);
      declare(mm, "TruthMET");
      declare(SmearedMET(mm, MET_SMEAR_ATLAS_RUN2), "RecoMET");

      PromptFinalState es(Cuts::abseta < 2.47 && Cuts::abspid == PID::ELECTRON, true, true);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2), "RecoElectrons");

      PromptFinalState mus(Cuts::abseta < 2.7 && Cuts::abspid == PID::MUON, true);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2), "RecoMuons");


      // Book histograms/counters
      _h_2jl = bookCounter("2jl");
      _h_2jm = bookCounter("2jm");
      _h_2jt = bookCounter("2jt");
      _h_4jt = bookCounter("4jt");
      _h_5j  = bookCounter("5j");
      _h_6jm = bookCounter("6jm");
      _h_6jt = bookCounter("6jt");


      // Book cut-flows
      const vector<string> cuts2j = {"Pre-sel+MET+pT1", "Njet", "Dphi_min(j,MET)", "pT2", "MET/sqrtHT", "m_eff(incl)"};
      _flows.addCutflow("2jl", cuts2j);
      _flows.addCutflow("2jm", cuts2j);
      _flows.addCutflow("2jt", cuts2j);
      const vector<string> cutsXj = {"Pre-sel+MET+pT1", "Njet", "Dphi_min(j,MET)", "pT2", "pT4", "Aplanarity", "MET/m_eff(Nj)", "m_eff(incl)"};
      _flows.addCutflow("4jt", cutsXj);
      _flows.addCutflow("5j",  cutsXj);
      _flows.addCutflow("6jm", cutsXj);
      _flows.addCutflow("6jt", cutsXj);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _flows.fillinit();

      // Same MET cut for all signal regions
      //const Vector3 vmet = -apply<MissingMomentum>(event, "TruthMET").vectorEt();
      const Vector3 vmet = -apply<SmearedMET>(event, "RecoMET").vectorEt();
      const double met = vmet.mod();
      if (met < 200*GeV) vetoEvent;

      // Get baseline electrons, muons, and jets
      Particles elecs = apply<ParticleFinder>(event, "RecoElectrons").particles(Cuts::pT > 10*GeV);
      Particles muons = apply<ParticleFinder>(event, "RecoMuons").particles(Cuts::pT > 10*GeV);
      Jets jets = apply<JetAlg>(event, "RecoJets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.8); ///< @todo Pile-up subtraction

      // Jet/electron/muons overlap removal and selection
      // Remove any |eta| < 2.8 jet within dR = 0.2 of a baseline electron
      for (const Particle& e : elecs)
        ifilter_discard(jets, deltaRLess(e, 0.2, RAPIDITY));
      // Remove any electron or muon with dR < 0.4 of a remaining (Nch > 3) jet
      for (const Jet& j : jets) {
        /// @todo Add track efficiency random filtering
        ifilter_discard(elecs, deltaRLess(j, 0.4, RAPIDITY));
        if (j.particles(Cuts::abscharge > 0 && Cuts::pT > 500*MeV).size() >= 3)
          ifilter_discard(muons, deltaRLess(j, 0.4, RAPIDITY));
      }
      // Discard the softer of any electrons within dR < 0.05
      for (size_t i = 0; i < elecs.size(); ++i) {
        const Particle& e1 = elecs[i];
        /// @todo Would be nice to pass a "tail view" for the filtering, but awkward without range API / iterator guts
        ifilter_discard(elecs, [&](const Particle& e2){ return e2.pT() < e1.pT() && deltaR(e1,e2) < 0.05; });
      }

      // Loose electron selection
      ifilter_select(elecs, ParticleEffFilter(ELECTRON_IDEFF_ATLAS_RUN2_LOOSE));

      // Veto the event if there are any remaining baseline leptons
      if (!elecs.empty()) vetoEvent;
      if (!muons.empty()) vetoEvent;

      // Signal jets have pT > 50 GeV
      const Jets jets50 = filter_select(jets, Cuts::pT > 50*GeV);
      if (jets50.size() < 2) vetoEvent;
      vector<double> jetpts; transform(jets, jetpts, pT);
      vector<double> jetpts50; transform(jets50, jetpts50, pT);
      const double j1pt = jetpts50[0];
      const double j2pt = jetpts50[1];
      if (j1pt < 200*GeV) vetoEvent;

      // Construct multi-jet observables
      const double ht = sum(jetpts, 0.0);
      const double met_sqrt_ht = met / sqrt(ht);
      const double meff_incl = sum(jetpts50, met);

      // Get dphis between MET and jets
      vector<double> dphimets50; transform(jets50, dphimets50, deltaPhiWRT(vmet));
      const double min_dphi_met_3 = min(head(dphimets50, 3));
      MSG_DEBUG(dphimets50 << ", " << min_dphi_met_3);

      // Jet aplanarity
      Sphericity sph; sph.calc(jets);
      const double aplanarity = sph.aplanarity();


      // Fill SR counters
      // 2-jet SRs
      if (_flows["2jl"].filltail({true, true, min_dphi_met_3 > 0.8, j2pt > 200*GeV,
              met_sqrt_ht > 15*sqrt(GeV), meff_incl > 1200*GeV})) _h_2jl->fill(event.weight());
      if (_flows["2jm"].filltail({j1pt > 300*GeV, true, min_dphi_met_3 > 0.4, j2pt > 50*GeV,
              met_sqrt_ht > 15*sqrt(GeV), meff_incl > 1600*GeV})) _h_2jm->fill(event.weight());
      if (_flows["2jt"].filltail({true, true, min_dphi_met_3 > 0.8, j2pt > 200*GeV,
              met_sqrt_ht > 20*sqrt(GeV), meff_incl > 2000*GeV})) _h_2jt->fill(event.weight());

      // Upper multiplicity SRs
      const double j4pt = jets50.size() > 3 ? jetpts50[3] : -1;
      const double j5pt = jets50.size() > 4 ? jetpts50[4] : -1;
      const double j6pt = jets50.size() > 5 ? jetpts50[5] : -1;
      const double meff_4 = jets50.size() > 3 ? sum(head(jetpts50, 4), met) : -1;
      const double meff_5 = jets50.size() > 4 ? meff_4 + jetpts50[4] : -1;
      const double meff_6 = jets50.size() > 5 ? meff_5 + jetpts50[5] : -1;
      const double met_meff_4 = met / meff_4;
      const double met_meff_5 = met / meff_5;
      const double met_meff_6 = met / meff_6;
      const double min_dphi_met_more = jets50.size() > 3 ? min(tail(dphimets50, -3)) : -1;

      if (_flows["4jt"].filltail({true, jets50.size() >= 4, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
              jetpts[1] > 100*GeV, j4pt > 100*GeV, aplanarity > 0.04, met_meff_4 > 0.20, meff_incl > 2200*GeV}))
        _h_4jt->fill(event.weight());
      if (_flows["5j"].filltail({true, jets50.size() >= 5, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
              jetpts[1] > 100*GeV, j4pt > 100*GeV && j5pt > 50*GeV, aplanarity > 0.04, met_meff_5 > 0.25, meff_incl > 1600*GeV}))
        _h_5j->fill(event.weight());
      if (_flows["6jm"].filltail({true, jets50.size() >= 6, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
              jetpts[1] > 100*GeV, j4pt > 100*GeV && j6pt > 50*GeV, aplanarity > 0.04, met_meff_6 > 0.25, meff_incl > 1600*GeV}))
        _h_6jm->fill(event.weight());
      if (_flows["6jt"].filltail({true, jets50.size() >= 6, min_dphi_met_3 > 0.4 && min_dphi_met_more > 0.2,
              jetpts[1] > 100*GeV, j4pt > 100*GeV && j6pt > 50*GeV, aplanarity > 0.04, met_meff_6 > 0.20, meff_incl > 2000*GeV}))
        _h_6jt->fill(event.weight());

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = 3.2*crossSection()/femtobarn/sumOfWeights();
      scale({_h_2jl, _h_2jl, _h_2jl}, sf);
      scale({_h_4jt, _h_5j}, sf);
      scale({_h_6jm, _h_6jt}, sf);

      MSG_INFO("CUTFLOWS:\n\n" << _flows);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    CounterPtr _h_2jl, _h_2jm, _h_2jt;
    CounterPtr _h_4jt, _h_5j;
    CounterPtr _h_6jm, _h_6jt;
    //@}

    /// Cut-flows
    Cutflows _flows;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1458270);


}
