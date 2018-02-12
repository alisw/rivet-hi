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


  /// @brief ATLAS 2016 0-lepton SUSY search, from 13/fb ICHEP'16 CONF note
  class ATLAS_2016_CONF_2016_078 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_CONF_2016_078);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 3.2);
      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_ATLAS_RUN2, //JET_BTAG_ATLAS_RUN2_MV2C10
                          [](const Jet& j) {
                            if (j.abseta() > 2.5) return 0.;
                            return j.bTagged(Cuts::pT > 5*GeV) ? 0.77 : j.cTagged(Cuts::pT > 5*GeV) ? 1/6. : 1/134.;
                          }), "RecoJets");

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
      _h_2j_0800 = bookCounter("2j-0800");
      _h_2j_1200 = bookCounter("2j-1200");
      _h_2j_1600 = bookCounter("2j-1600");
      _h_2j_2000 = bookCounter("2j-2000");
      _h_3j_1200 = bookCounter("2j-2000");
      _h_4j_1000 = bookCounter("4j-1000");
      _h_4j_1400 = bookCounter("4j-1400");
      _h_4j_1800 = bookCounter("4j-1800");
      _h_4j_2200 = bookCounter("4j-2200");
      _h_4j_2600 = bookCounter("4j-2600");
      _h_5j_1400 = bookCounter("5j-1400");
      _h_6j_1800 = bookCounter("6j-1800");
      _h_6j_2200 = bookCounter("6j-2200");


      // Book cut-flows
      const vector<string> cuts23j = {"Pre-sel+MET+pT1+meff", "Njet", "Dphi_min(j123,MET)", "Dphi_min(j4+,MET)", "pT2", "eta_j12", "MET/sqrtHT", "m_eff(incl)"};
      _flows.addCutflow("2j-0800", cuts23j);
      _flows.addCutflow("2j-1200", cuts23j);
      _flows.addCutflow("2j-1600", cuts23j);
      _flows.addCutflow("2j-2000", cuts23j);
      _flows.addCutflow("3j-1200", cuts23j);
      const vector<string> cuts456j = {"Pre-sel+MET+pT1+meff", "Njet", "Dphi_min(j123,MET)", "Dphi_min(j4+,MET)", "pT4", "eta_j1234", "Aplanarity", "MET/m_eff(Nj)", "m_eff(incl)"};
      _flows.addCutflow("4j-1000", cuts456j);
      _flows.addCutflow("4j-1400", cuts456j);
      _flows.addCutflow("4j-1800", cuts456j);
      _flows.addCutflow("4j-2200", cuts456j);
      _flows.addCutflow("4j-2600", cuts456j);
      _flows.addCutflow("5j-1400", cuts456j);
      _flows.addCutflow("6j-1800", cuts456j);
      _flows.addCutflow("6j-2200", cuts456j);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      _flows.fillinit();

      // Same MET cut for all signal regions
      const Vector3 vmet = -apply<SmearedMET>(event, "RecoMET").vectorEt();
      const double met = vmet.mod();
      if (met < 250*GeV) vetoEvent;

      // Get baseline electrons, muons, and jets
      Particles elecs = apply<ParticleFinder>(event, "RecoElectrons").particles(Cuts::pT > 10*GeV);
      Particles muons = apply<ParticleFinder>(event, "RecoMuons").particles(Cuts::pT > 10*GeV);
      Jets jets = apply<JetAlg>(event, "RecoJets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.8); ///< @todo Pile-up subtraction

      // Jet/electron/muons overlap removal and selection
      // Remove electrons within dR = 0.2 of a b-tagged jet
      for (const Jet& j : jets)
        if (j.abseta() < 2.5 && j.pT() > 50*GeV && j.bTagged(Cuts::pT > 5*GeV))
          ifilter_discard(elecs, deltaRLess(j, 0.2, RAPIDITY));
      // Remove any |eta| < 2.8 jet within dR = 0.2 of a remaining electron
      for (const Particle& e : elecs)
        ifilter_discard(jets, deltaRLess(e, 0.2, RAPIDITY));
      // Remove any electron with dR in [0.2, 0.4] of a remaining jet
      for (const Jet& j : jets)
        ifilter_discard(elecs, [&](const Particle& e) { return inRange(deltaR(e,j, RAPIDITY), 0.2, 0.4); });
      // Remove any muon with dR close to a remaining jet, via a functional form
      for (const Jet& j : jets)
        ifilter_discard(muons, [&](const Particle& m) { return deltaR(m,j, RAPIDITY) < min(0.4, 0.04 + 10*GeV/m.pT()); });
      // Remove any |eta| < 2.8 jet within dR = 0.2 of a remaining muon if track conditions are met
      for (const Particle& m : muons)
        /// @todo Add track efficiency random filtering
        ifilter_discard(jets, [&](const Jet& j) {
            if (deltaR(j,m, RAPIDITY) > 0.2) return false;
            const Particles trks = j.particles(Cuts::abscharge > 0 && Cuts::pT > 0.5*GeV);
            return trks.size() < 3 || (m.pT() > 2*j.pT() && m.pT() > 0.7*sum(trks, pT, 0.0));
          });
      // Loose electron selection
      ifilter_select(elecs, ParticleEffFilter(ELECTRON_IDEFF_ATLAS_RUN2_LOOSE));

      // Veto the event if there are any remaining baseline leptons
      if (!elecs.empty()) vetoEvent;
      if (!muons.empty()) vetoEvent;

      // Passed presel & MET
      _flows.fill(1);

      // Get jets and their pTs
      const Jets jets20 = jets;
      const Jets jets50 = filterBy(jets, Cuts::pT > 50*GeV);
      const size_t njets50 = jets50.size(), njets20 = jets20.size();
      if (jets50.size() < 2) vetoEvent;
      vector<double> jetpts20, jetpts50;
      transform(jets20, jetpts20, pT);
      transform(jets50, jetpts50, pT);

      // Construct multi-jet observables
      const double ht = sum(jetpts20, 0.0);
      const double met_sqrtHT = met / sqrt(ht);
      const double meff_incl = sum(jetpts50, met);
      const double meff_4 = (njets50 >= 4) ? sum(head(jetpts50, 4), met) : -1;
      const double meff_5 = (njets50 >= 5) ? sum(head(jetpts50, 5), met) : -1;
      const double meff_6 = (njets50 >= 6) ? sum(head(jetpts50, 6), met) : -1;
      const double met_meff_4 = met / meff_4;
      const double met_meff_5 = met / meff_5;
      const double met_meff_6 = met / meff_6;

      // Jet |eta|s
      vector<double> jetetas20; transform(jets20, jetetas20, abseta);
      const double etamax_2 = (njets20 >= 2) ? max(head(jetetas20, 2)) : -1;
      const double etamax_4 = (njets20 >= 4) ? max(head(jetetas20, 4)) : -1;
      const double etamax_6 = (njets20 >= 6) ? max(head(jetetas20, 6)) : -1;

      // Get dphis between MET and jets
      vector<double> dphimets50; transform(jets50, dphimets50, deltaPhiWRT(vmet));
      const vector<double> dphimets50_123 = head(dphimets50, 3);
      const vector<double> dphimets50_more = tail(dphimets50, -3);
      const double dphimin_123 = !dphimets50_123.empty() ? min(dphimets50_123) : -1;
      const double dphimin_more = !dphimets50_more.empty() ? min(dphimets50_more) : -1;

      // Jet aplanarity
      Sphericity sph; sph.calc(jets50);
      const double aplanarity = sph.aplanarity();


      //////////////////


      const double w = event.weight();

      // 2 jet regions
      if (dphimin_123 > 0.8 && dphimin_more > 0.4) {
        if (jetpts50[1] > 200*GeV && etamax_2 < 0.8) { //< implicit pT[0] cut
          if (met_sqrtHT > 14*sqrt(GeV) && meff_incl > 800*GeV) _h_2j_0800->fill(w);
        }
        if (jetpts50[1] > 250*GeV && etamax_2 < 1.2) { //< implicit pT[0] cut
          if (met_sqrtHT > 16*sqrt(GeV) && meff_incl > 1200*GeV) _h_2j_1200->fill(w);
          if (met_sqrtHT > 18*sqrt(GeV) && meff_incl > 1600*GeV) _h_2j_1600->fill(w);
          if (met_sqrtHT > 20*sqrt(GeV) && meff_incl > 2000*GeV) _h_2j_2000->fill(w);
        }
      }

      // 3 jet region
      if (njets50 >= 3 && dphimin_123 > 0.4 && dphimin_more > 0.2) {
        if (jetpts50[0] > 600*GeV && jetpts50[2] > 50*GeV) { //< implicit pT[1] cut
          if (met_sqrtHT > 16*sqrt(GeV) && meff_incl > 1200*GeV) _h_3j_1200->fill(w);
        }
      }

      // 4 jet regions (note implicit pT[1,2] cuts)
      if (njets50 >= 4 && dphimin_123 > 0.4 && dphimin_more > 0.4 && jetpts50[0] > 200*GeV && aplanarity > 0.04) {
        if (jetpts50[3] > 100*GeV && etamax_4 < 1.2 && met_meff_4 > 0.25*sqrt(GeV) && meff_incl > 1000*GeV) _h_4j_1000->fill(w);
        if (jetpts50[3] > 100*GeV && etamax_4 < 2.0 && met_meff_4 > 0.25*sqrt(GeV) && meff_incl > 1400*GeV) _h_4j_1400->fill(w);
        if (jetpts50[3] > 100*GeV && etamax_4 < 2.0 && met_meff_4 > 0.20*sqrt(GeV) && meff_incl > 1800*GeV) _h_4j_1800->fill(w);
        if (jetpts50[3] > 150*GeV && etamax_4 < 2.0 && met_meff_4 > 0.20*sqrt(GeV) && meff_incl > 2200*GeV) _h_4j_2200->fill(w);
        if (jetpts50[3] > 150*GeV &&                   met_meff_4 > 0.20*sqrt(GeV) && meff_incl > 2600*GeV) _h_4j_2600->fill(w);
      }

      // 5 jet region (note implicit pT[1,2,3] cuts)
      if (njets50 >= 5 && dphimin_123 > 0.4 && dphimin_more > 0.2 && jetpts50[0] > 500*GeV) {
        if (jetpts50[4] > 50*GeV && met_meff_5 > 0.3*sqrt(GeV) && meff_incl > 1400*GeV) _h_5j_1400->fill(w);
      }

      // 6 jet regions (note implicit pT[1,2,3,4] cuts)
      if (njets50 >= 6 && dphimin_123 > 0.4 && dphimin_more > 0.2 && jetpts50[0] > 200*GeV && aplanarity > 0.08) {
        if (jetpts50[5] >  50*GeV && etamax_6 < 2.0 && met_meff_6*sqrt(GeV) > 0.20 && meff_incl > 1800*GeV) _h_6j_1800->fill(w);
        if (jetpts50[5] > 100*GeV &&                   met_meff_6*sqrt(GeV) > 0.15 && meff_incl > 2200*GeV) _h_6j_2200->fill(w);
      }

      // Cutflows
      _flows["2j-0800"].filltail({true, dphimin_123 > 0.8, dphimin_more > 0.4, jetpts50[1] > 200*GeV, etamax_2 < 0.8, met_sqrtHT > 14*sqrt(GeV), meff_incl >  800*GeV});
      _flows["2j-1200"].filltail({true, dphimin_123 > 0.8, dphimin_more > 0.4, jetpts50[1] > 250*GeV, etamax_2 < 1.2, met_sqrtHT > 16*sqrt(GeV), meff_incl > 1200*GeV});
      _flows["2j-1600"].filltail({true, dphimin_123 > 0.8, dphimin_more > 0.4, jetpts50[1] > 250*GeV, etamax_2 < 1.2, met_sqrtHT > 18*sqrt(GeV), meff_incl > 1600*GeV});
      _flows["2j-2000"].filltail({true, dphimin_123 > 0.8, dphimin_more > 0.4, jetpts50[1] > 250*GeV, etamax_2 < 1.2, met_sqrtHT > 20*sqrt(GeV), meff_incl > 2000*GeV});
      _flows["3j-1200"].filltail({njets50 >= 3, dphimin_123 > 0.4, dphimin_more > 0.2, jetpts50[0] > 600*GeV && jetpts50[2] > 50*GeV, true, met_sqrtHT > 16*sqrt(GeV), meff_incl > 1200*GeV});
      _flows["4j-1000"].filltail({njets50 >= 4, dphimin_123 > 0.4, dphimin_more > 0.4, jetpts50[0] > 200*GeV && jetpts50[3] > 100*GeV, etamax_4 < 1.2, aplanarity > 0.04, met_meff_4 > 0.25*sqrt(GeV), meff_incl > 1000*GeV});
      _flows["4j-1400"].filltail({njets50 >= 4, dphimin_123 > 0.4, dphimin_more > 0.4, jetpts50[0] > 200*GeV && jetpts50[3] > 100*GeV, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.25*sqrt(GeV), meff_incl > 1400*GeV});
      _flows["4j-1800"].filltail({njets50 >= 4, dphimin_123 > 0.4, dphimin_more > 0.4, jetpts50[0] > 200*GeV && jetpts50[3] > 100*GeV, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.20*sqrt(GeV), meff_incl > 1800*GeV});
      _flows["4j-2200"].filltail({njets50 >= 4, dphimin_123 > 0.4, dphimin_more > 0.4, jetpts50[0] > 200*GeV && jetpts50[3] > 150*GeV, etamax_4 < 2.0, aplanarity > 0.04, met_meff_4 > 0.20*sqrt(GeV), meff_incl > 2200*GeV});
      _flows["4j-2600"].filltail({njets50 >= 4, dphimin_123 > 0.4, dphimin_more > 0.4, jetpts50[0] > 200*GeV && jetpts50[3] > 150*GeV, true,           aplanarity > 0.04, met_meff_4 > 0.20*sqrt(GeV), meff_incl > 2600*GeV});
      _flows["5j-1400"].filltail({njets50 >= 5, dphimin_123 > 0.4, dphimin_more > 0.2, jetpts50[0] > 500*GeV && jetpts50[4] > 50*GeV, true, true, met_meff_5 > 0.3*sqrt(GeV), meff_incl > 1400*GeV});
      _flows["6j-1800"].filltail({njets50 >= 6, dphimin_123 > 0.4, dphimin_more > 0.2, jetpts50[0] > 200*GeV && jetpts50[5] >  50*GeV, etamax_6 < 2.0, aplanarity > 0.08, met_meff_6 > 0.20*sqrt(GeV), meff_incl > 1800*GeV});
      _flows["6j-2200"].filltail({njets50 >= 6, dphimin_123 > 0.4, dphimin_more > 0.2, jetpts50[0] > 200*GeV && jetpts50[5] > 100*GeV, true,           aplanarity > 0.08, met_meff_6 > 0.15*sqrt(GeV), meff_incl > 2200*GeV});

    }


    /// Normalise counters after the run
    void finalize() {

      const double sf = 13.3*crossSection()/femtobarn/sumOfWeights();
      scale({_h_2j_0800, _h_2j_1200, _h_2j_1600, _h_2j_2000}, sf);
      scale( _h_3j_1200, sf);
      scale({_h_4j_1000, _h_4j_1400, _h_4j_1800, _h_4j_2200, _h_4j_2600}, sf);
      scale( _h_5j_1400, sf);
      scale({_h_6j_1800, _h_6j_2200}, sf);

      _flows.scale(sf);
      MSG_INFO("CUTFLOWS:\n\n" << _flows);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    CounterPtr _h_2j_0800, _h_2j_1200, _h_2j_1600, _h_2j_2000, _h_3j_1200;
    CounterPtr _h_4j_1000, _h_4j_1400, _h_4j_1800, _h_4j_2200, _h_4j_2600;
    CounterPtr _h_5j_1400, _h_6j_1800, _h_6j_2200;
    //@}

    /// Cut-flows
    Cutflows _flows;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_CONF_2016_078);


}
