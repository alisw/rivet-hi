// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalStates.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/Smearing.hh"
#include <bitset>

namespace Rivet {


  /// @brief Add a short analysis description here
  class CMS_2013_I1223519 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2013_I1223519);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 5.0);
      declare(calofs, "Clusters");

      MissingMomentum mm(calofs);
      declare(mm, "TruthMET");
      declare(SmearedMET(mm, MET_SMEAR_CMS_RUN2), "MET");

      FastJets fj(calofs, FastJets::ANTIKT, 0.5);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_CMS_RUN2, [](const Jet& j) {
            if (j.abseta() > 2.4) return 0.;
            return j.bTagged() ? 0.65 : 0.01; }), "Jets"); ///< @note Charm mistag and exact b-tag eff not given

      FinalState ys(Cuts::abspid == PID::PHOTON && Cuts::abseta < 5.0);
      declare(ys, "TruthPhotons");
      declare(SmearedParticles(ys, PHOTON_EFF_CMS_RUN2 /*, PHOTON_SMEAR_CMS_RUN2 */), "Photons");

      FinalState es(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_CMS_RUN2, ELECTRON_SMEAR_CMS_RUN2), "Electrons");

      FinalState mus(Cuts::abspid == PID::MUON && Cuts::abseta < 2.4);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_CMS_RUN2, MUON_SMEAR_CMS_RUN2), "Muons");

      ChargedFinalState cfs(Cuts::abseta < 2.5);
      declare(cfs, "TruthTracks");
      declare(SmearedParticles(cfs, TRK_EFF_CMS_RUN2), "Tracks");


      // Book histograms
      _h_alphaT23 = bookHisto1D("alphaT23", 15, 0, 3);
      _h_alphaT4 = bookHisto1D("alphaT4", 15, 0, 3);
      /// @todo Add HT histograms

      // Book counters
      _h_srcounters.resize(8*7 + 3);
      for (size_t inj = 0; inj < 2; ++inj) {
        const size_t njmax = inj + 3;
        for (size_t nb = 0; nb < njmax; ++nb) {
          for (size_t iht = 0; iht < 8; ++iht) {
            const size_t i = 8 * ((inj == 0 ? 0 : 3) + nb) + iht;
            _h_srcounters[i] = bookCounter("srcount_j" + toString(njmax) + "_b" + toString(nb) + "_ht" + toString(iht+1));
          }
        }
      }
      // Special nj >= 4, nb >= 4 bins
      for (size_t iht = 0; iht < 3; ++iht) {
        _h_srcounters[8*7 + iht] = bookCounter("srcount_j4_b4_ht" + toString(iht+1));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get baseline photons, electrons & muons
      Particles photons = apply<ParticleFinder>(event, "Photons").particles(Cuts::pT > 25*GeV);
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particles(Cuts::pT > 10*GeV);
      Particles muons = apply<ParticleFinder>(event, "Muons").particles(Cuts::pT > 10*GeV);

      // Electron/muon isolation (guesswork/copied from other CMS analysis -- paper is unspecific)
      const Particles calofs = apply<ParticleFinder>(event, "Clusters").particles();
      ifilter_discard(photons, [&](const Particle& y) {
          double ptsum = -y.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,y) < 0.3) ptsum += p.pT();
          return ptsum / y.pT() > 0.1;
        });
      ifilter_discard(elecs, [&](const Particle& e) {
          double ptsum = -e.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,e) < 0.3) ptsum += p.pT();
          return ptsum / e.pT() > 0.1;
        });
      ifilter_discard(muons, [&](const Particle& m) {
          double ptsum = -m.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,m) < 0.3) ptsum += p.pT();
          return ptsum / m.pT() > 0.2;
        });

      // Veto the event if there are any remaining baseline photons or leptons
      if (!photons.empty()) vetoEvent;
      if (!elecs.empty()) vetoEvent;
      if (!muons.empty()) vetoEvent;


      // Get jets and apply jet-based event-selection cuts
      const JetAlg& jetproj = apply<JetAlg>(event, "Jets");
      const Jets alljets = jetproj.jetsByPt(Cuts::abseta < 3.0 && Cuts::Et > 37*GeV); //< most inclusive jets requirement
      if (filter_select(alljets, Cuts::Et > 73*GeV).size() < 2) vetoEvent; //< most inclusive lead jets requirement

      // Filter jets into different Et requirements & compute corresponding HTs
      /// @note It's not clear if different HTs are used to choose the HT bins
      const Jets jets37 = filter_select(alljets, Cuts::Et > 37*GeV);
      const Jets jets43 = filter_select(jets37, Cuts::Et > 43*GeV);
      const Jets jets50 = filter_select(jets43, Cuts::Et > 50*GeV);
      const double ht37 = sum(jets37, Et, 0.0);
      const double ht43 = sum(jets43, Et, 0.0);
      const double ht50 = sum(jets50, Et, 0.0);

      // Find the relevant HT bin and apply leading jet event-selection cuts
      static const vector<double> htcuts = { /* 275., 325., */ 375., 475., 575., 675., 775., 875.}; //< comment to avoid jets50 "fall-down"
      const int iht = inRange(ht37, 275*GeV, 325*GeV) ? 0 : inRange(ht43, 325*GeV, 375*GeV) ? 1 : (2+binIndex(ht50, htcuts, true));
      MSG_TRACE("HT = {" << ht37 << ", " << ht43 << ", " << ht50 << "} => IHT = " << iht);
      if (iht < 0) vetoEvent;
      if (iht == 1 && filter_select(jets43, Cuts::Et > 78*GeV).size() < 2) vetoEvent;
      if (iht >= 2 && filter_select(jets50, Cuts::Et > 100*GeV).size() < 2) vetoEvent;

      // Create references for uniform access to relevant set of jets & HT
      const double etcut = iht == 0 ? 37. : iht == 1 ? 43. : 50.;
      const double& ht = iht == 0 ? ht37 : iht == 1 ? ht43 : ht50;
      const Jets& jets = iht == 0 ? jets37 : iht == 1 ? jets43 : jets50;
      if (!jetproj.jets(Cuts::abseta > 3 && Cuts::Et > etcut*GeV).empty()) vetoEvent;
      const size_t nj = jets.size();
      const size_t nb = count_if(jets.begin(), jets.end(), [](const Jet& j) { return j.bTagged(Cuts::pT > 5*GeV); });

      // Compute HTmiss = pT of 4-vector sum of jet momenta
      const FourMomentum jsum = sum(jets, mom, FourMomentum());
      const double htmiss = jsum.pT();

      // Require HTmiss / ETmiss < 1.25
      const double etmiss = apply<SmearedMET>(event, "MET").met();
      if (htmiss/etmiss > 1.25) vetoEvent;

      // Compute DeltaHT = minimum difference of "dijet" ETs, i.e. max(|1+2-3|, |1+3-2|, |2+3-1|)
      double deltaht = -1;
      vector<double> jetets; transform(jets, jetets, Et);
      for (int i = 1; i < (1 << (jetets.size()-1)); ++i) { // count from 1 to 2**N-1, i.e. through all heterogeneous bitmasks with MSB(2**N)==0
        const bitset<10> bits(i); /// @warning There'd better not be more than 10 jets...
        const double htdiff = partition_diff(bits, jetets);
        // MSG_INFO(bits.to_string() << " => " << htdiff);
        if (deltaht < 0 || htdiff < deltaht) deltaht = htdiff;
      }
      MSG_DEBUG("dHT_bitmask = " << deltaht);

      // Cross-check calculation in 2- and 3-jet cases
      // if (jets.size() == 2) {
      //   MSG_INFO("dHT2 = " << fabs(jets[0].Et() - jets[1].Et()));
      // } else if (jets.size() == 3) {
      //   double deltaht_01_2 = fabs(jets[0].Et()+jets[1].Et()-jets[2].Et());
      //   double deltaht_02_1 = fabs(jets[0].Et()+jets[2].Et()-jets[1].Et());
      //   double deltaht_12_0 = fabs(jets[1].Et()+jets[2].Et()-jets[0].Et());
      //   MSG_INFO("dHT3 = " << min({deltaht_01_2, deltaht_02_1, deltaht_12_0}));
      // }

      // Compute alphaT from the above
      double alphaT = fabs(0.5*((ht-deltaht)/(sqrt((ht*ht)-(htmiss*htmiss)))));
      if (alphaT < 0.55) vetoEvent;

      /// @todo Need to include trigger efficiency sampling or weighting?

      // Fill histograms
      const double weight = event.weight();
      const size_t inj = nj < 4 ? 0 : 1;
      const size_t inb = nb < 4 ? nb : 4;
      if (iht >= 2)
        (inj == 0 ? _h_alphaT23 : _h_alphaT4)->fill(alphaT, weight);

      // Fill the appropriate counter -- after working out the irregular SR bin index! *sigh*
      size_t i = 8 * ((inj == 0 ? 0 : 3) + inb) + iht;
      if (inj == 1 && inb == 4) i = 8*7 + (iht < 3 ? iht : 2);
      MSG_INFO("inj = " << inj << ", inb = " << inb << ", i = " << i);
      _h_srcounters[i]->fill(weight);

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection()/femtobarn*11.7/sumOfWeights();
      scale({_h_alphaT23,_h_alphaT4}, sf);
      for (size_t i = 0; i < 8*7+3; ++i)
        scale(_h_srcounters[i], sf);

    }

    //@}


    /// @name Utility functions for partitioning jet pTs into two groups and summing/diffing them
    //@{

    /// Sum the given values into two subsets according to the provided bitmask
    template <size_t N>
    pair<double, double> partition_sum(const bitset<N>& mask, const vector<double>& vals) const {
      pair<double, double> rtn(0., 0.);
      for (size_t i = 0; i < vals.size(); ++i) {
        (!mask[vals.size()-1-i] ? rtn.first : rtn.second) += vals[i];
      }
      return rtn;
    }

    /// Return the difference between summed subsets according to the provided bitmask
    template <size_t N>
    double partition_diff(const bitset<N>& mask, const vector<double>& vals) const {
      const pair<double, double> sums = partition_sum(mask, vals);
      const double diff = fabs(sums.first - sums.second);
      MSG_TRACE(mask.to_string() << ": " << sums.first << "/" << sums.second << " => " << diff);
      return diff;
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_alphaT23, _h_alphaT4;
    vector<CounterPtr> _h_srcounters;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1223519);


}
