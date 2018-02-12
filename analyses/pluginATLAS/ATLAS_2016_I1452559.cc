#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedMET.hh"

namespace Rivet {


  /// ATLAS 13 TeV monojet search with 3.2/fb of pp data
  class ATLAS_2016_I1452559 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1452559);

    void init() {

      FastJets jets(FinalState(Cuts::abseta < 4.9), FastJets::ANTIKT, 0.4);
      SmearedJets recojets(jets, JET_SMEAR_ATLAS_RUN1);
      declare(recojets, "Jets");

      FinalState electrons(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.47 && Cuts::pT > 20*GeV);
      SmearedParticles recoelectrons(electrons, ELECTRON_EFF_ATLAS_RUN1);
      declare(recoelectrons, "Electrons");

      FinalState muons(Cuts::abspid == PID::MUON && Cuts::abseta < 2.50 && Cuts::pT > 10*GeV);
      SmearedParticles recomuons(muons, MUON_EFF_ATLAS_RUN1);
      declare(recomuons, "Muons");

      VisibleFinalState calofs(Cuts::abseta < 4.9 && Cuts::abspid != PID::MUON);
      MissingMomentum met(calofs);
      SmearedMET recomet(met, MET_SMEAR_ATLAS_RUN1);
      declare(recomet, "MET");


      /// Book histograms
      for (size_t i = 0; i < 7; ++i)
        _count_IM[i] = bookCounter("count_IM" + toString(i+1));
      for (size_t i = 0; i < 6; ++i)
        _count_EM[i] = bookCounter("count_EM" + toString(i+1));

    }


    void analyze(const Event& event) {

      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.8);
      const Particles elecs = apply<ParticleFinder>(event, "Electrons").particlesByPt();
      const Particles mus = apply<ParticleFinder>(event, "Muons").particlesByPt();
      MSG_DEBUG("Number of raw jets, electrons, muons = "
                << jets.size() << ", " << elecs.size() << ", " << mus.size());

      // Discard jets very close to electrons, or with low track multiplicity and close to muons
      const Jets isojets = filter_discard(jets, [&](const Jet& j) {
          /// @todo Add track efficiency random filtering
          if (any(elecs, deltaRLess(j, 0.2))) return true;
          if (j.particles(Cuts::abscharge > 0 && Cuts::pT > 0.4*GeV).size() < 3 &&
              any(mus, deltaRLess(j, 0.4))) return true;
          return false;
        });

      // Discard electrons close to remaining jets
      const Particles isoelecs = filter_discard(elecs, [&](const Particle& e) {
          return any(isojets, deltaRLess(e, 0.4));
        });

      // Discard muons close to remaining jets
      const Particles isomus = filter_discard(mus, [&](const Particle& m) {
          for (const Jet& j : isojets) {
            if (deltaR(j,m) > 0.4) continue;
            if (j.particles(Cuts::abscharge > 0 && Cuts::pT > 0.4*GeV).size() > 3) return true;
          }
          return false;
        });

      // Calculate ETmiss
      //const Vector3& vet = apply<MissingMomentum>(event, "MET").vectorEt();
      const Vector3& vet = apply<SmearedMET>(event, "MET").vectorEt();
      const double etmiss = vet.perp();


      // Event selection cuts
      if (etmiss < 250*GeV) vetoEvent;
      // Require at least one jet with pT > 250 GeV and |eta| < 2.4
      if (filter_select(isojets, Cuts::pT > 250*GeV && Cuts::abseta < 2.4).empty()) vetoEvent;
      // Require at most 4 jets with pT > 30 GeV and |eta| < 2.8
      if (filter_select(isojets, Cuts::pT > 30*GeV).size() > 4) vetoEvent;
      // Require no isolated jets within |dphi| < 0.4 of the MET vector
      if (any(isojets, deltaPhiLess(-vet, 0.4))) vetoEvent;
      // Require no isolated electrons or muons
      if (!isoelecs.empty() || !isomus.empty()) vetoEvent;


      ////////////////////


      const double weight = event.weight();

      // Get ETmiss bin number and fill counters
      const int i_etmiss = binIndex(etmiss/GeV, ETMISS_CUTS);
      // Inclusive ETmiss bins
      for (int ibin = 0; ibin < 7; ++ibin)
        if (i_etmiss >= ibin) _count_IM[ibin]->fill(weight);
      // Exclusive ETmiss bins
      if (inRange(i_etmiss, 0, 6)) _count_EM[i_etmiss]->fill(weight);

    }


    void finalize() {
      const double norm = 3.2*crossSection()/femtobarn;
      scale(_count_IM, norm/sumOfWeights());
      scale(_count_EM, norm/sumOfWeights());
    }


  private:

    const vector<double> ETMISS_CUTS = { 250, 300, 350, 400, 500, 600, 700, 13000 };
    CounterPtr _count_IM[7], _count_EM[6];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1452559);

}
