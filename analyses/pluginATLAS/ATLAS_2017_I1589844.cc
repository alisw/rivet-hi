// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// kT splittings in Z events at 8 TeV
  class ATLAS_2017_I1589844 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructors
    ATLAS_2017_I1589844(string name="ATLAS_2017_I1589844") : Analysis(name) {
      _mode = 1; // pick electron channel by default
      setNeedsCrossSection(true);
    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;

      const Cut cuts = (_mode) ?
        Cuts::pT > 25*GeV && (Cuts::abseta <= 1.37 || (Cuts::abseta >= 1.52 && Cuts::abseta < 2.47)) : //< electron channel version
        (Cuts::pT > 25*GeV) && (Cuts::abseta < 2.4); //< muon channel version

      IdentifiedFinalState bareleptons(fs);
      bareleptons.acceptIdPair(_mode? PID::ELECTRON : PID::MUON);
      const DressedLeptons leptons(fs, bareleptons, 0.1, cuts, true);
      declare(leptons, "leptons");

      const ChargedFinalState cfs(Cuts::abseta < 2.5 && Cuts::pT > 0.4*GeV);
      VetoedFinalState jet_fs(cfs);
      jet_fs.addVetoOnThisFinalState(leptons);
      declare(FastJets(jet_fs, FastJets::KT, 0.4), "Kt04Jets");
      declare(FastJets(jet_fs, FastJets::KT, 1.0), "Kt10Jets");

      VetoedFinalState jet_fs_all(Cuts::abseta < 2.5 && Cuts::pT > 0.4*GeV);
      jet_fs_all.addVetoOnThisFinalState(leptons);
      FastJets jetpro04_all(jet_fs_all, FastJets::KT, 0.4);
      jetpro04_all.useInvisibles();
      declare(jetpro04_all, "Kt04Jets_all");
      FastJets jetpro10_all(jet_fs_all, FastJets::KT, 1.0);
      jetpro10_all.useInvisibles();
      declare(jetpro10_all, "Kt10Jets_all");

      // Histograms with data binning
      _ndij = 8;
      for (size_t i = 0; i < _ndij; ++i) {
        string label = "d" + to_str(i) + "_kT4";
        _h[label] = bookHisto1D(i + 1, 1, _mode + 1);
        _h[label + "_all"] = bookHisto1D(i + 1, 1, _mode + 5);
        label = "d" + to_str(i) + "_kT10";
        _h[label] = bookHisto1D(i + 1, 1, _mode + 3);
        _h[label + "_all"] = bookHisto1D(i + 1, 1, _mode + 7);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {

      // Check we have a Z candidate:
      const vector<DressedLepton>& leptons = apply<DressedLeptons>(e, "leptons").dressedLeptons();
      if (leptons.size() != 2)  vetoEvent;
      if (leptons[0].charge3()*leptons[1].charge3() > 0) vetoEvent;
      const double dilepton_mass = (leptons[0].momentum() + leptons[1].momentum()).mass();
      if (!inRange(dilepton_mass, 71*GeV, 111*GeV)) vetoEvent;

      const double weight = e.weight();

      // Get kT splitting scales (charged particles only)
      const FastJets& jetpro04 = applyProjection<FastJets>(e, "Kt04Jets");
      const shared_ptr<fastjet::ClusterSequence> seq04 = jetpro04.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq04->n_particles()); ++i) {
        const double dij = sqrt(seq04->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = "d" + to_str(i) + "_kT4";
        _h[label]->fill(dij, weight);
      }
      const FastJets& jetpro10 = applyProjection<FastJets>(e, "Kt10Jets");
      const shared_ptr<fastjet::ClusterSequence> seq10 = jetpro10.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq10->n_particles()); ++i) {
        const double dij = sqrt(seq10->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = "d" + to_str(i) + "_kT10";
        _h[label]->fill(dij, weight);
      }

      // Get kT splitting scales (all particles)
      const FastJets& jetpro04_all = applyProjection<FastJets>(e, "Kt04Jets_all");
      const shared_ptr<fastjet::ClusterSequence> seq04_all = jetpro04_all.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq04_all->n_particles()); ++i) {
        const double dij = sqrt(seq04_all->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = "d" + to_str(i) + "_kT4_all";
        _h[label]->fill(dij, weight);
      }
      const FastJets& jetpro10_all = applyProjection<FastJets>(e, "Kt10Jets_all");
      const shared_ptr<fastjet::ClusterSequence> seq10_all = jetpro10_all.clusterSeq();
      for (size_t i = 0; i < min(_ndij, (size_t)seq10_all->n_particles()); ++i) {
        const double dij = sqrt(seq10_all->exclusive_dmerge_max(i))/GeV;
        if (dij <= 0.0) continue;
        const string label = "d" + to_str(i) + "_kT10_all";
        _h[label]->fill(dij, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSectionPerEvent();
      for (auto& kv : _h) scale(kv.second, sf);
    }

    //@}


  protected:

    // Data members like post-cuts event weight counters go here
    size_t _mode, _ndij;


  private:

    // Histograms
    map<string, Histo1DPtr> _h;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1589844);



  /// kT splittings in Z events at 8 TeV (electron channel)
  struct ATLAS_2017_I1589844_EL : public ATLAS_2017_I1589844 {
    ATLAS_2017_I1589844_EL() : ATLAS_2017_I1589844("ATLAS_2017_I1589844_EL") { _mode = 1; }
  };
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1589844_EL);


  /// kT splittings in Z events at 8 TeV (muon channel)
  struct ATLAS_2017_I1589844_MU : public ATLAS_2017_I1589844 {
    ATLAS_2017_I1589844_MU() : ATLAS_2017_I1589844("ATLAS_2017_I1589844_MU") { _mode = 0; }
  };
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1589844_MU);


}
