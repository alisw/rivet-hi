// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Higgs differential cross-section combination between the ATLAS measurements in the yy and 4l channels
  ///
  /// Computes Higgs transverse momentum, rapidity, jet multiplicity and leading jet pT.
  ///
  /// @author Michaela Queitsch-Maitland <michaela.queitsch-maitland@cern.ch>
  /// @author Dag Gillberg <dag.gillberg@cern.ch>
  /// @author Florian Bernlochner <florian.bernlochner@cern.ch>
  /// @author Sarah Heim <sarah.heim@cern.ch>
  class ATLAS_2015_I1364361 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1364361);


    /// Book histograms and initialise projections before the run
    void init() {

      // All final state particles
      const FinalState fs;
      declare(fs, "FS");

      // Histograms with data bins
      _h_pTH_incl   = bookHisto1D(1,1,1);
      _h_yH_incl    = bookHisto1D(2,1,1);
      _h_Njets_incl = bookHisto1D(3,1,1);
      _h_pTj1_incl  = bookHisto1D(4,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the final state particles ordered by pT
      const Particles& fs = apply<FinalState>(event, "FS").particlesByPt();

      // Find a stable Higgs (mandatory)
      const auto higgsiter = std::find_if(fs.begin(), fs.end(), [](const Particle& p){ return p.pid() == PID::HIGGSBOSON; });
      if (higgsiter == fs.end()) vetoEvent;
      const Particle& higgs = *higgsiter;

      // bool stable_higgs = false;
      // const Particle* higgs = 0;
      // for (const Particle& p : fs) {
      //   if (p.pid() == PID::HIGGSBOSON) {
      //     stable_higgs = true;
      //     higgs = &p;
      //     break;
      //   }
      // }
      // if (!stable_higgs) {
      //   MSG_WARNING("FATAL: No stable Higgs found in event record.\n");
      //   vetoEvent;
      // }


      // Loop over final state particles and fill various particle vectors
      Particles leptons, photons, jet_ptcls;
      foreach ( const Particle& ptcl, fs ) {
        // Do not include the Higgs in jet finding!
        if ( ptcl.pid() == PID::HIGGSBOSON ) continue;
        // Neutrinos not from hadronisation
        if ( ptcl.isNeutrino() && !ptcl.fromHadron() ) continue;
        // Electrons and muons not from hadronisation
        if ( ( ptcl.abspid() == PID::ELECTRON || ptcl.abspid() == PID::MUON ) && !ptcl.fromHadron() ) {
          leptons.push_back(ptcl);
          continue;
        }
        // Photons not from hadronisation
        if ( ptcl.abspid() == PID::PHOTON && !ptcl.fromHadron() ) {
          photons.push_back(ptcl);
          continue;
        }
        // Add particle to jet inputs
        jet_ptcls.push_back(ptcl);
      }

      // Match FS photons to leptons within cone R=0.1
      // If they are not 'dressing' photons, add to jet particle vector
      foreach ( const Particle& ph, photons ) {
        bool fsr_photon = false;
        foreach ( const Particle& lep, leptons ) {
          if ( deltaR(ph.momentum(),lep.momentum()) < 0.1 ){
            fsr_photon=true;
            continue;
          }
        }
        if ( !fsr_photon ) jet_ptcls.push_back(ph);
      }

      // Let's build the jets! By hand...
      const PseudoJets pjs_in = mkPseudoJets(jet_ptcls);
      const fastjet::JetDefinition jdef(fastjet::antikt_algorithm, 0.4);
      const Jets alljets = mkJets(fastjet::ClusterSequence(pjs_in, jdef).inclusive_jets());
      const Jets jets = sortByPt(filterBy(alljets, Cuts::pT > 30*GeV && Cuts::absrap < 4.4));
      // FastJets jet_pro(FastJets::ANTIKT, 0.4);
      // jet_pro.calc(jet_ptcls);
      // Jets jets = jet_pro.jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);

      const double weight = event.weight();
      _h_pTH_incl->fill(higgs.pT(), weight);
      _h_yH_incl->fill(higgs.absrap(), weight);
      _h_Njets_incl->fill(jets.size() > 3 ? 3 : jets.size(), weight);
      _h_pTj1_incl->fill(jets.empty() ? 0 : jets[0].pT(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double xs = crossSectionPerEvent();
      scale(_h_pTH_incl, xs);
      scale(_h_yH_incl, xs);
      scale(_h_Njets_incl, xs);
      scale(_h_pTj1_incl, xs);
    }


  private:

    Histo1DPtr _h_pTH_incl, _h_yH_incl, _h_Njets_incl, _h_pTj1_incl;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1364361);


}
