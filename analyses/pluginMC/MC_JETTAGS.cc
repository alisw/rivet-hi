// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {



  /// @brief MC validation analysis for jet events
  class MC_JETTAGS : public Analysis {
  public:

    MC_JETTAGS()
      : Analysis("MC_JETTAGS")
    {    }


    void init() {
      FinalState fs;
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets04");
      declare(FastJets(fs, FastJets::ANTIKT, 0.6), "Jets06");

      _h_numBTagsPerJet[0] = bookHisto1D("numBTagsPer04Jet", 5, -0.5, 4.5);
      _h_numBTagsPerJet[1] = bookHisto1D("numBTagsPer06Jet", 5, -0.5, 4.5);
      _h_numCTagsPerJet[0] = bookHisto1D("numCTagsPer04Jet", 5, -0.5, 4.5);
      _h_numCTagsPerJet[1] = bookHisto1D("numCTagsPer06Jet", 5, -0.5, 4.5);
      _h_numTauTagsPerJet[0] = bookHisto1D("numTauTagsPer04Jet", 5, -0.5, 4.5);
      _h_numTauTagsPerJet[1] = bookHisto1D("numTauTagsPer06Jet", 5, -0.5, 4.5);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      const Jets jets04 = apply<FastJets>(event, "Jets04").jetsByPt(20*GeV);
      const Jets jets06 = apply<FastJets>(event, "Jets06").jetsByPt(20*GeV);

      for (const Jet& j : jets04) {
        _h_numBTagsPerJet[0]->fill(j.bTags().size(), weight);
        _h_numCTagsPerJet[0]->fill(j.cTags().size(), weight);
        _h_numTauTagsPerJet[0]->fill(j.tauTags().size(), weight);
      }
      for (const Jet& j : jets06) {
        _h_numBTagsPerJet[1]->fill(j.bTags().size(), weight);
        _h_numCTagsPerJet[1]->fill(j.cTags().size(), weight);
        _h_numTauTagsPerJet[1]->fill(j.tauTags().size(), weight);
      }
    }


    void finalize() {
      for (size_t i = 0; i < 2; ++i) {
        normalize(_h_numBTagsPerJet[i]);
        normalize(_h_numCTagsPerJet[i]);
        normalize(_h_numTauTagsPerJet[i]);
      }
    }


  private:

    Histo1DPtr _h_numBTagsPerJet[2], _h_numCTagsPerJet[2], _h_numTauTagsPerJet[2];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_JETTAGS);

}
