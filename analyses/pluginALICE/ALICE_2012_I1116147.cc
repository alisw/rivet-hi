//-*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  class ALICE_2012_I1116147 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALICE_2012_I1116147);


    /// Initialise projections and histograms
    void init() {

      const UnstableFinalState ufs(Cuts::absrap < RAPMAX);
      addProjection(ufs, "UFS");

      // Check if cm energy is 7 TeV or 0.9 TeV
      if (fuzzyEquals(sqrtS()/GeV, 900, 1E-3))       _cm_energy_case = 1;
      else if (fuzzyEquals(sqrtS()/GeV, 7000, 1E-3)) _cm_energy_case = 2;
      if (_cm_energy_case == 0)
        throw UserError("Center of mass energy of the given input is neither 900 nor 7000 GeV.");

      // Book histos
      if (_cm_energy_case == 1) {
        _h_pi0 = bookHisto1D(2,1,1);
      } else {
        _h_pi0 = bookHisto1D(1,1,1);
        _h_eta = bookHisto1D(3,1,1);
        _h_etaToPion = bookScatter2D(4,1,1);
      }

      // Temporary plots with the binning of _h_etaToPion to construct the eta/pi0 ratio
      _temp_h_pion = bookHisto1D("TMP/h_pion", refData(4,1,1));
      _temp_h_eta = bookHisto1D("TMP/h_eta", refData(4,1,1));
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const FinalState& ufs = apply<UnstableFinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const double normfactor = TWOPI*p.pT()/GeV*2*RAPMAX;
        if (p.pid() == 111) {
          // Neutral pion; ALICE corrects for pi0 feed-down from K_0_s and Lambda
          if (p.hasAncestor(310) || p.hasAncestor(3122) || p.hasAncestor(-3122)) continue; //< K_0_s, Lambda, Anti-Lambda
          _h_pi0->fill(p.pT()/GeV, weight/normfactor);
          _temp_h_pion->fill(p.pT()/GeV, weight);
        } else if (p.pid() == 221 && _cm_energy_case == 2) {
          // eta meson (only for 7 TeV)
          _h_eta->fill(p.pT()/GeV, weight/normfactor);
          _temp_h_eta->fill(p.pT()/GeV, weight);
        }
      }
    }


    /// Normalize histos and construct ratio
    void finalize() {
      scale(_h_pi0, crossSection()/microbarn/sumOfWeights());
      if (_cm_energy_case == 2) {
        divide(_temp_h_eta, _temp_h_pion, _h_etaToPion);
        scale(_h_eta, crossSection()/microbarn/sumOfWeights());
      }
    }


  private:

    const double RAPMAX = 0.8;
    int _cm_energy_case = 0;

    Histo1DPtr _h_pi0, _h_eta;
    Histo1DPtr _temp_h_pion, _temp_h_eta;
    Scatter2DPtr _h_etaToPion;

  };


  DECLARE_RIVET_PLUGIN(ALICE_2012_I1116147);

}
