//-*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  class ALICE_2017_I1512110 : public Analysis {
  public:

    /// Constructor
    ALICE_2017_I1512110()
      : Analysis("ALICE_2017_I1512110"),
        _rapmax(0.8)
    {    }


    void init() {

      const UnstableFinalState ufs(Cuts::absrap < _rapmax);
      addProjection(ufs, "UFS");

      _h_pi0 = bookHisto1D(3,1,1);
      _h_eta = bookHisto1D(4,1,1);
      _h_etaToPion = bookScatter2D(5,1,1);

      // temporary plots with the binning of _h_etaToPion
      // to construct the eta/pi0 ratio in the end
      _temp_h_pion = bookHisto1D("TMP/h_pion",refData(5,1,1));
      _temp_h_eta = bookHisto1D("TMP/h_eta",refData(5,1,1));
    }


    void analyze(const Event& event) {

      const double weight = event.weight();
      const UnstableFinalState& ufs = applyProjection<UnstableFinalState>(event, "UFS");

      for (const Particle& p : ufs.particles()) {

        if (p.pid() == 111) {
          // neutral pion; ALICE corrects for pi0 feed-down
          if ( !(p.hasAncestor(310)  || p.hasAncestor(130)   || // K0_s, K0_l
                 p.hasAncestor(321)  || p.hasAncestor(-321)  || // K+,K-
                 p.hasAncestor(3122) || p.hasAncestor(-3122) || // Lambda, Anti-Lambda
                 p.hasAncestor(3212) || p.hasAncestor(-3212) || // Sigma0
                 p.hasAncestor(3222) || p.hasAncestor(-3222) || // Sigmas
                 p.hasAncestor(3112) || p.hasAncestor(-3112) || // Sigmas
                 p.hasAncestor(3322) || p.hasAncestor(-3322) || // Cascades
                 p.hasAncestor(3312) || p.hasAncestor(-3312) )) // Cascades
            {
              _h_pi0->fill(p.pT()/GeV, weight /(TWOPI*p.pT()/GeV*2*_rapmax));
              _temp_h_pion->fill(p.pT()/GeV, weight);
            }
        }
        else if (p.pid() == 221){
          // eta meson
          _h_eta->fill(p.pT()/GeV, weight /(TWOPI*p.pT()/GeV*2*_rapmax));
          _temp_h_eta->fill(p.pT()/GeV, weight);
        }
      }
    }


    void finalize() {

      scale(_h_pi0, crossSection()/picobarn/sumOfWeights());
      scale(_h_eta, crossSection()/picobarn/sumOfWeights());
      divide(_temp_h_eta, _temp_h_pion, _h_etaToPion);

    }


  private:

    double _rapmax;
    Histo1DPtr _h_pi0;
    Histo1DPtr _h_eta;
    Histo1DPtr _temp_h_pion;
    Histo1DPtr _temp_h_eta;
    Scatter2DPtr _h_etaToPion;

  };


  DECLARE_RIVET_PLUGIN(ALICE_2017_I1512110);

}
