// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Jet fragmentation at 7 TeV
  class ATLAS_2011_I929691 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2011_I929691);


    /// Initialisation
    void init() {
      const FinalState fs(Cuts::abseta < 2.0);

      FastJets antikt_06_jets(fs, FastJets::ANTIKT, 0.6, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      declare(antikt_06_jets, "jets");

      ChargedFinalState tracks(Cuts::pT > 0.5*GeV && Cuts::abseta < 2.0);
      declare(tracks, "tracks");

      // Set up the histograms (each element is a binning in jet pT)
      for (size_t i = 0; i < 10; i++) {
        _p_F_z[i]     = bookProfile1D(i+1, 1, 1);
        _p_rho_r[i]   = bookProfile1D(i+11, 1, 1);
        _p_f_pTrel[i] = bookProfile1D(i+21, 1, 1);
      }

    }


    // Per-event analysis
    void analyze(const Event& event) {

      const Jets alljets = apply<FastJets>(event, "jets").jetsByPt(Cuts::absrap < 1.2);
      const Particles& tracks = apply<ChargedFinalState>(event, "tracks").particlesByPt();

      for (size_t i = 0; i < 10; ++i) {

        const Jets jets = filter_select(alljets, Cuts::pT > bedges[i] && Cuts::pT < bedges[i+1]);
        const int n_jets = jets.size();
        if (n_jets == 0) continue;

        // First... count the tracks
        Histo1D h_ntracks_z(*_p_F_z[i]), h_ntracks_r(*_p_rho_r[i]), h_ntracks_pTrel(*_p_f_pTrel[i]);

        for (const Jet& j : jets) {
          for (const Particle& p : tracks) {
            const double dr = deltaR(j, p, RAPIDITY);
            if (dr > 0.6) continue; // The paper uses pseudorapidity, but this is a requirement for filling the histogram
            h_ntracks_z.fill(z(j, p), 1.0/n_jets);
            h_ntracks_r.fill(dr, 1.0/n_jets);
            h_ntracks_pTrel.fill(pTrel(j, p), 1.0/n_jets);
          }
        }

        // Then... calculate the observable and fill the profiles
        const double weight = event.weight();
        for (const HistoBin1D& b : h_ntracks_z.bins())
          _p_F_z[i]->fill(b.xMid(), b.height(), weight);
        for (const HistoBin1D& b : h_ntracks_r.bins())
          _p_rho_r[i]->fill(b.xMid(), b.area()/annulus_area(b.xMin(), b.xMax()), weight);
        for (const HistoBin1D& b : h_ntracks_pTrel.bins())
          _p_f_pTrel[i]->fill(b.xMid(), b.height(), weight);

      }

    }


    double z (const Jet& jet, const Particle& ch) {
      return dot(jet.p3(), ch.p3()) / jet.p3().mod2();
    }

    double pTrel (const Jet& jet, const Particle& ch) {
      return (ch.p3().cross(jet.p3())).mod()/(jet.p3().mod());
    }

    // To calculate the area of the annulus in an r bin
    double annulus_area(double r1, double r2) {
      return M_PI*(sqr(r2) - sqr(r1));
    }


  private:

    Profile1DPtr _p_F_z[10], _p_rho_r[10], _p_f_pTrel[10];
    const vector<double> bedges = { 25., 40., 60., 80., 110., 160., 210., 260., 310., 400., 500. };

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I929691);


}
