// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// CMS cross-section and angular correlations in Z boson + b-hadrons events at 7 TeV
  class CMS_2013_I1256943 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2013_I1256943);


    /// Add projections and book histograms
    void init() {
      _sumW = 0;
      _sumW50 = 0;
      _sumWpT = 0;

      FinalState fs(Cuts::abseta < 2.4 && Cuts::pT > 20*GeV);
      declare(fs, "FS");

      UnstableFinalState ufs(Cuts::abseta < 2 && Cuts::pT > 15*GeV);
      declare(ufs, "UFS");

      Cut zetacut = Cuts::abseta < 2.4;

      ZFinder zfindermu(fs, zetacut, PID::MUON, 81.0*GeV, 101.0*GeV, 0.1, ZFinder::NOCLUSTER, ZFinder::TRACK, 91.2*GeV);
      declare(zfindermu, "ZFinderMu");

      ZFinder zfinderel(fs, zetacut, PID::ELECTRON, 81.0*GeV, 101.0*GeV, 0.1, ZFinder::NOCLUSTER, ZFinder::TRACK, 91.2*GeV);
      declare(zfinderel, "ZFinderEl");


      // Histograms in non-boosted region of Z pT
      _h_dR_BB = bookHisto1D(1, 1, 1);
      _h_dphi_BB = bookHisto1D(2, 1, 1);
      _h_min_dR_ZB = bookHisto1D(3, 1, 1);
      _h_A_ZBB = bookHisto1D(4, 1, 1);

      // Histograms in boosted region of Z pT (pT > 50 GeV)
      _h_dR_BB_boost = bookHisto1D(5, 1, 1);
      _h_dphi_BB_boost = bookHisto1D(6, 1, 1);
      _h_min_dR_ZB_boost = bookHisto1D(7, 1, 1);
      _h_A_ZBB_boost = bookHisto1D(8, 1, 1);

      _h_min_ZpT = bookHisto1D(9,1,1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      vector<FourMomentum> Bmom;

      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      const ZFinder& zfindermu = apply<ZFinder>(e, "ZFinderMu");
      const ZFinder& zfinderel = apply<ZFinder>(e, "ZFinderEl");

      // Look for a Z --> mu+ mu- event in the final state
      if (zfindermu.empty() && zfinderel.empty()) vetoEvent;

      const Particles& z = !zfindermu.empty() ? zfindermu.bosons() : zfinderel.bosons();
      const bool is_boosted = ( z[0].pT() > 50*GeV );

      // Loop over the unstable particles
      for (const Particle& p : ufs.particles()) {
        const PdgId pid = p.pid();

        // Look for particles with a bottom quark
        if (PID::hasBottom(pid)) {

          bool good_B = false;
          const GenParticle* pgen = p.genParticle();
          const GenVertex* vgen = pgen -> end_vertex();

          // Loop over the decay products of each unstable particle, looking for a b-hadron pair
          /// @todo Avoid HepMC API
          for (GenVertex::particles_out_const_iterator it = vgen->particles_out_const_begin(); it !=  vgen->particles_out_const_end(); ++it) {
            // If the particle produced has a bottom quark do not count it and go to the next loop cycle.
            if (!( PID::hasBottom( (*it)->pdg_id() ) ) ) {
              good_B = true;
              continue;
            } else {
              good_B = false;
              break;
            }
          }
          if (good_B ) Bmom.push_back( p.momentum() );
        }
        else continue;
      }

      // If there are more than two B's in the final state veto the event
      if (Bmom.size() != 2 ) vetoEvent;

      // Calculate the observables
      double dphiBB = deltaPhi(Bmom[0], Bmom[1]);
      double dRBB = deltaR(Bmom[0], Bmom[1]);

      const FourMomentum& pZ = z[0].momentum();
      const bool closest_B = ( deltaR(pZ, Bmom[0]) < deltaR(pZ, Bmom[1]) );
      const double mindR_ZB = closest_B ? deltaR(pZ, Bmom[0]) : deltaR(pZ, Bmom[1]);
      const double maxdR_ZB = closest_B ? deltaR(pZ, Bmom[1]) : deltaR(pZ, Bmom[0]);
      const double AZBB = ( maxdR_ZB - mindR_ZB ) / ( maxdR_ZB + mindR_ZB );

      // Get event weight for histogramming
      const double weight = e.weight();

      // Fill the histograms in the non-boosted region
      _h_dphi_BB->fill(dphiBB, weight);
      _h_dR_BB->fill(dRBB, weight);
      _h_min_dR_ZB->fill(mindR_ZB, weight);
      _h_A_ZBB->fill(AZBB, weight);
      _sumW += weight;
      _sumWpT += weight;

      // Fill the histograms in the boosted region
      if (is_boosted) {
        _sumW50 += weight;
        _h_dphi_BB_boost->fill(dphiBB, weight);
        _h_dR_BB_boost->fill(dRBB, weight);
        _h_min_dR_ZB_boost->fill(mindR_ZB, weight);
        _h_A_ZBB_boost->fill(AZBB, weight);
      }

      // Fill Z pT (cumulative) histogram
      _h_min_ZpT->fill(0, weight);
      if (pZ.pT() > 40*GeV ) {
        _sumWpT += weight;
        _h_min_ZpT->fill(40, weight);
      }
      if (pZ.pT() > 80*GeV ) {
        _sumWpT += weight;
        _h_min_ZpT->fill(80, weight);
      }
      if (pZ.pT() > 120*GeV ) {
        _sumWpT += weight;
        _h_min_ZpT->fill(120, weight);
      }

      Bmom.clear();
    }


    /// Finalize
    void finalize() {

      // Normalize excluding overflow bins (d'oh)
      normalize(_h_dR_BB, 0.7*crossSection()*_sumW/sumOfWeights(), false);  // d01-x01-y01
      normalize(_h_dphi_BB, 0.53*crossSection()*_sumW/sumOfWeights(), false);   // d02-x01-y01
      normalize(_h_min_dR_ZB, 0.84*crossSection()*_sumW/sumOfWeights(), false); // d03-x01-y01
      normalize(_h_A_ZBB, 0.2*crossSection()*_sumW/sumOfWeights(), false);  // d04-x01-y01

      normalize(_h_dR_BB_boost, 0.84*crossSection()*_sumW50/sumOfWeights(), false); // d05-x01-y01
      normalize(_h_dphi_BB_boost, 0.63*crossSection()*_sumW50/sumOfWeights(), false);   // d06-x01-y01
      normalize(_h_min_dR_ZB_boost, 1*crossSection()*_sumW50/sumOfWeights(), false);    // d07-x01-y01
      normalize(_h_A_ZBB_boost, 0.25*crossSection()*_sumW50/sumOfWeights(), false); // d08-x01-y01

      normalize(_h_min_ZpT, 40*crossSection()*_sumWpT/sumOfWeights(), false);   // d09-x01-y01
    }


  private:

    /// @name Weight counters
    //@{
    double _sumW, _sumW50, _sumWpT;
    //@}

    /// @name Histograms
    //@{
    Histo1DPtr _h_dphi_BB, _h_dR_BB, _h_min_dR_ZB, _h_A_ZBB;
    Histo1DPtr _h_dphi_BB_boost, _h_dR_BB_boost, _h_min_dR_ZB_boost, _h_A_ZBB_boost, _h_min_ZpT;
    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2013_I1256943);

}
