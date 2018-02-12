// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief ATLAS 13 TeV minimum bias analysis for low-pT tracks
  class ATLAS_2016_I1467230 : public Analysis {
  public:

    /// Particle types included
    enum PartTypes {
      k_NoStrange,
      k_AllCharged,
      kNPartTypes
    };

    /// Phase space regions
    enum regionID {
      k_pt100_nch2_eta25,
      kNregions
    };


    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1467230);


    /// Initialization, called once before running
    void init() {

      for (int iT = 0; iT < kNPartTypes; ++iT) {
        for (int iR = 0; iR < kNregions; ++iR) {
          _sumW[iT][iR] = 0.;
        }
      }

      // Initialize and register projections
      declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 100*MeV), "CFS100_25");

      for (int iT = 0; iT < kNPartTypes; ++iT) {
        for (int iR = 0; iR < kNregions; ++iR) {
          _hist_nch  [iT][iR] = bookHisto1D  ( 1, iR + 1, iT + 1);
          _hist_pt   [iT][iR] = bookHisto1D  ( 2, iR + 1, iT + 1);
          _hist_eta  [iT][iR] = bookHisto1D  ( 3, iR + 1, iT + 1);
          _hist_ptnch[iT][iR] = bookProfile1D( 4, iR + 1, iT + 1);
        }
      }

    }


    /// Fill histograms for the given particle selection and phase-space region
    void fillPtEtaNch(const Particles& particles, int ptype, int iRegion, double weight) {

      // Skip if event fails multiplicity cut
      const size_t nch = particles.size();
      if (nch < 2) return;

      // Fill event weight info
      _sumW[ptype][iRegion] += weight;

      // Fill nch
      _hist_nch[ptype][iRegion]->fill(nch, weight);

      // Loop over particles, fill pT, eta and ptnch
      for (const Particle& p : particles)  {
        const double pt  = p.pT()/GeV;
        const double eta = p.eta();
        _hist_pt   [ptype][iRegion]->fill(pt , weight/pt);
        _hist_eta  [ptype][iRegion]->fill(eta, weight);
        _hist_ptnch[ptype][iRegion]->fill(nch, pt, weight);
      }
    }


    /// Per-event analysis
    void analyze(const Event& event) {

      // Get all charged particles
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS100_25");
      const Particles& pall = cfs.particles();

      // Get charged particles, filtered to omit charged strange baryons
      const Cut& pcut = Cuts::abspid != PID::SIGMAMINUS && Cuts::abspid != PID::SIGMAPLUS && Cuts::abspid != PID::XIMINUS && Cuts::abspid != PID::OMEGAMINUS;
      const Particles& pnostrange = cfs.particles(pcut);

      // Fill all histograms
      for (int iR = 0; iR < kNregions; ++iR)  {
        fillPtEtaNch(pall,       k_AllCharged, iR, event.weight());
        fillPtEtaNch(pnostrange, k_NoStrange,  iR, event.weight());
      }

    }


    /// Post-run data manipulation
    void finalize() {

      // Scale all histograms
      for (int iT = 0; iT < kNPartTypes; ++iT) {
        for (int iR = 0; iR < kNregions; ++iR) {
          if (_sumW[iT][iR] > 0) {
            scale(_hist_nch[iT][iR], 1.0/_sumW[iT][iR]);
            scale(_hist_pt [iT][iR], 1.0/_sumW[iT][iR]/TWOPI/5.);
            scale(_hist_eta[iT][iR], 1.0/_sumW[iT][iR]);
          }
        }
      }

    }


  private:

    /// Weight sums
    double _sumW[kNPartTypes][kNregions];

    /// @name Histogram arrays
    //@{
    Histo1DPtr   _hist_nch    [kNPartTypes][kNregions];
    Histo1DPtr   _hist_pt     [kNPartTypes][kNregions];
    Histo1DPtr   _hist_eta    [kNPartTypes][kNregions];
    Profile1DPtr _hist_ptnch  [kNPartTypes][kNregions];
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1467230);

}
