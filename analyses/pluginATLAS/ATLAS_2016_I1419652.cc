// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ATLAS_2016_I1419652 : public Analysis {
  public:

    /// Particle types included
    enum PartTypes {
      k_NoStrange,
      k_AllCharged,
      kNPartTypes
    };

    /// Phase space regions
    enum RegionID {
      k_pt500_nch1_eta25,
      k_pt500_nch1_eta08,
      kNregions
    };

    /// Nch cut for each region
    const static int nchCut[kNregions];


    /// Default constructor
    ATLAS_2016_I1419652() : Analysis("ATLAS_2016_I1419652") {
      for (int iT = 0; iT < kNPartTypes; ++iT)  {
        for (int iR = 0; iR < kNregions; ++iR)  {
          _sumW[iT][iR]  = 0.;
        }
      }
    }


    /// Initialization, called once before running
    void init() {

      // Projections
      const ChargedFinalState cfs500_25(-2.5, 2.5, 500.0*MeV);
      declare(cfs500_25, "CFS500_25");

      const ChargedFinalState cfs500_08(-0.8, 0.8, 500.0*MeV);
      declare(cfs500_08, "CFS500_08");

      for (int iT = 0; iT < kNPartTypes; ++iT)  {
        for (int iR = 0; iR < kNregions; ++iR)  {
          _hist_nch  [iT][iR] = bookHisto1D  ( 1, iR + 1, iT + 1);
          _hist_pt   [iT][iR] = bookHisto1D  ( 2, iR + 1, iT + 1);
          _hist_eta  [iT][iR] = bookHisto1D  ( 3, iR + 1, iT + 1);
          _hist_ptnch[iT][iR] = bookProfile1D( 4, iR + 1, iT + 1);
        }
      }
    }


    void analyze(const Event& event) {
      string fsName;
      for (int iR = 0; iR < kNregions; ++iR)  {
        switch (iR) {
        case k_pt500_nch1_eta25:  fsName = "CFS500_25"; break;
        case k_pt500_nch1_eta08:  fsName = "CFS500_08"; break;
        }

        const ChargedFinalState& cfs = apply<ChargedFinalState>(event, fsName);

        /// What's the benefit in separating this code which is only called from one place?!
        fillPtEtaNch(cfs, iR, event.weight());
      }
    }



    void finalize() {
      // Standard histograms
      for (int iT = 0; iT < kNPartTypes; ++iT)  {
        for (int iR = 0; iR < kNregions; ++iR)  {

          double etaRangeSize = -999.0; //intentionally crazy
          switch (iR) {
            case k_pt500_nch1_eta25  : etaRangeSize = 5.0 ;  break;
            case k_pt500_nch1_eta08  : etaRangeSize = 1.6 ;  break;
            default: etaRangeSize = -999.0; break; //intentionally crazy
          }

          if (_sumW[iT][iR] > 0) {
            scale(_hist_nch[iT][iR], 1.0/_sumW[iT][iR]);
            scale(_hist_pt [iT][iR], 1.0/_sumW[iT][iR]/TWOPI/etaRangeSize);
            scale(_hist_eta[iT][iR], 1.0/_sumW[iT][iR]);
          } else {
            MSG_WARNING("Sum of weights is zero (!) in type/region: " << iT << " " << iR);
          }
        }
      }
    }


    /// Helper for collectively filling Nch, pT, eta, and pT vs. Nch histograms
    void fillPtEtaNch(const ChargedFinalState& cfs, int iRegion, double weight) {

      // Get number of particles
      int nch[kNPartTypes];
      int nch_noStrange = 0;
      foreach (const Particle& p, cfs.particles()) {
        PdgId pdg = p.abspid ();
        if ( pdg == 3112 || // Sigma-
             pdg == 3222 || // Sigma+
             pdg == 3312 || // Xi-
             pdg == 3334 )  // Omega-
	        continue;
	      nch_noStrange++;
      }
      nch[k_AllCharged] = cfs.size();
      nch[k_NoStrange ] = nch_noStrange;

      // Skip if event fails cut for all charged (noStrange will always be less)
      if (nch[k_AllCharged] < nchCut[iRegion]) return;

      // Fill event weight info
      _sumW[k_AllCharged][iRegion] += weight;
      if (nch[k_NoStrange ] >= nchCut[iRegion])	{
        _sumW[k_NoStrange][iRegion] += weight;
      }

      // Fill nch
      _hist_nch[k_AllCharged][iRegion]->fill(nch[k_AllCharged], weight);
      if (nch[k_NoStrange ] >= nchCut[iRegion])	{
        _hist_nch [k_NoStrange][iRegion]->fill(nch[k_NoStrange ], weight);
      }

      // Loop over particles, fill pT, eta and ptnch
      foreach (const Particle& p, cfs.particles())  {
        const double pt  = p.pT()/GeV;
        const double eta = p.eta();
        _hist_pt     [k_AllCharged][iRegion]->fill(pt , weight/pt);
        _hist_eta    [k_AllCharged][iRegion]->fill(eta, weight);
        _hist_ptnch  [k_AllCharged][iRegion]->fill(nch[k_AllCharged], pt, weight);

        // Make sure nch cut is passed for nonStrange particles!
        if (nch[k_NoStrange ] >= nchCut[iRegion])  {
          PdgId pdg = p.abspid ();
          if ( pdg == 3112 || // Sigma-
               pdg == 3222 || // Sigma+
               pdg == 3312 || // Xi-
               pdg == 3334 )  // Omega-
	          continue;
          // Here we don't have strange particles anymore
          _hist_pt   [k_NoStrange][iRegion]->fill(pt , weight/pt);
          _hist_eta  [k_NoStrange][iRegion]->fill(eta, weight);
          _hist_ptnch[k_NoStrange][iRegion]->fill(nch[k_NoStrange], pt, weight);
        }
      }
    }


  private:

    double _sumW[kNPartTypes][kNregions];

    Histo1DPtr   _hist_nch  [kNPartTypes][kNregions];
    Histo1DPtr   _hist_pt   [kNPartTypes][kNregions];
    Histo1DPtr   _hist_eta  [kNPartTypes][kNregions];
    Profile1DPtr _hist_ptnch[kNPartTypes][kNregions];

  };


  // Constants: pT & eta regions
  const int ATLAS_2016_I1419652::nchCut[] = {1, 1};


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1419652);

}
