// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ATLAS_2016_I1426695 : public Analysis {
  public:

    //phase space regions
    enum regionID {
      k_pt100_nch2 = 0,
      k_pt500_nch1,
      k_pt500_nch6,
      k_pt500_nch20,
      k_pt500_nch50,
      kNregions
    };
    
    
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1426695);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      for (int iR=0; iR < kNregions; ++iR)  {
        _sumW[iR]  = 0.;
      }

      // Initialise and register projections
      declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 100*MeV), "CFS_100");
      declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV), "CFS_500");

      // Book histograms
      for (int iR=0; iR < kNregions; ++iR)  {
        if (iR == k_pt100_nch2 || iR == k_pt500_nch1) {
          _hist_nch  [iR] = bookHisto1D  (1, iR + 1, 1);
          _hist_ptnch[iR] = bookProfile1D(4, iR + 1, 1);
        }
        _hist_pt [iR] = bookHisto1D(2, iR + 1, 1);
        _hist_eta[iR] = bookHisto1D(3, iR + 1, 1);
      }
    }

    void fillPtEtaNch(const Particles particles, int nMin, int iRegion, double weight) {

      //skip if event fails multiplicity cut
      int nch =particles.size();
      if (nch < nMin)  return;

      //fill event weight info
      _sumW[iRegion] += weight;

      // Fill nch
      if (iRegion == k_pt100_nch2 || iRegion == k_pt500_nch1) {
        _hist_nch[iRegion]->fill(nch, weight);
      }
     
      for (const Particle&p : particles) { 
      // Loop over particles, fill pT, eta and ptnch
        const double pt  = p.pT()/GeV;
        const double eta = p.eta();

        _hist_pt [iRegion]->fill(pt , weight/pt);
        _hist_eta[iRegion]->fill(eta, weight);

        if (iRegion == k_pt100_nch2 || iRegion == k_pt500_nch1) {
          _hist_ptnch[iRegion]->fill(nch, pt, weight);
        }
      } //end loop over particles
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get charged particles, omitting some strange heavies
      const Cut& pcut = (
          (Cuts::abspid!=PID::SIGMAMINUS) && (Cuts::abspid!=PID::SIGMAPLUS) &&
          (Cuts::abspid!=PID::XIMINUS)    && (Cuts::abspid!=PID::OMEGAMINUS));
      const Particles& p_100 = apply<ChargedFinalState>(event, "CFS_100").particles(pcut);
      const Particles& p_500 = apply<ChargedFinalState>(event, "CFS_500").particles(pcut);

      fillPtEtaNch(p_100,  2, 0, event.weight());
      fillPtEtaNch(p_500,  1, 1, event.weight());
      fillPtEtaNch(p_500,  6, 2, event.weight());
      fillPtEtaNch(p_500, 20, 3, event.weight());
      fillPtEtaNch(p_500, 50, 4, event.weight());
    }


    void finalize() {

      for (int iR = 0; iR < kNregions; ++iR)  {
        if (_sumW[iR] > 0) {
          if (iR == k_pt100_nch2 || iR == k_pt500_nch1) {
            scale(_hist_nch[iR], 1.0/_sumW[iR]);
          }
          scale(_hist_pt [iR], 1.0/_sumW[iR]/TWOPI/5.);
          scale(_hist_eta[iR], 1.0/_sumW[iR]);
        } 
      }
    }

    //@}


  private:

    double _sumW[kNregions];

    /// @name Histograms
    Histo1DPtr   _hist_nch    [kNregions];
    Histo1DPtr   _hist_pt     [kNregions];
    Histo1DPtr   _hist_eta    [kNregions];
    Profile1DPtr _hist_ptnch  [kNregions];

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1426695);


}
