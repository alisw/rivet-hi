// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Per-event yield of the highest transverse momentum charged particle and charged-particle jet
  class CMS_2015_I1380605 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1380605);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(-7., 7., 0.0*GeV);
      declare(cfs, "CFS");
      declare(FastJets(cfs, FastJets::ANTIKT, 0.5), "Jets");

      _h_tracks = bookHisto1D(1, 1, 1);
      _h_jets = bookHisto1D(2, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Veto events without forward activity on both sides
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      const size_t count_plus  = cfs.particles(Cuts::eta > 5.3  && Cuts::eta < 6.5).size();
      const size_t count_minus = cfs.particles(Cuts::eta < -5.3 && Cuts::eta > -6.5).size();
      if (!(count_plus > 0 || count_minus > 0)) vetoEvent; //< @todo Should this really also veto the jet analysis?
      /// @warning Needs migration to an AO Counter
      /// @note This isn't the number of tracks, it's the sum of event weights passing the veto
      const double weight = event.weight();
      _ntracks += weight;

      // Do track analysis here
      // Find pttrackmax
      double track_ptmax = 0;
      for (const Particle& p : cfs.particles(Cuts::abseta < 2.4)) track_ptmax = max(track_ptmax, p.pT());
      // Fill track analysis histograms
      for (size_t i = 0; i < _h_tracks->numBins(); ++i) {
        const double binlimitlow_t = _h_tracks->bin(i).xMin();
        const double weightbw_t = weight * _h_tracks->bin(i).xWidth();
        const double xbin_t = _h_tracks->bin(i).xMid();
        if (track_ptmax > binlimitlow_t) _h_tracks -> fill(xbin_t, weightbw_t);
      }

      // Do jet analysis here
      const Jets jetsdeta = apply<FastJets>(event,"Jets").jets(Cuts::pT > 1*GeV && Cuts::pT < 60*GeV && Cuts::abseta < 1.9);
      // Find ptjetmax
      double jet_ptmax = 0;
      for (const Jet& j : jetsdeta) jet_ptmax = max(jet_ptmax, j.pT());
      // Fill jet analysis histograms
      for (size_t i = 0; i < _h_jets->numBins(); ++i) {
        const double binlimitlow_j = _h_jets->bin(i).xMin();
        const double weightbw_j = weight * _h_jets->bin(i).xWidth();
        const double xbin_j = _h_jets->bin(i).xMid();
        if (jet_ptmax > binlimitlow_j) _h_jets -> fill(xbin_j, weightbw_j);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double norm_t0 = _h_tracks->bin(7).height()/2.056170e-03;
      const double norm_t1 = _h_tracks->bin(7).sumW()/2.056170e-03;
      const double norm_j0 = _h_jets->bin(13).height()/3.575290e-03;
      const double norm_j1 = _h_jets->bin(13).sumW()/3.575290e-03;
      MSG_DEBUG("Norm track " << norm_t0 << " " << norm_t1);
      MSG_DEBUG("Norm  jets " << norm_j0 << " " << norm_j1);
      if (norm_t0 > 0 ) scale(_h_tracks, 1./ norm_t0);
      if (norm_j0 > 0 ) scale(_h_jets, 1./ norm_j0);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_tracks, _h_jets;
    double _ntracks = 0;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2015_I1380605);

}
