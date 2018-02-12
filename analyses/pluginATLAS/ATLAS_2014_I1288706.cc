// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Particle.fhh"

namespace Rivet {


  class ATLAS_2014_I1288706 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1288706()
      : Analysis("ATLAS_2014_I1288706")
    {
       _sumw_ext_mu_dressed = 0;
       _sumw_mu_dressed     = 0;
       _sumw_el_dressed     = 0;
     }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

       // Set up projections
       FinalState fs;

       ZFinder zfinder_ext_dressed_mu(fs, Cuts::abseta<2.4 && Cuts::pT>6.0*GeV, PID::MUON, 12.0*GeV, 66.0*GeV, 0.1);
       declare(zfinder_ext_dressed_mu, "ZFinder_ext_dressed_mu");

       ZFinder zfinder_dressed_mu(fs, Cuts::abseta<2.4 && Cuts::pT>12*GeV, PID::MUON, 26.0*GeV, 66.0*GeV, 0.1);
       declare(zfinder_dressed_mu, "ZFinder_dressed_mu");

       ZFinder zfinder_dressed_el(fs, Cuts::abseta<2.4 && Cuts::pT>12*GeV, PID::ELECTRON, 26.0*GeV, 66.0*GeV, 0.1);
       declare(zfinder_dressed_el, "ZFinder_dressed_el");

       _hist_ext_mu_dressed = bookHisto1D(1, 1, 1); // muon, dressed level, extended phase space
       _hist_mu_dressed	    = bookHisto1D(2, 1, 1); // muon, dressed level, nominal phase space
       _hist_el_dressed	    = bookHisto1D(2, 1, 2); // electron, dressed level, nominal phase space
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      const ZFinder& zfinder_ext_dressed_mu = apply<ZFinder>(event, "ZFinder_ext_dressed_mu");
      const ZFinder& zfinder_dressed_mu     = apply<ZFinder>(event, "ZFinder_dressed_mu"    );
      const ZFinder& zfinder_dressed_el     = apply<ZFinder>(event, "ZFinder_dressed_el"    );

      fillPlots(zfinder_ext_dressed_mu, _hist_ext_mu_dressed, 9*GeV, weight);
      fillPlots(zfinder_dressed_mu,     _hist_mu_dressed,    15*GeV, weight);
      fillPlots(zfinder_dressed_el,     _hist_el_dressed,    15*GeV, weight);
    }


    /// Helper function to fill the histogram if a Z is found with the required lepton cuts
    void fillPlots(const ZFinder& zfinder, Histo1DPtr hist, double leading_pT, double weight) {
      if (zfinder.bosons().size() != 1) return;
      const FourMomentum el1 = zfinder.rawParticles()[0];
      const FourMomentum el2 = zfinder.rawParticles()[1];
      if (el1.pT() < leading_pT && el2.pT() < leading_pT) return;
      const double mass = zfinder.bosons()[0].mass();
      hist->fill(mass/GeV, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hist_ext_mu_dressed, crossSection()/sumOfWeights());
      scale(_hist_mu_dressed,     crossSection()/sumOfWeights());
      scale(_hist_el_dressed,     crossSection()/sumOfWeights());
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here
    double _sumw_ext_mu_dressed, _sumw_mu_dressed, _sumw_el_dressed;

    /// Histograms
    Histo1DPtr _hist_ext_mu_dressed, _hist_mu_dressed, _hist_el_dressed;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1288706);

}
