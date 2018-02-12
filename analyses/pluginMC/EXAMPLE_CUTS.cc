// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Just measures a few random things as an example.
  class EXAMPLE_CUTS : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(EXAMPLE_CUTS);


    /// @name Analysis methods
    //@{

    /// Set up projections and book histograms
    void init() {

      // Projections
      const FinalState cnfs(Cuts::abseta < 4);
      declare(cnfs, "FS");

      // Histograms
      _histPt         = bookHisto1D("pT", 30, 0, 30);
      _histMass       = bookHisto1D("Mass", 20, 0, 1);

    }


    /// Do the analysis
    void analyze(const Event& event) {
      // Make sure to always include the event weight in histogram fills!
      const double weight = event.weight();

      const Particles ps = apply<FinalState>(event, "FS").particlesByPt();

      Cut ptcut = Cuts::range(Cuts::pT, 5, 20);
      Cut masscut = Cuts::range(Cuts::mass, 0, 0.2);
      Cut combine = ptcut && masscut; //< Possible to combine cuts

      for (const Particle& p : ps) {
        if ( ptcut->accept(p) )
          _histPt->fill(p.pT(), weight);
        if ( combine->accept(p) )
          _histMass->fill(p.mass(), weight);
      }
    }


    /// Finalize
    void finalize() {
      normalize({_histPt, _histMass});
    }

    //@}


    //@{
    /// Histograms
    Histo1DPtr _histPt, _histMass;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EXAMPLE_CUTS);

}
