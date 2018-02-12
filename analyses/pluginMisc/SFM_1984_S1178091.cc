// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief SFM charged multiplicities in NSD and inelastic minbias events
  class SFM_1984_S1178091 : public Analysis {
  public:

    /// Constructor
    SFM_1984_S1178091() : Analysis("SFM_1984_S1178091") {}


    /// @name Analysis methods
    //@{

    void init() {
      // Projections
      // 
      declare(ChargedFinalState(Cuts::absrap<5 && Cuts::pT>250*MeV && Cuts::pT<3*GeV), "FS");

      // Histograms
      if (fuzzyEquals(sqrtS()/GeV, 30.4, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 1);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 1);
      } else if (fuzzyEquals(sqrtS(), 44.5, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 2);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 2);
      } else if (fuzzyEquals(sqrtS(), 52.2, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 3);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 3);
      } else if (fuzzyEquals(sqrtS(), 62.2, 1E-1)) {
        _hist_multiplicity_inel = bookHisto1D(1, 1, 4);
        _hist_multiplicity_nsd = bookHisto1D(2, 1, 4);
      }

    }


    // Analyse each event
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");

      // Trigger
      if (fs.particles().size() <1 ) vetoEvent;

      // Event classification: 
      int n_left(0), n_right(0), n_large_x(0);
      foreach (const Particle& p, fs.particles()) {
        // Calculate the particles' Feynman x
        const double x_feyn = 2.0 * fabs(p.pz())/sqrtS();
        if (x_feyn > 0.8 ) n_large_x += 1;

        // Pseudorapidity
        const double eta = p.eta();
        if (eta > 0.0) n_right += 1;
        else if (eta < 0.0) n_left += 1;
      }
      MSG_DEBUG("N_left: " << n_left << ", "
                << "N_right: " << n_right << ", "
                << "N_large_x: " << n_large_x);

      
      // Single diffractive: either one large x particle or 0 particles in the one hemisphere but more than 7 in the other hemisphere
      bool isDiffractive = (n_large_x == 1) ||  ( ((n_left==0) && (fs.particles().size() < 7)) || ((n_right==0) && (fs.particles().size() < 7)) );


      _hist_multiplicity_inel->fill(fs.particles().size(), weight);
      if (!isDiffractive) _hist_multiplicity_nsd->fill(fs.particles().size(), weight);
    }


    void finalize() {
      normalize(_hist_multiplicity_inel);
      normalize(_hist_multiplicity_nsd);
    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _hist_multiplicity_inel;
    Histo1DPtr _hist_multiplicity_nsd;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SFM_1984_S1178091);

}
