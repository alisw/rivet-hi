// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  class CMSTOTEM_2014_I1294140 : public Analysis {
  public:

    CMSTOTEM_2014_I1294140()
      : Analysis("CMSTOTEM_2014_I1294140")
    {     }


    void init() {
      ChargedFinalState cfs(-7.0, 7.0, 0.0*GeV);
      declare(cfs, "CFS");

      _Nevt_after_cuts_or = 0;
      _Nevt_after_cuts_and = 0;
      _Nevt_after_cuts_xor = 0;

      if (fuzzyEquals(sqrtS(), 8000*GeV, 1E-3)) {
        _h_dNch_dEta_OR = bookHisto1D(1, 1, 1);
        _h_dNch_dEta_AND = bookHisto1D(2, 1, 1);
        _h_dNch_dEta_XOR = bookHisto1D(3, 1, 1);
      }
    }


    void analyze(const Event& event) {
      // Count forward and backward charged particles
      const ChargedFinalState& charged = apply<ChargedFinalState>(event, "CFS");
      int count_plus = 0, count_minus = 0;
      foreach (const Particle& p, charged.particles()) {
        if (inRange(p.eta(),  5.3,  6.5)) count_plus++;
        if (inRange(p.eta(), -6.5, -5.3)) count_minus++;
      }

      // Cut combinations
      const bool cutsor  = (count_plus > 0 || count_minus > 0);
      const bool cutsand = (count_plus > 0 && count_minus > 0);
      const bool cutsxor = ( (count_plus > 0 && count_minus == 0) || (count_plus == 0 && count_minus > 0) );

      // Increment counters and fill histos
      const double weight = event.weight();
      if (cutsor)  _Nevt_after_cuts_or  += weight;
      if (cutsand) _Nevt_after_cuts_and += weight;
      if (cutsxor) _Nevt_after_cuts_xor += weight;
      foreach (const Particle& p, charged.particles()) {
        if (cutsor)  _h_dNch_dEta_OR ->fill(p.abseta(), weight);
        if (cutsand) _h_dNch_dEta_AND->fill(p.abseta(), weight);
        if (cutsxor) _h_dNch_dEta_XOR->fill(p.abseta(), weight);
      }

    }


    void finalize() {
      scale(_h_dNch_dEta_OR,  0.5/_Nevt_after_cuts_or);
      scale(_h_dNch_dEta_AND, 0.5/_Nevt_after_cuts_and);
      scale(_h_dNch_dEta_XOR, 0.5/_Nevt_after_cuts_xor);
    }


  private:

    Histo1DPtr _h_dNch_dEta_OR, _h_dNch_dEta_AND, _h_dNch_dEta_XOR;
    double _Nevt_after_cuts_or, _Nevt_after_cuts_and, _Nevt_after_cuts_xor;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMSTOTEM_2014_I1294140);

}
