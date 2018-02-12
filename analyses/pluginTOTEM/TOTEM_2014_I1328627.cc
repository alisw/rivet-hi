// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class TOTEM_2014_I1328627 : public Analysis {
  public:


    TOTEM_2014_I1328627()
      : Analysis("TOTEM_2014_I1328627")
    {    }



    void init() {
      ChargedFinalState cfsm(-7.0, -6.0, 0.0*GeV);
      ChargedFinalState cfsp( 3.7,  4.8, 0.0*GeV);
      declare(cfsm, "CFSM");
      declare(cfsp, "CFSP");

      _h_eta = bookHisto1D(1, 1, 1);
      _sumofweights = 0.;
    }


    void analyze(const Event& event) {
      const ChargedFinalState cfsm = apply<ChargedFinalState>(event, "CFSM");
      const ChargedFinalState cfsp = apply<ChargedFinalState>(event, "CFSP");
      if (cfsm.size() == 0 && cfsp.size() == 0) vetoEvent;

      _sumofweights += event.weight();
      foreach (const Particle& p, cfsm.particles() + cfsp.particles()) {
        _h_eta->fill(p.abseta(), event.weight());
      }
    }


    void finalize() {
      scale(_h_eta, 1./_sumofweights);
    }


  private:

    double _sumofweights;
    Histo1DPtr _h_eta;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TOTEM_2014_I1328627);


}
