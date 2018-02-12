// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"

namespace Rivet {


  /// Rivet analysis class for ATLAS min bias event shapes
  class ATLAS_2012_I1124167 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1124167()
      : Analysis("ATLAS_2012_I1124167") {  }


    /// Initialization, called once before running
    void init() {
      // Projections
      ChargedFinalState cfs(-2.5, 2.5, 0.5*GeV);
      declare(cfs, "CFS");
      declare(Sphericity(cfs), "Sphericity");

      // Book histograms
      _hist_T_05_25 = bookHisto1D(1,1,1);
      _hist_T_05    = bookHisto1D(2,1,1);
      _hist_T_25_50 = bookHisto1D(1,1,2);
      _hist_T_25    = bookHisto1D(2,1,2);
      _hist_T_50_75 = bookHisto1D(1,1,3);
      _hist_T_50    = bookHisto1D(2,1,3);
      _hist_T_75_100= bookHisto1D(1,1,4);
      _hist_T_75    = bookHisto1D(2,1,4);
      _hist_T_100   = bookHisto1D(2,1,5);

      _hist_TM_05_25 = bookHisto1D(3,1,1);
      _hist_TM_05    = bookHisto1D(4,1,1);
      _hist_TM_25_50 = bookHisto1D(3,1,2);
      _hist_TM_25    = bookHisto1D(4,1,2);
      _hist_TM_50_75 = bookHisto1D(3,1,3);
      _hist_TM_50    = bookHisto1D(4,1,3);
      _hist_TM_75_100= bookHisto1D(3,1,4);
      _hist_TM_75    = bookHisto1D(4,1,4);
      _hist_TM_100   = bookHisto1D(4,1,5);

      _hist_S_05_25 = bookHisto1D(5,1,1);
      _hist_S_05    = bookHisto1D(6,1,1);
      _hist_S_25_50 = bookHisto1D(5,1,2);
      _hist_S_25    = bookHisto1D(6,1,2);
      _hist_S_50_75 = bookHisto1D(5,1,3);
      _hist_S_50    = bookHisto1D(6,1,3);
      _hist_S_75_100= bookHisto1D(5,1,4);
      _hist_S_75    = bookHisto1D(6,1,4);
      _hist_S_100   = bookHisto1D(6,1,5);


      _hist_T_N  = bookProfile1D(7,1,1);
      _hist_TM_N = bookProfile1D(7,1,2);
      _hist_S_N  = bookProfile1D(7,1,3);

      _hist_T_S  = bookProfile1D(8,1,1);
      _hist_TM_S = bookProfile1D(8,1,2);
      _hist_S_S  = bookProfile1D(8,1,3);
    }


    void analyze(const Event& event) {
      const double weight = event.weight();

      // CFS projection and particles
      const ChargedFinalState& cfs500 = apply<ChargedFinalState>(event, "CFS");
      ParticleVector particles500 = cfs500.particlesByPt();

      // Require at least 6 charged particles
      if (cfs500.size() < 6) vetoEvent;

      // Preparation for Thrust calculation
      vector<Vector3> momenta;

      // Counters
      double num500 = 0;
      double ptSum500 = 0;

      double pTlead = particles500[0].pT()/GeV;

      // Loop over particles
      foreach (const Particle& p, particles500) {
        num500 += 1;
        ptSum500 += p.pT()/GeV;

        // Transverse Thrust calculation requires p_z to be set to 0
        Vector3 mom = p.p3();
        mom.setZ(0.0);
        momenta.push_back(mom);
      }

      // If only 2 particles, we need to use a ghost so that Thrust.calc() doesn't return 1.
      if (momenta.size() == 2) {
        momenta.push_back(Vector3(1e-10*MeV, 0., 0.));
      }

      // Actual thrust calculation
      Thrust thrust;
      thrust.calc(momenta);

      const double T  = 1.0 - thrust.thrust();
      const double TM = thrust.thrustMajor();

      Sphericity sphericity;
      sphericity.calc(momenta);

      const double S = sphericity.transSphericity();

      // Fill histos, most inclusive first

      // pTlead > 0.5
      _hist_T_05->fill(T , weight);
      _hist_TM_05->fill(TM, weight);
      _hist_S_05->fill(S , weight);

      // pTlead 0.5 - 2.5
      if (pTlead <= 2.5) {
        _hist_T_05_25->fill(T , weight);
        _hist_TM_05_25->fill(TM, weight);
        _hist_S_05_25->fill(S , weight);
      }

      // pTlead > 2.5
      if (pTlead > 2.5) {
        _hist_T_25->fill(T , weight);
        _hist_TM_25->fill(TM, weight);
        _hist_S_25->fill(S , weight);
      }

      // pTlead 2.5 - .5
      if (inRange(pTlead, 2.5, 5.0)) {
        _hist_T_25_50->fill(T , weight);
        _hist_TM_25_50->fill(TM, weight);
        _hist_S_25_50->fill(S , weight);
      }

      // pTlead > 5
      if (pTlead > 5) {
        _hist_T_50->fill(T , weight);
        _hist_TM_50->fill(TM, weight);
        _hist_S_50->fill(S , weight);
      }

      // pTlead 5 - 7.5
      if (inRange(pTlead, 5.0, 7.5)) {
        _hist_T_50_75->fill(T , weight);
        _hist_TM_50_75->fill(TM, weight);
        _hist_S_50_75->fill(S , weight);
      }

      // pTlead > 7.5
      if (pTlead > 7.5) {
        _hist_T_75->fill(T , weight);
        _hist_TM_75->fill(TM, weight);
        _hist_S_75->fill(S , weight);
      }

      // pTlead 7.5 - 10
      if (inRange(pTlead, 7.5, 10)) {
        _hist_T_75_100->fill(T , weight);
        _hist_TM_75_100->fill(TM, weight);
        _hist_S_75_100->fill(S , weight);
      }

      // pTlead > 10
      if (pTlead > 10) {
        _hist_T_100->fill(T , weight);
        _hist_TM_100->fill(TM, weight);
        _hist_S_100->fill(S , weight);
      }


      // Profiles Nch vs. ES
      _hist_T_N->fill(num500, T, weight);
      _hist_TM_N->fill(num500, TM, weight);
      _hist_S_N->fill(num500, S, weight);

      // Profiles pTsum vs. ES
      _hist_T_S->fill(ptSum500, T, weight);
      _hist_TM_S->fill(ptSum500, TM, weight);
      _hist_S_S->fill(ptSum500, S, weight);
    }


    void finalize() {
      normalize(_hist_T_05_25);
      normalize(_hist_T_05);
      normalize(_hist_T_25_50);
      normalize(_hist_T_25);
      normalize(_hist_T_50_75);
      normalize(_hist_T_50);
      normalize(_hist_T_75_100);
      normalize(_hist_T_75);
      normalize(_hist_T_100);

      normalize(_hist_TM_05_25);
      normalize(_hist_TM_05);
      normalize(_hist_TM_25_50);
      normalize(_hist_TM_25);
      normalize(_hist_TM_50_75);
      normalize(_hist_TM_50);
      normalize(_hist_TM_75_100);
      normalize(_hist_TM_75);
      normalize(_hist_TM_100);

      normalize(_hist_S_05_25);
      normalize(_hist_S_05);
      normalize(_hist_S_25_50);
      normalize(_hist_S_25);
      normalize(_hist_S_50_75);
      normalize(_hist_S_50);
      normalize(_hist_S_75_100);
      normalize(_hist_S_75);
      normalize(_hist_S_100);
    }


  private:

    Histo1DPtr _hist_T_05_25, _hist_T_05, _hist_T_25_50, _hist_T_25, _hist_T_50_75, _hist_T_50, _hist_T_75_100, _hist_T_75, _hist_T_100;
    Histo1DPtr _hist_TM_05_25, _hist_TM_05, _hist_TM_25_50, _hist_TM_25, _hist_TM_50_75, _hist_TM_50, _hist_TM_75_100, _hist_TM_75, _hist_TM_100;
    Histo1DPtr _hist_S_05_25, _hist_S_05, _hist_S_25_50, _hist_S_25, _hist_S_50_75, _hist_S_50, _hist_S_75_100, _hist_S_75, _hist_S_100;
    Profile1DPtr _hist_T_N, _hist_TM_N, _hist_S_N;
    Profile1DPtr _hist_T_S, _hist_TM_S, _hist_S_S;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1124167);

}
