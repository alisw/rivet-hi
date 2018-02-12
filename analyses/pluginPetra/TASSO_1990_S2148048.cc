// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class TASSO_1990_S2148048 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TASSO_1990_S2148048);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(-MAXDOUBLE, MAXDOUBLE, 0.1/GeV);
      declare(cfs, "CFS");

      // Thrust and sphericity
      declare(Thrust(cfs), "Thrust");
      declare(Sphericity(cfs), "Sphericity");

      // Histos
      int offset = 0;
      switch (int(sqrtS()/GeV)) {
        case 14:
          offset = 0;
          break;
        case 22:
          offset = 1;
          break;
        case 35:
          offset = 2;
          break;
        case 44:
          offset = 3;
          break;
      }
      //_h_xp         = bookHisto1D( 2, 1, 1+offset);
      _h_sphericity = bookHisto1D( 6, 1, 1+offset);
      _h_aplanarity = bookHisto1D( 7, 1, 1+offset);
      _h_thrust     = bookHisto1D( 8, 1, 1+offset);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      //// Get beams and average beam momentum
      //const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      //const double meanBeamMom = ( beams.first.p3().mod() +
                                   //beams.second.p3().mod() ) / 2.0;

      // TASSO hadronic event selection TODO: move this into a trigger definition
      // See page 2 in publication
      // Condition 1)  --- require at least 5 (4) 'good' tracks
      int nch = cfs.particles().size();
      if ( (int(sqrtS()/GeV) > 27 && nch < 5) || (int(sqrtS()/GeV) <= 27 && nch < 4 ) ) {
        MSG_DEBUG("Failed # good tracks cut: " << nch);
        vetoEvent;
      }
      // Condition 2) ---
      // Condition 5) --- scalar momentum (not pT!!!) sum >= 0.265*s
      double momsum = 0.0;
      for (const Particle& p : cfs.particles()) {
        const double mom = p.p3().mod();
        momsum += mom;
      }
      if (momsum <=0.265 * sqrtS()/GeV) {
        MSG_DEBUG("Failed pTsum cut: " << momsum << " < " << 0.265 * sqrtS()/GeV);
        vetoEvent;
      }

      // Raise counter for events that pass trigger conditions
      //_sumWPassed += event.weight();

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      //const Vector3 & thrustAxis = thrust.thrustAxis ();
      //double theta = thrustAxis.theta();
      //if ( fabs(cos(theta)) >= 0.8 ) {
        //MSG_DEBUG("Failed thrust angle cut: " << fabs(cos(theta)));
        //vetoEvent;
      //}

      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      //// Fill histograms in order of appearance in paper
      //foreach (const Particle& p, cfs.particles()) {
        //// Get momentum and energy of each particle.
        //const Vector3 mom3 = p.p3();
        //// Scaled momenta.
        //const double mom = mom3.mod();
        //const double scaledMom = mom/meanBeamMom;
        //_h_xp->fill(scaledMom, weight);
      //}
      //
      _h_sphericity->fill(sphericity.sphericity(), weight);
      _h_aplanarity->fill(sphericity.aplanarity(), weight);
      _h_thrust->fill(thrust.thrust(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      //scale(_h_xp, _sumWPassed/(crossSection()*sumOfWeights()));
      normalize(_h_sphericity);
      normalize(_h_aplanarity);
      normalize(_h_thrust    );
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_xp, _h_sphericity, _h_aplanarity, _h_thrust;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TASSO_1990_S2148048);

}
