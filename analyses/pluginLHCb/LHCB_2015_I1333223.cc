// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Math/Units.hh"
#include <vector>

using namespace std;

namespace Rivet {


  class  LHCB_2015_I1333223 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
     LHCB_2015_I1333223()
      : Analysis("LHCB_2015_I1333223")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Charged particles
      declare(ChargedFinalState(Cuts::eta> 2.0 && Cuts::eta <4.5 && Cuts::pT >0.2*GeV), "CFS");
      // Reproducing only measurement for prompt charged particles
      _hInelasticXs = bookHisto1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double             weight  = event.weight();
      const ChargedFinalState   &cfs    = apply<ChargedFinalState> (event, "CFS");

      // eliminate non-inelastic events and empty events in LHCb
      if (cfs.particles().size() == 0) vetoEvent;

      // See if this event has at least one prompt particle
      foreach (const Particle &myp, cfs.particles()){
          double dPV = getPVDCA(myp);
          // if IP > 200 microns the particle is not considered prompt
          if ((dPV < 0.) || (dPV > 0.2 * millimeter)) {
            MSG_DEBUG(" Vetoing " << myp.pdgId() << " at " << dPV);
            continue;
          }
          // histo gets filled only for inelastic events (at least one prompt charged particle)
          _hInelasticXs->fill(sqrtS(), weight);
          break;
      } //end loop on particles

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_hInelasticXs, crossSection()/sumOfWeights()/millibarn);
    }

    //@}


  private:

    /*
     * Compute Distance of Closest Approach in z range for one particle. 
     * Assuming length unit is mm.
     * Returns -1. if unable to compute the DCA to PV.
     */
    double getPVDCA(const Particle& p) {
      const HepMC::GenVertex* vtx = p.genParticle()->production_vertex();
      if ( 0 == vtx ) return -1.;
      
      // Unit vector of particle's MOMENTUM three vector
      const Vector3 u = p.momentum().p3().unit();
      
      // The interaction point is always at (0, 0,0,0) hence the
      // vector pointing from the PV to the particle production vertex is:
      Vector3 d(vtx->position().x(), vtx->position().y(), vtx->position().z());

      // Subtract projection of d onto u from d
      double proj = d.dot(u);
      d -= (u * proj);

      // d should be orthogonal to u and it's length give the distance of 
      // closest approach
      return d.mod();
    }


    /// @name Histograms
    //@{
    Histo1DPtr _hInelasticXs;
    //@}
    //
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2015_I1333223);

}
