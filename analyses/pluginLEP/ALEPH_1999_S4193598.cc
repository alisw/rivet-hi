// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {

  class ALEPH_1999_S4193598 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ALEPH_1999_S4193598()
      : Analysis("ALEPH_1999_S4193598")
    { }

    //@}


  public:

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(UnstableFinalState(), "UFS");
      declare(ChargedFinalState(), "CFS");

      _h_Xe_Ds = bookHisto1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Trigger condition
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      if (cfs.size() < 5) vetoEvent;

      const UnstableFinalState& ufs = apply<UnstableFinalState>(event, "UFS");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0/GeV;

      // Accept all D*+- decays. Normalisation to data in finalize
      for (const Particle& p : filter_select(ufs.particles(), Cuts::abspid==PID::DSTARPLUS)) {
          // Scaled energy.
          const double energy = p.E()/GeV;
          const double scaledEnergy = energy/meanBeamMom;
          _h_Xe_Ds->fill(scaledEnergy, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      // Normalize to data integral
      normalize(_h_Xe_Ds, 0.00498);
    }

  private:

    Histo1DPtr _h_Xe_Ds;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALEPH_1999_S4193598);

}
