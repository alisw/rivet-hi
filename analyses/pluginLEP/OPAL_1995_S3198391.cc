// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief OPAL Delta++ fragmentation function paper
  /// @author Peter Richardson
  class OPAL_1995_S3198391 : public Analysis {
  public:

    /// Constructor
    OPAL_1995_S3198391()
      : Analysis("OPAL_1995_S3198391")
    {}


    /// @name Analysis methods
    //@{

    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableFinalState(), "UFS");
      _histXpDelta   = bookHisto1D( 1, 1, 1);
    }


    void analyze(const Event& e) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get event weight for histo filling
      const double weight = e.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");

      foreach (const Particle& p, ufs.particles()) {
        if(p.abspid()==2224) {
          double xp = p.p3().mod()/meanBeamMom;
          _histXpDelta->fill(xp, weight);
        }
      }
    }


    /// Finalize
    void finalize() {
      scale(_histXpDelta, 1./sumOfWeights());
    }

    //@}


  private:

      Histo1DPtr _histXpDelta;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_1995_S3198391);

}
