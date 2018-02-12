// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"


namespace Rivet {


  /// @brief SLD b-fragmentation measurement
  /// @author Peter Richardson
  class SLD_2002_S4869273 : public Analysis {
  public:

    /// Constructor
    SLD_2002_S4869273()
      : Analysis("SLD_2002_S4869273")
    {    }


    /// @name Helper functions
    /// @note The PID:: namespace functions would be preferable, but don't have exactly the same behaviour. Preserving the original form.
    //@{
    // bool isParton(int id) { return abs(id) <= 100 && abs(id) != 22 && (abs(id) < 11 || abs(id) > 18); }
    // bool isBHadron(int id) { return ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999); }
    //@}


    /// @name Analysis methods
    //@{

    /// Book projections and histograms
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      _histXbweak     = bookHisto1D(1, 1, 1);
    }


    void analyze(const Event& e) {
      const FinalState& fs = apply<FinalState>(e, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed ncharged cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed ncharged cut");

      // Get event weight for histo filling
      const double weight = e.weight();

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(e, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);


      for (const GenParticle* p : particles(e.genEvent())) {
        const GenVertex* dv = p->end_vertex();
        if (PID::isBottomHadron(p->pdg_id())) {
          const double xp = p->momentum().e()/meanBeamMom;

          // If the B-hadron has no B-hadron as a child, it decayed weakly:
          if (dv) {
            bool is_weak = true;
            for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin() ;
                 pp != dv->particles_out_const_end() ; ++pp) {
              if (PID::isBottomHadron((*pp)->pdg_id())) {
                is_weak = false;
              }
            }
            if (is_weak) {
              _histXbweak->fill(xp, weight);
            }
          }

        }
      }
    }


    // Finalize
    void finalize() {
      normalize(_histXbweak);
    }


  private:

    /// Store the weighted sums of numbers of charged / charged+neutral
    /// particles - used to calculate average number of particles for the
    /// inclusive single particle distributions' normalisations.

    Histo1DPtr _histXbweak;

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SLD_2002_S4869273);

}
