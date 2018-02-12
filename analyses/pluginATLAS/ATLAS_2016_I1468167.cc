// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// Measurement of the inelastic proton-proton cross-section at \sqrt{s} = 13 TeV
  class ATLAS_2016_I1468167 : public Analysis {
  public:


    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1468167);


    /// Initialisation
    void init() {
      declare(FinalState(), "FS");
      _h_sigma = bookHisto1D(1, 1, 1);
    }


    /// Per-event analysis
    void analyze(const Event& event) {

      // Get all particles, sorted from minus to plus in eta
      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.size() < 2) vetoEvent; // need at least two particles to calculate gaps
      const Particles particles = fs.particles(cmpMomByEta);

      // Find this event's largest gap size and center
      double etapre = particles.front().eta();
      double gapcenter = 0.;
      double gapsize = -1;
      for (const Particle& p : particles) {
        const double gap = fabs(p.eta() - etapre);
        if (gap > gapsize) { // new largest gap
          gapsize = gap;
          gapcenter = (p.eta() + etapre)/2.;
        }
        etapre = p.eta();
      }

      // Calculate xi variable of the more massive side of the event, and apply xi cut
      FourMomentum mxFourVector, myFourVector;
      for (const Particle& p : particles) {
        ((p.eta() > gapcenter) ? mxFourVector : myFourVector) += p;
      }
      const double M2 = max(mxFourVector.mass2(), myFourVector.mass2());
      const double xi = M2/sqr(sqrtS()); // sqrt(s)=7000 GeV, note that units cancel
      if (xi < 1e-6) vetoEvent;

      // Fill the histogram
      _h_sigma->fill(sqrtS()/GeV, event.weight());
    }


    /// Scale the acceptance histogram to inelastic cross-section
    void finalize() {
      scale(_h_sigma, crossSection()/millibarn/sumOfWeights());
    }


    /// Histogram
    Histo1DPtr _h_sigma;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1468167);

}
