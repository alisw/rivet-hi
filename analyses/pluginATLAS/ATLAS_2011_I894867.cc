// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class ATLAS_2011_I894867 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2011_I894867);

    void init() {
      declare(FinalState(), "FS");
      _h_sigma = bookHisto1D(1, 1, 1);
    }


    void analyze(const Event& event) {

      const FinalState& fs = apply<FinalState>(event, "FS");
      if (fs.size() < 2) vetoEvent; // need at least two particles to calculate gaps

      const Particles particles = fs.particles(cmpMomByEta);
      double etaprev = particles.front().eta();
      double gapcenter = etaprev;
      double detamax = -1;
      for (const Particle& p : particles) { // sorted from minus to plus
        const double deta = p.eta() - etaprev; // guaranteed positive
        if (deta > detamax) { // largest gap so far
          detamax = deta;
          gapcenter = (p.eta() + etaprev)/2.; // find the center of the gap to separate the X and Y systems.
        }
        etaprev = p.eta();
      }

      FourMomentum mxFourVector, myFourVector;
      for (const Particle& p : particles)
        (p.eta() > gapcenter ? mxFourVector : myFourVector) += p.momentum();
      const double m2 = max(mxFourVector.mass2(), myFourVector.mass2());
      const double xi = m2/sqr(sqrtS()); // sqrt(s) = 7000 GeV
      if (xi < 5e-6) vetoEvent;

      _h_sigma->fill(sqrtS()/GeV, event.weight());
    }


    void finalize() {
      scale(_h_sigma, crossSection()/millibarn/sumOfWeights());
    }


    Histo1DPtr _h_sigma;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I894867);

}
