// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BABAR_2013_I1116411 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BABAR_2013_I1116411);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableFinalState(), "UFS");

      // Book histograms
      _h_q2 = bookHisto1D(1, 1, 1);

    }
    
    // Calculate the Q2 using mother and daughter charged lepton
    double q2(const Particle& B) {
      const Particle chlept = filter_select(B.children(), Cuts::pid==PID::POSITRON || Cuts::pid==PID::ANTIMUON)[0];
      FourMomentum q = B.mom() - chlept.mom();
      return q*q;
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get B+ Mesons
      foreach(const Particle& p, apply<UnstableFinalState>(event, "UFS").particles(Cuts::pid==PID::BPLUS)) {
        if (isSemileptonicDecay(p, {PID::OMEGA, PID::POSITRON, PID::NU_E}) ||
            isSemileptonicDecay(p, {PID::OMEGA, PID::ANTIMUON, PID::NU_MU})) {
            _h_q2->fill(q2(p), event.weight());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_q2, 1.21); // normalize to BF

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_q2;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BABAR_2013_I1116411);


}
