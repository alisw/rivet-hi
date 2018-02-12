// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2013_I1238273 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2013_I1238273);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableFinalState(), "UFS");

      // Book histograms
      _h_q2_B0bar_pi     = bookHisto1D(1, 1, 1);
      _h_q2_B0bar_rho    = bookHisto1D(3, 1, 1);
      _h_q2_Bminus_pi    = bookHisto1D(2, 1, 1);
      _h_q2_Bminus_rho   = bookHisto1D(4, 1, 1);
      _h_q2_Bminus_omega = bookHisto1D(5, 1, 1);

    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - filter_select(B.children(), Cuts::pid==mesonID)[0];
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
      // Loop over B0bar Mesons
      foreach(const Particle& p, apply<UnstableFinalState>(event, "UFS").particles(Cuts::pid==PID::B0BAR)) {
        if (isSemileptonicDecay(p, {PID::PIPLUS, PID::ELECTRON, PID::NU_EBAR}) ||
            isSemileptonicDecay(p, {PID::PIPLUS, PID::MUON,     PID::NU_MUBAR})) {
            _h_q2_B0bar_pi->fill(q2(p, PID::PIPLUS), event.weight());
        }
        if (isSemileptonicDecay(p, {PID::RHOPLUS, PID::ELECTRON, PID::NU_EBAR}) ||
            isSemileptonicDecay(p, {PID::RHOPLUS, PID::MUON,     PID::NU_MUBAR})) {
            _h_q2_B0bar_rho->fill(q2(p, PID::RHOPLUS), event.weight());
        }
      }
      // Loop over B- Mesons
      foreach(const Particle& p, apply<UnstableFinalState>(event, "UFS").particles(Cuts::pid==PID::BMINUS)) {
        if (isSemileptonicDecay(p, {PID::PI0, PID::ELECTRON, PID::NU_EBAR}) ||
            isSemileptonicDecay(p, {PID::PI0, PID::MUON,     PID::NU_MUBAR})) {
            _h_q2_Bminus_pi->fill(q2(p, PID::PI0), event.weight());
        }
        if (isSemileptonicDecay(p, {PID::RHO0, PID::ELECTRON, PID::NU_EBAR}) ||
            isSemileptonicDecay(p, {PID::RHO0, PID::MUON,    PID::NU_MUBAR})) {
            _h_q2_Bminus_rho->fill(q2(p,PID::RHO0), event.weight());
        }
        if (isSemileptonicDecay(p, {PID::OMEGA, PID::ELECTRON, PID::NU_EBAR}) ||
            isSemileptonicDecay(p, {PID::OMEGA, PID::MUON,     PID::NU_MUBAR})) {
            _h_q2_Bminus_omega->fill(q2(p, PID::OMEGA), event.weight());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_q2_B0bar_pi    , 298.8);  // normalize to BF*dQ2
      normalize(_h_q2_B0bar_rho   , 1304.8); // normalize to BF*dQ2
      normalize(_h_q2_Bminus_pi   , 324.8);  // normalize to BF*dQ2
      normalize(_h_q2_Bminus_rho  , 367.0);  // normalize to BF*dQ2
      normalize(_h_q2_Bminus_omega, 793.1);  // normalize to BF*dQ2

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_q2_B0bar_pi    ;
    Histo1DPtr _h_q2_B0bar_rho   ;
    Histo1DPtr _h_q2_Bminus_pi   ;
    Histo1DPtr _h_q2_Bminus_rho  ;
    Histo1DPtr _h_q2_Bminus_omega;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2013_I1238273);


}
