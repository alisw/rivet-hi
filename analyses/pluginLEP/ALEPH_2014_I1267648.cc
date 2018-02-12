// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ALEPH_2014_I1267648 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALEPH_2014_I1267648);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableFinalState(), "UFS");

      // Book histograms
      _h_pip0  = bookHisto1D(1, 1, 1);
      _h_pi2p0 = bookHisto1D(2, 1, 1);
      _h_pi3p0 = bookHisto1D(3, 1, 1);
      _h_3pi   = bookHisto1D(4, 1, 1);
      _h_3pip0 = bookHisto1D(5, 1, 1);

    }

    // Helper function to look for specific decays
    bool isSpecificDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      
      // Specific bits for tau -> pi decays
      unsigned int n_pi0(0), n_piplus(0), n_piminus(0), n_nutau(0), n_nutaubar(0);
      for (int id : ids) {
        if      (id == PID::PI0)        n_pi0++;
        else if (id == PID::PIPLUS)     n_piplus++;
        else if (id == PID::PIMINUS)    n_piminus++;
        else if (id == PID::NU_TAU)     n_nutau++;
        else if (id == PID::NU_TAUBAR)  n_nutaubar++;
      }
 
      // Check for the explicit decay -- easy as we only deal with pi0 and pi+/-
      if ( count(children, hasPID(PID::PI0))       != n_pi0      ) return false;
      if ( count(children, hasPID(PID::PIPLUS))    != n_piplus   ) return false;
      if ( count(children, hasPID(PID::PIMINUS))   != n_piminus  ) return false;
      if ( count(children, hasPID(PID::NU_TAU))    != n_nutau    ) return false;
      if ( count(children, hasPID(PID::NU_TAUBAR)) != n_nutaubar ) return false;

      return true;
 
    }
    

    // Conveniece function to get m2 of sum of all hadronic tau decay product 4-vectors
    double hadronicm2(const Particle& mother) {
      FourMomentum p_tot(0,0,0,0);
      // Iterate over all children that are mesons
      for (const Particle & meson : filter_select(mother.children(), isMeson)) {
        // Add this mesons 4-momentum to total 4-momentum
        p_tot += meson.momentum();
      }
      return p_tot.mass2();
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Loop over taus
      for(const Particle& tau : apply<UnstableFinalState>(event, "UFS").particles(Cuts::abspid==PID::TAU)) {
        // tau -> pi pi0 nu_tau (both charges)
        if (isSpecificDecay(tau,  {PID::PIPLUS, PID::PI0, PID::NU_TAUBAR}) ||  
            isSpecificDecay(tau,  {PID::PIMINUS, PID::PI0, PID::NU_TAU}) ) {
          _h_pip0->fill(hadronicm2(tau), event.weight());
        }
        // tau -> pi pi0 pi0 nu_tau (both charges)
        else if (isSpecificDecay(tau,  {PID::PIPLUS, PID::PI0, PID::PI0, PID::NU_TAUBAR}) ||  
                 isSpecificDecay(tau,  {PID::PIMINUS, PID::PI0, PID::PI0, PID::NU_TAU}) ) {
          _h_pi2p0->fill(hadronicm2(tau), event.weight());
        }
        //    tau -> pi pi0 pi0 pi0         (3,1,1)
        else if (isSpecificDecay(tau,  {PID::PIPLUS,  PID::PI0, PID::PI0, PID::PI0, PID::NU_TAUBAR}) ||
                 isSpecificDecay(tau,  {PID::PIMINUS, PID::PI0, PID::PI0, PID::PI0, PID::NU_TAU}) ) {
          _h_pi3p0->fill(hadronicm2(tau), event.weight());
        }
        //    tau -> 3 charged pions        (4,1,1)
        else if (isSpecificDecay(tau,  {PID::PIPLUS,  PID::PIPLUS,  PID::PIMINUS, PID::NU_TAUBAR}) ||
                 isSpecificDecay(tau,  {PID::PIMINUS, PID::PIMINUS, PID::PIPLUS, PID::NU_TAU}) ) {
          _h_3pi->fill(hadronicm2(tau), event.weight());
        }
        //    tau -> 3 charged pions + pi0  (5,1,1)
        else if (isSpecificDecay(tau,  {PID::PIPLUS,  PID::PIPLUS,  PID::PIMINUS, PID::PI0, PID::NU_TAUBAR}) ||
                 isSpecificDecay(tau,  {PID::PIMINUS, PID::PIMINUS, PID::PIPLUS,  PID::PI0, PID::NU_TAU}) ) {
          _h_3pip0->fill(hadronicm2(tau), event.weight());
        }
        //
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_pip0);  // normalize to unity
      normalize(_h_pi2p0); // normalize to unity
      normalize(_h_pi3p0); // nor\pi^0malize to unity
      normalize(_h_3pi);   // normalize to unity
      normalize(_h_3pip0); // normalize to unity

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_pip0;
    Histo1DPtr _h_pi2p0;
    Histo1DPtr _h_pi3p0;
    Histo1DPtr _h_3pi;
    Histo1DPtr _h_3pip0;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALEPH_2014_I1267648);


}
