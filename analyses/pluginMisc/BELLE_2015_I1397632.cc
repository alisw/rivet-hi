// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2015_I1397632 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2015_I1397632);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      declare(UnstableFinalState(), "UFS");

      // Book histograms
      _h_B_Denu      = bookHisto1D(1, 1, 1);
      _h_B_Dmunu     = bookHisto1D(1, 1, 2);
      _h_B_Deplusnu  = bookHisto1D(1, 1, 3);
      _h_B_Dmuplusnu = bookHisto1D(1, 1, 4);
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }
    
    // Calculate the recoil w using mother and daugher meson
    double recoilW(const Particle& B, int mesonID) {
      // TODO why does that not work with const?
      Particle D = filter_select(B.children(), Cuts::pid==mesonID)[0];
      FourMomentum q = B.mom() - D.mom();
      return (B.mom()*B.mom() + D.mom()*D.mom() - q*q )/ (2. * sqrt(B.mom()*B.mom()) * sqrt(D.mom()*D.mom()) );
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Get B0 Mesons
      foreach(const Particle& p, apply<UnstableFinalState>(event, "UFS").particles(Cuts::pid==PID::B0)) {
        if (isSemileptonicDecay(p, {PID::DMINUS,PID::POSITRON,PID::NU_E})) _h_B_Denu->fill( recoilW(p, PID::DMINUS), event.weight());
        if (isSemileptonicDecay(p, {PID::DMINUS,PID::ANTIMUON,PID::NU_MU})) _h_B_Dmunu->fill(recoilW(p, PID::DMINUS), event.weight());
      }
      // Get B+ Mesons
      foreach(const Particle& p, apply<UnstableFinalState>(event, "UFS").particles(Cuts::pid==PID::BPLUS)) {
        if (isSemileptonicDecay(p, {PID::D0BAR,PID::POSITRON,PID::NU_E})) _h_B_Deplusnu->fill( recoilW(p, PID::D0BAR), event.weight());
        if (isSemileptonicDecay(p, {PID::D0BAR,PID::ANTIMUON,PID::NU_MU})) _h_B_Dmuplusnu->fill(recoilW(p, PID::D0BAR), event.weight());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_B_Denu);     
      normalize(_h_B_Dmunu);    
      normalize(_h_B_Deplusnu); 
      normalize(_h_B_Dmuplusnu);

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_B_Denu;
    Histo1DPtr _h_B_Dmunu;
    Histo1DPtr _h_B_Deplusnu;
    Histo1DPtr _h_B_Dmuplusnu;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2015_I1397632);


}
