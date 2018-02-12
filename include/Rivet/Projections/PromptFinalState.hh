// -*- C++ -*-
#ifndef RIVET_PromptFinalState_HH
#define RIVET_PromptFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Find final state particles directly connected to the hard process.
  ///
  /// The definition of "prompt" used in Rivet is that from high-scale physics, i.e.
  /// particles directly connected to the hard process in an interaction, regardless
  /// of realistic reconstructibility of displaced vertices, etc. By construction
  /// hadrons cannot be considered prompt as they will be colour connected to other
  /// parts of the event through non-perturbative effects: this projection can
  /// return electrons, muons, photons, and exotic particles which do not have a
  /// hadron in their post-hadronization ancestor chain. Flags exist to choose
  /// whether intermediate tau or muon decays invalidate a particle's promptness.
  ///
  /// @todo Decide how to treat brem photons off prompt leptons -- are they also prompt? "Decay" does not change the lepton PID...
  class PromptFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    /// Constructor without cuts
    PromptFinalState(bool accepttaudecays=false, bool acceptmudecays=false);

    /// Constructor from a Cut
    PromptFinalState(const Cut& c, bool accepttaudecays=false, bool acceptmudecays=false);

    // Constructor from a FinalState
    PromptFinalState(const FinalState& fsp, bool accepttaudecays=false, bool acceptmudecays=false);

    // /// Constructor from a Cut and optional FinalState.
    // PromptFinalState(const Cut& c, const FinalState& fsp=FinalState(), bool accepttaudecays, bool acceptmudecays);

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PromptFinalState);

    //@}

    /// Accept leptons from decays of prompt muons as themselves being prompt?
    void acceptMuonDecays(bool acc=true) { _acceptMuDecays = acc; }
    /// Accept leptons from decays of prompt taus as themselves being prompt?
    void acceptTauDecays(bool acc=true) { _acceptTauDecays = acc; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  private:

    bool _acceptMuDecays, _acceptTauDecays;

  };


  /// Alias with a more correct name
  using DirectFinalState = PromptFinalState;

}


#endif
