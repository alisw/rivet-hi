// -*- C++ -*-
#ifndef RIVET_NonPromptFinalState_HH
#define RIVET_NonPromptFinalState_HH

#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Find final state particles NOT directly connected to the hard process.
  ///
  /// See PromptFinalState for details.
  class NonPromptFinalState : public FinalState {
  public:

    /// @name Constructors
    //@{

    // Constructor from a final state.
    NonPromptFinalState(const FinalState& fsp, bool accepttaudecays=false, bool acceptmudecays=false);

    /// Constructor from a Cut (and implicit general FS).
    NonPromptFinalState(const Cut& c, bool accepttaudecays=false, bool acceptmudecays=false);

    // /// Constructor from a Cut and optional FinalState.
    // NonPromptFinalState(const Cut& c, const FinalState& fsp=FinalState(),
    //                     bool accepttaudecays=false, bool acceptmudecays=false);

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(NonPromptFinalState);

    //@}

    /// Treat particles from decays of prompt muons as non-prompt?
    void acceptMuonDecays(bool acc=true) { _acceptMuDecays = acc; }
    /// Treat particles from decays of prompt taus as non-prompt?
    void acceptTauDecays(bool acc=true) { _acceptTauDecays = acc; }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;

  private:

    bool _acceptMuDecays, _acceptTauDecays;

  };

}


#endif
