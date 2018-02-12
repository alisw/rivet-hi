// -*- C++ -*-
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  PromptFinalState::PromptFinalState(bool accepttaudecays, bool acceptmudecays)
    : _acceptMuDecays(acceptmudecays), _acceptTauDecays(accepttaudecays)
  {
    setName("PromptFinalState");
    addProjection(FinalState(), "FS");
  }

  PromptFinalState::PromptFinalState(const Cut& c, bool accepttaudecays, bool acceptmudecays)
    : _acceptMuDecays(acceptmudecays), _acceptTauDecays(accepttaudecays)
  {
    setName("PromptFinalState");
    addProjection(FinalState(c), "FS");
  }

  PromptFinalState::PromptFinalState(const FinalState& fsp, bool accepttaudecays, bool acceptmudecays)
    : _acceptMuDecays(acceptmudecays), _acceptTauDecays(accepttaudecays)
  {
    setName("PromptFinalState");
    addProjection(fsp, "FS");
  }



  int PromptFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;
    const PromptFinalState& other = dynamic_cast<const PromptFinalState&>(p);
    return cmp(_acceptMuDecays, other._acceptMuDecays) || cmp(_acceptTauDecays, other._acceptTauDecays);
  }


  void PromptFinalState::project(const Event& e) {
    _theParticles.clear();

    const Particles& particles = applyProjection<FinalState>(e, "FS").particles();
    for (const Particle& p : particles)
      if (isPrompt(p, _acceptTauDecays, _acceptMuDecays)) _theParticles.push_back(p);
    MSG_DEBUG("Number of final state particles not from hadron decays = " << _theParticles.size());

    if (getLog().isActive(Log::TRACE)) {
      for (const Particle& p : _theParticles)
        MSG_TRACE("Selected: " << p.pid() << ", charge = " << p.charge());
    }
  }


}
