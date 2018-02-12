// -*- C++ -*-
#include "Rivet/Projections/NonPromptFinalState.hh"

namespace Rivet {


  NonPromptFinalState::NonPromptFinalState(const FinalState& fsp, bool accepttaudecays, bool acceptmudecays)
    : _acceptMuDecays(acceptmudecays), _acceptTauDecays(accepttaudecays)
  {
    setName("NonPromptFinalState");
    addProjection(fsp, "FS");
  }


  NonPromptFinalState::NonPromptFinalState(const Cut& c, bool accepttaudecays, bool acceptmudecays)
    : _acceptMuDecays(acceptmudecays), _acceptTauDecays(accepttaudecays)
  {
    setName("NonPromptFinalState");
    addProjection(FinalState(c), "FS");
  }


  int NonPromptFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;
    const NonPromptFinalState& other = dynamic_cast<const NonPromptFinalState&>(p);
    return cmp(_acceptMuDecays, other._acceptMuDecays) || cmp(_acceptTauDecays, other._acceptTauDecays);
  }


  void NonPromptFinalState::project(const Event& e) {
    _theParticles.clear();

    const Particles& particles = applyProjection<FinalState>(e, "FS").particles();
    for (const Particle& p : particles)
      if (!isPrompt(p, !_acceptTauDecays, !_acceptMuDecays)) _theParticles.push_back(p);
    MSG_DEBUG("Number of final state particles from hadron decays = " << _theParticles.size());

    if (getLog().isActive(Log::TRACE)) {
      for (const Particle& p : _theParticles)
        MSG_TRACE("Selected: " << p.pid() << ", charge = " << p.charge());
    }
  }


}
