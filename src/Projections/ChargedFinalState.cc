// -*- C++ -*-
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  ChargedFinalState::ChargedFinalState(const FinalState& fsp) {
    setName("ChargedFinalState");
    addProjection(fsp, "FS");
  }

  ChargedFinalState::ChargedFinalState(const Cut& c) {
    setName("ChargedFinalState");
    addProjection(FinalState(c), "FS");
  }

  ChargedFinalState::ChargedFinalState(double mineta, double maxeta, double minpt) {
    setName("ChargedFinalState");
    addProjection(FinalState(mineta, maxeta, minpt), "FS");
  }

  int ChargedFinalState::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }
}

namespace {
  inline bool chargedParticleFilter(const Rivet::Particle& p) {
    return Rivet::PID::threeCharge(p.pdgId()) == 0;
  }
}

namespace Rivet {
  void ChargedFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(),
                        std::back_inserter(_theParticles), chargedParticleFilter);
    MSG_DEBUG("Number of charged final-state particles = " << _theParticles.size());
    if (getLog().isActive(Log::TRACE)) {
      for (vector<Particle>::iterator p = _theParticles.begin(); p != _theParticles.end(); ++p) {
        MSG_TRACE("Selected: " << p->pdgId() << ", charge = " << PID::threeCharge(p->pdgId())/3.0);
      }
    }
  }


}
