// -*- C++ -*-
#include "Rivet/Projections/VisibleFinalState.hh"

namespace Rivet {


  int VisibleFinalState::compare(const Projection& p) const {
    return mkNamedPCmp(p, "FS");
  }


  // Since we remove invisibles from the FinalState in project(),
  // we need a filter where invisible --> true
  namespace {
    bool isInvisible(const Particle& p) { return !p.isVisible(); }
  }


  void VisibleFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    std::remove_copy_if(fs.particles().begin(), fs.particles().end(),
                        std::back_inserter(_theParticles), isInvisible);
    MSG_DEBUG("Number of visible final-state particles = " << _theParticles.size());
  }


}
