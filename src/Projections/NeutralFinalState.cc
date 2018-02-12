// -*- C++ -*-
#include "Rivet/Projections/NeutralFinalState.hh"

namespace Rivet {


  int NeutralFinalState::compare(const Projection& p) const {
    const NeutralFinalState& other = dynamic_cast<const NeutralFinalState&>(p);
    return mkNamedPCmp(other, "FS") || cmp(_Etmin, other._Etmin);
  }


  void NeutralFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    for (const Particle& p : fs.particles()) {
      if (p.charge3() == 0 && p.Et() > _Etmin) {
        _theParticles.push_back(p);
        MSG_TRACE("Selected: ID = " << p.pid()
                  << ", Et = " << p.Et()
                  << ", eta = " << p.eta()
                  << ", charge = " << p.charge());
      }
    }
    MSG_DEBUG("Number of neutral final-state particles = " << _theParticles.size());
  }


}
