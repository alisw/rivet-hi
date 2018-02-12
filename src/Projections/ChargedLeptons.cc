// -*- C++ -*-
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  int ChargedLeptons::compare(const Projection& other) const {
    return mkNamedPCmp(other, "ChFS");
  }


  void ChargedLeptons::project(const Event& evt) {
    // Reset result
    _theParticles.clear();

    // Loop over charged particles and fill vector with leptons
    const FinalState& fs = applyProjection<FinalState>(evt, "ChFS");
    for (const Particle& p : fs.particles()) {
      if (PID::isChargedLepton(p.pid()))
        _theParticles += Particle(p);
    }
    sortBy(_theParticles, cmpMomByPt);
  }


}
