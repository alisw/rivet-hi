#include "Rivet/Tools/ParticleUtils.hh"
#include "Rivet/Tools/Cuts.hh"

namespace Rivet {


  FirstParticleWith::FirstParticleWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  FirstParticleWithout::FirstParticleWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  LastParticleWith::LastParticleWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  LastParticleWithout::LastParticleWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }


  HasParticleAncestorWith::HasParticleAncestorWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleAncestorWithout::HasParticleAncestorWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleParentWith::HasParticleParentWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleParentWithout::HasParticleParentWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleChildWith::HasParticleChildWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleChildWithout::HasParticleChildWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleDescendantWith::HasParticleDescendantWith(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }

  HasParticleDescendantWithout::HasParticleDescendantWithout(const Cut& c)
    : fn([&](const Particle& p){ return c->accept(p); }) { }


  Particles& ifilter_select(Particles& particles, const Cut& c) {
    if (c == Cuts::OPEN) return particles;
    // return ifilter_select(particles, *c);
    return ifilter_select(particles, [&](const Particle& p){return c->accept(p);});
  }

  Particles& ifilter_discard(Particles& particles, const Cut& c) {
    if (c == Cuts::OPEN) { particles.clear(); return particles; }
    // return ifilter_discard(particles, *c);
    return ifilter_discard(particles, [&](const Particle& p){return c->accept(p);});
  }


}
