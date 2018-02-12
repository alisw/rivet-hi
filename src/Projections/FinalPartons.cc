// -*- C++ -*-

#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Tools/RivetHepMC.hh"

namespace Rivet {


  bool FinalPartons::accept(const Particle& p) const {

    // Reject if *not* a parton
    if (!isParton(p))
      return false;

    // Accept partons if they end on a standard hadronization vertex
    if (p.genParticle()->end_vertex() != nullptr && p.genParticle()->end_vertex()->id() == 5)
      return true;

    // Reject if p has a parton child.
    for (const Particle& c : p.children())
      if (isParton(c))
        return false;

    // Reject if from a hadron or tau decay.
    if (p.fromDecay())
      return false;

    return _cuts->accept(p);
  }


  void FinalPartons::project(const Event& e) {
    _theParticles.clear();
    for (const GenParticle* gp : Rivet::particles(e.genEvent())) {
      if (!gp) continue;
      const Particle p(gp);
      if (accept(p)) _theParticles.push_back(p);
    }
  }


}
