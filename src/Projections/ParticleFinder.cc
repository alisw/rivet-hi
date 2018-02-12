// -*- C++ -*-
#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {

  /// @todo HOW DO WE COMPARE CUTS OBJECTS?
  int ParticleFinder::compare(const Projection& p) const {
    const ParticleFinder& other = dynamic_cast<const ParticleFinder&>(p);

    //MSG_TRACE("FS::compare: " << 1 << " " << this << " " << &p);
    //    std::vector<std::pair<double, double> > eta1(_etaRanges);
    //std::vector<std::pair<double, double> > eta2(other._etaRanges);
    //std::sort(eta1.begin(), eta1.end());
    //std::sort(eta2.begin(), eta2.end());

    //MSG_TRACE("FS::compare: " << 2 << " " << this << " " << &p);
    //if (eta1 < eta2) return ORDERED;
    //else if (eta2 < eta1) return UNORDERED;

    //MSG_TRACE("FS::compare: " << 3 << " " << this << " " << &p);
    //return cmp(_ptmin, other._ptmin);
    return _cuts == other._cuts ? EQUIVALENT : UNDEFINED;
  }

  // void ParticleFinder::project(const Event& e) {
  //   _theParticles.clear();
  // }


}
