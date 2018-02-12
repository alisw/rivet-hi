// -*- C++ -*-
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  IdentifiedFinalState::IdentifiedFinalState(const FinalState& fsp, const vector<PdgId>& pids) {
    setName("IdentifiedFinalState");
    addProjection(fsp, "FS");
    acceptIds(pids);
  }

  IdentifiedFinalState::IdentifiedFinalState(const vector<PdgId>& pids, const FinalState& fsp) {
    setName("IdentifiedFinalState");
    addProjection(fsp, "FS");
    acceptIds(pids);
  }

  IdentifiedFinalState::IdentifiedFinalState(const FinalState& fsp, PdgId pid) {
    setName("IdentifiedFinalState");
    addProjection(fsp, "FS");
    acceptId(pid);
  }

  IdentifiedFinalState::IdentifiedFinalState(PdgId pid, const FinalState& fsp) {
    setName("IdentifiedFinalState");
    addProjection(fsp, "FS");
    acceptId(pid);
  }


  IdentifiedFinalState::IdentifiedFinalState(const Cut& c, const vector<PdgId>& pids) {
    setName("IdentifiedFinalState");
    addProjection(FinalState(c), "FS");
    acceptIds(pids);
  }

  IdentifiedFinalState::IdentifiedFinalState(const vector<PdgId>& pids, const Cut& c) {
    setName("IdentifiedFinalState");
    addProjection(FinalState(c), "FS");
    acceptIds(pids);
  }

  IdentifiedFinalState::IdentifiedFinalState(const Cut& c, PdgId pid) {
    setName("IdentifiedFinalState");
    addProjection(FinalState(c), "FS");
    acceptId(pid);
  }

  IdentifiedFinalState::IdentifiedFinalState(PdgId pid, const Cut& c) {
    setName("IdentifiedFinalState");
    addProjection(FinalState(c), "FS");
    acceptId(pid);
  }



  int IdentifiedFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;

    const IdentifiedFinalState& other = dynamic_cast<const IdentifiedFinalState&>(p);
    int pidssize = cmp(_pids.size(), other._pids.size());
    if (pidssize != EQUIVALENT) return pidssize;
    return cmp(_pids, other._pids);
  }


  void IdentifiedFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());
    _remainingParticles.clear();
    _remainingParticles.reserve(fs.particles().size());
    for (const Particle& p : fs.particles()) {
      if (acceptedIds().find(p.pid()) != acceptedIds().end()) {
        _theParticles.push_back(p);       // Identified
      } else {
        _remainingParticles.push_back(p); // Remaining
      }
    }
  }


}
