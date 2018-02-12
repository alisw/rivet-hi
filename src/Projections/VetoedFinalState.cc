// -*- C++ -*-
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  int VetoedFinalState::compare(const Projection& p) const {
    const PCmp fscmp = mkNamedPCmp(p, "FS");
    if (fscmp != EQUIVALENT) return fscmp;
    /// @todo We can do better than this...
    if (_vetofsnames.size() != 0) return UNDEFINED;
    const VetoedFinalState& other = dynamic_cast<const VetoedFinalState&>(p);
    return \
      cmp(_vetoCodes, other._vetoCodes) ||
      cmp(_compositeVetoes, other._compositeVetoes) ||
      cmp(_parentVetoes, other._parentVetoes);
  }


  void VetoedFinalState::project(const Event& e) {
    const FinalState& fs = applyProjection<FinalState>(e, "FS");
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size());


    // Veto by PID code
    if (getLog().isActive(Log::TRACE)) {
      /// @todo Should be PdgId, but _vetoCodes is currently a long
      vector<long> codes;
      for (auto& code : _vetoCodes) codes += code.first;
      MSG_TRACE("Veto codes = " << codes << " (" << codes.size() << ")");
    }
    if (_vetoCodes.empty()) {
      _theParticles = fs.particles();
    } else {
      // Test every particle against the codes
      for (const Particle& p : fs.particles()) {
        VetoDetails::iterator iter = _vetoCodes.find(p.pid());
        if (iter == _vetoCodes.end()) {
          // MSG_TRACE("Storing with PDG code = " << p.pid() << ", pT = " << p.pT());
          _theParticles.push_back(p);
        } else {
          // This particle code is listed as a possible veto... check pT.
          // Make sure that the pT range is sensible:
          BinaryCut ptrange = iter->second;
          assert(ptrange.first <= ptrange.second);
          stringstream rangess;
          if (ptrange.first < numeric_limits<double>::max()) rangess << ptrange.second;
          rangess << " - ";
          if (ptrange.second < numeric_limits<double>::max()) rangess << ptrange.second;
          MSG_TRACE("ID = " << p.pid() << ", pT range = " << rangess.str());
          stringstream debugline;
          debugline << "with PDG code = " << p.pid() << " pT = " << p.pT();
          if (p.pT() < ptrange.first || p.pT() > ptrange.second) {
            MSG_TRACE("Storing " << debugline.str());
            _theParticles.push_back(p);
          } else {
            MSG_TRACE("Vetoing " << debugline.str());
          }
        }
      }
    }

    /// @todo What is this block? Mass vetoing?
    set<Particles::iterator> toErase;
    for (set<int>::iterator nIt = _nCompositeDecays.begin(); nIt != _nCompositeDecays.end() && !_theParticles.empty(); ++nIt) {
      map<set<Particles::iterator>, FourMomentum> oldMasses;
      map<set<Particles::iterator>, FourMomentum> newMasses;
      set<Particles::iterator> start;
      start.insert(_theParticles.begin());
      oldMasses.insert(pair<set<Particles::iterator>, FourMomentum>(start, _theParticles.begin()->momentum()));
      for (int nParts = 1; nParts != *nIt; ++nParts) {
        for (map<set<Particles::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
             mIt != oldMasses.end(); ++mIt) {
          Particles::iterator pStart = *(mIt->first.rbegin());
          for (Particles::iterator pIt = pStart + 1; pIt != _theParticles.end(); ++pIt) {
            FourMomentum cMom = mIt->second + pIt->momentum();
            set<Particles::iterator> pList(mIt->first);
            pList.insert(pIt);
            newMasses[pList] = cMom;
          }
        }
        oldMasses = newMasses;
        newMasses.clear();
      }
      for (map<set<Particles::iterator>, FourMomentum>::iterator mIt = oldMasses.begin();
           mIt != oldMasses.end(); ++mIt) {
        double mass2 = mIt->second.mass2();
        if (mass2 >= 0.0) {
          double mass = sqrt(mass2);
          for (CompositeVeto::iterator cIt = _compositeVetoes.lower_bound(*nIt);
               cIt != _compositeVetoes.upper_bound(*nIt); ++cIt) {
            BinaryCut massRange = cIt->second;
            if (mass < massRange.second && mass > massRange.first) {
              for (set<Particles::iterator>::iterator lIt = mIt->first.begin();
                   lIt != mIt->first.end(); ++lIt) {
                toErase.insert(*lIt);
              }
            }
          }
        }
      }
    }
    for (set<Particles::iterator>::reverse_iterator p = toErase.rbegin(); p != toErase.rend(); ++p) {
      _theParticles.erase(*p);
    }


    // Remove particles whose parents match entries in the parent veto PDG ID codes list
    /// @todo There must be a nice way to do this -- an STL algorithm (or we provide a nicer wrapper)
    for (PdgId vetoid : _parentVetoes) {
      for (Particles::iterator ip = _theParticles.begin(); ip != _theParticles.end(); ++ip) {
        const GenVertex* startVtx = ip->genParticle()->production_vertex();
        if (startVtx == NULL) continue;
        // Loop over parents and test their IDs
        /// @todo Could use any() here?
        for (const GenParticle* parent : Rivet::particles(startVtx, HepMC::ancestors)) {
          if (vetoid == parent->pdg_id()) {
            ip = _theParticles.erase(ip); --ip; //< Erase this _theParticles entry
            break;
          }
        }
      }
    }

    // Finally veto on the registered FSes
    for (const string& ifs : _vetofsnames) {
      const ParticleFinder& vfs = applyProjection<ParticleFinder>(e, ifs);
      const Particles& pvetos = vfs.rawParticles();
      ifilter_discard(_theParticles, [&](const Particle& pcheck) {
          if (pcheck.genParticle() == nullptr) return false;
          for (const Particle& pveto : pvetos) {
            if (pveto.genParticle() == nullptr) continue;
            if (pveto.genParticle() == pcheck.genParticle()) { MSG_TRACE("Vetoing: " << pcheck); return true; }
          }
          return false;
        });
    }

    MSG_DEBUG("FS vetoing from #particles = " << fs.size() << " -> " << _theParticles.size());
  }


}
