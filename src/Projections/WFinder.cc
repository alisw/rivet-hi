// -*- C++ -*-
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  WFinder::WFinder(const FinalState& inputfs,
                   const Cut& leptoncuts,
                   PdgId pid,
                   double minmass, double maxmass,
                   double missingET,
                   double dRmax,
                   ChargedLeptons chLeptons,
                   ClusterPhotons clusterPhotons,
                   PhotonTracking trackPhotons,
                   MassWindow masstype,
                   double masstarget) {
    setName("WFinder");

    _etMissMin = missingET;
    _minmass = minmass;
    _maxmass = maxmass;
    _masstarget = masstarget;
    _pid = abs(pid);
    _trackPhotons = trackPhotons;
    _useTransverseMass = (masstype == TRANSMASS);

    // Check that the arguments are legal
    if (_pid != PID::ELECTRON && _pid != PID::MUON)
      throw Error("Invalid charged lepton PID given to WFinder");

    // Identify bare leptons for dressing
    // Bit of a code nightmare -- FS projection copy constructors don't work?
    /// @todo Fix FS copy constructors!!
    if (chLeptons == PROMPTCHLEPTONS) {
      PromptFinalState inputfs_prompt(inputfs);
      IdentifiedFinalState bareleptons(inputfs_prompt);
      bareleptons.acceptIdPair(_pid);
      declare(bareleptons, "BareLeptons");
    } else {
      IdentifiedFinalState bareleptons(inputfs);
      bareleptons.acceptIdPair(_pid);
      declare(bareleptons, "BareLeptons");
    }

    // Dress the bare leptons
    const bool doClustering = (clusterPhotons != NOCLUSTER);
    const bool useDecayPhotons = (clusterPhotons == CLUSTERALL);
    DressedLeptons leptons(inputfs, get<FinalState>("BareLeptons"), (doClustering ? dRmax : -1.), leptoncuts, useDecayPhotons);
    addProjection(leptons, "DressedLeptons");

    // Add MissingMomentum proj to calc MET
    MissingMomentum vismom(inputfs);
    addProjection(vismom, "MissingET");

    // Identify the non-Z part of the event
    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    addProjection(remainingFS, "RFS");
  }



  /////////////////////////////////////////////////////



  const Particles WFinder::constituentLeptons() const {
    if (empty()) return Particles();
    return boson().constituents(isChargedLepton);
  }


  const Particles WFinder::constituentNeutrinos() const {
    if (empty()) return Particles();
    return boson().constituents(isNeutrino);
  }


  const VetoedFinalState& WFinder::remainingFinalState() const {
    return getProjection<VetoedFinalState>("RFS");
  }


  const MissingMomentum& WFinder::missingMom() const {
    return getProjection<MissingMomentum>("MissingET");
  }


  int WFinder::compare(const Projection& p) const {
    PCmp dlcmp = mkNamedPCmp(p, "DressedLeptons");
    if (dlcmp != EQUIVALENT) return dlcmp;

    const WFinder& other = dynamic_cast<const WFinder&>(p);
    return (cmp(_minmass, other._minmass) ||
            cmp(_maxmass, other._maxmass) ||
            cmp(_useTransverseMass, other._useTransverseMass) ||
            cmp(_etMissMin, other._etMissMin) ||
            cmp(_pid, other._pid) ||
            cmp(_trackPhotons, other._trackPhotons));
  }


  void WFinder::project(const Event& e) {
    clear();

    // Check missing ET
    const MissingMomentum& missmom = applyProjection<MissingMomentum>(e, "MissingET");
    const double met = missmom.vectorEt().mod();
    MSG_TRACE("MET = " << met/GeV << " GeV vs. required > " << _etMissMin/GeV << " GeV");
    if (met < _etMissMin) {
      MSG_DEBUG("Not enough missing ET: " << met/GeV << " GeV vs. required > " << _etMissMin/GeV << " GeV");
      return;
    }

    // Get lepton
    const DressedLeptons& leptons = applyProjection<DressedLeptons>(e, "DressedLeptons");
    if ( leptons.dressedLeptons().empty() ) {
      MSG_DEBUG("No dressed leptons");
      return;
    }
    MSG_DEBUG("Found at least one dressed lepton: " << leptons.dressedLeptons().front().momentum() );

    // Get missing momentum 4-vector, assuming a massless invisible particle
    const FourMomentum pmiss = missmom.missingMomentum(0*GeV);
    MSG_DEBUG("Found missing 4-momentum: " << pmiss);

    // Compute an invariant mass final state for the W decay leptons (using pseudo-neutrinos from ETmiss)
    PdgId _nu_pid = _pid + 1;
    assert(_nu_pid == PID::NU_E || _nu_pid == PID::NU_MU);
    vector<pair<PdgId, PdgId> > l_nu_ids;
    l_nu_ids += make_pair(_pid, -_nu_pid);
    l_nu_ids += make_pair(-_pid, _nu_pid);
    InvMassFinalState imfs(l_nu_ids, _minmass, _maxmass, _masstarget);
    imfs.useTransverseMass(_useTransverseMass);
    Particles tmp = leptons.particles();
    tmp += { Particle( _nu_pid, pmiss), Particle(-_nu_pid, pmiss) }; // fake (anti)neutrinos from ETmiss vector
    imfs.calc(tmp);
    if (imfs.particlePairs().size() < 1) return;

    // Assemble a pseudo-W particle
    const ParticlePair Wconstituents = imfs.particlePairs().front();
    const Particle& p1(Wconstituents.first), p2(Wconstituents.second);
    const FourMomentum pW = p1.momentum() + p2.momentum();
    const int wcharge3 = p1.charge3() + p2.charge3();
    assert(abs(wcharge3) == 3);
    const int wcharge = wcharge3/3;
    const PdgId wpid = (wcharge == 1) ? PID::WPLUSBOSON : PID::WMINUSBOSON;
    Particle w(wpid, pW);
    MSG_DEBUG(w << " reconstructed from: " << p1 << " + " << p2);

    // Add (dressed) lepton constituents to the W (skipping photons if requested)
    /// @todo Do we need to add all used invisibles to _theParticles ?
    const Particle l = p1.isChargedLepton() ? p1 : p2;
    w.addConstituent(_trackPhotons == TRACK ? l : l.constituents().front());
    const Particle nu = p1.isNeutrino() ? p1 : p2;
    w.addConstituent(nu);

    // Register the completed W
    _theParticles.push_back(w);
  }


}
