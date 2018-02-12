// -*- C++ -*-
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  ZFinder::ZFinder(const FinalState& inputfs,
		   const Cut & fsCut,
                   PdgId pid,
                   double minmass, double maxmass,
                   double dRmax,
                   ChargedLeptons chLeptons,
                   ClusterPhotons clusterPhotons,
                   PhotonTracking trackPhotons,
                   double masstarget)
  {
    setName("ZFinder");

    _minmass = minmass;
    _maxmass = maxmass;
    _masstarget = masstarget;
    _pid = abs(pid);
    _trackPhotons = trackPhotons;

    // Identify bare leptons for dressing
    // Bit of a code nightmare -- FS projection copy constructors don't work?
    /// @todo Fix FS copy constructors!!
    if (chLeptons == PROMPTCHLEPTONS) {
      PromptFinalState inputfs_prompt(inputfs);
      IdentifiedFinalState bareleptons = IdentifiedFinalState(inputfs_prompt);
      bareleptons.acceptIdPair(_pid);
      declare(bareleptons, "BareLeptons");
    } else {
      IdentifiedFinalState bareleptons = IdentifiedFinalState(inputfs);
      bareleptons.acceptIdPair(_pid);
      declare(bareleptons, "BareLeptons");
    }

    // Dress the bare leptons
    const bool doClustering = (clusterPhotons != NOCLUSTER);
    const bool useDecayPhotons = (clusterPhotons == CLUSTERALL);
    DressedLeptons leptons(inputfs, get<FinalState>("BareLeptons"), (doClustering ? dRmax : -1.0), fsCut, useDecayPhotons);
    addProjection(leptons, "DressedLeptons");

    // Identify the non-Z part of the event
    VetoedFinalState remainingFS;
    remainingFS.addVetoOnThisFinalState(*this);
    addProjection(remainingFS, "RFS");
  }



  /////////////////////////////////////////////////////



  Particles ZFinder::constituentLeptons() const {
    if (empty()) return Particles();
    return boson().constituents();
    // return boson().constituents(isChargedLepton);
  }


  Particles ZFinder::constituentLeptons(const ParticleSorter& cmp) const {
    return sortBy(constituentLeptons(), cmp);
  }


  const VetoedFinalState& ZFinder::remainingFinalState() const {
    return getProjection<VetoedFinalState>("RFS");
  }


  int ZFinder::compare(const Projection& p) const {
    PCmp LCcmp = mkNamedPCmp(p, "DressedLeptons");
    if (LCcmp != EQUIVALENT) return LCcmp;

    const ZFinder& other = dynamic_cast<const ZFinder&>(p);
    return (cmp(_minmass, other._minmass) ||
            cmp(_maxmass, other._maxmass) ||
            cmp(_pid, other._pid) ||
            cmp(_trackPhotons, other._trackPhotons));
  }


  void ZFinder::project(const Event& e) {
    clear();

    // Get leptons and find an acceptable invariant mass OSSF pair
    const DressedLeptons& leptons = applyProjection<DressedLeptons>(e, "DressedLeptons");
    InvMassFinalState imfs({_pid, -_pid}, _minmass, _maxmass, _masstarget);
    imfs.calc(leptons.particles());
    if (imfs.particlePairs().empty()) {
      MSG_TRACE("No acceptable inv-mass lepton/antilepton pairs found");
      return;
    }

    // Assemble a pseudo-Z particle
    const ParticlePair& Zconstituents = imfs.particlePairs().front();
    const Particle& p1(Zconstituents.first), p2(Zconstituents.second);
    const FourMomentum pZ = p1.momentum() + p2.momentum();
    assert(p1.charge3() + p2.charge3() == 0);
    Particle z(PID::Z0BOSON, pZ);
    MSG_DEBUG(z << " reconstructed from: " << p1 << " + " << p2);

    // Add (dressed) lepton constituents to the Z (skipping photons if requested)
    // Keep the DressedLeptons found by the ZFinder
    const Particle& l1 = p1.charge() > 0 ? p1 : p2;
    const Particle& l2 = p2.charge() < 0 ? p2 : p1;
    MSG_TRACE("l1 = " << l1.constituents());
    MSG_TRACE("l2 = " << l2.constituents());
    z.addConstituent(_trackPhotons == TRACK ? l1 : l1.constituents().front());
    z.addConstituent(_trackPhotons == TRACK ? l2 : l2.constituents().front());
    MSG_DEBUG("Number of stored raw Z constituents = " << z.rawConstituents().size() << "  " << z.rawConstituents());

    // Register the completed Z
    _theParticles.push_back(z);
  }


}
