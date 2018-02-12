// -*- C++ -*-
#ifndef RIVET_ZFinder_HH
#define RIVET_ZFinder_HH

#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Convenience finder of leptonically decaying Zs
  ///
  /// Chain together different projections as convenience for finding Z's
  /// from two leptons in the final state, including photon clustering.
  ///
  /// @todo Alias then rename as Dileptons
  class ZFinder : public ParticleFinder {
  public:

    enum ChargedLeptons { PROMPTCHLEPTONS=0, ALLCHLEPTONS };
    enum ClusterPhotons { NOCLUSTER=0, CLUSTERNODECAY=1, CLUSTERALL };
    enum PhotonTracking { NOTRACK=0, TRACK=1 };

    /// @name Constructors
    //@{

    /// Constructor taking cuts object
    /// @param inputfs Input final state
    /// @param cuts lepton cuts
    /// @param pid type of the leptons
    /// @param minmass,maxmass mass window
    /// @param dRmax maximum dR of photons around leptons to take into account
    ///  for Z reconstruction (only relevant if one of the following are true)
    /// @param clusterPhotons whether such photons are supposed to be
    ///  clustered to the lepton objects and thus Z mom
    /// @param trackPhotons whether such photons should be added to _theParticles
    ///  (cf. _trackPhotons)
    ZFinder(const FinalState& inputfs,
            const Cut& cuts,
            PdgId pid,
            double minmass, double maxmass,
            double dRmax=0.1,
            ChargedLeptons chLeptons=PROMPTCHLEPTONS,
            ClusterPhotons clusterPhotons=CLUSTERNODECAY,
            PhotonTracking trackPhotons=NOTRACK,
            double masstarget=91.2*GeV);

    /// Backward-compatible constructor with implicit chLeptons mode = PROMPTCHLEPTONS
    /// @deprecated Remove this and always use the constructor with chLeptons argument.
    ZFinder(const FinalState& inputfs,
            const Cut& cuts,
            PdgId pid,
            double minmass, double maxmass,
            double dRmax,
            ClusterPhotons clusterPhotons,
            PhotonTracking trackPhotons=NOTRACK,
            double masstarget=91.2*GeV)
      : ZFinder(inputfs, cuts, pid, minmass, maxmass,
                dRmax, PROMPTCHLEPTONS, clusterPhotons, trackPhotons, masstarget)
    {   }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(ZFinder);

    //@}


    /// Access to the found bosons
    ///
    /// @note Currently either 0 or 1 boson can be found.
    const Particles& bosons() const { return particles(); }
    /// Access to the found boson (assuming it exists).
    const Particle boson() const { return bosons().front(); }


    /// Access to the Z constituent clustered leptons
    ///
    /// For example, to make more fine-grained cuts on the clustered leptons.
    /// The positive charge constituent is first in the list (if not empty), and
    /// the negative one second.
    Particles constituentLeptons() const;
    Particles constituents() const { return constituentLeptons(); }


    /// Access to the Z constituent clustered leptons, sorted by a comparison functor
    ///
    /// Unlike the no-arg version, this returns by value (i.e. is less efficient)
    Particles constituentLeptons(const ParticleSorter& cmp) const;
    Particles constituents(const ParticleSorter& cmp) const { return constituentLeptons(cmp); }


    // /// Access to all DressedLeptons in the fiducial region.
    // ///
    // /// This includes those DressedLeptons that could not
    // /// be paired up with any other DressedLepton to form a Z candidate.
    // const vector<DressedLepton>& allLeptons() const { return _allLeptons; }

    /// Access to the particles other than the Z leptons and clustered photons
    ///
    /// Useful for e.g. input to a jet finder
    const VetoedFinalState& remainingFinalState() const;


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& e);

    /// Compare projections.
    int compare(const Projection& p) const;


  public:

    /// Clear the projection
    void clear() { _theParticles.clear(); }


  private:

    /// Mass cuts to apply to clustered leptons (cf. InvMassFinalState)
    double _minmass, _maxmass, _masstarget;

    /// Switch for tracking of photons (whether to include them in the Z particle)
    /// This is relevant when the clustered photons need to be excluded from e.g. a jet finder
    PhotonTracking _trackPhotons;

    /// Lepton flavour
    PdgId _pid;

  };


}

#endif
