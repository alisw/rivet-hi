// -*- C++ -*-
#ifndef RIVET_PartonicTops_HH
#define RIVET_PartonicTops_HH

#include "Rivet/Projections/ParticleFinder.hh"

namespace Rivet {


  /// @brief Convenience finder of partonic top quarks
  ///
  /// @warning Requires there to be tops in the event record. A fiducial pseudo-top
  /// analysis approach is strongly recommended instead of this.
  class PartonicTops : public ParticleFinder {
  public:


    /// @brief Enum for categorising top quark decay modes
    ///
    /// More specifically, the decay mode of the W from the top. We presume top decay to a W and b quark.
    enum DecayMode { ELECTRON, MUON, TAU, E_MU, E_MU_TAU, HADRONIC, ALL };


    /// @name Constructors
    //@{

    /// Constructor optionally taking cuts object
    PartonicTops(const Cut& c=Cuts::OPEN)
      : ParticleFinder(c), _decaymode(ALL), _emu_from_prompt_tau(true), _include_hadronic_taus(false)
    {  }

    /// Constructor taking decay mode details (and an optional cuts object)
    PartonicTops(DecayMode decaymode, bool emu_from_prompt_tau=true, bool include_hadronic_taus=false, const Cut& c=Cuts::OPEN)
      : ParticleFinder(c), _decaymode(decaymode),
        _emu_from_prompt_tau(emu_from_prompt_tau), _include_hadronic_taus(include_hadronic_taus)
    {  }

    /// Constructor taking decay mode details (and an optional cuts object)
    PartonicTops(DecayMode decaymode, const Cut& c, bool emu_from_prompt_tau=true, bool include_hadronic_taus=false)
      : ParticleFinder(c), _decaymode(decaymode),
        _emu_from_prompt_tau(emu_from_prompt_tau), _include_hadronic_taus(include_hadronic_taus)
    {  }


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(PartonicTops);

    //@}


    /// Access to the found partonic tops
    const Particles& tops() const { return _theParticles; }


    /// Clear the projection
    void clear() {
      _theParticles.clear();
    }


  protected:

    /// Apply the projection on the supplied event.
    void project(const Event& event) {
      // Find partonic tops
      _theParticles = filter_select(event.allParticles(_cuts), lastParticleWith(isTop));
      // Filtering by decay mode
      if (_decaymode != ALL) {
        const auto fn = [&](const Particle& t) {
          const Particles descendants = t.allDescendants();
          const bool prompt_e = any(descendants, [&](const Particle& p){ return p.abspid() == PID::ELECTRON && p.isPrompt(_emu_from_prompt_tau) && !p.hasAncestor(PID::PHOTON); });
          const bool prompt_mu = any(descendants, [&](const Particle& p){ return p.abspid() == PID::MUON && p.isPrompt(_emu_from_prompt_tau) && !p.hasAncestor(PID::PHOTON); });
          if (prompt_e && (_decaymode == ELECTRON || _decaymode == E_MU || _decaymode == E_MU_TAU)) return true;
          if (prompt_mu && (_decaymode == MUON || _decaymode == E_MU || _decaymode == E_MU_TAU)) return true;
          const bool prompt_tau = any(descendants, [&](const Particle& p){ return p.abspid() == PID::TAU && p.isPrompt(); });
          const bool prompt_hadronic_tau = any(descendants, [&](const Particle& p){ return p.abspid() == PID::TAU && p.isPrompt() && none(p.children(), isChargedLepton); });
          if (prompt_tau && (_decaymode == TAU || _decaymode == E_MU_TAU)) return (_include_hadronic_taus || !prompt_hadronic_tau);
          if (_decaymode == HADRONIC && (!prompt_e && !prompt_mu && (!prompt_tau || (_include_hadronic_taus && prompt_hadronic_tau)))) return true; //< logical hairiness...
          return false;
        };
        ifilter_select(_theParticles, fn);
      }
    }


    /// Compare projections.
    int compare(const Projection& p) const {
      const PartonicTops& other = dynamic_cast<const PartonicTops&>(p);
      return cmp(_cuts, other._cuts) || cmp(_decaymode, other._decaymode) ||
        cmp(_emu_from_prompt_tau, other._emu_from_prompt_tau) ||
        cmp(_include_hadronic_taus, other._include_hadronic_taus);
    }


  private:

    DecayMode _decaymode;

    bool _emu_from_prompt_tau, _include_hadronic_taus;

  };


}


#endif
