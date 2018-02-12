// -*- C++ -*-
#ifndef RIVET_SmearedJets_HH
#define RIVET_SmearedJets_HH

#include "Rivet/Jet.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Projections/JetAlg.hh"
#include "Rivet/Tools/SmearingFunctions.hh"
#include <functional>

namespace Rivet {


  // // Recursive variadic template arg decoding
  // namespace {
  //   template<typename T>
  //   vector<JetEffSmearFn>& toEffSmearFns(vector<JetEffSmearFn>& v, const T& t) {
  //     v.push_back(JetEffSmearFn(t));
  //     return v;
  //   }
  //   template<typename T, typename... ARGS>
  //   vector<JetEffSmearFn>& toEffSmearFns(vector<JetEffSmearFn>& v, const T& first, ARGS... args) {
  //     v.push_back(JetEffSmearFn(first));
  //     toEffSmearFns(v, args...);
  //     return v;
  //   }
  // }



  /// Wrapper projection for smearing {@link Jet}s with detector resolutions and efficiencies
  class SmearedJets : public JetAlg {
  public:

    /// @name Constructors etc.
    //@{

    /// @brief Constructor with a reco efficiency and optional tagging efficiencies
    SmearedJets(const JetAlg& ja,
                const JetSmearFn& smearFn,
                const JetEffFn& bTagEffFn=JET_BTAG_PERFECT,
                const JetEffFn& cTagEffFn=JET_CTAG_PERFECT)
      : SmearedJets(ja, vector<JetEffSmearFn>{smearFn}, bTagEffFn, cTagEffFn)
    {    }


    /// @brief Constructor with tagging efficiencies, plus an ordered init-list of efficiency and smearing functions
    SmearedJets(const JetAlg& ja,
                const JetEffFn& bTagEffFn=JET_BTAG_PERFECT,
                const JetEffFn& cTagEffFn=JET_CTAG_PERFECT,
                const initializer_list<JetEffSmearFn>& effSmearFns={})
      : SmearedJets(ja, vector<JetEffSmearFn>{effSmearFns}, bTagEffFn, cTagEffFn)
    {    }

    /// @brief Constructor with tagging efficiencies, plus an ordered vector of efficiency and smearing functions
    SmearedJets(const JetAlg& ja,
                const JetEffFn& bTagEffFn=JET_BTAG_PERFECT,
                const JetEffFn& cTagEffFn=JET_CTAG_PERFECT,
                const vector<JetEffSmearFn>& effSmearFns={})
      : SmearedJets(ja, effSmearFns, bTagEffFn, cTagEffFn)
    {    }


    /// @brief Constructor with an ordered init-list of efficiency and smearing functions, plus optional tagging efficiencies
    SmearedJets(const JetAlg& ja,
                const initializer_list<JetEffSmearFn>& effSmearFns,
                const JetEffFn& bTagEffFn=JET_BTAG_PERFECT,
                const JetEffFn& cTagEffFn=JET_CTAG_PERFECT)
      : SmearedJets(ja, vector<JetEffSmearFn>{effSmearFns}, bTagEffFn, cTagEffFn)
    {    }

    /// @brief Constructor with an ordered vector of efficiency and smearing functions, plus optional tagging efficiencies
    SmearedJets(const JetAlg& ja,
                const vector<JetEffSmearFn>& effSmearFns,
                const JetEffFn& bTagEffFn=JET_BTAG_PERFECT,
                const JetEffFn& cTagEffFn=JET_CTAG_PERFECT)
      : _detFns(effSmearFns), _bTagEffFn(bTagEffFn), _cTagEffFn(cTagEffFn)
    {
      setName("SmearedJets");
      addProjection(ja, "TruthJets");
    }


    /// @brief Constructor with trailing efficiency arg
    /// @deprecated Use the version with pair-smearing list as 2nd argument
    SmearedJets(const JetAlg& ja,
                const JetSmearFn& smearFn,
                const JetEffFn& bTagEffFn,
                const JetEffFn& cTagEffFn,
                const JetEffFn& jetEffFn)
      : SmearedJets(ja, {jetEffFn,smearFn}, bTagEffFn, cTagEffFn)
    {    }


    /// @todo How to include tagging effs?
    /// @todo Variadic eff/smear fn list?
    /// @todo Add a trailing Cut arg cf. SmearedParticles? -- wrap into an eff function


    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SmearedJets);

    //@}


    /// Compare to another SmearedJets
    int compare(const Projection& p) const {
      // Compare truth jets definitions
      const int teq = mkPCmp(p, "TruthJets");
      if (teq != EQUIVALENT) return UNEQUAL;

      // Compare lists of detector functions
      const SmearedJets& other = dynamic_cast<const SmearedJets&>(p);
      if (_detFns.size() != other._detFns.size()) return UNEQUAL;
      for (size_t i = 0; i < _detFns.size(); ++i) {
        const int feq = _detFns[i].cmp(other._detFns[i]);
        if (feq != EQUIVALENT) return UNEQUAL;
      }

      // If we got this far, we're equal
      return EQUIVALENT;
    }


    /// Perform the jet finding & smearing calculation
    void project(const Event& e) {
      // Copying and filtering
      const Jets& truthjets = apply<JetAlg>(e, "TruthJets").jetsByPt();
      _recojets.clear(); _recojets.reserve(truthjets.size());
      // Apply jet smearing and efficiency transforms
      for (const Jet& j : truthjets) {
        Jet jdet = j;
        bool keep = true;
        MSG_DEBUG("Truth jet: " << "mom=" << jdet.mom()/GeV << " GeV, pT=" << jdet.pT()/GeV << ", eta=" << jdet.eta());
        for (const JetEffSmearFn& fn : _detFns) {
          double jeff = -1;
          tie(jdet, jeff) = fn(jdet); // smear & eff
          // Re-add constituents & tags if (we assume accidentally) they were lost by the smearing function
          if (jdet.particles().empty() && !j.particles().empty()) jdet.particles() = j.particles();
          if (jdet.tags().empty() && !j.tags().empty()) jdet.tags() = j.tags();
          MSG_DEBUG("         ->" << "mom=" << jdet.mom()/GeV << " GeV, pT=" << jdet.pT()/GeV << ", eta=" << jdet.eta());
          // MSG_DEBUG("New det jet: "
          //           << "mom=" << jdet.mom()/GeV << " GeV, pT=" << jdet.pT()/GeV << ", eta=" << jdet.eta()
          //           << ", b-tag=" << boolalpha << jdet.bTagged()
          //           << ", c-tag=" << boolalpha << jdet.cTagged()
          //           << " : eff=" << 100*jeff << "%");
          if (jeff <= 0) { keep = false; break; } //< no need to roll expensive dice (and we deal with -ve probabilities, just in case)
          if (jeff < 1 && rand01() > jeff)  { keep = false; break; } //< roll dice (and deal with >1 probabilities, just in case)
        }
        if (keep) _recojets.push_back(jdet);
      }
      // Apply tagging efficiencies, using smeared kinematics as input to the tag eff functions
      for (Jet& j : _recojets) {
        // Decide whether or not there should be a b-tag on this jet
        const double beff = _bTagEffFn ? _bTagEffFn(j) : j.bTagged();
        const bool btag = beff == 1 || (beff != 0 && rand01() < beff);
        // Remove b-tags if needed, and add a dummy one if needed
        if (!btag && j.bTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasBottom), j.tags().end());
        if (btag && !j.bTagged()) j.tags().push_back(Particle(PID::BQUARK, j.mom())); ///< @todo Or could use the/an actual clustered b-quark momentum?
        // Decide whether or not there should be a c-tag on this jet
        const double ceff = _cTagEffFn ? _cTagEffFn(j) : j.cTagged();
        const bool ctag = ceff == 1 || (ceff != 0 && rand01() < beff);
        // Remove c-tags if needed, and add a dummy one if needed
        if (!ctag && j.cTagged()) j.tags().erase(std::remove_if(j.tags().begin(), j.tags().end(), hasCharm), j.tags().end());
        if (ctag && !j.cTagged()) j.tags().push_back(Particle(PID::CQUARK, j.mom())); ///< @todo As above... ?
      }
    }


    /// Return the full jet list for the JetAlg methods to use
    Jets _jets() const { return _recojets; }


    /// Reset the projection. Smearing functions will be unchanged.
    void reset() { _recojets.clear(); }


  private:

    /// Smeared jets
    Jets _recojets;

    /// Stored efficiency & smearing functions
    vector<JetEffSmearFn> _detFns;

    /// Stored efficiency functions
    JetEffFn _bTagEffFn, _cTagEffFn;

  };


}

#endif
