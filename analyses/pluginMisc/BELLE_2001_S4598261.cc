// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief BELLE pi0 spectrum at Upsilon(4S)
  /// @author Peter Richardson
  class BELLE_2001_S4598261 : public Analysis {
  public:

    BELLE_2001_S4598261()
      : Analysis("BELLE_2001_S4598261"), _weightSum(0.)
    { }


    void init() {
      declare(UnstableFinalState(), "UFS");
      _histdSigDp = bookHisto1D(1, 1, 1); // spectrum
      _histMult   = bookHisto1D(2, 1, 1); // multiplicity
    }


    void analyze(const Event& e) {
      const double weight = e.weight();

      // Find the upsilons
      Particles upsilons;
      // First in unstable final state
      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      foreach (const Particle& p, ufs.particles())
        if (p.pid()==300553) upsilons.push_back(p);
      // Then in whole event if fails
      if (upsilons.empty()) {
        foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
          if (p->pdg_id() != 300553) continue;
          const GenVertex* pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            /// @todo Use better looping
            for (GenVertex::particles_in_const_iterator pp = pv->particles_in_const_begin() ; pp != pv->particles_in_const_end() ; ++pp) {
              if ( p->pdg_id() == (*pp)->pdg_id() ) {
                passed = false;
                break;
              }
            }
          }
          if (passed) upsilons.push_back(Particle(p));
        }
      }

      // Find upsilons
      foreach (const Particle& p, upsilons) {
        _weightSum += weight;
        // Find the neutral pions from the decay
        vector<GenParticle *> pions;
        findDecayProducts(p.genParticle(), pions);
        const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        for (size_t ix=0; ix<pions.size(); ++ix) {
          const double pcm = cms_boost.transform(FourMomentum(pions[ix]->momentum())).p();
          _histdSigDp->fill(pcm,weight);
        }
        _histMult->fill(0., pions.size()*weight);
      }
    }


    void finalize() {
      scale(_histdSigDp, 1./_weightSum);
      scale(_histMult  , 1./_weightSum);
    }


  private:

    //@{
    // count of weights
    double _weightSum;
    /// Histograms
    Histo1DPtr _histdSigDp;
    Histo1DPtr _histMult;
    //@}


    void findDecayProducts(const GenParticle* p, vector<GenParticle*>& pions) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        const int id = (*pp)->pdg_id();
        if (id == 111) {
          pions.push_back(*pp);
        } else if ((*pp)->end_vertex())
          findDecayProducts(*pp, pions);
      }
    }


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2001_S4598261);

}
