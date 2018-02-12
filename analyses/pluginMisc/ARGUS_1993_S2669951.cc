// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Production of the $\eta'(958)$ and $f_0(980)$ in $e^+e^-$ annihilation in the Upsilon region
  /// @author Peter Richardson
  class ARGUS_1993_S2669951 : public Analysis {
  public:

    ARGUS_1993_S2669951()
      : Analysis("ARGUS_1993_S2669951"),
        _count_etaPrime_highZ(2, 0.),
        _count_etaPrime_allZ(3, 0.),
        _count_f0(3, 0.),
        _weightSum_cont(0.),
        _weightSum_Ups1(0.),
        _weightSum_Ups2(0.)
    {   }


    void init() {
      declare(UnstableFinalState(), "UFS");

      _hist_cont_f0 = bookHisto1D(2, 1, 1);
      _hist_Ups1_f0 = bookHisto1D(3, 1, 1);
      _hist_Ups2_f0 = bookHisto1D(4, 1, 1);
    }


    void analyze(const Event& e) {

      // Find the Upsilons among the unstables
      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      Particles upsilons;

      // First in unstable final state
      foreach (const Particle& p, ufs.particles())
        if (p.pid() == 553 || p.pid() == 100553)
          upsilons.push_back(p);
      // Then in whole event if fails
      if (upsilons.empty()) {
        /// @todo Replace HepMC digging with Particle::descendents etc. calls
        foreach (const GenParticle* p, Rivet::particles(e.genEvent())) {
          if ( p->pdg_id() != 553 && p->pdg_id() != 100553 ) continue;
          // Discard it if its parent has the same PDG ID code (avoid duplicates)
          const GenVertex* pv = p->production_vertex();
          bool passed = true;
          if (pv) {
            foreach (const GenParticle* pp, particles_in(pv)) {
              if ( p->pdg_id() == pp->pdg_id() ) {
                passed = false;
                break;
              }
            }
          }
          if (passed) upsilons.push_back(Particle(*p));
        }
      }


      // Finding done, now fill counters
      const double weight = e.weight();
      if (upsilons.empty()) { // Continuum
        MSG_DEBUG("No Upsilons found => continuum event");

        _weightSum_cont += weight;
        unsigned int nEtaA(0), nEtaB(0), nf0(0);
        foreach (const Particle& p, ufs.particles()) {
          const int id = p.abspid();
          const double xp = 2.*p.E()/sqrtS();
          const double beta = p.p3().mod() / p.E();
          if (id == 9010221) {
            _hist_cont_f0->fill(xp, weight/beta);
            nf0 += 1;
          } else if (id == 331) {
            if (xp > 0.35) nEtaA += 1;
            nEtaB += 1;
          }
        }
        _count_f0[2]             += nf0*weight;
        _count_etaPrime_highZ[1] += nEtaA*weight;
        _count_etaPrime_allZ[2]  += nEtaB*weight;

      } else { // Upsilon(s) found
        MSG_DEBUG("Upsilons found => resonance event");

        foreach (const Particle& ups, upsilons) {
          const int parentId = ups.pid();
          ((parentId == 553) ? _weightSum_Ups1 : _weightSum_Ups2) += weight;
          Particles unstable;
          // Find the decay products we want
          findDecayProducts(ups.genParticle(), unstable);
          LorentzTransform cms_boost;
          if (ups.p3().mod() > 1*MeV)
            cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.momentum().betaVec());
          const double mass = ups.mass();
          unsigned int nEtaA(0), nEtaB(0), nf0(0);
          foreach(const Particle& p, unstable) {
            const int id = p.abspid();
            const FourMomentum p2 = cms_boost.transform(p.momentum());
            const double xp = 2.*p2.E()/mass;
            const double beta = p2.p3().mod()/p2.E();
            if (id == 9010221) { //< ?
              ((parentId == 553) ? _hist_Ups1_f0 : _hist_Ups2_f0)->fill(xp, weight/beta);
              nf0 += 1;
            } else if (id == 331) { //< ?
              if (xp > 0.35) nEtaA += 1;
              nEtaB += 1;
            }
          }
          if (parentId == 553) {
            _count_f0[0]             +=   nf0*weight;
            _count_etaPrime_highZ[0] += nEtaA*weight;
            _count_etaPrime_allZ[0]  += nEtaB*weight;
          } else {
            _count_f0[1] += nf0*weight;
            _count_etaPrime_allZ[1]  += nEtaB*weight;
          }
        }
      }
    }


    void finalize() {

      // High-Z eta' multiplicity
      Scatter2DPtr s111 = bookScatter2D(1, 1, 1, true);
      if (_weightSum_Ups1 > 0) // Point at 9.460
        s111->point(0).setY(_count_etaPrime_highZ[0] / _weightSum_Ups1, 0);
      if (_weightSum_cont > 0) // Point at 9.905
        s111->point(1).setY(_count_etaPrime_highZ[1] / _weightSum_cont, 0);

      // All-Z eta' multiplicity
      Scatter2DPtr s112 = bookScatter2D(1, 1, 2, true);
      if (_weightSum_Ups1 > 0) // Point at 9.460
        s112->point(0).setY(_count_etaPrime_allZ[0] / _weightSum_Ups1, 0);
      if (_weightSum_cont > 0) // Point at 9.905
        s112->point(1).setY(_count_etaPrime_allZ[2] / _weightSum_cont, 0);
      if (_weightSum_Ups2 > 0) // Point at 10.02
        s112->point(2).setY(_count_etaPrime_allZ[1] / _weightSum_Ups2, 0);

      // f0 multiplicity
      Scatter2DPtr s511 = bookScatter2D(5, 1, 1, true);
      if (_weightSum_Ups1 > 0) // Point at 9.46
        s511->point(0).setY(_count_f0[0] / _weightSum_Ups1, 0);
      if (_weightSum_Ups2 > 0) // Point at 10.02
        s511->point(1).setY(_count_f0[1] / _weightSum_Ups2, 0);
      if (_weightSum_cont > 0) // Point at 10.45
        s511->point(2).setY(_count_f0[2] / _weightSum_cont, 0);

      // Scale histos
      if (_weightSum_cont > 0.) scale(_hist_cont_f0, 1./_weightSum_cont);
      if (_weightSum_Ups1 > 0.) scale(_hist_Ups1_f0, 1./_weightSum_Ups1);
      if (_weightSum_Ups2 > 0.) scale(_hist_Ups2_f0, 1./_weightSum_Ups2);
    }


  private:

    /// @name Counters
    //@{
    vector<double> _count_etaPrime_highZ, _count_etaPrime_allZ, _count_f0;
    double _weightSum_cont,_weightSum_Ups1,_weightSum_Ups2;
    //@}


    /// Histos
    Histo1DPtr _hist_cont_f0, _hist_Ups1_f0, _hist_Ups2_f0;


    /// Recursively walk the HepMC tree to find decay products of @a p
    void findDecayProducts(const GenParticle* p, Particles& unstable) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        const int id = abs((*pp)->pdg_id());
        if (id == 331 || id == 9010221) unstable.push_back(Particle(*pp));
        else if ((*pp)->end_vertex()) findDecayProducts(*pp, unstable);
      }
    }


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ARGUS_1993_S2669951);

}
