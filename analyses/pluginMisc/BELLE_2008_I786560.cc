// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief BELLE tau lepton to pi pi
  /// @author Peter Richardson
  class BELLE_2008_I786560 : public Analysis {
  public:

    BELLE_2008_I786560()
      : Analysis("BELLE_2008_I786560"),
        _weight_total(0),
        _weight_pipi(0)
    {   }


    void init() {
      declare(UnstableFinalState(), "UFS");
      _hist_pipi = bookHisto1D( 1, 1, 1);
    }


    void analyze(const Event& e) {
      // Find the taus
      Particles taus;
      const UnstableFinalState& ufs = apply<UnstableFinalState>(e, "UFS");
      foreach (const Particle& p, ufs.particles()) {
        if (p.abspid() != PID::TAU) continue;
        _weight_total += 1.;
        Particles pip, pim, pi0;
        unsigned int nstable = 0;
        // get the boost to the rest frame
        LorentzTransform cms_boost;
        if (p.p3().mod() > 1*MeV)
          cms_boost = LorentzTransform::mkFrameTransformFromBeta(p.momentum().betaVec());
        // find the decay products we want
        findDecayProducts(p.genParticle(), nstable, pip, pim, pi0);
        if (p.pid() < 0) {
          swap(pip, pim);
        }
        if (nstable != 3) continue;
        // pipi
        if (pim.size() == 1 && pi0.size() == 1) {
          _weight_pipi += 1.;
          _hist_pipi->fill((pi0[0].momentum()+pim[0].momentum()).mass2(),1.);
        }
      }
    }


    void finalize() {
      if (_weight_pipi > 0.) scale(_hist_pipi, 1./_weight_pipi);
    }


  private:

    //@{

    // Histograms
    Histo1DPtr _hist_pipi;

    // Weights counters
    double _weight_total, _weight_pipi;

    //@}

    void findDecayProducts(const GenParticle* p,
                           unsigned int & nstable,
                           Particles& pip, Particles& pim,
                           Particles& pi0) {
      const GenVertex* dv = p->end_vertex();
      /// @todo Use better looping
      for (GenVertex::particles_out_const_iterator pp = dv->particles_out_const_begin(); pp != dv->particles_out_const_end(); ++pp) {
        int id = (*pp)->pdg_id();
        if (id == PID::PI0 ) {
          pi0.push_back(Particle(**pp));
          ++nstable;
	}
        else if (id == PID::K0S)
          ++nstable;
        else if (id == PID::PIPLUS) {
          pip.push_back(Particle(**pp));
          ++nstable;
        }
        else if (id == PID::PIMINUS) {
          pim.push_back(Particle(**pp));
          ++nstable;
        }
        else if (id == PID::KPLUS) {
          ++nstable;
        }
        else if (id == PID::KMINUS) {
          ++nstable;
        }
        else if ((*pp)->end_vertex()) {
          findDecayProducts(*pp, nstable, pip, pim, pi0);
        }
        else
          ++nstable;
      }
    }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2008_I786560);

}
