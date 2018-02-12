// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  class ATLAS_2014_I1282441 : public Analysis {
  public:

    ATLAS_2014_I1282441()
      : Analysis("ATLAS_2014_I1282441")
    { }


    void init() {

      // Use a large eta range such that we can discriminate on y
      /// @todo Convert to use a y-cut directly
      UnstableFinalState ufs(Cuts::abseta < 10 && Cuts::pT > 500*MeV);
      IdentifiedFinalState phis(ufs);
      phis.acceptIdPair(PID::PHI);
      declare(phis, "Phis");

      IdentifiedFinalState kpms(Cuts::abseta < 2.0 && Cuts::pT > 230*MeV);
      kpms.acceptIdPair(PID::KPLUS);
      declare(kpms, "Kpms");

      _h_phi_rapidity = bookHisto1D(1,1,1);
      _h_phi_pT       = bookHisto1D(2,1,1);
    }


    void analyze(const Event& event) {
      const Particles& ks_all = apply<IdentifiedFinalState>(event, "Kpms").particles();
      Particles kp, km;
      foreach (const Particle& p, ks_all) {
        if (!p.hasAncestor(PID::PHI)) { MSG_DEBUG("-- K not from phi."); continue; }
        if (p.p3().mod() > 800*MeV) { MSG_DEBUG("-- p K too high."); continue; }
        (p.charge() > 0 ? kp : km).push_back(p);
      }

      const Particles& phis_all = apply<FinalState>(event, "Phis").particles();
      Particles phis;
      /// @todo Use particles(Cuts&) instead
      foreach (const Particle& p, phis_all) {
        if ( p.absrap() > 0.8 ) { MSG_DEBUG("-- phi Y too high."); continue; }
        if ( p.pT() > 1.2*GeV ) { MSG_DEBUG("-- phi pT too high."); continue; }
        phis.push_back(p);
      }

      // Find Phi -> KK decays through matching of the kinematics
      if (!kp.empty() && !km.empty() && !phis.empty()) {
        const double w = event.weight();
        MSG_DEBUG("Numbers of particles:  #phi=" << phis.size() << ", #K+=" << kp.size() << ", #K-=" << km.size());
        for (size_t ip = 0; ip < phis.size(); ++ip) {
          const Particle& phi = phis[ip];
          for (size_t ikm = 0; ikm < km.size(); ++ikm) {
            for (size_t ikp = 0; ikp < kp.size(); ++ikp) {
              const FourMomentum mom = kp[ikp].mom() + km[ikm].mom();
              if ( fuzzyEquals(mom.mass(), phi.mass(), 1e-5) ) {
                MSG_DEBUG("Accepted combinatoric: phi#:" << ip << " K+#:" << ikp << " K-#:" << ikm);
                _h_phi_rapidity->fill(phi.absrap(), w);
                _h_phi_pT->fill(phi.pT()/MeV, w);
              } else {
                MSG_DEBUG("Rejected combinatoric: phi#:" << ip << " K+#:" << ikp << " K-#:" << ikm << " Mass difference is " << mom.mass()-phi.mass());
              }
            }
          }
        }
      }

    }


    void finalize() {
      scale(_h_phi_rapidity, crossSection()/microbarn/sumOfWeights());
      scale(_h_phi_pT, crossSection()/microbarn/sumOfWeights());
    }


  private:

    Histo1DPtr _h_phi_rapidity, _h_phi_pT;

  };



  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1282441);

}
