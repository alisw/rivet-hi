// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief CDF <what is this analysis doing?>
  class CDF_1993_S2742446 : public Analysis {
  public:

    CDF_1993_S2742446()
      : Analysis("CDF_1993_S2742446")
    {    }


  public:

    void init() {

      // The photon selection has been corrected to pTmin=22 GeV (vs. 23 in the trigger)
      LeadingParticlesFinalState photonfs(FinalState(-0.9, 0.9, 22.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // FS excluding the leading photon
      VetoedFinalState vfs(FinalState(-4.2, 4.2));
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "VFS");

      // Jets
      declare(FastJets(vfs, FastJets::CDFJETCLU, 0.7), "Jets");

      _h_costheta = bookHisto1D(1, 1, 1);

    }


    void analyze(const Event& event) {

      const double weight = event.weight();

      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size()!=1 || photons[0].pT()>45.0*GeV) {
        vetoEvent;
      }
      FourMomentum leadingPhoton = photons[0].momentum();
      double eta_P = leadingPhoton.eta();
      double phi_P = leadingPhoton.phi();

      // photon isolation: less than 2 GeV EM E_T
      double Etsum=0.0;
      foreach (const Particle& p, apply<VetoedFinalState>(event, "VFS").particles()) {
        if (p.charge() != 0 && deltaR(eta_P, phi_P, p.eta(), p.phi()) < 0.7) Etsum += p.Et();
      }
      if (Etsum > 2*GeV) vetoEvent;

      FourMomentum jetsum;
      Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::pT > 10*GeV, cmpMomByPt);
  
      // Require at least one jet with pT>10 GeV
      if (jets.size()==0) vetoEvent;

      // Require the leading jet to be in the opposite (phi) hemisphere w.r.t. the photon
      if (jets[0].phi() - phi_P <= M_PI) vetoEvent;
      
      // sum all jets in the opposite hemisphere in phi from the photon
      foreach (const Jet& jet, jets) {
        if (fabs(jet.phi()-phi_P) > M_PI) jetsum+=jet.momentum();
      }

      // c.m. cuts, see Table 1
      double etaboost = (jetsum.eta()+eta_P)/2.;
      if (!inRange(etaboost, -1.2, 1.2)) vetoEvent;
      
      double etastar = (jetsum.eta()-eta_P)/2.;
      if (!inRange(etastar, -1.1, 1.1))  vetoEvent;
      
      double pstar = photons[0].pT()*cosh(etastar);
      if (!inRange(pstar, 27.8, 47.0)) vetoEvent;

      const double costheta = fabs(tanh((eta_P-jetsum.eta())/2.0));
      if (!inRange(costheta, 0., 0.8)) vetoEvent;
      
      // Fill histo
      _h_costheta->fill(costheta, weight);
    }


    void finalize() {
      /// @todo Take fixed norm direct from ref histo
      normalize(_h_costheta, 1.4271); // fixed norm ok
    }


  private:

    Histo1DPtr _h_costheta;

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CDF_1993_S2742446);

}
