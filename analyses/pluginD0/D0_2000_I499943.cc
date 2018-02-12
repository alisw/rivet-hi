// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  class D0_2000_I499943 : public Analysis {
  public:

    /// Constructor
    D0_2000_I499943()
      : Analysis("D0_2000_I499943")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      IdentifiedFinalState muons(Cuts::abseta < 0.8 && Cuts::pT > 4.0*GeV);
      muons.acceptIdPair(PID::MUON);
      declare(muons, "Muons");

      FastJets jetproj(fs, FastJets::D0ILCONE, 0.7);
      jetproj.useInvisibles();
      declare(jetproj, "Jets");

      // Book histograms
      _h_pt_leading_mu = bookHisto1D(1, 1, 1);
      _h_dphi_mumu = bookHisto1D(3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(12*GeV);
      if (jets.size() < 2) vetoEvent;

      const Particles& muons = apply<IdentifiedFinalState>(event, "Muons").particlesByPt();
      if (muons.size() < 2) vetoEvent;

      // Muon selection: require the muons to be *close* to jets, not the usual overlap vetoing!
      Particles cand_mu;
      foreach (const Particle& mu, muons) {
        // Ignore muons in "bad" region 80 < phi < 110 degrees
        /// @todo Is this really not corrected for?!
        if (inRange(mu.phi(), 1.4, 1.92)) continue;

        // A muon is a good candidate if within R = 0.8 of a jet
        foreach (const Jet& jet, jets) {
          if (deltaR(mu, jet) < 0.8) {
            cand_mu.push_back(mu);
            break;
          }
        }
      }

      // Must find at least two jet-matched muons in the event
      if (cand_mu.size() < 2) vetoEvent;

      /// @todo Is this cut needed? Does space angle mean dR or 3D opening angle in lab frame?
      // Remove muon pairs closer than 165 deg in space angle (cosmic veto)
      // double dR_mumu = deltaR(cand_mu[0].momentum(), cand_mu[1].momentum());
      // if (dR_mumu < 165*degree) vetoEvent;

      // Selecting muon pairs with 6 < mass < 35 GeV (we use the two with highest pT)
      double m_mumu = (cand_mu[0].momentum() + cand_mu[1].momentum()).mass();
      if (!inRange(m_mumu, 6*GeV, 35*GeV)) vetoEvent;

      // Get phi angle between muons in degrees
      double dphi_mumu = deltaPhi(cand_mu[0], cand_mu[1]) * 180/M_PI;

      // Fill histos
      _h_pt_leading_mu->fill(cand_mu[0].pt()/GeV, event.weight());
      _h_dphi_mumu->fill(dphi_mumu, event.weight());
    }


    // Normalise histograms to cross-section
    void finalize() {
      scale(_h_pt_leading_mu, crossSection()/sumOfWeights()/nanobarn);
      scale(_h_dphi_mumu, crossSection()/sumOfWeights()/nanobarn);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_pt_leading_mu, _h_dphi_mumu;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2000_I499943);

}
