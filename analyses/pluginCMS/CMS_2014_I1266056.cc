// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Measurement of gamma + jets + X triple differential cross-sections
  ///
  /// @author David Grellscheid
  class CMS_2014_I1266056 : public Analysis {
  public:

    // Constructor
    CMS_2014_I1266056()
      : Analysis("CMS_2014_I1266056")
    {    }


    // Book histograms and initialise projections before the run
    void init() {
      // Final state
      FinalState fs(-3, 3);
      addProjection(fs, "FS");

      // Leading photon
      LeadingParticlesFinalState photonfs(FinalState(-2.5, 2.5, 40.0*GeV));
      photonfs.addParticleId(PID::PHOTON);
      addProjection(photonfs, "LeadingPhoton");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      addProjection(vfs, "JetFS");

      // Jets
      FastJets jetpro(vfs, FastJets::ANTIKT, 0.5);
      //jetpro.useInvisibles();
      addProjection(jetpro, "Jets");

      _h_phverycentral_jetcentral = bookHisto1D(1, 1, 1);
      _h_phcentral_jetcentral     = bookHisto1D(2, 1, 1);
      _h_phforward_jetcentral     = bookHisto1D(3, 1, 1);
      _h_phveryforward_jetcentral = bookHisto1D(4, 1, 1);

      _h_phverycentral_jetforward = bookHisto1D(1, 1, 2);
      _h_phcentral_jetforward     = bookHisto1D(2, 1, 2);
      _h_phforward_jetforward     = bookHisto1D(3, 1, 2);
      _h_phveryforward_jetforward = bookHisto1D(4, 1, 2);

    }

    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the photon
      const FinalState& photonfs = applyProjection<FinalState>(event, "LeadingPhoton");
      if (photonfs.particles().empty()) vetoEvent;
      const FourMomentum photon = photonfs.particles().front().momentum();

      // Get the jet
      Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      if (jets.empty()) vetoEvent;
      FourMomentum leadingJet;
      for ( const Jet & j : jets ) {
        leadingJet = j.momentum();
        // keep the first separated jet
        if (deltaR(photon, leadingJet) > 0.5)
          break;
      }
      if (deltaR(photon, leadingJet) < 0.5)
          vetoEvent;

      // Veto if leading jet is outside plotted rapidity regions
      if (leadingJet.abseta() > 2.5) vetoEvent;

      // TODO: photon isolation 'IsoGamma' needed?

      // Fill histos
      const double abs_jet_eta = leadingJet.abseta();
      const double photon_pt = photon.pT()/GeV;
      const double abs_photon_eta = photon.abseta();
      const double weight = event.weight();

      if (abs_jet_eta < 1.5) {
        if      (abs_photon_eta < 0.9)  _h_phverycentral_jetcentral->fill(photon_pt, weight);
        else if (abs_photon_eta < 1.44) _h_phcentral_jetcentral->fill(    photon_pt, weight);
        else if (abs_photon_eta < 1.57) {}
        else if (abs_photon_eta < 2.1)  _h_phforward_jetcentral->fill(    photon_pt, weight);
        else if (abs_photon_eta < 2.5)  _h_phveryforward_jetcentral->fill(photon_pt, weight);
      }
      else if (abs_jet_eta < 2.5) {
        if      (abs_photon_eta < 0.9)  _h_phverycentral_jetforward->fill(photon_pt, weight);
        else if (abs_photon_eta < 1.44) _h_phcentral_jetforward->fill(    photon_pt, weight);
        else if (abs_photon_eta < 1.57) {}
        else if (abs_photon_eta < 2.1)  _h_phforward_jetforward->fill(    photon_pt, weight);
        else if (abs_photon_eta < 2.5)  _h_phveryforward_jetforward->fill(photon_pt, weight);
      }
    }
    


    /// Normalise histograms etc., after the run
    void finalize() {
      const double scale_jetcentral = crossSection()/sumOfWeights(); // *3 (jet eta < 1.5)
      scale(_h_phverycentral_jetcentral, scale_jetcentral); // * 1.8 (photon eta < 0.9)
      scale(_h_phcentral_jetcentral    , scale_jetcentral); // * 1.08 (0.9 .. 1.44)
      scale(_h_phforward_jetcentral    , scale_jetcentral); // * 1.06 (1.57 .. 2.1)
      scale(_h_phveryforward_jetcentral, scale_jetcentral); // * 0.8  (2.1 .. 2.5)

      const double scale_jetforward = crossSection()/sumOfWeights(); // *2 (1.5 < eta < 2.5)
      scale(_h_phverycentral_jetforward, scale_jetforward); // .. as above ..
      scale(_h_phcentral_jetforward    , scale_jetforward); // .. as above ..
      scale(_h_phforward_jetforward    , scale_jetforward); // .. as above ..
      scale(_h_phveryforward_jetforward, scale_jetforward); // .. as above ..

    }


  private:

    Histo1DPtr _h_phverycentral_jetcentral;
    Histo1DPtr _h_phcentral_jetcentral    ;
    Histo1DPtr _h_phforward_jetcentral    ;
    Histo1DPtr _h_phveryforward_jetcentral;

    Histo1DPtr _h_phverycentral_jetforward;
    Histo1DPtr _h_phcentral_jetforward    ;
    Histo1DPtr _h_phforward_jetforward    ;
    Histo1DPtr _h_phveryforward_jetforward;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2014_I1266056);


}
