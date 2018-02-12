// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  /// CMS 8 TeV dilepton channel ttbar charge asymmetry analysis
  class CMS_2016_I1430892 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1430892);


    /// Book histograms and initialise projections
    void init() {

      // Complete final state
      FinalState fs(-MAXDOUBLE, MAXDOUBLE, 0*GeV);

      // Projection for dressed electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressed_electrons(photons, electrons, 0.1);
      addProjection(dressed_electrons, "DressedElectrons");

      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressed_muons(photons, muons, 0.1);
      addProjection(dressed_muons, "DressedMuons");

      // Parton-level top quarks
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicPartonTops");


      // Booking of histograms

      // This histogram is independent of the parton-level information, and is an
      // addition to the original analysis. It is compared to the same data as
      // the parton-level delta_abseta histogram d05-x01-y01.
      _h_dabsetadressedleptons = bookHisto1D("d00-x01-y01", _bins_dabseta);

      // The remaining histos use parton-level information
      _h_dabseta = bookHisto1D("d05-x01-y01", _bins_dabseta);
      _h_dabsrapidity = bookHisto1D("d02-x01-y01", _bins_dabsrapidity);

      // 2D histos
      _h_dabsrapidity_var[0] = bookHisto2D("d11-x01-y01", _bins_dabsrapidity, _bins_tt_mass);
      _h_dabseta_var[0] = bookHisto2D("d17-x01-y01", _bins_dabseta, _bins_tt_mass);

      _h_dabsrapidity_var[1] = bookHisto2D("d23-x01-y01", _bins_dabsrapidity, _bins_tt_pT);
      _h_dabseta_var[1] = bookHisto2D("d29-x01-y01", _bins_dabseta, _bins_tt_pT);

      _h_dabsrapidity_var[2] = bookHisto2D("d35-x01-y01", _bins_dabsrapidity, _bins_tt_absrapidity);
      _h_dabseta_var[2] = bookHisto2D("d41-x01-y01", _bins_dabseta, _bins_tt_absrapidity);

      // Profile histos for asymmetries
      _h_dabsrapidity_profile[0] = bookProfile1D("d08-x01-y01", _bins_tt_mass);
      _h_dabseta_profile[0] = bookProfile1D("d14-x01-y01", _bins_tt_mass);

      _h_dabsrapidity_profile[1] = bookProfile1D("d20-x01-y01", _bins_tt_pT);
      _h_dabseta_profile[1] = bookProfile1D("d26-x01-y01", _bins_tt_pT);

      _h_dabsrapidity_profile[2] = bookProfile1D("d32-x01-y01", _bins_tt_absrapidity);
      _h_dabseta_profile[2] = bookProfile1D("d38-x01-y01", _bins_tt_absrapidity);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Use particle-level leptons for the first histogram
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");

      const vector<DressedLepton> dressedels = dressed_electrons.dressedLeptons();
      const vector<DressedLepton> dressedmus = dressed_muons.dressedLeptons();

      const size_t ndressedel = dressedels.size();
      const size_t ndressedmu = dressedmus.size();

      // For the particle-level histogram, require exactly one electron and exactly one muon, to select the ttbar->emu channel.
      // Note this means ttbar->emu events with additional PromptFinalState dilepton pairs from the shower are vetoed - for PYTHIA8,
      // this affects ~0.5% of events, so the effect is well below the level of sensitivity of the measured distribution.
      if ( ndressedel == 1 && ndressedmu == 1 ) {

        const int electrontouse = 0, muontouse = 0;

        // Opposite-charge leptons only
        if ( sameSign(dressedels[electrontouse], dressedmus[muontouse]) ) {
          MSG_INFO("Error, e and mu have same charge, skipping event");
        } else {
          // Get the four-momenta of the positively- and negatively-charged leptons
          FourMomentum lepPlus = dressedels[electrontouse].charge() > 0 ? dressedels[electrontouse] : dressedmus[muontouse];
          FourMomentum lepMinus = dressedels[electrontouse].charge() > 0 ? dressedmus[muontouse] : dressedels[electrontouse];

          // Now calculate the variable
          double dabseta_temp = lepPlus.abseta() - lepMinus.abseta();

          fillWithUFOF( _h_dabsetadressedleptons, dabseta_temp, weight );
        }

      }


      // The remaining variables use parton-level information.

      // Get the leptonically decaying tops
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      Particles chargedleptons;

      unsigned int ntrueleptonictops = 0;
      bool oppositesign = false;

      if ( leptonicpartontops.size() == 2 ) {
        for (size_t k = 0; k < leptonicpartontops.size(); ++k) {

          // Get the lepton
          const Particle lepTop = leptonicpartontops[k];
          const auto isPromptChargedLepton = [](const Particle& p){return (isChargedLepton(p) && isPrompt(p, false, false));};
          Particles lepton_candidates = lepTop.allDescendants(firstParticleWith(isPromptChargedLepton), false);
          if ( lepton_candidates.size() < 1 ) MSG_WARNING("error, PartonicTops::E_MU top quark had no daughter lepton candidate, skipping event.");

          // In some cases there is no lepton from the W decay but only leptons from the decay of a radiated gamma.
          // These hadronic PartonicTops are currently being mistakenly selected by PartonicTops::E_MU (as of April 2017), and need to be rejected.
          // PartonicTops::E_MU is being fixed in Rivet, and when it is the veto below should do nothing.
          /// @todo Should no longer be necessary -- remove
          bool istrueleptonictop = false;
          for (size_t i = 0; i < lepton_candidates.size(); ++i) {
            const Particle& lepton_candidate = lepton_candidates[i];
            if ( lepton_candidate.hasParent(PID::PHOTON) ) {
              MSG_DEBUG("Found gamma parent, top: " << k+1 << " of " << leptonicpartontops.size() << " , lepton: " << i+1 << " of " << lepton_candidates.size());
              continue;
            }
            if ( !istrueleptonictop && sameSign(lepTop,lepton_candidate) ) {
              chargedleptons.push_back(lepton_candidate);
              istrueleptonictop = true;
            }
            else MSG_WARNING("Error, found extra prompt charged lepton from top decay (and without gamma parent), ignoring it.");
          }
          if ( istrueleptonictop ) ++ntrueleptonictops;
        }
      }

      if ( ntrueleptonictops == 2 ) {
        oppositesign = !( sameSign(chargedleptons[0],chargedleptons[1]) );
        if ( !oppositesign ) MSG_WARNING("error, same charge tops, skipping event.");
      }

      if ( ntrueleptonictops == 2 && oppositesign ) {

        // Get the four-momenta of the positively- and negatively-charged leptons
        const FourMomentum lepPlus = chargedleptons[0].charge() > 0 ? chargedleptons[0] : chargedleptons[1];
        const FourMomentum lepMinus = chargedleptons[0].charge() > 0 ? chargedleptons[1] : chargedleptons[0];

        const double dabseta_temp = lepPlus.abseta() - lepMinus.abseta();

        // Get the four-momenta of the positively- and negatively-charged tops
        const FourMomentum topPlus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[0] : leptonicpartontops[1];
        const FourMomentum topMinus_p4 = leptonicpartontops[0].pdgId() > 0 ? leptonicpartontops[1] : leptonicpartontops[0];

        const FourMomentum ttbar_p4 = topPlus_p4 + topMinus_p4;

        const double tt_mass_temp = ttbar_p4.mass();
        const double tt_absrapidity_temp = ttbar_p4.absrapidity();
        const double tt_pT_temp = ttbar_p4.pT();
        const double dabsrapidity_temp = topPlus_p4.absrapidity() - topMinus_p4.absrapidity();

        // Fill parton-level histos
        fillWithUFOF( _h_dabseta, dabseta_temp, weight );
        fillWithUFOF( _h_dabsrapidity, dabsrapidity_temp, weight );

        // Now fill the same variables in the 2D and profile histos vs ttbar invariant mass, pT, and absolute rapidity
        for (int i_var = 0; i_var < 3; ++i_var) {
          double var;
          if ( i_var == 0 ) {
            var = tt_mass_temp;
          } else if ( i_var == 1 ) {
            var = tt_pT_temp;
          } else {
            var = tt_absrapidity_temp;
          }

          fillWithUFOF( _h_dabsrapidity_var[i_var], dabsrapidity_temp, var, weight );
          fillWithUFOF( _h_dabseta_var[i_var], dabseta_temp, var, weight );

          fillWithUFOF( _h_dabsrapidity_profile[i_var], dabsrapidity_temp, var, weight, (_h_dabsrapidity->xMax() + _h_dabsrapidity->xMin())/2. );
          fillWithUFOF( _h_dabseta_profile[i_var], dabseta_temp, var, weight, (_h_dabseta->xMax() + _h_dabseta->xMin())/2. );
        }

      }

    }


    /// Normalise histograms to unit area
    void finalize() {

      normalize(_h_dabsetadressedleptons);

      normalize(_h_dabseta);
      normalize(_h_dabsrapidity);

      for (int i_var = 0; i_var < 3; ++i_var) {
        normalize(_h_dabsrapidity_var[i_var]);
        normalize(_h_dabseta_var[i_var]);
      }

    }


  private:

    Histo1DPtr _h_dabsetadressedleptons, _h_dabseta, _h_dabsrapidity;
    Histo2DPtr _h_dabseta_var[3], _h_dabsrapidity_var[3];
    Profile1DPtr _h_dabseta_profile[3], _h_dabsrapidity_profile[3];

    const vector<double> _bins_tt_mass = {300., 430., 530., 1200.};
    const vector<double> _bins_tt_pT = {0., 41., 92., 300.};
    const vector<double> _bins_tt_absrapidity = {0., 0.34, 0.75, 1.5};
    const vector<double> _bins_dabseta = { -2., -68./60., -48./60., -32./60., -20./60., -8./60., 0., 8./60., 20./60., 32./60., 48./60., 68./60., 2.};
    const vector<double> _bins_dabsrapidity = {-2., -44./60., -20./60., 0., 20./60., 44./60., 2.};

    void fillWithUFOF(Histo1DPtr h, double x, double w) {
      h->fill(std::max(std::min(x, h->xMax()-1e-9),h->xMin()+1e-9), w);
    }

    void fillWithUFOF(Histo2DPtr h, double x, double y, double w) {
      h->fill(std::max(std::min(x, h->xMax()-1e-9),h->xMin()+1e-9), std::max(std::min(y, h->yMax()-1e-9),h->yMin()+1e-9), w);
    }

    void fillWithUFOF(Profile1DPtr h, double x, double y, double w, double c) {
      h->fill(std::max(std::min(y, h->xMax()-1e-9),h->xMin()+1e-9), float(x > c) - float(x < c), w);
    }


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1430892);


}
