// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ATLAS 7 TeV pseudo-top analysis
  ///
  /// @author K .Finelli <kevin.finelli@cern.ch>
  /// @author A. Saavedra <a.saavedra@physics.usyd.edu.au>
  /// @author L. Lan <llan@physics.usyd.edu.au>
  class ATLAS_2015_I1345452 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1345452);


    void init() {
      // Eta ranges
      Cut eta_full = (Cuts::abseta < 5.0) & (Cuts::pT >= 1.0*MeV);
      Cut eta_lep = (Cuts::abseta < 2.5);

      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);

      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      declare(electrons, "electrons");

      DressedLeptons dressedelectrons(photons, electrons, 0.1, eta_lep && Cuts::pT > 25*GeV, true);
      declare(dressedelectrons, "dressedelectrons");

      DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full, true);
      declare(ewdressedelectrons, "ewdressedelectrons");

      DressedLeptons vetodressedelectrons(photons, electrons, 0.1, eta_lep && Cuts::pT > 15*GeV, true);
      declare(vetodressedelectrons, "vetodressedelectrons");

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      declare(muons, "muons");
      DressedLeptons dressedmuons(photons, muons, 0.1, eta_lep && Cuts::pT > 25*GeV, true);
      declare(dressedmuons, "dressedmuons");
      DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full, true);
      declare(ewdressedmuons, "ewdressedmuons");
      DressedLeptons vetodressedmuons(photons, muons, 0.1, eta_lep && Cuts::pT > 15*GeV, true);
      declare(vetodressedmuons, "vetodressedmuons");

      // Projection to find neutrinos and produce MET
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      declare(neutrinos, "neutrinos");

      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(ewdressedelectrons);
      vfs.addVetoOnThisFinalState(ewdressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");

      //pseudotop leptons and hadrons
      _h["ptpseudotophadron_mu"]     = bookHisto1D( 1, 1, 2);
      _h["ptpseudotophadron_el"]     = bookHisto1D( 2, 1, 2);
      _h["absrappseudotophadron_mu"] = bookHisto1D( 3, 1, 2);
      _h["absrappseudotophadron_el"] = bookHisto1D( 4, 1, 2);
      _h["ptpseudotoplepton_mu"]     = bookHisto1D( 5, 1, 2);
      _h["ptpseudotoplepton_el"]     = bookHisto1D( 6, 1, 2);
      _h["absrappseudotoplepton_mu"] = bookHisto1D( 7, 1, 2);
      _h["absrappseudotoplepton_el"] = bookHisto1D( 8, 1, 2);
      _h["ptttbar_mu"]               = bookHisto1D( 9, 1, 2);
      _h["ptttbar_el"]               = bookHisto1D(10, 1, 2);
      _h["absrapttbar_mu"]           = bookHisto1D(11, 1, 2);
      _h["absrapttbar_el"]           = bookHisto1D(12, 1, 2);
      _h["ttbarmass_mu"]             = bookHisto1D(13, 1, 2);
      _h["ttbarmass_el"]             = bookHisto1D(14, 1, 2);
      _h["ptpseudotophadron"]        = bookHisto1D(15, 1, 2);
      _h["absrappseudotophadron"]    = bookHisto1D(16, 1, 2);
      _h["ptpseudotoplepton"]        = bookHisto1D(17, 1, 2);
      _h["absrappseudotoplepton"]    = bookHisto1D(18, 1, 2);
      _h["ptttbar"]                  = bookHisto1D(19, 1, 2);
      _h["absrapttbar"]              = bookHisto1D(20, 1, 2);
      _h["ttbarmass"]                = bookHisto1D(21, 1, 2);

    }

    void analyze(const Event& event) {

      // Get the selected objects, using the projections.
      _dressedelectrons     = apply<DressedLeptons>(  event, "dressedelectrons").dressedLeptons();
      _vetodressedelectrons = apply<DressedLeptons>(  event, "vetodressedelectrons").dressedLeptons();
      _dressedmuons         = apply<DressedLeptons>(  event, "dressedmuons").dressedLeptons();
      _vetodressedmuons     = apply<DressedLeptons>(  event, "vetodressedmuons").dressedLeptons();
      _neutrinos            = apply<PromptFinalState>(event, "neutrinos").particlesByPt();
      const Jets& all_jets  = apply<FastJets>(        event, "jets").jetsByPt(Cuts::pT > 25.0*GeV && Cuts::abseta < 2.5);

      //get true l+jets events by removing events with more than 1 electron||muon neutrino
      unsigned int n_elmu_neutrinos = 0;
      foreach (const Particle p, _neutrinos) {
        if (p.abspid() == 12 || p.abspid() == 14)  ++n_elmu_neutrinos;
      }
      if (n_elmu_neutrinos != 1)  vetoEvent;

      DressedLepton *lepton;
      if ( _dressedelectrons.size())  lepton = &_dressedelectrons[0];
      else if (_dressedmuons.size())  lepton = &_dressedmuons[0];
      else vetoEvent;

      // Calculate the missing ET, using the prompt neutrinos only (really?)
      /// @todo Why not use MissingMomentum?
      FourMomentum met;
      foreach (const Particle& p, _neutrinos)  met += p.momentum();

      //remove jets if they are within dR < 0.2 of lepton
      Jets jets;
      foreach(const Jet& jet, all_jets) {
        bool keep = true;
        foreach (const DressedLepton& el, _vetodressedelectrons) {
          keep &= deltaR(jet, el) >= 0.2;
        }
        if (keep)  jets += jet;
      }

      bool overlap = false;
      Jets bjets, lightjets;
      for (unsigned int i = 0; i < jets.size(); ++i) {
        const Jet& jet = jets[i];
        foreach (const DressedLepton& el, _dressedelectrons)  overlap |= deltaR(jet, el) < 0.4;
        foreach (const DressedLepton& mu, _dressedmuons)      overlap |= deltaR(jet, mu) < 0.4;
        for (unsigned int j = i + 1; j < jets.size(); ++j) {
          overlap |= deltaR(jet, jets[j]) < 0.5;
        }
        //// Count the number of b-tags
        bool b_tagged = false;           //  This is closer to the
        Particles bTags = jet.bTags();   //  analysis. Something
        foreach ( Particle b, bTags ) {  //  about ghost-associated
          b_tagged |= b.pT() > 5*GeV;    //  B-hadrons
        }                                //
        if ( b_tagged )  bjets += jet;
        else lightjets += jet;
      }

      // remove events with object overlap
      if (overlap) vetoEvent;

      if (bjets.size() < 2 || lightjets.size() < 2)  vetoEvent;

      FourMomentum pbjet1; //Momentum of bjet1
      FourMomentum pbjet2; //Momentum of bjet2
      if ( deltaR(bjets[0], *lepton) <= deltaR(bjets[1], *lepton) ) {
        pbjet1 = bjets[0].momentum();
        pbjet2 = bjets[1].momentum();
      } else {
        pbjet1 = bjets[1].momentum();
        pbjet2 = bjets[0].momentum();
      }

      FourMomentum pjet1; // Momentum of jet1
      if (lightjets.size())  pjet1 = lightjets[0].momentum();

      FourMomentum pjet2; // Momentum of jet 2
      if (lightjets.size() > 1)  pjet2 = lightjets[1].momentum();

      double pz = computeneutrinoz(lepton->momentum(), met);
      FourMomentum ppseudoneutrino( sqrt(sqr(met.px()) + sqr(met.py()) + sqr(pz)), met.px(), met.py(), pz);

      //compute leptonic, hadronic, combined pseudo-top
      FourMomentum ppseudotoplepton = lepton->momentum() + ppseudoneutrino + pbjet1;
      FourMomentum ppseudotophadron = pbjet2 + pjet1 + pjet2;
      FourMomentum pttbar = ppseudotoplepton + ppseudotophadron;

      // Evaluate basic event selection
      bool pass_eljets = (_dressedelectrons.size() == 1) &&
                                                                                           (_vetodressedelectrons.size() < 2) &&
        (_vetodressedmuons.empty()) &&
        (met.pT() > 30*GeV) &&
                                                                                           (_mT(_dressedelectrons[0].momentum(), met) > 35*GeV) &&
                                                                                           (jets.size() >= 4);
      bool pass_mujets = (_dressedmuons.size() == 1) &&
        (_vetodressedmuons.size() < 2) &&
        (_vetodressedelectrons.empty()) &&
        (met.pT() > 30*GeV) &&
        (_mT(_dressedmuons[0].momentum(), met) > 35*GeV) &&
        (jets.size() >= 4);

      // basic event selection requirements
      if (!pass_eljets && !pass_mujets) vetoEvent;

      // Fill histograms
      const double weight = event.weight();

      //pseudotop hadrons and leptons fill histogram
      _h["ptpseudotoplepton"]->fill(    ppseudotoplepton.pt(),     weight); //pT of pseudo top lepton
      _h["absrappseudotoplepton"]->fill(ppseudotoplepton.absrap(), weight); //absolute rapidity of pseudo top lepton
      _h["ptpseudotophadron"]->fill(    ppseudotophadron.pt(),     weight); //pT of pseudo top hadron
      _h["absrappseudotophadron"]->fill(ppseudotophadron.absrap(), weight); //absolute rapidity of pseudo top hadron
      _h["absrapttbar"]->fill(          pttbar.absrap(),           weight); //absolute rapidity of ttbar
      _h["ttbarmass"]->fill(            pttbar.mass(),             weight); //mass of ttbar
      _h["ptttbar"]->fill(              pttbar.pt(),               weight); //fill pT of ttbar in combined channel

      if (pass_eljets) { // electron channel fill histogram
        _h["ptpseudotoplepton_el"]->fill(    ppseudotoplepton.pt(),     weight); //pT of pseudo top lepton
        _h["absrappseudotoplepton_el"]->fill(ppseudotoplepton.absrap(), weight); //absolute rapidity of pseudo top lepton
        _h["ptpseudotophadron_el"]->fill(    ppseudotophadron.pt(),     weight); //pT of pseudo top hadron
        _h["absrappseudotophadron_el"]->fill(ppseudotophadron.absrap(), weight); //absolute rapidity of pseudo top hadron
        _h["absrapttbar_el"]->fill(          pttbar.absrap(),           weight); //absolute rapidity of ttbar
        _h["ttbarmass_el"]->fill(            pttbar.mass(),             weight); // fill electron channel ttbar mass
        _h["ptttbar_el"]->fill(              pttbar.pt(),               weight); //fill pT of ttbar in electron channel
      }
      else { // muon channel fill histogram
        _h["ptpseudotoplepton_mu"]->fill(    ppseudotoplepton.pt(),     weight); //pT of pseudo top lepton
        _h["absrappseudotoplepton_mu"]->fill(ppseudotoplepton.absrap(), weight); //absolute rapidity of pseudo top lepton
        _h["ptpseudotophadron_mu"]->fill(    ppseudotophadron.pt(),     weight); //pT of pseudo top hadron
        _h["absrappseudotophadron_mu"]->fill(ppseudotophadron.absrap(), weight); //absolute rapidity of pseudo top hadron
        _h["absrapttbar_mu"]->fill(          pttbar.absrap(),           weight); //absolute rapidity of ttbar
        _h["ttbarmass_mu"]->fill(            pttbar.mass(),             weight); //fill muon channel histograms
        _h["ptttbar_mu"]->fill(              pttbar.pt(),               weight); //fill pT of ttbar in electron channel
      }
    }

    void finalize() {
      // Normalize to cross-section
      const double scalefactor(crossSection() / sumOfWeights());
      for (map<string, Histo1DPtr>::iterator hit = _h.begin(); hit != _h.end(); ++hit) {
        double sf = scalefactor;
        if ( (hit->first).find("_") == std::string::npos )  sf *= 0.5;
        scale(hit->second, sf);
      }
    }

  private:


    double computeneutrinoz(const FourMomentum& lepton, FourMomentum& met) const {
      //computing z component of neutrino momentum given lepton and met
      double pzneutrino;
      double m_W = 80.399; // in GeV, given in the paper
      double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.px() + lepton.py() * met.py());
      double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
      double b = -2*k*lepton.pz();
      double c = sqr( lepton.E() ) * sqr( met.pT() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns
      if (discriminant < 0)  pzneutrino = - b / (2 * a); //if the discriminant is negative
      else { //if the discriminant is greater than or equal to zero, take the soln with smallest absolute value
        double absquad[2];
        for (int n=0; n<2; ++n)  absquad[n] = fabs(quad[n]);
        if (absquad[0] < absquad[1])  pzneutrino = quad[0];
        else                          pzneutrino = quad[1];
      }
      if ( !std::isfinite(pzneutrino) )  std::cout << "Found non-finite value" << std::endl;
      return pzneutrino;
    }

    double _mT(const FourMomentum &l, FourMomentum &nu) const {
      return sqrt( 2 * l.pT() * nu.pT() * (1 - cos(deltaPhi(l, nu))) );
    }

    /// @name Objects that are used by the event selection decisions
    vector<DressedLepton> _dressedelectrons, _vetodressedelectrons, _dressedmuons, _vetodressedmuons;
    Particles _neutrinos;
    map<string, Histo1DPtr> _h;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1345452);

}
