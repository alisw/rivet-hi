// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VisibleFinalState.hh"

namespace Rivet {


  class ATLAS_2015_I1404878 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1404878);


    void init() {
      // Eta ranges
      Cut eta_full = (Cuts::abseta < 4.2) & (Cuts::pT >= 1.0*MeV);
      Cut lep_cuts = (Cuts::abseta < 2.5) && (Cuts::pT > 25*GeV);

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

      DressedLeptons dressedelectrons(photons, electrons, 0.1, lep_cuts, true);
      declare(dressedelectrons, "dressedelectrons");

      DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full, true);
      declare(ewdressedelectrons, "ewdressedelectrons");

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);

      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      declare(muons, "muons");

      DressedLeptons dressedmuons(photons, muons, 0.1, lep_cuts, true);
      declare(dressedmuons, "dressedmuons");

      DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full, true);
      declare(ewdressedmuons, "ewdressedmuons");

      // Projection to find neutrinos
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);

      // get MET from generic invisibles
      VetoedFinalState inv_fs(fs);
      inv_fs.addVetoOnThisFinalState(VisibleFinalState(fs));
      declare(inv_fs, "InvisibleFS");

      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(ewdressedelectrons);
      vfs.addVetoOnThisFinalState(ewdressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "jets");


      // Histogram booking
      _h["massttbar"]                  = bookHisto1D( 1, 1, 1);
      _h["massttbar_norm"]             = bookHisto1D( 2, 1, 1);
      _h["ptttbar"]                    = bookHisto1D( 3, 1, 1);
      _h["ptttbar_norm"]               = bookHisto1D( 4, 1, 1);
      _h["absrapttbar"]                = bookHisto1D( 5, 1, 1);
      _h["absrapttbar_norm"]           = bookHisto1D( 6, 1, 1);
      _h["ptpseudotophadron"]          = bookHisto1D( 7, 1, 1);
      _h["ptpseudotophadron_norm"]     = bookHisto1D( 8, 1, 1);
      _h["absrappseudotophadron"]      = bookHisto1D( 9, 1, 1);
      _h["absrappseudotophadron_norm"] = bookHisto1D(10, 1, 1);
      _h["absPout"]                    = bookHisto1D(11, 1, 1);
      _h["absPout_norm"]               = bookHisto1D(12, 1, 1);
      _h["dPhittbar"]                  = bookHisto1D(13, 1, 1);
      _h["dPhittbar_norm"]             = bookHisto1D(14, 1, 1);
      _h["HTttbar"]                    = bookHisto1D(15, 1, 1);
      _h["HTttbar_norm"]               = bookHisto1D(16, 1, 1);
      _h["Yboost"]                     = bookHisto1D(17, 1, 1);
      _h["Yboost_norm"]                = bookHisto1D(18, 1, 1);
      _h["chittbar"]                   = bookHisto1D(19, 1, 1);
      _h["chittbar_norm"]              = bookHisto1D(20, 1, 1);
      _h["RWt"]                        = bookHisto1D(21, 1, 1);
      _h["RWt_norm"]                   = bookHisto1D(22, 1, 1);

    }

    void analyze(const Event& event) {

      const double weight = event.weight();

      // Get the selected objects, using the projections.
      vector<DressedLepton> electrons = applyProjection<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      vector<DressedLepton> muons     = applyProjection<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      const Jets& jets = applyProjection<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      const FinalState& ifs = applyProjection<FinalState>(event, "InvisibleFS");

      // Calculate MET
      FourMomentum met;
      for (const Particle& p : ifs.particles())  met += p.momentum();

      // Count the number of b-tags
      Jets bjets, lightjets;
      for (const Jet& jet : jets){
        bool b_tagged = jet.bTags(Cuts::pT > 5*GeV).size();
        if ( b_tagged && bjets.size() < 2 )  bjets += jet;
        else lightjets += jet;
      }

      bool single_electron = (electrons.size() == 1) && (muons.empty());
      bool single_muon     = (muons.size() == 1) && (electrons.empty());

      DressedLepton* lepton = nullptr;
      if (single_electron)   lepton = &electrons[0];
      else if (single_muon)  lepton = &muons[0];

      if(!single_electron && !single_muon)  vetoEvent;
      if (jets.size()  < 4 || bjets.size() < 2)  vetoEvent;

      FourMomentum pbjet1; // Momentum of bjet1
      FourMomentum pbjet2; // Momentum of bjet2
      if ( deltaR(bjets[0], *lepton) <= deltaR(bjets[1], *lepton) ) {
        pbjet1 = bjets[0].momentum();
        pbjet2 = bjets[1].momentum();
      } else {
        pbjet1 = bjets[1].momentum();
        pbjet2 = bjets[0].momentum();
      }


      double bestWmass = 1000.0*TeV;
      double mWPDG = 80.399*GeV;
      int Wj1index = -1, Wj2index = -1;
      for (unsigned int i = 0; i < (lightjets.size() - 1); ++i) {
        for (unsigned int j = i + 1; j < lightjets.size(); ++j) {
          double wmass = (lightjets[i].momentum() + lightjets[j].momentum()).mass();
          if (fabs(wmass - mWPDG) < fabs(bestWmass - mWPDG)) {
            bestWmass = wmass;
            Wj1index = i;
            Wj2index = j;
          }
        }
      }

      // Compute hadronic W boson
      FourMomentum pWhadron = lightjets[Wj1index].momentum() + lightjets[Wj2index].momentum();
      double pz = _computeneutrinoz(lepton->momentum(), met);
      FourMomentum ppseudoneutrino( sqrt(sqr(met.px()) + sqr(met.py()) + sqr(pz)), met.px(), met.py(), pz);


      // Compute leptonic, hadronic, combined pseudo-top
      FourMomentum ppseudotoplepton = lepton->momentum() + ppseudoneutrino + pbjet1;
      FourMomentum ppseudotophadron = pbjet2 + pWhadron;
      FourMomentum pttbar = ppseudotoplepton + ppseudotophadron;

      Vector3 z_versor(0,0,1);
      Vector3 vpseudotophadron = ppseudotophadron.vector3();
      Vector3 vpseudotoplepton = ppseudotoplepton.vector3();
      // Observables
      double ystar = 0.5 * deltaRap(ppseudotophadron, ppseudotoplepton);
      double chi_ttbar = exp(2 * fabs(ystar));
      double deltaPhi_ttbar = deltaPhi(ppseudotoplepton,ppseudotophadron);
      double HT_ttbar = ppseudotophadron.pt() + ppseudotoplepton.pt();
      double Yboost = 0.5 * fabs(ppseudotophadron.rapidity() + ppseudotoplepton.rapidity());
      double R_Wt = pWhadron.pt() / ppseudotophadron.pt();
      double absPout = fabs(vpseudotophadron.dot((vpseudotoplepton.cross(z_versor))/(vpseudotoplepton.cross(z_versor).mod())));

      // absolute cross sections
      _h["ptpseudotophadron"]->fill(    ppseudotophadron.pt(),     weight); //pT of pseudo top hadron
      _h["ptttbar"]->fill(              pttbar.pt(),               weight); //fill pT of ttbar in combined channel
      _h["absrappseudotophadron"]->fill(ppseudotophadron.absrap(), weight);
      _h["absrapttbar"]->fill(          pttbar.absrap(),           weight);
      _h["massttbar"]->fill(            pttbar.mass(),             weight);
      _h["absPout"]->fill(              absPout,                   weight);
      _h["chittbar"]->fill(             chi_ttbar,                 weight);
      _h["dPhittbar"]->fill(            deltaPhi_ttbar,            weight);
      _h["HTttbar"]->fill(              HT_ttbar,                  weight);
      _h["Yboost"]->fill(               Yboost,                    weight);
      _h["RWt"]->fill(                  R_Wt,                      weight);
      // normalised cross sections
      _h["ptpseudotophadron_norm"]->fill(    ppseudotophadron.pt(),     weight); //pT of pseudo top hadron
      _h["ptttbar_norm"]->fill(              pttbar.pt(),               weight); //fill pT of ttbar in combined channel
      _h["absrappseudotophadron_norm"]->fill(ppseudotophadron.absrap(), weight);
      _h["absrapttbar_norm"]->fill(          pttbar.absrap(),           weight);
      _h["massttbar_norm"]->fill(            pttbar.mass(),             weight);
      _h["absPout_norm"]->fill(              absPout,                   weight);
      _h["chittbar_norm"]->fill(             chi_ttbar,                 weight);
      _h["dPhittbar_norm"]->fill(            deltaPhi_ttbar,            weight);
      _h["HTttbar_norm"]->fill(              HT_ttbar,                  weight);
      _h["Yboost_norm"]->fill(               Yboost,                    weight);
      _h["RWt_norm"]->fill(                  R_Wt,                      weight);

    }

    void finalize() {
      // Normalize to cross-section
      const double sf = crossSection() / sumOfWeights();
      for (auto& k_h : _h) {
        scale(k_h.second, sf);
        if (k_h.first.find("_norm") != string::npos) normalize(k_h.second);
      }
    }


  private:

    // Compute z component of neutrino momentum given lepton and met
    double _computeneutrinoz(const FourMomentum& lepton, FourMomentum& met) const {
      double m_W = 80.399; // in GeV, given in the paper
      double k = (( sqr( m_W ) - sqr( lepton.mass() ) ) / 2 ) + (lepton.px() * met.px() + lepton.py() * met.py());
      double a = sqr ( lepton.E() )- sqr ( lepton.pz() );
      double b = -2*k*lepton.pz();
      double c = sqr( lepton.E() ) * sqr( met.pT() ) - sqr( k );
      double discriminant = sqr(b) - 4 * a * c;
      double quad[2] = { (- b - sqrt(discriminant)) / (2 * a), (- b + sqrt(discriminant)) / (2 * a) }; //two possible quadratic solns

      double pzneutrino;
      if (discriminant < 0) { // if the discriminant is negative:
        pzneutrino = - b / (2 * a);
      } else { // if the discriminant is positive, take the soln with smallest absolute value
        pzneutrino = (fabs(quad[0]) < fabs(quad[1])) ? quad[0] : quad[1];
      }
      return pzneutrino;
    }

    /// @todo Replace with central version
    double _mT(const FourMomentum &l, FourMomentum &nu) const {
      return sqrt( 2 * l.pT() * nu.pT() * (1 - cos(deltaPhi(l, nu))) );
    }


    /// @name Objects that are used by the event selection decisions
    map<string, Histo1DPtr> _h;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1404878);


}
