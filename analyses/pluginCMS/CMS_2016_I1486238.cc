// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Studies of 2 b-jet + 2 jet production in proton-proton collisions at 7 TeV
  class CMS_2016_I1486238 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1486238);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FastJets akt(FinalState(), FastJets::ANTIKT, 0.5);
      addProjection(akt, "antikT");

      _h_Deltaphi_newway = bookHisto1D(1,1,1);
      _h_deltaphiafterlight = bookHisto1D(9,1,1);
      _h_SumPLight = bookHisto1D(5,1,1);

      _h_LeadingBJetpt = bookHisto1D(11,1,1);
      _h_SubleadingBJetpt = bookHisto1D(15,1,1);
      _h_LeadingLightJetpt = bookHisto1D(13,1,1);
      _h_SubleadingLightJetpt = bookHisto1D(17,1,1);

      _h_LeadingBJeteta = bookHisto1D(10,1,1);
      _h_SubleadingBJeteta = bookHisto1D(14,1,1);
      _h_LeadingLightJeteta = bookHisto1D(12,1,1);
      _h_SubleadingLightJeteta = bookHisto1D(16,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<JetAlg>(event, "antikT").jetsByPt(Cuts::absrap < 4.7 && Cuts::pT > 20*GeV);
      if (jets.size() < 4) vetoEvent;

      // Initial quarks
      /// @note Quark-level tagging...
      Particles bquarks;
      for (const GenParticle* p : particles(event.genEvent())) {
        if (abs(p->pdg_id()) == PID::BQUARK) bquarks += Particle(p);
      }
      Jets bjets, ljets;
      for (const Jet& j : jets) {
        const bool btag = any(bquarks, deltaRLess(j, 0.3));
        // for (const Particle& b : bquarks) if (deltaR(j, b) < 0.3) btag = true;
        (btag && j.abseta() < 2.4 ? bjets : ljets).push_back(j);
      }

      // Fill histograms
      const double weight = event.weight();
      if (bjets.size() >= 2 && ljets.size() >= 2) {
        _h_LeadingBJetpt->fill(bjets[0].pT()/GeV, weight);
        _h_SubleadingBJetpt->fill(bjets[1].pT()/GeV, weight);
        _h_LeadingLightJetpt->fill(ljets[0].pT()/GeV, weight);
        _h_SubleadingLightJetpt->fill(ljets[1].pT()/GeV, weight);
        //
        _h_LeadingBJeteta->fill(bjets[0].eta(), weight);
        _h_SubleadingBJeteta->fill(bjets[1].eta(), weight);
        _h_LeadingLightJeteta->fill(ljets[0].eta(), weight);
        _h_SubleadingLightJeteta->fill(ljets[1].eta(), weight);

        const double lightdphi = deltaPhi(ljets[0], ljets[1]);
        _h_deltaphiafterlight->fill(lightdphi, weight);

        const double vecsumlightjets = sqrt(sqr(ljets[0].px()+ljets[1].px()) + sqr(ljets[0].py()+ljets[1].py())); //< @todo Just (lj0+lj1).pT()? Or use add_quad
        const double term2 = vecsumlightjets/(sqrt(sqr(ljets[0].px()) + sqr(ljets[0].py())) + sqrt(sqr(ljets[1].px()) + sqr(ljets[1].py()))); //< @todo lj0.pT() + lj1.pT()? Or add_quad
        _h_SumPLight->fill(term2, weight);

        const double pxBsyst2 = bjets[0].px()+bjets[1].px(); // @todo (bj0+bj1).px()
        const double pyBsyst2 = bjets[0].py()+bjets[1].py(); // @todo (bj0+bj1).py()
        const double pxJetssyst2 = ljets[0].px()+ljets[1].px(); // @todo (lj0+lj1).px()
        const double pyJetssyst2 = ljets[0].py()+ljets[1].py(); // @todo (lj0+lj1).py()
        const double modulusB2 = sqrt(sqr(pxBsyst2)+sqr(pyBsyst2)); //< @todo add_quad
        const double modulusJets2 = sqrt(sqr(pxJetssyst2)+sqr(pyJetssyst2)); //< @todo add_quad
        const double cosphiBsyst2 = pxBsyst2/modulusB2;
        const double cosphiJetssyst2 = pxJetssyst2/modulusJets2;
        const double phiBsyst2 = ((pyBsyst2 > 0) ? 1 : -1) * acos(cosphiBsyst2); //< @todo sign(pyBsyst2)
        const double phiJetssyst2 = sign(pyJetssyst2) * acos(cosphiJetssyst2);
        const double Dphi2 = deltaPhi(phiBsyst2, phiJetssyst2);
        _h_Deltaphi_newway->fill(Dphi2,weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double invlumi = crossSection()/picobarn/sumOfWeights();
      normalize({_h_SumPLight, _h_deltaphiafterlight, _h_Deltaphi_newway});
      scale({_h_LeadingLightJetpt, _h_SubleadingLightJetpt, _h_LeadingBJetpt, _h_SubleadingBJetpt}, invlumi);
      scale({_h_LeadingLightJeteta, _h_SubleadingLightJeteta, _h_LeadingBJeteta, _h_SubleadingBJeteta}, invlumi);
    }

    //@}


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_deltaphiafterlight, _h_Deltaphi_newway, _h_SumPLight;
    Histo1DPtr _h_LeadingBJetpt, _h_SubleadingBJetpt, _h_LeadingLightJetpt, _h_SubleadingLightJetpt;
    Histo1DPtr _h_LeadingBJeteta, _h_SubleadingBJeteta, _h_LeadingLightJeteta, _h_SubleadingLightJeteta;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1486238);

}
