// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Hard double-parton scattering in four-jet events at 7 TeV
  class ATLAS_2016_I1479760 : public Analysis {
  public:


    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1479760);


    /// Book histograms and initialise projections before the run
    void init() {

      /// Declare AntiKt 0.6 jets without muons and neutrinos
      FastJets fastJets(FinalState(), FastJets::ANTIKT, 0.6);
      fastJets.useInvisibles(JetAlg::NO_INVISIBLES);
      fastJets.useMuons(JetAlg::NO_MUONS);
      addProjection(fastJets, "AntiKt6Jets");

      _hists["deltaPt34"]        = bookHisto1D( 1, 1, 1);
      _hists["deltaPhi34"]       = bookHisto1D( 2, 1, 1);
      _hists["deltaPt12"]        = bookHisto1D( 3, 1, 1);
      _hists["deltaPt13"]        = bookHisto1D( 4, 1, 1);
      _hists["deltaPt23"]        = bookHisto1D( 5, 1, 1);
      _hists["deltaPt14"]        = bookHisto1D( 6, 1, 1);
      _hists["deltaPt24"]        = bookHisto1D( 7, 1, 1);
      _hists["deltaPhi12"]       = bookHisto1D( 8, 1, 1);
      _hists["deltaPhi13"]       = bookHisto1D( 9, 1, 1);
      _hists["deltaPhi23"]       = bookHisto1D(10, 1, 1);
      _hists["deltaPhi14"]       = bookHisto1D(11, 1, 1);
      _hists["deltaPhi24"]       = bookHisto1D(12, 1, 1);
      _hists["deltaY12"]         = bookHisto1D(13, 1, 1);
      _hists["deltaY34"]         = bookHisto1D(14, 1, 1);
      _hists["deltaY13"]         = bookHisto1D(15, 1, 1);
      _hists["deltaY23"]         = bookHisto1D(16, 1, 1);
      _hists["deltaY14"]         = bookHisto1D(17, 1, 1);
      _hists["deltaY24"]         = bookHisto1D(18, 1, 1);
      _hists["deltaPhiPlanes12"] = bookHisto1D(19, 1, 1);
      _hists["deltaPhiPlanes13"] = bookHisto1D(20, 1, 1);
      _hists["deltaPhiPlanes14"] = bookHisto1D(21, 1, 1);
    }


    /// Calculate the DeltaPt variable
    double calcDeltaPt(const Jet& j1, const Jet& j2) {
      return  (j1.momentum() + j2.momentum()).pT() / (j1.pT() + j2.pT());
    }

    /// Calculate the DeltaPhi variable between event planes
    double calcDeltaPhiPlanes(const Jet& j1, const Jet& j2, const Jet& j3, const Jet& j4) {
      const FourMomentum sumVec1 = j1.momentum() + j2.momentum();
      const FourMomentum sumVec2 = j3.momentum() + j4.momentum();
      return  deltaPhi(sumVec1, sumVec2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve all anti-kt R=0.6 jets with pT above 20 GeV and eta < 4.4
      const Jets jets = applyProjection<JetAlg>(event, "AntiKt6Jets").jetsByPt(Cuts::pT >= 20*GeV && Cuts::abseta <= 4.4);

      // Require at least 4 jets, with the leading jet pT above 42.5 GeV
      if (jets.size() < 4) vetoEvent;
      if (jets[0].pT() < 42.5*GeV) vetoEvent;

      /// Fill histograms
      const double weight = event.weight();
      _hists["deltaPt12"]->fill( calcDeltaPt( jets[0], jets[1] ), weight );
      _hists["deltaPt34"]->fill( calcDeltaPt( jets[2], jets[3] ), weight );
      _hists["deltaPt13"]->fill( calcDeltaPt( jets[0], jets[2] ), weight );
      _hists["deltaPt23"]->fill( calcDeltaPt( jets[1], jets[2] ), weight );
      _hists["deltaPt14"]->fill( calcDeltaPt( jets[0], jets[3] ), weight );
      _hists["deltaPt24"]->fill( calcDeltaPt( jets[1], jets[3] ), weight );
      //
      _hists["deltaPhi12"]->fill( deltaPhi( jets[0],jets[1] ), weight );
      _hists["deltaPhi34"]->fill( deltaPhi( jets[2],jets[3] ), weight );
      _hists["deltaPhi13"]->fill( deltaPhi( jets[0],jets[2] ), weight );
      _hists["deltaPhi23"]->fill( deltaPhi( jets[1],jets[2] ), weight );
      _hists["deltaPhi14"]->fill( deltaPhi( jets[0],jets[3] ), weight );
      _hists["deltaPhi24"]->fill( deltaPhi( jets[1],jets[3] ), weight );
      //
      _hists["deltaY12"]->fill( deltaRap( jets[0], jets[1] ), weight );
      _hists["deltaY34"]->fill( deltaRap( jets[2], jets[3] ), weight );
      _hists["deltaY13"]->fill( deltaRap( jets[0], jets[2] ), weight );
      _hists["deltaY23"]->fill( deltaRap( jets[1], jets[2] ), weight );
      _hists["deltaY14"]->fill( deltaRap( jets[0], jets[3] ), weight );
      _hists["deltaY24"]->fill( deltaRap( jets[1], jets[3] ), weight );
      //
      _hists["deltaPhiPlanes12"]->fill( calcDeltaPhiPlanes(jets[0], jets[1], jets[2], jets[3] ), weight );
      _hists["deltaPhiPlanes13"]->fill( calcDeltaPhiPlanes(jets[0], jets[2], jets[1], jets[3] ), weight );
      _hists["deltaPhiPlanes14"]->fill( calcDeltaPhiPlanes(jets[0], jets[3], jets[1], jets[2] ), weight );
    }


    /// Post-run processing
    void finalize() {
      for (auto& key_hist : _hists)
        normalize(key_hist.second);
    }

    //@}


    /// Histograms
    map<string, Histo1DPtr> _hists;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1479760);

}
