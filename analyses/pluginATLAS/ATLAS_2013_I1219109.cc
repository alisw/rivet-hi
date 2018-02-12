// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

namespace Rivet {


  /// @brief ATLAS W+b measurement
  class ATLAS_2013_I1219109: public Analysis {
  public:

    ATLAS_2013_I1219109(string name = "ATLAS_2013_I1219109")
      : Analysis(name)
    {
      // the electron mode is used by default
      _mode = 1;
    }


    void init() {
      FinalState fs;
      declare(fs, "FinalState");

      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT >= 25*GeV;

      // W finder for electrons and muons
      WFinder wf(fs, cuts, _mode==3? PID::MUON : PID::ELECTRON, 0.0*GeV, MAXDOUBLE, 0.0, 0.1,
                 WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(wf, "WF");

      // jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("WF"));
      FastJets fj(jet_fs, FastJets::ANTIKT, 0.4);
      fj.useInvisibles();
      declare(fj, "Jets");
      declare(HeavyHadrons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV), "BHadrons");


      // book histograms
      _njet     = bookHisto1D(1, 1, _mode); // dSigma / dNjet
      _jet1_bPt = bookHisto1D(2, 1, _mode); // dSigma / dBjetPt for Njet = 1
      _jet2_bPt = bookHisto1D(2, 2, _mode); // dSigma / dBjetPt for Njet = 2

    }


    void analyze(const Event& event) {

      const double weight = event.weight();

      //  retrieve W boson candidate
      const WFinder& wf = apply<WFinder>(event, "WF");
      if( wf.bosons().size() != 1 )  vetoEvent; // only one W boson candidate
      if( !(wf.mT() > 60.0*GeV) )    vetoEvent;
      //const Particle& Wboson  = wf.boson();


      // retrieve constituent neutrino
      const Particle& neutrino = wf.constituentNeutrino();
      if( !(neutrino.pT() > 25.0*GeV) )  vetoEvent;

      // retrieve constituent lepton
      const Particle& lepton = wf.constituentLepton();

      // count good jets, check if good jet contains B hadron
      const Particles& bHadrons = apply<HeavyHadrons>(event, "BHadrons").bHadrons();
      const Jets& jets = apply<JetAlg>(event, "Jets").jetsByPt(25*GeV);
      int goodjets = 0, bjets = 0;
      double bPt = 0.;
      foreach(const Jet& j, jets) {
        if( (j.abseta() < 2.1) && (deltaR(lepton, j) > 0.5) ) {
          // this jet passes the selection!
          ++goodjets;
          // j.bTagged() uses ghost association which is
          // more elegant, but not what has been used in
          // this analysis originally, will match B had-
          // rons in eta-phi space instead
          foreach(const Particle& b, bHadrons) {
            if( deltaR(j, b) < 0.3 ) {
              // jet matched to B hadron!
              if(!bPt)  bPt = j.pT() * GeV; // leading b-jet pT
              ++bjets; // count number of b-jets
              break;
            }
          }
        }
      }
      if( goodjets > 2 )  vetoEvent; // at most two jets
      if( !bjets )  vetoEvent; // at least one of them b-tagged

      double njets = double(goodjets);
      double ncomb = 3.0;
      _njet->fill(njets, weight);
      _njet->fill(ncomb, weight);

      if(     goodjets == 1)  _jet1_bPt->fill(bPt, weight);
      else if(goodjets == 2)  _jet2_bPt->fill(bPt, weight);
    }


    void finalize() {

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double sumw(sumOfWeights());
      MSG_INFO("Cross-Section/pb: " << xs_pb      );
      MSG_INFO("Sum of weights  : " << sumw       );
      MSG_INFO("nEvents         : " << numEvents());

      const double sf(xs_pb / sumw);

      scale(_njet,     sf);
      scale(_jet1_bPt, sf);
      scale(_jet2_bPt, sf);
    }

  protected:

    size_t _mode;

  private:

    Histo1DPtr _njet;
    Histo1DPtr _jet1_bPt;
    Histo1DPtr _jet2_bPt;

    //bool _isMuon;

  };

  class ATLAS_2013_I1219109_EL : public ATLAS_2013_I1219109 {
  public:
    ATLAS_2013_I1219109_EL()
      : ATLAS_2013_I1219109("ATLAS_2013_I1219109_EL")
    {
      _mode = 2;
    }
  };

  class ATLAS_2013_I1219109_MU : public ATLAS_2013_I1219109 {
  public:
    ATLAS_2013_I1219109_MU()
      : ATLAS_2013_I1219109("ATLAS_2013_I1219109_MU")
    {
      _mode = 3;
    }
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1219109);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1219109_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1219109_MU);

}
