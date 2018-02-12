// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  class ATLAS_2013_I1217863_W : public Analysis {
  public:

    /// Constructor
    ATLAS_2013_I1217863_W(string name="ATLAS_2013_I1217863_W")
      : Analysis(name)
    {
      // the electron mode is used by default
      _mode = 1;
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      declare(fs, "FS");

      Cut cuts = Cuts::abseta < 2.47 && Cuts::pT > 25*GeV;

      // W finder for electrons and muons
      WFinder wf(fs, cuts, _mode==3? PID::MUON : PID::ELECTRON, 0.0*GeV, 1000.0*GeV, 35.0*GeV, 0.1,
                                     WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(wf, "WF");

      // leading photon
      LeadingParticlesFinalState photonfs(FinalState(Cuts::abseta < 2.37 && Cuts::pT > 15*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("WF"));
      jet_fs.addVetoOnThisFinalState(getProjection<LeadingParticlesFinalState>("LeadingPhoton"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "Jets");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "isolatedFS");


      // Book histograms
      _hist_EgammaT_incl   = bookHisto1D( 7, 1, _mode); // dSigma / dE^gamma_T for Njet >= 0
      _hist_EgammaT_excl   = bookHisto1D( 8, 1, _mode); // dSigma / dE^gamma_T for Njet = 0
      _hist_Njet_EgammaT15 = bookHisto1D(15, 1, _mode); // dSigma / dNjet for E^gamma_T > 15
      _hist_Njet_EgammaT60 = bookHisto1D(16, 1, _mode); // dSigma / dNjet for E^gamma_T > 60
      _hist_mWgammaT       = bookHisto1D(19, 1, _mode); // dSigma / dm^{Wgamma}_T

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // retrieve leading photon
      Particles photons = apply<LeadingParticlesFinalState>(event, "LeadingPhoton").particles();
      if (photons.size() != 1)  vetoEvent;
      const Particle& leadingPhoton = photons[0];
      if (leadingPhoton.Et() < 15.0*GeV) vetoEvent;
      if (leadingPhoton.abseta() > 2.37) vetoEvent;

      // check photon isolation
      double coneEnergy(0.0);
      Particles fs = apply<VetoedFinalState>(event, "isolatedFS").particles();
      for (const Particle& p : fs) {
        if ( deltaR(leadingPhoton, p) < 0.4 )  coneEnergy += p.E();
      }
      if ( coneEnergy / leadingPhoton.E() >= 0.5 )  vetoEvent;

      // retrieve W boson candidate
      const WFinder& wf = apply<WFinder>(event, "WF");
      if ( wf.bosons().size() != 1 )  vetoEvent; // only one W boson candidate
      //const Particle& Wboson  = wf.boson();

      // retrieve constituent neutrino
      const Particle& neutrino = wf.constituentNeutrino();
      if ( !(neutrino.pT() > 35.0*GeV) )  vetoEvent;

      // retrieve constituent lepton
      const Particle& lepton = wf.constituentLepton();
      if ( !(lepton.pT() > 25.0*GeV && lepton.abseta() < 2.47) )  vetoEvent;

      // check photon-lepton overlap
      if ( !(deltaR(leadingPhoton, lepton) > 0.7) )  vetoEvent;

      // count jets
      const FastJets& jetfs = apply<FastJets>(event, "Jets");
      Jets jets = jetfs.jets(cmpMomByEt);
      int goodJets = 0;
      for (const Jet& j : jets) {
        if ( !(j.Et() > 30.0*GeV) )  break;
        if ( (j.abseta() < 4.4) && \
             (deltaR(leadingPhoton, j) > 0.3) &&            \
             (deltaR(lepton,        j) > 0.3) )  ++goodJets;
      }

      double Njets = double(goodJets) + 0.5;
      double photonEt = leadingPhoton.Et()*GeV;

      const FourMomentum& lep_gamma = lepton.momentum() + leadingPhoton.momentum();
      double term1 = sqrt(lep_gamma.mass2() + lep_gamma.pT2()) + neutrino.Et();
      double term2 = (lep_gamma + neutrino.momentum()).pT2();
      double mWgammaT = sqrt(term1 * term1 - term2) * GeV;

      _hist_EgammaT_incl->fill(photonEt, weight);

      _hist_Njet_EgammaT15->fill(Njets, weight);

      if ( !goodJets )  _hist_EgammaT_excl->fill(photonEt, weight);

      if (photonEt > 40.0*GeV) {
        _hist_mWgammaT->fill(mWgammaT, weight);
        if (photonEt > 60.0*GeV)  _hist_Njet_EgammaT60->fill(Njets, weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xs_fb = crossSection()/femtobarn;
      const double sumw = sumOfWeights();
      const double sf = xs_fb / sumw;

      scale(_hist_EgammaT_excl, sf);
      scale(_hist_EgammaT_incl, sf);

      normalize(_hist_Njet_EgammaT15);
      normalize(_hist_Njet_EgammaT60);
      normalize(_hist_mWgammaT);

    }

    //@}

  protected:

    size_t _mode;

  private:

    /// @name Histograms
    //@{

    Histo1DPtr _hist_EgammaT_incl;
    Histo1DPtr _hist_EgammaT_excl;
    Histo1DPtr _hist_Njet_EgammaT15;
    Histo1DPtr _hist_Njet_EgammaT60;
    Histo1DPtr _hist_mWgammaT;

    //@}

  };


  class ATLAS_2013_I1217863_W_EL : public ATLAS_2013_I1217863_W {
  public:
    ATLAS_2013_I1217863_W_EL()
      : ATLAS_2013_I1217863_W("ATLAS_2013_I1217863_W_EL")
    {
      _mode = 2;
    }
  };


  class ATLAS_2013_I1217863_W_MU : public ATLAS_2013_I1217863_W {
  public:
    ATLAS_2013_I1217863_W_MU()
      : ATLAS_2013_I1217863_W("ATLAS_2013_I1217863_W_MU")
    {
      _mode = 3;
    }
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1217863_W);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1217863_W_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1217863_W_MU);

}
