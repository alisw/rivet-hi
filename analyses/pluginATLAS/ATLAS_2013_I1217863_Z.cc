// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"

namespace Rivet {


  class ATLAS_2013_I1217863_Z : public Analysis {
  public:

    /// Constructor
    ATLAS_2013_I1217863_Z(string name="ATLAS_2013_I1217863_Z")
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

      // Z finder
      ZFinder zf(fs, cuts, _mode==3? PID::MUON : PID::ELECTRON, 40.0*GeV, 1000.0*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zf, "ZF");

      // leading photon
      LeadingParticlesFinalState photonfs(FinalState(Cuts::abseta < 2.37 && Cuts::pT > 15*GeV));
      photonfs.addParticleId(PID::PHOTON);
      declare(photonfs, "LeadingPhoton");

      // jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<ZFinder>("ZF"));
      jet_fs.addVetoOnThisFinalState(getProjection<LeadingParticlesFinalState>("LeadingPhoton"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "Jets");

      // FS excluding the leading photon
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(photonfs);
      declare(vfs, "isolatedFS");


      // Book histograms
      _hist_EgammaT_incl   = bookHisto1D(11, 1, _mode); // dSigma / dE^gamma_T for Njet >= 0
      _hist_EgammaT_excl   = bookHisto1D(12, 1, _mode); // dSigma / dE^gamma_T for Njet = 0
      _hist_Njet_EgammaT15 = bookHisto1D(17, 1, _mode); // dSigma / dNjet for E^gamma_T >= 15
      _hist_Njet_EgammaT60 = bookHisto1D(18, 1, _mode); // dSigma / dNjet for E^gamma_T >= 60
      _hist_mZgamma        = bookHisto1D(20, 1, _mode); // dSigma / dm^{Zgamma}

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
      if (coneEnergy / leadingPhoton.E() >= 0.5 )  vetoEvent;

      // retrieve Z boson candidate
      const ZFinder& zf = apply<ZFinder>(event, "ZF");
      if ( zf.bosons().size() != 1 )  vetoEvent; // only one Z boson candidate
      const Particle& Zboson  = zf.boson();
      if ( !(Zboson.mass() > 40.0*GeV) )  vetoEvent;

      // check charge of constituent leptons
      const ParticleVector& leptons = zf.constituents();
      if (leptons.size() != 2 || leptons[0].charge() * leptons[1].charge() > 0.)  vetoEvent;

      // check photon-lepton overlap
      for (const Particle& p : leptons) {
        if ( !(p.pT() > 25.0*GeV && p.abseta() < 2.47 && deltaR(leadingPhoton, p) > 0.7) )  vetoEvent;
      }

      // count jets
      const FastJets& jetfs = apply<FastJets>(event, "Jets");
      Jets jets = jetfs.jets(cmpMomByEt);
      int goodJets = 0;
      foreach (const Jet& j, jets) {
        if ( !(j.Et() > 30.0*GeV) )  break;
        if ( (j.abseta() < 4.4) && \
             (deltaR(leadingPhoton, j) > 0.3) &&    \
             (deltaR(leptons[0],    j) > 0.3) &&            \
             (deltaR(leptons[1],    j) > 0.3) )  ++goodJets;
      }

      double Njets = double(goodJets) + 0.5;
      double photonEt = leadingPhoton.Et()*GeV;
      double mZgamma = (Zboson.momentum() + leadingPhoton.momentum()).mass() * GeV;

      _hist_EgammaT_incl->fill(photonEt, weight);

      _hist_Njet_EgammaT15->fill(Njets, weight);

      if ( !goodJets )   _hist_EgammaT_excl->fill(photonEt, weight);

      if (photonEt >= 40.0*GeV) {
        _hist_mZgamma->fill(mZgamma, weight);
        if (photonEt >= 60.0*GeV)  _hist_Njet_EgammaT60->fill(Njets, weight);
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
      normalize(_hist_mZgamma);

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
    Histo1DPtr _hist_mZgamma;

    //@}

  };


  class ATLAS_2013_I1217863_Z_EL : public ATLAS_2013_I1217863_Z {
  public:
    ATLAS_2013_I1217863_Z_EL()
      : ATLAS_2013_I1217863_Z("ATLAS_2013_I1217863_Z_EL")
    {
      _mode = 2;
    }
  };


  class ATLAS_2013_I1217863_Z_MU : public ATLAS_2013_I1217863_Z {
  public:
    ATLAS_2013_I1217863_Z_MU()
      : ATLAS_2013_I1217863_Z("ATLAS_2013_I1217863_Z_MU")
    {
      _mode = 3;
    }
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1217863_Z);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1217863_Z_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1217863_Z_MU);

}
