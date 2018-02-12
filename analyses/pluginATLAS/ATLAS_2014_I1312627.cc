// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Measurement of V+jets distributions, taken as ratios between W and Z channels
  class ATLAS_2014_I1312627 : public Analysis {
  public:

    /// @name Plotting helper types
    //@{

    struct Plots {
      string ref;
      Histo1DPtr comp[2]; // (de)nominator components
      Scatter2DPtr ratio; // Rjets plot
    };

    typedef map<string, Plots> PlotMap;
    typedef PlotMap::value_type PlotMapPair;

    //@}


    /// Constructor
    ATLAS_2014_I1312627(std::string name="ATLAS_2014_I1312627")
      : Analysis(name)
    {
      _mode = 0; // using electron channel for combined data by default
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Set up cuts
      Cut cuts;
      if (_mode == 2) { // muon channel
        cuts = (Cuts::pT > 25.0*GeV) & Cuts::etaIn(-2.4, 2.4);
      } else if (_mode) { // electron channel
        cuts = (Cuts::pT > 25.0*GeV) & ( Cuts::etaIn(-2.47, -1.52) | Cuts::etaIn(-1.37, 1.37) | Cuts::etaIn(1.52, 2.47) );
      } else { // combined data extrapolated to common phase space
        cuts = (Cuts::pT > 25.0*GeV) & Cuts::etaIn(-2.5, 2.5);
      }

      // Boson finders
      FinalState fs;
      WFinder wfinder(fs, cuts, _mode > 1? PID::MUON : PID::ELECTRON, 40.0*GeV, MAXDOUBLE, 0.0*GeV, 0.1,
                      WFinder::CLUSTERNODECAY, WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(wfinder, "WF");

      ZFinder zfinder(fs, cuts, _mode > 1? PID::MUON : PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder, "ZF");

      // Jets
      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("WF"));
      jet_fs.addVetoOnThisFinalState(getProjection<ZFinder>("ZF"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "Jets");


      // Book Rjets plots
      _suff = string("-y0") + to_str(_mode + 1);
      hInit("Njets_incl",  "d01"); // inclusive number of jets
      hInit("Njets_excl",  "d04"); // exclusive number of jets
      hInit("Pt1_N1incl",  "d05"); // leading jet pT, at least 1 jet
      hInit("Pt1_N1excl",  "d06"); // leading jet pT, exactly 1 jet
      hInit("Pt1_N2incl",  "d07"); // leading jet pT, at least 2 jets
      hInit("Pt1_N3incl",  "d08"); // leading jet pT, at least 3 jets
      hInit("Pt2_N2incl",  "d09"); // subleading jet pT, at least 2 jets
      hInit("Pt3_N3incl",  "d10"); // sub-subleading jet pT, at least 3 jets
      hInit("ST_N2incl",   "d11"); // scalar jet pT sum, at least 2 jets
      hInit("ST_N2excl",   "d12"); // scalar jet pT sum, exactly 2 jets
      hInit("ST_N3incl",   "d13"); // scalar jet pT sum, at least 3 jets
      hInit("ST_N3excl",   "d14"); // scalar jet pT sum, exactly 3 jets
      hInit("DR_N2incl",   "d15"); // deltaR(j1, j2), at least 2 jets
      hInit("DPhi_N2incl", "d16"); // deltaPhi(j1, j2), at least 2 jets
      hInit("Mjj_N2incl",  "d17"); // mjj, at least 2 jets
      hInit("Rap1_N1incl", "d18"); // leading jet rapidity, at least 1 jet
      hInit("Rap2_N2incl", "d19"); // subleading jet rapidity, at least 2 jets
      hInit("Rap3_N3incl", "d20"); // sub-subleading jet rapidity, at least 3 jets

      // Also book numerator and denominator for Rjets plots
      foreach (PlotMapPair& str_plot, _plots) {
        str_plot.second.comp[0] = bookHisto1D( str_plot.second.ref + "2" + _suff, *(str_plot.second.ratio) );
        str_plot.second.comp[1] = bookHisto1D( str_plot.second.ref + "3" + _suff, *(str_plot.second.ratio) );
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Retrieve boson candidate
      const WFinder& wf = apply<WFinder>(event, "WF");
      const ZFinder& zf = apply<ZFinder>(event, "ZF");
      if (wf.empty() && zf.empty())  vetoEvent;

      // Retrieve jets
      const JetAlg& jetfs = apply<JetAlg>(event, "Jets");
      Jets all_jets = jetfs.jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);

      // Apply boson cuts and fill histograms
      const double weight = event.weight();
      if (zf.size() == 2) {
        const Particles& leptons = zf.constituents();
        if (oppSign(leptons[0], leptons[1]) && deltaR(leptons[0], leptons[1]) > 0.2)
          fillPlots(leptons, all_jets, weight, 1);
      }
      if (!wf.empty()) {
        const Particles& leptons = wf.constituentLeptons();
        if (wf.constituentNeutrino().pT() > 25*GeV && wf.mT() > 40*GeV )
          fillPlots(leptons, all_jets, weight, 0);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      ///  Normalise, scale and otherwise manipulate histograms here
      const double sf( crossSection() / sumOfWeights() );
      foreach (PlotMapPair& str_plot, _plots) {
        scale(str_plot.second.comp[0], sf);
        scale(str_plot.second.comp[1], sf);
        divide(str_plot.second.comp[0], str_plot.second.comp[1], str_plot.second.ratio);
      }
    }
    //@}


    /// Analysis helper functions
    //@{

    /// Histogram filling function, to avoid duplication
    void fillPlots(const Particles& leptons, Jets& all_jets, const double& weight, int isZ) {
      // Do jet-lepton overlap removal
      Jets jets;
      foreach (const Jet& j, all_jets) {
        /// @todo A nice place for a lambda and any() logic when C++11 is available
        bool keep = true;
        foreach (const Particle& l, leptons) keep &= deltaR(j, l) > 0.5;
        if (keep)  jets.push_back(j);
      }

      // Calculate jet ST
      const size_t njets = jets.size();
      double ST = 0.0; // scalar pT sum of all selected jets
      for (size_t i = 0; i < njets; ++i)
        ST += jets[i].pT() * GeV;

      // Fill jet histos
      _plots["Njets_excl"].comp[isZ]->fill(njets + 0.5, weight);
      for (size_t i = 0; i <= njets; ++i)
        _plots["Njets_incl"].comp[isZ]->fill(i + 0.5, weight);

      if (njets > 0) {
        const double pT1  = jets[0].pT()/GeV;
        const double rap1 = jets[0].absrap();
        _plots["Pt1_N1incl" ].comp[isZ]->fill(pT1,  weight);
        _plots["Rap1_N1incl"].comp[isZ]->fill(rap1, weight);

        if (njets == 1) {
          _plots["Pt1_N1excl"].comp[isZ]->fill(pT1, weight);
        } else if (njets > 1) {
          const double pT2  = jets[1].pT()/GeV;
          const double rap2 = jets[1].absrap();
          const double dR   = deltaR(jets[0], jets[1]);
          const double dPhi = deltaPhi(jets[0], jets[1]);
          const double mjj  = (jets[0].momentum() + jets[1].momentum()).mass()/GeV;
          _plots["Pt1_N2incl" ].comp[isZ]->fill(pT1,  weight);
          _plots["Pt2_N2incl" ].comp[isZ]->fill(pT2,  weight);
          _plots["Rap2_N2incl"].comp[isZ]->fill(rap2, weight);
          _plots["DR_N2incl"  ].comp[isZ]->fill(dR,   weight);
          _plots["DPhi_N2incl"].comp[isZ]->fill(dPhi, weight);
          _plots["Mjj_N2incl" ].comp[isZ]->fill(mjj,  weight);
          _plots["ST_N2incl"  ].comp[isZ]->fill(ST,   weight);

          if (njets == 2) {
            _plots["ST_N2excl"].comp[isZ]->fill(ST, weight);
          } else if (njets > 2) {
            const double pT3  = jets[2].pT()/GeV;
            const double rap3 = jets[2].absrap();
            _plots["Pt1_N3incl" ].comp[isZ]->fill(pT1,  weight);
            _plots["Pt3_N3incl" ].comp[isZ]->fill(pT3,  weight);
            _plots["Rap3_N3incl"].comp[isZ]->fill(rap3, weight);
            _plots["ST_N3incl"  ].comp[isZ]->fill(ST,   weight);

            if (njets == 3)
              _plots["ST_N3excl"].comp[isZ]->fill(ST, weight);
          }
        }
      }
    }


    /// Helper for histogram initialisation
    void hInit(string label, string ident) {
      string pre = ident + "-x0";
      _plots[label].ref = pre;
      _plots[label].ratio = bookScatter2D(pre + "1" + _suff, true);
    }

    //@}


  protected:

    // Data members
    size_t _mode;
    string _suff;


  private:

    /// @name Map of histograms
    PlotMap _plots;

  };



  /// Electron-specific version of the ATLAS_2014_I1312627 R-jets analysis
  class ATLAS_2014_I1312627_EL : public ATLAS_2014_I1312627 {
  public:
    ATLAS_2014_I1312627_EL()
      : ATLAS_2014_I1312627("ATLAS_2014_I1312627_EL")
    { _mode = 1; }
  };


  /// Muon-specific version of the ATLAS_2014_I1312627 R-jets analysis
  class ATLAS_2014_I1312627_MU : public ATLAS_2014_I1312627 {
  public:
    ATLAS_2014_I1312627_MU()
      : ATLAS_2014_I1312627("ATLAS_2014_I1312627_MU")
    { _mode = 2; }
  };


  // Hooks for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1312627);
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1312627_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1312627_MU);


}
