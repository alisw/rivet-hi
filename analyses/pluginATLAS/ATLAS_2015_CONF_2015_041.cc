// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// Z + jets in pp at 13 TeV 
  /// @note This base class contains a "mode" variable for combined, e, and mu channel derived classes
  class ATLAS_2015_CONF_2015_041 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2015_CONF_2015_041(string name="ATLAS_2015_CONF_2015_041")
      : Analysis(name),
        _weights(5, 0.0)
    {
      // This class uses the combined e+mu mode
      _mode = 0;
    }

    //@}


    /// Book histograms and initialise projections before the run
    void init() {
      const FinalState fs;

      Cut cuts = (Cuts::pT > 25*GeV) & (Cuts::abseta < 2.5);
      ZFinder zfinder(fs, cuts, _mode? PID::MUON : PID::ELECTRON, 66*GeV, 116*GeV);
      declare(zfinder, "zfinder");

      // Define veto FS in order to prevent Z-decay products entering the jet algorithm
      VetoedFinalState had_fs;
      had_fs.addVetoOnThisFinalState(zfinder);
      FastJets jets(had_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "jets");

      // individual channels
      _hNjets      = bookHisto1D(1, 1, _mode + 1);
      _hNjetsRatio = bookScatter2D(2, 1, _mode + 1, true);
      // combination
      _hNjets_comb      = bookHisto1D(1, 2, _mode + 1);
      _hNjetsRatio_comb = bookScatter2D(2, 2, _mode + 1, true);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      const ZFinder& zfinder = apply<ZFinder>(event, "zfinder");
      const Particles& leptons = zfinder.constituents();
      if (leptons.size() != 2)  vetoEvent;

      Jets jets;
      foreach (Jet j, apply<JetAlg>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 2.5)) {
        bool keep = true;
        foreach(const Particle& l, leptons)  keep &= deltaR(j, l) > 0.4;
        if (keep)  jets += j;
      }

      size_t njets = jets.size();

      for(size_t i = 0; i <= njets; ++i) {
        _hNjets->fill(i + 0.5, weight);
        _hNjets_comb->fill(i + 0.5, weight);
      }

      for (size_t i = 0; i < 5; ++i) {
        if (njets >= i) _weights[i] += weight;
      }

    }

    /// @name Ratio calculator util functions
    //@{

    /// Calculate the ratio, being careful about div-by-zero
    double ratio(double a, double b) {
      return (b != 0) ? a/b : 0;
    }

    /// Calculate the ratio error, being careful about div-by-zero
    double ratio_err(double a, double b) {
      return (b != 0) ? sqrt(a/b*(1-a/b)/b) : 0;
    }

    //@}

    void finalize() {
      for (size_t i = 0; i < 4; ++i) {
        double  n = _hNjets->bin(i + 1).sumW();
        double dN = _hNjets->bin(i + 1).sumW2();
        double  d = _hNjets->bin(i).sumW();
        double dD = _hNjets->bin(i).sumW2();
        double r = safediv(n, d);
        double e = sqrt( safediv(r * (1 - r), d) );
        if ( _hNjets->effNumEntries() != _hNjets->numEntries() ) {
          // use F. James's approximation for weighted events:
          e = sqrt( safediv((1 - 2 * r) * dN + r * r * dD, d * d) );
        }
        _hNjetsRatio->point(i).setY(r, e);
        _hNjetsRatio_comb->point(i).setY(r, e);
      }

      scale(_hNjets,      crossSectionPerEvent() );
      scale(_hNjets_comb, crossSectionPerEvent() );
    }

    //@}


  protected:

    size_t _mode;


  private:

    vector<double> _weights;
    Scatter2DPtr _hNjetsRatio, _hNjetsRatio_comb;
    Histo1DPtr _hNjets, _hNjets_comb;
  };



  class ATLAS_2015_CONF_2015_041_EL : public ATLAS_2015_CONF_2015_041 {
  public:
    ATLAS_2015_CONF_2015_041_EL()
      : ATLAS_2015_CONF_2015_041("ATLAS_2015_CONF_2015_041_EL")
    {
      _mode = 0;
    }
  };



  class ATLAS_2015_CONF_2015_041_MU : public ATLAS_2015_CONF_2015_041 {
  public:
    ATLAS_2015_CONF_2015_041_MU()
      : ATLAS_2015_CONF_2015_041("ATLAS_2015_CONF_2015_041_MU")
    {
      _mode = 1;
    }
  };



  DECLARE_RIVET_PLUGIN(ATLAS_2015_CONF_2015_041);
  DECLARE_RIVET_PLUGIN(ATLAS_2015_CONF_2015_041_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2015_CONF_2015_041_MU);
}
