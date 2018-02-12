#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class ATLAS_2015_I1408516 : public Analysis {
  public:

    /// Constructor
    ATLAS_2015_I1408516(string name="ATLAS_2015_I1408516", size_t mode=0)
      : Analysis(name), _mode(mode) // using electron channel for combined data
    { }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Configure projections
      FinalState fs;
      Cut cuts = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;
      ZFinder zfinder_el(fs, cuts, (_mode ? PID::MUON : PID::ELECTRON),
                         12*GeV, 150*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::NOTRACK);
      declare(zfinder_el, "ZFinder");

      // Book histograms
      const size_t offset = _mode ? 4 : 1;

      _h["phistar_lo_00_08"] = bookHisto1D( 2, 1, offset);
      _h["phistar_lo_08_16"] = bookHisto1D( 3, 1, offset);
      _h["phistar_lo_16_24"] = bookHisto1D( 4, 1, offset);

      _h["phistar_me_00_04"] = bookHisto1D( 5, 1, offset);
      _h["phistar_me_04_08"] = bookHisto1D( 6, 1, offset);
      _h["phistar_me_08_12"] = bookHisto1D( 7, 1, offset);
      _h["phistar_me_12_16"] = bookHisto1D( 8, 1, offset);
      _h["phistar_me_16_20"] = bookHisto1D( 9, 1, offset);
      _h["phistar_me_20_24"] = bookHisto1D(10, 1, offset);

      _h["phistar_hi_00_08"] = bookHisto1D(11, 1, offset);
      _h["phistar_hi_08_16"] = bookHisto1D(12, 1, offset);
      _h["phistar_hi_16_24"] = bookHisto1D(13, 1, offset);

      _h["phistar_mll_46_66"  ] = bookHisto1D(14, 1, offset);
      _h["phistar_mll_66_116" ] = bookHisto1D(15, 1, offset);
      _h["phistar_mll_116_150"] = bookHisto1D(16, 1, offset);

      _h["zpt_00_04"] = bookHisto1D(17, 1, offset);
      _h["zpt_04_08"] = bookHisto1D(18, 1, offset);
      _h["zpt_08_12"] = bookHisto1D(19, 1, offset);
      _h["zpt_12_16"] = bookHisto1D(20, 1, offset);
      _h["zpt_16_20"] = bookHisto1D(21, 1, offset);
      _h["zpt_20_24"] = bookHisto1D(22, 1, offset);

      _h["zpt_mll_12_20"  ] = bookHisto1D(23, 1, offset);
      _h["zpt_mll_20_30"  ] = bookHisto1D(24, 1, offset);
      _h["zpt_mll_30_46"  ] = bookHisto1D(25, 1, offset);
      _h["zpt_mll_46_66"  ] = bookHisto1D(26, 1, offset);
      _h["zpt_mll_66_116" ] = bookHisto1D(27, 1, offset);
      _h["zpt_mll_116_150"] = bookHisto1D(28, 1, offset);

      _h["zpt_00_04_xsec"] = bookHisto1D(29, 1, offset);
      _h["zpt_04_08_xsec"] = bookHisto1D(30, 1, offset);
      _h["zpt_08_12_xsec"] = bookHisto1D(31, 1, offset);
      _h["zpt_12_16_xsec"] = bookHisto1D(32, 1, offset);
      _h["zpt_16_20_xsec"] = bookHisto1D(33, 1, offset);
      _h["zpt_20_24_xsec"] = bookHisto1D(34, 1, offset);

      _h["zpt_mll_12_20_xsec"  ] = bookHisto1D(35, 1, offset);
      _h["zpt_mll_20_30_xsec"  ] = bookHisto1D(36, 1, offset);
      _h["zpt_mll_30_46_xsec"  ] = bookHisto1D(37, 1, offset);
      _h["zpt_mll_46_66_xsec"  ] = bookHisto1D(38, 1, offset);
      _h["zpt_mll_66_116_xsec" ] = bookHisto1D(39, 1, offset);
      _h["zpt_mll_116_150_xsec"] = bookHisto1D(40, 1, offset);

      _h["mll_xsec"] = bookHisto1D(41, 1, 1 + _mode);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get leptonic Z boson
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");
      if (zfinder.bosons().size() != 1 ) vetoEvent;
      const Particle& Zboson = zfinder.boson();

      // Get/cut on heavily used Z boson properties
      const double zpt   = Zboson.pT();
      const double zrap  = Zboson.absrap();
      const double zmass = Zboson.mass();
      if (zrap > 2.4) vetoEvent;

      // Get/cut on Z boson leptons
      const ParticleVector& leptons = zfinder.constituents();
      if (leptons.size() != 2 || leptons[0].threeCharge() * leptons[1].threeCharge() > 0) vetoEvent;
      const Particle& lminus = leptons[0].charge() < 0 ? leptons[0] : leptons[1];
      const Particle& lplus  = leptons[0].charge() < 0 ? leptons[1] : leptons[0];

      // Compute phi*
      const double phi_acop = M_PI - deltaPhi(lminus, lplus);
      const double costhetastar = tanh( 0.5 * (lminus.eta() - lplus.eta()) );
      const double sin2thetastar = (costhetastar > 1) ? 0.0 : (1.0 - sqr(costhetastar));
      const double phistar = tan(0.5 * phi_acop) * sqrt(sin2thetastar);

      // Event weight for histogramming
      const double weight = event.weight();

      // Inclusive mll
      if (zmass > 46*GeV || zpt > 45*GeV) {
        // 46 GeV < mll < 150 GeV OR (12 GeV < mll < 46 GeV AND ZpT >45 GeV)
        _h["mll_xsec"]->fill(zmass, weight);
      }

      // 12 GeV < mll < 150 GeV observables
      if (zmass < 20*GeV) {
        // 12 GeV < mll < 20 GeV
        if (zpt > 45*GeV) { // ZpT cut only for low-mass regions
          _h["zpt_mll_12_20_xsec"]->fill(zpt, weight);
          _h["zpt_mll_12_20"     ]->fill(zpt, weight);
        }
      } else if (zmass < 30*GeV) {
        // 20 GeV < mll < 30 GeV
        if (zpt > 45*GeV) { // ZpT cut only for low-mass regions
          _h["zpt_mll_20_30_xsec"]->fill(zpt, weight);
          _h["zpt_mll_20_30"     ]->fill(zpt, weight);
        }
      } else if (zmass <  46*GeV) {
        // 30 GeV < mll < 46 GeV
        if (zpt > 45*GeV) { // ZpT cut only for low-mass regions
          _h["zpt_mll_30_46_xsec"]->fill(zpt, weight);
          _h["zpt_mll_30_46"     ]->fill(zpt, weight);
        }
      } else if (zmass <  66*GeV) {
        // 46 GeV < mll < 66 GeV
        _h["zpt_mll_46_66_xsec"]->fill(zpt, weight);
        _h["zpt_mll_46_66"     ]->fill(zpt, weight);

        _h["phistar_mll_46_66"]->fill(phistar, weight);
        if      (zrap < 0.8)  _h["phistar_lo_00_08"]->fill(phistar, weight);
        else if (zrap < 1.6)  _h["phistar_lo_08_16"]->fill(phistar, weight);
        else                  _h["phistar_lo_16_24"]->fill(phistar, weight);

      } else if (zmass < 116*GeV) {
        // 66 GeV < mll < 116 GeV
        _h["zpt_mll_66_116_xsec"]->fill(zpt, weight);
        _h["zpt_mll_66_116"     ]->fill(zpt, weight);

        if (zrap < 0.4) {
          _h["zpt_00_04_xsec"]->fill(zpt, weight);
          _h["zpt_00_04"]->fill(zpt, weight);
        } else if (zrap < 0.8) {
          _h["zpt_04_08_xsec"]->fill(zpt, weight);
          _h["zpt_04_08"]->fill(zpt, weight);
        } else if (zrap < 1.2) {
          _h["zpt_08_12_xsec"]->fill(zpt, weight);
          _h["zpt_08_12"]->fill(zpt, weight);
        } else if (zrap < 1.6) {
          _h["zpt_12_16_xsec"]->fill(zpt, weight);
          _h["zpt_12_16"]->fill(zpt, weight);
        } else if (zrap < 2.0) {
          _h["zpt_16_20_xsec"]->fill(zpt, weight);
          _h["zpt_16_20"]->fill(zpt, weight);
        } else {
          _h["zpt_20_24_xsec"]->fill(zpt, weight);
          _h["zpt_20_24"]->fill(zpt, weight);
        }

        _h["phistar_mll_66_116"]->fill(phistar, weight);
        if      (zrap < 0.4)  _h["phistar_me_00_04"]->fill(phistar, weight);
        else if (zrap < 0.8)  _h["phistar_me_04_08"]->fill(phistar, weight);
        else if (zrap < 1.2)  _h["phistar_me_08_12"]->fill(phistar, weight);
        else if (zrap < 1.6)  _h["phistar_me_12_16"]->fill(phistar, weight);
        else if (zrap < 2.0)  _h["phistar_me_16_20"]->fill(phistar, weight);
        else                  _h["phistar_me_20_24"]->fill(phistar, weight);

      } else {

        // 116 GeV < mll < 150 GeV
        _h["zpt_mll_116_150_xsec"]->fill(zpt, weight);
        _h["zpt_mll_116_150"     ]->fill(zpt, weight);

        _h["phistar_mll_116_150"]->fill(phistar, weight);
        if      (zrap < 0.8)  _h["phistar_hi_00_08"]->fill(phistar, weight);
        else if (zrap < 1.6)  _h["phistar_hi_08_16"]->fill(phistar, weight);
        else                  _h["phistar_hi_16_24"]->fill(phistar, weight);

      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // Scale non-xsec plots to cross-section
      const double sf = crossSection() / picobarn / sumOfWeights();
      for (const auto& key_hist : _h) {
        scale(key_hist.second, sf);
        if (!contains(key_hist.first, "_xsec")) normalize(key_hist.second);
      }

      // M(ll) plot isn't a differential cross section so shouldn't be divided by bin width
      for (size_t i = 0; i < 6; ++i) {
        double bw = _h["mll_xsec"]->bin(i).xWidth();
        _h["mll_xsec"]->bin(i).scaleW(bw);
      }
    }
    //@}


  protected:

    size_t _mode;


  private:

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    //@}

  };



  class ATLAS_2015_I1408516_EL : public ATLAS_2015_I1408516 {
  public:
    ATLAS_2015_I1408516_EL() : ATLAS_2015_I1408516("ATLAS_2015_I1408516_EL", 0) { }
  };

  class ATLAS_2015_I1408516_MU : public ATLAS_2015_I1408516 {
  public:
    ATLAS_2015_I1408516_MU() : ATLAS_2015_I1408516("ATLAS_2015_I1408516_MU", 1) { }
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1408516);
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1408516_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1408516_MU);

}
