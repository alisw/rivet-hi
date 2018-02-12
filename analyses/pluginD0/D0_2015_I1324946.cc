// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {

  class D0_2015_I1324946 : public Analysis {
  public:

    /// Constructor
    D0_2015_I1324946()
      : Analysis("D0_2015_I1324946")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs;
      ZFinder zfinder_mm(fs, Cuts::abseta < 2 && Cuts::pT > 15*GeV, PID::MUON, 30*GeV, 500*GeV, 0.0, ZFinder::NOCLUSTER, ZFinder::NOTRACK);
      declare(zfinder_mm, "zfinder_mm");

      _h_phistar_mm_peak_central = bookHisto1D(1, 1, 1);
      _h_phistar_mm_peak_forward = bookHisto1D(1, 1, 2);
      _h_phistar_mm_low_central = bookHisto1D(2, 1, 1);
      _h_phistar_mm_low_forward = bookHisto1D(2, 1, 2);
      _h_phistar_mm_high1 = bookHisto1D(3, 1, 1);
      _h_phistar_mm_high2 = bookHisto1D(4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();


      //70<Mmm<105
      const ZFinder& zfinder_mm = apply<ZFinder>(event, "zfinder_mm");
      if (zfinder_mm.bosons().size() == 1) {
	Particles mm = zfinder_mm.constituents();
	std::sort(mm.begin(), mm.end(), cmpMomByPt);
	const FourMomentum& mminus = PID::threeCharge(mm[0].pid()) < 0 ? mm[0].momentum() : mm[1].momentum();
	const FourMomentum& mplus  = PID::threeCharge(mm[0].pid()) < 0 ? mm[1].momentum() : mm[0].momentum();
	double phi_acop = M_PI - mapAngle0ToPi(mminus.phi() - mplus.phi());
	double costhetastar = tanh((mminus.eta() - mplus.eta())/2);
	double sin2thetastar = 1 - sqr(costhetastar);
	if (sin2thetastar < 0) sin2thetastar = 0;
	const double phistar = tan(phi_acop/2) * sqrt(sin2thetastar);
	const FourMomentum& zmom = zfinder_mm.bosons()[0].momentum();
        if (zmom.mass()<30*GeV || zmom.mass() >500*GeV) vetoEvent;
	
        if( zmom.mass()>70 && zmom.mass()<100 && zmom.absrap()<1.0) _h_phistar_mm_peak_central->fill(phistar, weight);
	if( zmom.mass()>70 && zmom.mass()<100 && zmom.absrap()>1.0  && zmom.absrap()<2.0) _h_phistar_mm_peak_forward->fill(phistar, weight);
	if( zmom.mass()>30 && zmom.mass()<60  && zmom.absrap()<1.0) _h_phistar_mm_low_central->fill(phistar, weight);
	if( zmom.mass()>30 && zmom.mass()<60 && zmom.absrap()>1.0 && zmom.absrap()<2.0) _h_phistar_mm_low_forward->fill(phistar, weight);
	if( zmom.mass()>160 && zmom.mass()<300) _h_phistar_mm_high1->fill(phistar, weight);
	if( zmom.mass()>300 && zmom.mass()<500) _h_phistar_mm_high2->fill(phistar, weight);
	
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(_h_phistar_mm_low_central);
      normalize(_h_phistar_mm_low_forward);
      normalize(_h_phistar_mm_peak_central);
      normalize(_h_phistar_mm_peak_forward);
      normalize(_h_phistar_mm_high1);
      normalize(_h_phistar_mm_high2);

    }


    //}

    //@}


  private:
    /// @name Histograms
    //@{

    Histo1DPtr _h_phistar_mm_low_central;
    Histo1DPtr _h_phistar_mm_low_forward;
    Histo1DPtr _h_phistar_mm_peak_central;
    Histo1DPtr _h_phistar_mm_peak_forward;
    Histo1DPtr _h_phistar_mm_high1;
    Histo1DPtr _h_phistar_mm_high2;
    //@}                                        
  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2015_I1324946);

}
