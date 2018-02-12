// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class ATLAS_2012_I1091481 : public Analysis {
  public:

    /// Constructor
    ATLAS_2012_I1091481()
      : Analysis("ATLAS_2012_I1091481")
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      ChargedFinalState cfs100(Cuts::abseta < 2.5 && Cuts::pT > 0.1*GeV);
      declare(cfs100,"CFS100");
      ChargedFinalState cfs500(Cuts::abseta < 2.5 && Cuts::pT > 0.5*GeV);
      declare(cfs500,"CFS500");

      // collision energy
      int isqrts = -1;
      if (fuzzyEquals(sqrtS(), 900*GeV)) isqrts = 2;
      if (fuzzyEquals(sqrtS(),   7*TeV)) isqrts = 1;
      assert(isqrts > 0);

      _sE_10_100   = bookHisto1D(isqrts, 1, 1);
      _sE_1_100    = bookHisto1D(isqrts, 1, 2);
      _sE_10_500   = bookHisto1D(isqrts, 1, 3);

      _sEta_10_100 = bookHisto1D(isqrts, 2, 1);
      _sEta_1_100  = bookHisto1D(isqrts, 2, 2);
      _sEta_10_500 = bookHisto1D(isqrts, 2, 3);

      norm_inclusive = 0.;
      norm_lowPt = 0.;
      norm_pt500 = 0.;
    }


    // Recalculate particle energy assuming pion mass
    double getPionEnergy(const Particle& p) {
      double m_pi = 0.1396*GeV;
      double p2 = p.p3().mod2()/(GeV*GeV);
      return sqrt(sqr(m_pi) + p2);
    }


    // S_eta core for one event
    //
    //  -1 + 1/Nch * |sum_j^Nch exp[i*(xi eta_j - Phi_j)]|^2
    //
    double getSeta(const Particles& part, double xi) {
      std::complex<double> c_eta (0.0, 0.0);
      foreach (const Particle& p, part) {
        double eta = p.eta();
        double phi = p.phi();
        double arg = xi*eta-phi;
         std::complex<double> temp(cos(arg), sin(arg));
         c_eta += temp;
      }
      return std::norm(c_eta)/part.size() - 1.0;
    }


    // S_E core for one event
    //
    //  -1 + 1/Nch * |sum_j^Nch exp[i*(omega X_j - Phi_j)]|^2
    //
    double getSE(const Particles& part, double omega) {
      double Xj = 0.0;
      std::complex<double> c_E (0.0, 0.0);
      for (unsigned int i=0; i < part.size(); ++i) {
        Xj += 0.5*getPionEnergy(part[i]);
        double phi = part[i].phi();
        double arg = omega*Xj - phi;
        std::complex<double> temp(cos(arg), sin(arg));
        c_E += temp;
        Xj += 0.5*getPionEnergy(part[i]);
      }
      return std::norm(c_E)/part.size() - 1.0;
    }


    // Convenient fill function
    void fillS(Histo1DPtr h, const Particles& part, double weight, bool SE=true) {
      // Loop over bins, take bin centers as parameter values
      for(size_t i=0; i < h->numBins(); ++i) {
        double x = h->bin(i).xMid();
        double width = h->bin(i).xMax() - h->bin(i).xMin();
        double y;
        if(SE)  y = getSE(part,   x);
        else    y = getSeta(part, x);
        h->fill(x, y * width * weight);
        // Histo1D objects will be converted to Scatter2D objects for plotting
        // As part of this conversion, Rivet will divide by bin width
        // However, we want the (x,y) of the Scatter2D to be the (binCenter, sumW) of
        // the current Histo1D. This is why in the above line we multiply by bin width,
        // so as to undo later division by bin width.
        //
        // Could have used Scatter2D objects in the first place, but they cannot be merged
        // as easily as Histo1Ds can using yodamerge (missing ScaledBy attribute)
        }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();

      // Charged fs
      const ChargedFinalState& cfs100  = apply<ChargedFinalState>(event, "CFS100");
      const Particles          part100 = cfs100.particles(cmpMomByEta);
      const ChargedFinalState& cfs500  = apply<ChargedFinalState>(event, "CFS500");
      const Particles&         part500 = cfs500.particles(cmpMomByEta);

      // Veto event if the most inclusive phase space has less than 10 particles and the max pT is > 10 GeV
      if (part100.size() < 11) vetoEvent;
      double ptmax = cfs100.particlesByPt()[0].pT()/GeV;
      if (ptmax > 10.0) vetoEvent;

      // Fill the pt>100, pTmax<10 GeV histos
      fillS(_sE_10_100, part100, weight, true);
      fillS(_sEta_10_100, part100, weight, false);
      norm_inclusive += weight;

      // Fill the pt>100, pTmax<1 GeV histos
      if (ptmax < 1.0) {
        fillS(_sE_1_100,   part100, weight, true);
        fillS(_sEta_1_100, part100, weight, false);
        norm_lowPt += weight;
      }

      // Fill the pt>500, pTmax<10 GeV histos
      if (part500.size() > 10) {
        fillS(_sE_10_500,   part500, weight, true );
        fillS(_sEta_10_500, part500, weight, false);
        norm_pt500 += weight;
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // The scaling takes the multiple fills per event into account
      scale(_sE_10_100, 1.0/norm_inclusive);
      scale(_sE_1_100 , 1.0/norm_lowPt);
      scale(_sE_10_500, 1.0/norm_pt500);

      scale(_sEta_10_100, 1.0/norm_inclusive);
      scale(_sEta_1_100 , 1.0/norm_lowPt);
      scale(_sEta_10_500, 1.0/norm_pt500);
    }

    //@}


  private:

    Histo1DPtr _sE_10_100;
    Histo1DPtr _sE_1_100;
    Histo1DPtr _sE_10_500;

    Histo1DPtr _sEta_10_100;
    Histo1DPtr _sEta_1_100;
    Histo1DPtr _sEta_10_500;

    double norm_inclusive;
    double norm_lowPt;
    double norm_pt500;
  };


  DECLARE_RIVET_PLUGIN(ATLAS_2012_I1091481);

}
