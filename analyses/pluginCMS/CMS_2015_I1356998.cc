// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  class CMS_2015_I1356998 : public Analysis {
  public:

    CMS_2015_I1356998()
      : Analysis("CMS_2015_I1356998"), edge(4.7)
    {    }


    void init() {

      declare(FinalState(),"FS");

      _h_noCASTORtag = bookHisto1D(1, 1, 1);
      _h_CASTORtag   = bookHisto1D(2, 1, 1);
      _h_centralGap  = bookHisto1D(3, 1, 1);
      _h_sigmaVis    = bookHisto1D(4, 1, 1);
      _h_maxFwdGap   = bookHisto1D(5, 1, 1);

    }


    void analyze(const Event& event) {

      const double weight = event.weight();
      const FinalState& fs = apply<FinalState>(event, "FS");

      // A vector containing a lot of eta values
      vector<double> detparticles;
      detparticles.push_back(-edge);
      foreach (const Particle& p, fs.particles(Cuts::pT > 0.2*GeV && Cuts::abseta<edge, cmpMomByEta) ) {
        detparticles.push_back(p.momentum().eta());
      }
      detparticles.push_back(edge);

      // Find maximum gap size
      vector <double>::iterator iter;
      vector<double> detgaps;
      for (iter = detparticles.begin()+1; iter != detparticles.end(); ++iter) {
        const double detgap = *iter - *(iter-1);
        detgaps.push_back(detgap);
      }
      double detgapbwd = detgaps.front();
      double detgapfwd = detgaps.back();
      double detfmax = max(detgapbwd, detgapfwd);

      // Fill rapidity gap histo
      if (detfmax != 2*edge ) {
        _h_maxFwdGap->fill(detfmax, weight);
      }
      // Everything that follows has to do with the cross-section measurements

      if (fs.size() < 2) vetoEvent;

      // Gap center calculations
      const ParticleVector particlesByRapidity = fs.particles(cmpMomByRap); //ByRapidity();

      vector<double> gaps;
      vector<double> midpoints;
      for (size_t ip = 1; ip < particlesByRapidity.size(); ++ip) {
        const Particle& p1 = particlesByRapidity[ip-1];
        const Particle& p2 = particlesByRapidity[ip];
        const double gap = p2.momentum().rapidity()  - p1.momentum().rapidity();
        const double mid = (p2.momentum().rapidity() + p1.momentum().rapidity()) / 2.;
        gaps.push_back(gap);
        midpoints.push_back(mid);
      }

      int imid = std::distance(gaps.begin(), max_element(gaps.begin(), gaps.end()));
      double gapcenter = midpoints[imid];

      // Calculations for cross-sections
      FourMomentum MxFourVector(0.,0.,0.,0.);
      FourMomentum MyFourVector(0.,0.,0.,0.);

      foreach(const Particle& p, fs.particles(cmpMomByEta)) {
        if (p.momentum().rapidity() > gapcenter) {
          MxFourVector += p.momentum();
        }
        else {
          MyFourVector += p.momentum();
        }
      }

      double Mx = MxFourVector.mass();
      double My = MyFourVector.mass();

      const double xix = (Mx*Mx)/(sqrtS()/GeV * sqrtS()/GeV);

      if (log10(My) < 0.5) {
        _h_noCASTORtag->fill(log10(xix), weight);
        if (log10(xix) > -5.5 && log10(xix) < -2.5) _h_sigmaVis->fill(0.5, weight);
      }
      else if (log10(My) < 1.1) {
        _h_CASTORtag->fill(log10(xix), weight);
        if (log10(xix) > -5.5 && log10(xix) < -2.5) _h_sigmaVis->fill(1.5, weight);
      }

      // Central gap x-section
      double xigen = (Mx*Mx) * (My*My) / (sqrtS()/GeV * sqrtS()/GeV * 0.93827 * 0.93827); // Proton masses...
      double dy0 = -log(xigen);

      if (dy0 > 3.) {
        if (log10(My) > 1.1 && log10(Mx) > 1.1) {
          _h_centralGap->fill(dy0, weight);
          _h_sigmaVis->fill(2.5, weight);
        }
      }

    }

    void finalize() {

      double xs = crossSection()/millibarn/sumOfWeights();
      scale(_h_noCASTORtag, xs);
      scale(_h_CASTORtag  , xs);
      scale(_h_centralGap , xs);
      scale(_h_sigmaVis   , xs);
      scale(_h_maxFwdGap  , xs);

    }

  private:

    Histo1DPtr _h_noCASTORtag;
    Histo1DPtr _h_CASTORtag;
    Histo1DPtr _h_centralGap;
    Histo1DPtr _h_sigmaVis;
    Histo1DPtr _h_maxFwdGap;
    double edge;

  };


  DECLARE_RIVET_PLUGIN(CMS_2015_I1356998);

}
