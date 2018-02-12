#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// Measurement of the inclusive energy spectrum in the very forward direction in pp collisions at 13 TeV
  class CMS_2017_I1511284 : public Analysis {
  public:

    /// Constructor
    CMS_2017_I1511284()
      : Analysis("CMS_2017_I1511284")
    {    }


    /// Book histograms and initialise projections
    void init() {
      declare(FinalState(), "FS");

      _h_totEnergy = bookHisto1D(1, 1, 1);
      _h_emEnergy = bookHisto1D(2, 1, 1);
      _h_hadEnergy = bookHisto1D(3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles fsparticles = apply<FinalState>(event, "FS").particles(cmpMomByRap);
      if (fsparticles.size() < 2) vetoEvent; // need at least two particles to calculate gaps

      double gapCenter = 0,largestGap = 0;
      double previousRapidity = fsparticles.front().rapidity();
      for (const Particle& p : fsparticles) {
        const double gap = fabs(p.rapidity() - previousRapidity);
        if (gap > largestGap) {
          largestGap = gap; // largest gap
          gapCenter = (p.rapidity()+previousRapidity)/2.; // find the center of the gap to separate the X and Y systems.
        }
        previousRapidity = p.rapidity();
      }

      FourMomentum mxFourVector, myFourVector;
      for (const Particle& p : fsparticles)
        (p.rapidity() > gapCenter ? mxFourVector : myFourVector) += p.momentum();
      const double xiX = mxFourVector.mass2()/sqr(sqrtS());
      const double xiY = myFourVector.mass2()/sqr(sqrtS());
      const double xi = max(xiX, xiY);
      if (xi < 1e-6) vetoEvent;

      double totEnergy = 0, emEnergy = 0, hadEnergy = 0;
      for (const Particle& p : fsparticles) {
        if (!inRange(p.eta(), -6.6, -5.2)) continue; //< @todo Should be abseta()?
        if (!p.isVisible() || p.abspid() == PID::MUON) continue;
        totEnergy += p.energy();
        if (p.abspid() == PID::ELECTRON || p.abspid() == PID::PHOTON || p.abspid() == 111) emEnergy += p.energy();
        if (p.abspid() != PID::ELECTRON && p.abspid() != PID::PHOTON && p.abspid() != 111) hadEnergy += p.energy();
      }

      const double weight = event.weight();
      _h_totEnergy->fill(totEnergy/GeV, weight);
      _h_emEnergy->fill(emEnergy/GeV, weight);
      _h_hadEnergy->fill(hadEnergy/GeV, weight);
    }


    /// Normalize histograms
    void finalize() {
      scale({_h_totEnergy, _h_emEnergy, _h_hadEnergy}, crossSection()/microbarn/sumOfWeights());
    }


  private:

    Histo1DPtr _h_totEnergy, _h_emEnergy, _h_hadEnergy;

  };


  DECLARE_RIVET_PLUGIN(CMS_2017_I1511284);

}
