// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class OPAL_2004_I631361 : public Analysis {
  public:

    /// Constructor
    OPAL_2004_I631361()
      : Analysis("OPAL_2004_I631361"), _sumW(0.0)
    {    }


    /// @name Analysis methods
    //@{
    void init() {
      declare(FinalState(), "FS");
      declare(ChargedFinalState(), "CFS");
      int ih(0), iy(0);
      if (inRange(0.5*sqrtS()/GeV, 5.0, 5.5)) {
        ih = 1;
	iy = 1;
      } else if (inRange(0.5*sqrtS()/GeV, 5.5, 6.5)) {
        ih = 1;
	iy = 2;
      } else if (inRange(0.5*sqrtS()/GeV, 6.5, 7.5)) {
        ih = 1;
	iy = 3;
      } else if (inRange(0.5*sqrtS()/GeV, 7.5, 9.5)) {
        ih = 2;
	iy = 1;
      } else if (inRange(0.5*sqrtS()/GeV, 9.5, 13.0)) {
        ih = 2;
	iy = 2;
      } else if (inRange(0.5*sqrtS()/GeV, 13.0, 16.0)) {
        ih = 3;
	iy = 1;
      } else if (inRange(0.5*sqrtS()/GeV, 16.0, 20.0)) {
        ih = 3;
	iy = 2;
      }
      assert(ih>0);
      _h_chMult     = bookHisto1D(ih,1,iy);
      if(ih==3)
	_h_chFragFunc = bookHisto1D(5,1,iy);
      else
	_h_chFragFunc = NULL;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      // find the initial gluons
      ParticleVector initial;
      for (const GenParticle* p : Rivet::particles(event.genEvent())) {
	const GenVertex* pv = p->production_vertex();
	const PdgId pid = p->pdg_id();
	if(pid!=21) continue;
	bool passed = false;
	for (const GenParticle* pp : particles_in(pv)) {
	  const PdgId ppid = abs(pp->pdg_id());
	  passed = (ppid == PID::ELECTRON || ppid == PID::HIGGS || 
		    ppid == PID::ZBOSON   || ppid == PID::GAMMA);
	  if(passed) break;
	}
	if(passed) initial.push_back(Particle(*p));
      }
      if(initial.size()!=2) vetoEvent;
      // use the direction for the event axis
      Vector3 axis = initial[0].momentum().p3().unit();
      // fill histograms
      const Particles& chps = applyProjection<FinalState>(event, "CFS").particles();
      unsigned int nMult[2] = {0,0};
      _sumW += 2.*weight;
      // distribution
      foreach(const Particle& p, chps) {
        double xE = 2.*p.E()/sqrtS();
	if(_h_chFragFunc) _h_chFragFunc->fill(xE, weight);
	if(p.momentum().p3().dot(axis)>0.)
	  ++nMult[0];
	else
	  ++nMult[1];
      }
      // multiplcities in jet
      _h_chMult->fill(nMult[0],weight);
      _h_chMult->fill(nMult[1],weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_chMult);
      if(_h_chFragFunc) scale(_h_chFragFunc, 1./_sumW);
    }

    //@}


  private:

    double _sumW;

    /// @name Histograms
    //@{
    Histo1DPtr _h_chMult;
    Histo1DPtr _h_chFragFunc;
    //@}
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2004_I631361);


}
