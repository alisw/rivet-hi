// -*- C++ -*-
#include "pluginALICE/HeavyIonAnalysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/Cuts.hh"
#include <fstream>

#define _USE_MATH_DEFINES
#include <math.h>

namespace Rivet {

  class ALICE_2012_I930312 : public HeavyIonAnalysis {

  public:

    ALICE_2012_I930312() :
      HeavyIonAnalysis("ALICE_2012_I930312")
    {  }


    void init() {
      HeavyIonAnalysis::init();

      // Set centrality method
      addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 10000, "method1");

      // Charged final states with |eta| < 1.0 and 8 < pT < 15 GeV/c for trigger particles
      const Cut& cutTrigger = Cuts::abseta < 1.0 && Cuts::pT > 8*GeV && Cuts::pT < 15*GeV;
      const ChargedFinalState cfsTrigger(cutTrigger);
      addProjection(cfsTrigger,"CFSTrigger");

      // Set limit values of pT bins
      pt_limits[0] = 3;
      pt_limits[1] = 4;
      pt_limits[2] = 6;
      pt_limits[3] = 8;
      pt_limits[4] = 10;

      // Chaged final states with |eta| < 1.0 and different pT bins for associated particles
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	Cut mycut = Cuts::abseta < 1.0 && Cuts::pT > pt_limits[ipt]*GeV && Cuts::pT < pt_limits[ipt + 1]*GeV;
	addProjection(ChargedFinalState(mycut), "CFSAssoc" + std::to_string(ipt));
      }

      // Create event strings
      event_string[0] = "pp";
      event_string[1] = "central";
      event_string[2] = "peripheral";
      event_string[3] = "other";

      // For each event type
      for (int itype = 0; itype < EVENT_TYPES; itype++) {
	// For each pT range
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  // Initialize yield histograms
	  _histYield[itype][ipt] = bookHisto1D("Yield_" + event_string[itype] + "_" + std::to_string(ipt),
					       36, -0.5 * M_PI, 1.5 * M_PI,
					       "Associated particle per trigger particle yield",
					       "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	  _histYieldBkgRemoved[itype][ipt] = bookHisto1D("Yield_" + event_string[itype] + "_nobkg_" + std::to_string(ipt),
							 36, -0.5 * M_PI, 1.5 * M_PI,
							 "Associated particle per trigger particle yield no bkg",
							 "#Delta#eta (rad)", "1 / N_trig dN_assoc / d#Delta#eta (rad^-1)");
	}
      }

      // Histogram for counting trigger particles for each event type
      _histTriggerCounter = bookHisto1D("Trigger", EVENT_TYPES, 0.0, EVENT_TYPES, "Trigger counter", "event type", "N");

      // Initialize IAA and ICP histograms
      _histIAA[0] = bookScatter2D(1, 1, 1);
      _histIAA[1] = bookScatter2D(2, 1, 1);
      _histIAA[2] = bookScatter2D(5, 1, 1);

      _histIAA[3] = bookScatter2D(3, 1, 1);
      _histIAA[4] = bookScatter2D(4, 1, 1);
      _histIAA[5] = bookScatter2D(6, 1, 1);

    }

    void analyze(const Event& event) {

      // Check if heavy-ion event
      if (HeavyIonAnalysis::is_heavy_ion(event)) {
	std::cout << "HI EVENT: ";
      } else {
	std::cout << "PP EVENT: ";
      }
      std::cout << event.genEvent()->heavy_ion();

      // Create charged final state for trigger particle
      const ChargedFinalState& triggerFinalState = applyProjection<ChargedFinalState>(event, "CFSTrigger");
      Particles triggerParticles = triggerFinalState.particlesByPt();
      std::cout << "trig: " << triggerParticles.size();

      // Create charged final state for associated particle
      ChargedFinalState associatedFinalState[PT_BINS];
      Particles associatedParticles[PT_BINS];
      std::cout << ", assoc: ";
      for (int ipt = 0; ipt < PT_BINS; ipt++) {
	associatedFinalState[ipt] = applyProjection<ChargedFinalState>(event, "CFSAssoc" + std::to_string(ipt));
	associatedParticles[ipt] = associatedFinalState[ipt].particlesByPt();
	std::cout << associatedParticles[ipt].size() << " ";
      }
      std::cout << std::endl;

      // Check event type
      if (HeavyIonAnalysis::is_heavy_ion(event)) {
	double centr = centrality(event, "method1");
	std::cout << "centrality: " << centr << std::endl;
	if (centr > 0.0 && centr < 5.0)
	  event_type = 1; // PbPb, central
	else if (centr > 60.0 && centr < 90.0)
	  event_type = 2; // PbPb, peripherial
	else
	  event_type = 3; // PbPb, other
      } else {
	event_type = 0; // pp
      }
      std::cout << "ev type: " << event_string[event_type] << " (" << event_type << ")" << std::endl;

      // Veto event if not valid event type
      if (event_type == 3)
	vetoEvent;

      // Fill trigger histogram for a proper event type
      _histTriggerCounter->fill(event_type, triggerParticles.size());

      // Loop over trigger particles
      foreach (const Particle& triggerParticle, triggerParticles) {
	std::cout << "Trigger particle" << std::endl;
	// For each pt bin
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  std::cout << "pt bin " << ipt << std::endl;
	  // Loop over associated particles
	  foreach (const Particle& associatedParticle, associatedParticles[ipt]) {
	    std::cout << "Associated particle" << std::endl;
	    // If associated and trigger particle are not the same particles
	    if (associatedParticle != triggerParticle) {
	      // If pt of trigger particle is greater than pt of associated particle
	      if (triggerParticle.pt() > associatedParticle.pt()) {
		// Calculate delta phi in range (-0.5*PI, 1.5*PI)
		double dPhi = triggerParticle.phi() - associatedParticle.phi();
		while (dPhi > 1.5 * M_PI)  { dPhi -= 2 * M_PI; }
		while (dPhi < -0.5 * M_PI) { dPhi += 2 * M_PI; }
		// Fill yield histogram for calculated delta phi
		std::cout << "Filling yield histogram for pt = " << pt_limits[ipt] << "-" << pt_limits[ipt + 1] << " GeV/c for delta phi = " << dPhi << std::endl;
		_histYield[event_type][ipt]->fill(dPhi, 1);
	      }
	    }
	  }
	}
      }
      std::cout << std::endl;

    }

    void finalize() {

      std::cout << "Finalizing..." << std::endl;

    }

    void post() {

      // Variables for near and away side peak calculation
      double nearSide[EVENT_TYPES][PT_BINS] = { {0.0} };
      double awaySide[EVENT_TYPES][PT_BINS] = { {0.0} };

      // Variables for background error calculation
      double background[EVENT_TYPES][PT_BINS] = { {0.0} };
      double backgroundError[EVENT_TYPES][PT_BINS] = { {0.0} };

      // Variables for integration error calculation
      double scalingFactor[EVENT_TYPES] = {0.0};
      double numberOfEntries[EVENT_TYPES][PT_BINS][2] = { { {0.0} } };
      int numberOfBins[EVENT_TYPES][PT_BINS][2] = { { {0} } };

      std::cout << "Trigger counter histogram entries: " <<
	_histTriggerCounter->bin(0).sumW() << " " <<
	_histTriggerCounter->bin(1).sumW() << " " <<
	_histTriggerCounter->bin(2).sumW() << std::endl;
      std::cout << "pp yield pt bin 0, entries: " << _histYield[0][0]->numEntries() << std::endl;
      std::cout << "pp yield pt bin 1, entries: " << _histYield[0][1]->numEntries() << std::endl;
      std::cout << "pp yield pt bin 2, entries: " << _histYield[0][2]->numEntries() << std::endl;
      std::cout << "pp yield pt bin 3, entries: " << _histYield[0][3]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 0, entries: " << _histYield[1][0]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 1, entries: " << _histYield[1][1]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 2, entries: " << _histYield[1][2]->numEntries() << std::endl;
      std::cout << "Central yield pt bin 3, entries: " << _histYield[1][3]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 0, entries: " << _histYield[2][0]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 1, entries: " << _histYield[2][1]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 2, entries: " << _histYield[2][2]->numEntries() << std::endl;
      std::cout << "Peripheral yield pt bin 3, entries: " << _histYield[2][3]->numEntries() << std::endl;

      // For each event type
      for (int itype = 0; itype < EVENT_TYPES; itype++) {
	// For each pT range
	for (int ipt = 0; ipt < PT_BINS; ipt++) {

	  // Check if histograms are fine
	  if (_histTriggerCounter->numEntries() == 0 || _histYield[itype][ipt]->numEntries() == 0) {
	    std::cout << _histTriggerCounter->numEntries() << std::endl;
	    std::cout << _histYield[itype][ipt]->numEntries() << std::endl;
	    std::cout << "There are no entries in one of the histograms" << std::endl;
	    continue;
	  }

	  // Scale yield histogram
	  std::cout << "Scaling yield histograms..." << std::endl;
	  if ((_histTriggerCounter->bin(itype).sumW() != 0)) {
	    std::cout << "Scaling histogram type " << itype << " pt bin " << ipt << "..." << std::endl;
	    std::cout << "Scaling factor: " << _histTriggerCounter->bin(itype).sumW() << std::endl;
	    scalingFactor[itype] = 1. / _histTriggerCounter->bin(itype).sumW();
	    std::cout << "Bin 1 before scaling: " << _histYield[itype][ipt]->bin(1).sumW() << std::endl;
	    scale(_histYield[itype][ipt], (1. / _histTriggerCounter->bin(itype).sumW()));
	    std::cout << "Bin 1 after scaling: " << _histYield[itype][ipt]->bin(1).sumW() << std::endl;
	  }

	  // Calculate background
	  std::cout << "Calculating background" << std::endl;
	  double sum = 0.0;
	  int nbins = 0;
	  for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	    std::cout << "Bin " << ibin << std::endl;
	    if ((_histYield[itype][ipt]->bin(ibin).xMid() > (-0.5 * M_PI) &&
		 _histYield[itype][ipt]->bin(ibin).xMid() < (-0.5 * M_PI + 0.4)) ||
		(_histYield[itype][ipt]->bin(ibin).xMid() > (0.5 * M_PI - 0.4) &&
		 _histYield[itype][ipt]->bin(ibin).xMid() < (0.5 * M_PI + 0.4)) ||
		(_histYield[itype][ipt]->bin(ibin).xMid() > (1.5 * M_PI - 0.4) &&
		 _histYield[itype][ipt]->bin(ibin).xMid() < (1.5 * M_PI))) {
	      std::cout << "Adding " << _histYield[itype][ipt]->bin(ibin).sumW() << std::endl;
	      sum += _histYield[itype][ipt]->bin(ibin).sumW();
	      nbins++;
	    }
	  }
	  if (nbins == 0) {
	    std::cout << "Failed to estimate background!" << std::endl;
	    continue;
	  }
	  std::cout << "Dividing " << sum << " / " << nbins << std::endl;
	  background[itype][ipt] = sum / nbins;
	  std::cout << "background: " << background[itype][ipt] << std::endl;

	  // Calculate background error
	  sum = 0.0;
	  nbins = 0;
	  for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	    if (_histYield[itype][ipt]->bin(ibin).xMid() > (0.5 * M_PI - 0.4) &&
		_histYield[itype][ipt]->bin(ibin).xMid() < (0.5 * M_PI + 0.4)) {
	      sum += (_histYield[itype][ipt]->bin(ibin).sumW() - background[itype][ipt]) *
		(_histYield[itype][ipt]->bin(ibin).sumW() - background[itype][ipt]);
	      nbins++;
	    }
	  }
	  backgroundError[itype][ipt] = sqrt(sum / (nbins - 1));
	  std::cout << "backgroundError: " << backgroundError[itype][ipt] << std::endl;

	  // Fill histograms with removed background
	  std::cout << "Filling histogram with removed background..." << std::endl;
	  for (unsigned int ibin = 0; ibin < _histYield[itype][ipt]->numBins(); ibin++) {
	    std::cout << "Bin " << ibin << ", value " << _histYield[itype][ipt]->bin(ibin).sumW() - background[itype][ipt] << std::endl;
	    _histYieldBkgRemoved[itype][ipt]->fillBin(ibin, _histYield[itype][ipt]->bin(ibin).sumW() - background[itype][ipt]);
	  }

	  std::cout << "Integrating the whole histogram: " << _histYieldBkgRemoved[itype][ipt]->integral() << std::endl;

	  // Integrate near-side yield
	  std::cout << "Integrating near-side yield..." << std::endl;
	  unsigned int lowerBin = _histYield[itype][ipt]->binIndexAt(-0.7 + 0.02);
	  unsigned int upperBin = _histYield[itype][ipt]->binIndexAt( 0.7 - 0.02) + 1;
	  nbins = upperBin - lowerBin;
	  numberOfBins[itype][ipt][0] = nbins;
	  std::cout << _histYield[itype][ipt]->integralRange(1, 1) << " " << _histYield[itype][ipt]->bin(1).sumW() << std::endl;
	  std::cout << _histYield[itype][ipt]->integralRange(1, 2) << " " << _histYield[itype][ipt]->bin(1).sumW() + _histYield[itype][ipt]->bin(2).sumW() << std::endl;
	  nearSide[itype][ipt] = _histYield[itype][ipt]->integralRange(lowerBin, upperBin) - nbins * background[itype][ipt];
	  numberOfEntries[itype][ipt][0] = _histYield[itype][ipt]->integralRange(lowerBin, upperBin) * _histTriggerCounter->bin(itype).sumW();
	  std::cout << "_histYield[" << itype << "][" << ipt << "]->integralRange(" << lowerBin << ", " << upperBin << ") - nbins * bkg = " <<
	    _histYield[itype][ipt]->integralRange(lowerBin, upperBin) << " - " << nbins << " * " << background[itype][ipt] << " = " <<
	    nearSide[itype][ipt] << std::endl;

	  // Integrate away-side yield
	  std::cout << "Integrating away-side yield..." << std::endl;
	  lowerBin = _histYield[itype][ipt]->binIndexAt(M_PI - 0.7 + 0.02);
	  upperBin = _histYield[itype][ipt]->binIndexAt(M_PI + 0.7 - 0.02) + 1;
	  nbins = upperBin - lowerBin;
	  numberOfBins[itype][ipt][1] = nbins;
	  awaySide[itype][ipt] = _histYield[itype][ipt]->integralRange(lowerBin, upperBin) - nbins * background[itype][ipt];
	  numberOfEntries[itype][ipt][1] = _histYield[itype][ipt]->integralRange(lowerBin, upperBin) * _histTriggerCounter->bin(itype).sumW();
	  std::cout << "_histYield[" << itype << "][" << ipt << "]->integralRange(" << lowerBin << ", " << upperBin << ") - nbins * bkg = " <<
	    _histYield[itype][ipt]->integralRange(lowerBin, upperBin) << " - " << nbins << " * " << background[itype][ipt] << " = " <<
	    awaySide[itype][ipt] << std::endl;

	}
      }

      // Variables for IAA/ICP plots
      double dI = 0.0;
      int near = 0;
      int away = 1;
      double xval[PT_BINS] = { 3.5, 5.0, 7.0, 9.0 };
      double xerr[PT_BINS] = { 0.5, 1.0, 1.0, 1.0 };

      int types1[3] = {1, 2, 1};
      int types2[3] = {0, 0, 2};
      string type_string[3] = {"I_AA 0-5%/pp", "I_AA 60-90%/pp", "I_CP 0-5%/60-90%"};

      // Fill IAA/ICP plots for near side peak
      std::cout << "Near side" << std::endl;
      for (int ihist = 0; ihist < 3; ihist++) {
	int type1 = types1[ihist];
	int type2 = types2[ihist];
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  std::cout << type_string[ihist] << ", pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << nearSide[type1][ipt] / nearSide[type2][ipt];
	  dI = scalingFactor[type1] * scalingFactor[type1] * numberOfEntries[type1][ipt][near] +
	    scalingFactor[type2] * scalingFactor[type2] * numberOfEntries[type2][ipt][near] * nearSide[type1][ipt] * nearSide[type1][ipt] / (nearSide[type2][ipt] * nearSide[type2][ipt]) +
	    numberOfBins[type1][ipt][near] * numberOfBins[type1][ipt][near] * backgroundError[type1][ipt] * backgroundError[type1][ipt] +
	    numberOfBins[type2][ipt][near] * numberOfBins[type2][ipt][near] * backgroundError[type2][ipt] * backgroundError[type2][ipt] * nearSide[type1][ipt] * nearSide[type1][ipt] / (nearSide[type2][ipt] * nearSide[type2][ipt]);
	  dI = sqrt(dI)/nearSide[type2][ipt];
	  std::cout << " +- " << dI << std::endl;
	  _histIAA[ihist]->addPoint(xval[ipt], nearSide[type1][ipt] / nearSide[type2][ipt], xerr[ipt], dI);
	}
      }

      // Fill IAA/ICP plots for away side peak
      std::cout << "Away side" << std::endl;
      for (int ihist = 0; ihist < 3; ihist++) {
	int type1 = types1[ihist];
	int type2 = types2[ihist];
	for (int ipt = 0; ipt < PT_BINS; ipt++) {
	  std::cout << type_string[ihist] << ", pt=(" << pt_limits[ipt] << ", " << pt_limits[ipt+1] << "): " << awaySide[type1][ipt] / awaySide[type2][ipt];
	  dI = scalingFactor[type1] * scalingFactor[type1] * numberOfEntries[type1][ipt][away] +
	    scalingFactor[type2] * scalingFactor[type2] * numberOfEntries[type2][ipt][away] * awaySide[type1][ipt] * awaySide[type1][ipt] / (awaySide[type2][ipt] * awaySide[type2][ipt]) +
	    numberOfBins[type1][ipt][away] * numberOfBins[type1][ipt][away] * backgroundError[type1][ipt] * backgroundError[type1][ipt] +
	    numberOfBins[type2][ipt][away] * numberOfBins[type2][ipt][away] * backgroundError[type2][ipt] * backgroundError[type2][ipt] * awaySide[type1][ipt] * awaySide[type1][ipt] / (awaySide[type2][ipt] * awaySide[type2][ipt]);
	  dI = sqrt(dI)/awaySide[type2][ipt];
	  std::cout << " +- " << dI << std::endl;
	  _histIAA[ihist + 3]->addPoint(xval[ipt], awaySide[type1][ipt] / awaySide[type2][ipt], xerr[ipt], dI);
	}
      }

    }

  private:

    static const int PT_BINS = 4;
    static const int EVENT_TYPES = 3;

    Histo1DPtr _histYield[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histYieldBkgRemoved[EVENT_TYPES][PT_BINS];
    Histo1DPtr _histTriggerCounter;
    Scatter2DPtr _histIAA[6];
    int pt_limits[5];
    int event_type;
    string event_string[EVENT_TYPES + 1];

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I930312);
}
