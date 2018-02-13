// -*- C++ -*-
#include "pluginALICE/HeavyIonAnalysis.hh"

namespace Rivet {
  HeavyIonAnalysis::HeavyIonAnalysis(const std::string &name) :
      Analysis(name)
    { }

  void HeavyIonAnalysis::addCentralityMethod(const CentralityMethod method,
					     const size_t numEventsRequired,
					     const string methodID,
					     const FinalState *fs) {

    std::cout << "Adding centrality method" << std::endl;

    switch (method) {

    case ImpactParameter:
      _methodNameVector.push_back("ImpactParameter");
      std::cout << "Method: ImpactParameter" << std::endl;
      _histCentralityCalibrationVector.push_back(bookHisto1D("calib_impactpar", 200, 0, 20.0, "Calibration histogram", "xlabel", "ylabel"));
      _histCentralityControlVector.push_back(bookHisto1D("control_impactpar", 200, 0, 20.0, "Control histogram", "xlabel", "ylabel"));
      break;

    case Multiplicity:
      if (!fs)
	throw Exception("something bad happened ...");
      addProjection(*fs, "FS" + methodID);
      _methodNameVector.push_back("Multiplicity");
      std::cout << "Method: Multiplicity" << std::endl;
      _histCentralityCalibrationVector.push_back(bookHisto1D("calib_multiplicity_" + trimString(fs->getCuts()->description()), 10000, 0, 10000.0, "Calibration histogram", "xlabel", "ylabel"));
      _histCentralityControlVector.push_back(bookHisto1D("control_multiplicity_" + trimString(fs->getCuts()->description()), 10000, 0, 10000.0, "Control histogram", "xlabel", "ylabel"));
      break;

    default:
      centrality_method_not_found();
    }

    _numEventsRequiredVector.push_back(numEventsRequired);
    _centmethodVector.push_back(method);
    _methodIDVector.push_back(methodID);
  }

  void HeavyIonAnalysis::centrality_method_not_found() const {
    throw Exception("Unimplemented method of centrality calibration.");
  }

  double HeavyIonAnalysis::centrality(const Event& event, const string methodID) const {

    // Check if there is at least one registered method
    if (_methodIDVector.size() < 1)
      throw Exception("No centrality methods registered!");

    int index = 0;
    // If there is no methodID specified
    if (methodID.empty()) {
      // Get the first method as default one
      index = std::distance(_methodIDVector.begin(), _methodIDVector.begin());
      MSG_WARNING("Centrality method ID not specified. Using the first declared method as default one...");
    } else {
      // methodID is provided. Check if it is registered
      auto iterID = std::find(_methodIDVector.begin(), _methodIDVector.end(), methodID);
      if (iterID == _methodIDVector.end())
	throw Exception("Method not found!");
      index = std::distance(_methodIDVector.begin(), iterID);
    }

    // Initialize value of the event observable to be filled to the histograms
    if (_centmethodVector.at(index) == Undefined)
      throw Exception("no centrality method defined");

    const float observable = calculateObservable(event, index);
    float centrality = -1.;

    // Check if there are enough events in the calibration histogram
    if (_histCentralityCalibrationVector.at(index)->numEntries() >= _numEventsRequiredVector.at(index)) {
      // Calculate centrality as percentile using observable calculated with the selected method
      centrality = quantile(observable, _histCentralityCalibrationVector.at(index)) * 100.;
      if (_centmethodVector.at(index) == Multiplicity)
	centrality = 100. - centrality;
      // Fill the control histogram with the impact parameter value and event weight
      MSG_INFO("Adding control point nr " << _histCentralityControlVector.at(index)->numEntries() << ": (" << observable << ", " << event.weight() << ")");
      _histCentralityControlVector.at(index)->fill(observable);//, event.weight());
    }
    // Otherwise, if there are not enough events in the calibration histogram
    else {
      // Fill the calibration histogram with the impact parameter value and event weight
      MSG_INFO("Adding calibration point nr " << _histCentralityCalibrationVector.at(index)->numEntries() << ": (" << observable << ", " << event.weight() << ")");
      _histCentralityCalibrationVector.at(index)->fill(observable);//, event.weight());
    }
    return centrality;
  }

  float HeavyIonAnalysis::calculateObservable(const Event& e, const int index) const {

    // Initialize observable
    float observable = -1.;

    // Calculate its value according to the selected method
    switch (_centmethodVector.at(index)) {

    case ImpactParameter:
	// Get impact parameter of the event
	observable = e.genEvent()->heavy_ion() ? e.genEvent()->heavy_ion()->impact_parameter() : -1.;
	break;

    case Multiplicity:
      {
	// Get multiplicity of the event
	const FinalState& fs = applyProjection<FinalState>(e, "FS" + _methodIDVector.at(index));
	observable = e.genEvent()->heavy_ion() ? fs.particles().size() : -1.;
	break;
      }

    default:
      centrality_method_not_found();
    }

    // Check if observable is correct. If not, skip this event
    checkObservable(observable, index);

    return observable;
  }


  void HeavyIonAnalysis::checkObservable(const float observable, const int index) const {

    try {
      // Check if the observable is greater than 0
      if (observable < 0.)
	throw UserError("Calculated observable is lower than 0!");
      // Check if the observable is inside histogram ranges
      if ((observable < _histCentralityCalibrationVector.at(index)->xMin()) ||
	  (observable > _histCentralityCalibrationVector.at(index)->xMax()) ||
	  (_histCentralityCalibrationVector.at(index)->binIndexAt(observable) < 0))
	throw RangeError("Calculated observable is out of bounds!");
    } catch (...) {
      MSG_WARNING("Skipping this event...");
    }

  }

  bool HeavyIonAnalysis::is_heavy_ion(const Event& event) const {
    // Get pair of beam particles from HepMC event
    std::pair<HepMC::GenParticle*,HepMC::GenParticle*> beam_pair = event.genEvent()->beam_particles();
    std::cout << "PDG ID 1: " << (std::get<0>(beam_pair))->pdg_id() << std::endl;
    std::cout << "PDG ID 2: " << (std::get<1>(beam_pair))->pdg_id() << std::endl;
    // Check if beam particles are heavy ions
    // @note improve this function to make it work for all kind of particles
    bool is_heavy_ion = ((std::get<0>(beam_pair))->pdg_id() == 1000822080) && ((std::get<1>(beam_pair))->pdg_id() == 1000822080);
    return is_heavy_ion;
  }

  string HeavyIonAnalysis::trimString(const string& str) const {
    // @note Which characters should be removed?
    string trimmedString = str;

    const char characters[] = "\\ ";
    for (const auto character : characters) {
      trimmedString.erase (std::remove (trimmedString.begin(), trimmedString.end(), character), trimmedString.end());
    }
    return trimmedString;
  }

  // Function for calculating the quantile from histogram at selected value
  // (below this value, including overflow)
  // @note this function is here temporarily, it will be moved and modified
  // to create a general function for calculating the quantile of a histogram
  double HeavyIonAnalysis::quantile(double value, const Histo1DPtr hist) const {

    // Check if there are entries in the histogram
    if (hist->numEntries() == 0) {
      throw WeightError("There are no entires in the histogram!");
    }

    // Integration ranges
    size_t upperBin = hist->binIndexAt(value);

    // Calculate centrality as percentile
    std::cout << "quantile: " << hist->integralTo(upperBin, true) << " / " << hist->integral(true) << std::endl;
    return hist->integralTo(upperBin, true) / hist->integral(true);

  }
}
