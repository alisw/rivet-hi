// -*- C++ -*-
#ifndef RIVET_HeavyIonAnalysis_HH
#define RIVET_HeavyIonAnalysis_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {

  class HeavyIonAnalysis : public Analysis {

  public:
    HeavyIonAnalysis(const std::string &name = "");

    // methods specific to heavy-ion analyses
    enum CentralityMethod {
      ImpactParameter,
      Multiplicity,
      Undefined
    };

    double quantile(double value, const Histo1DPtr hist) const;

    void addCentralityMethod(const CentralityMethod method,
			     const size_t numEventsRequired,
			     const string methodID,
			     const FinalState *fs = 0x0);

    double centrality(const Event& event, const string methodID = "") const;

    bool is_heavy_ion(const Event& event) const;

  private:

    void centrality_method_not_found() const;

    float calculateObservable(const Event &event, const int index) const;

    void checkObservable(const float observable, const int index) const;

    string trimString(const string& str) const;

    /// Histogram for centrality calibration
    std::vector<Histo1DPtr> _histCentralityCalibrationVector;

    /// Histogram for centrality control. It may be used to compare distribution
    /// in the current run to the one provided in calibration histogram.
    std::vector<Histo1DPtr> _histCentralityControlVector;

    /// String with the cuts
    std::vector<string> _cutStringVector;

    /// Method of the centrality calibration
    std::vector<CentralityMethod> _centmethodVector;

    /// Number of events required for a selected method of centrality calibration
    std::vector<size_t> _numEventsRequiredVector;

    /// Name of the centrality calibration method
    std::vector<string> _methodNameVector;

    /// ID of the centrality method
    std::vector<string> _methodIDVector;

  };
}

#endif
