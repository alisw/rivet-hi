// -*- C++ -*-
#ifndef RIVET_BINNEDHISTOGRAM_HH
#define RIVET_BINNEDHISTOGRAM_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {

  class Analysis;


  /// A set of booked Histo1DPtr, each in a bin of a second variable.
  ///
  /// BinnedHistogram contains a series of histograms of the same quantity
  /// each in a different region of a second quantity.  For example, a
  /// BinnedHistogram may contain histograms of the cross-section differential
  /// in \f$ p_T \f$ in different \f$ \eta \f$  regions.
  template<typename T>
  class BinnedHistogram {
  public:

    /// Create a new empty BinnedHistogram
    BinnedHistogram() {    }

    /// Create a new BinnedHistogram with the given bin edges and contents
    BinnedHistogram(const vector<T>& edges, const vector<Histo1DPtr>& histos) {
      assert(edges.size() == histos.size()+1);
      for (size_t i = 0; i < histos.size(); ++i)
        addHistogram(edges[i], edges[i+1], histos[i]);
    }

    /// @todo Can we have an "emplace constructor", passing tuples of bookHisto1D args?


    ///  Add a histogram in the @c T bin between @a binMin and @a binMax
    const BinnedHistogram<T>& add(const T& binMin, const T& binMax, Histo1DPtr histo);
    /// Clumsier alias
    /// @deprecated Prefer add()
    const BinnedHistogram<T>& addHistogram(const T& binMin, const T& binMax, Histo1DPtr histo) {
      return add(binMin, binMax, histo);
    }


    /// Fill the histogram in the same bin as @a binval with value @a val and weight @a weight
    void fill(const T& binval, double val, double weight);


    /// @brief Get the histogram in the same bin as @a binval (const)
    /// @note Throws a RangeError if @a binval doesn't fall in a declared bin
    const Histo1DPtr histo(const T& binval) const;
    /// @brief Get the histogram in the same bin as @a binval
    /// @note Throws a RangeError if @a binval doesn't fall in a declared bin
    Histo1DPtr histo(const T& binval);

    /// Get the contained histograms (const)
    const vector<Histo1DPtr>& histos() const { return _histos; }
    /// Get the contained histograms
    vector<Histo1DPtr>& histos() { return _histos; }
    /// @deprecated Prefer histos()
    const vector<Histo1DPtr>& getHistograms() const { return _histos; }
    /// @deprecated Prefer histos()
    vector<Histo1DPtr>& getHistograms() { return _histos; }


    /// Scale histograms taking into account its "external" binwidth, i.e. by scale/binWidth
    /// @note The Analysis pointer is passed in order to call the analysis' scale(h) method: can we avoid that?
    void scale(const T& scale, Analysis* ana);


  private:

    map<T, Histo1DPtr> _histosByUpperBound, _histosByLowerBound;
    vector<Histo1DPtr> _histos;
    map<Histo1DPtr, T> _binWidths;

  };


}

#endif
