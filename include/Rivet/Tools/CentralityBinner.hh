// -*- C++ -*-
#ifndef RIVET_CENTRALITYBINNER_HH
#define RIVET_CENTRALITYBINNER_HH
#include <tuple>
#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {

/**
   @brief Base class for projections giving the value of an
   observable sensitive to the centrality of a collision.

   @author Leif Lönnblad

   The centrality of a collision is not really an observable, but the
   concept is anyway often used in the heavy ion community as if it
   were just that.

   This base class can be used to provide a an estimator for the
   centrality by projecting down to a single number which then can be
   used by a CentralityBinner object to select a histogram to be
   filled with another observable depending on centrality percentile.

   The estimate() should be a non-negative number with large values
   indicating a higher overlap than small ones. A negative value
   indicates that the centrality estimate could not be calculated.

   In the best of all worlds the centrality estimator should be a
   proper hadron-level observable corrected for detector effects,
   however, this base class only returns the inverse of the
   impact_parameter member of the GenHeavyIon object in an GenEvent
   if present and zero otherwise.
*/
class CentralityEstimator : public Projection {
public:

  /// Constructor.
  CentralityEstimator(): _estimate(-1.0) {}

  /// Clone on the heap.
  DEFAULT_RIVET_PROJ_CLONE(CentralityEstimator);

protected:

  /// Perform the projection on the Event
  void project(const Event& e) {
    _estimate = -1.0;
    const HepMC::HeavyIon * hi = e.genEvent()->heavy_ion();
    if ( hi ) _estimate = hi->impact_parameter() > 0.0?
                1.0/hi->impact_parameter(): numeric_limits<double>::max();
  }

  /// Compare projections
  int compare(const Projection& p) const {
    return mkNamedPCmp(p, "CentEst");
  }


public:

  /// The value of the centrality estimate.
  double estimate() const { return _estimate; }


protected:

  /// The value of the centrality estimate.
  double _estimate;

};


/// This is a traits class describing how to handle object handles by
/// CentralityBinner. The default implementation basically describes
/// what to do with Histo1DPtr.
template <typename T>
struct CentralityBinTraits {

  /// Make a clone of the given object.
  static T clone(const T & t) {
    return T(t->newclone());
  }

  /// Add the contents of @a o to @a t.
  static void add(T & t, const T & o) {
    *t += *o;
  }

  /// Scale the contents of a given object.
  static void scale(T & t, double f) {
    t->scaleW(f);
  }

  /// Normalize the AnalysisObject to the sum of weights in a
  /// centrality bin.
  static void normalize(T & t, double sumw) {
    if ( t->sumW() > 0.0 ) t->normalize(t->sumW()/sumw);
  }

  /// Return the path of an AnalysisObject.
  static string path(T t) {
    return t->path();
  }

};

/// The sole purpose of the MergeDistance class is to provide a
/// "distance" for a potential merging of two neighboring bins in a
/// CentralityBinner.
struct MergeDistance {

  /// This function should return a generalized distance between two
  /// adjecent centrality bins to be merged. CentralityBinner will
  /// always try to merge bins with the smallest distance. @a cestLo
  /// and @cestHi are the lower and upper edges of resulting bin. @a
  /// weight is the resulting sum of event weights in the bin. @a
  /// centLo and @a centHi are the lower and upper prcentile limits
  /// where the two bins currently resides. The two last arguments are
  /// the total number of events in the two bins and the total number
  /// of previous mergers repectively.
  static double dist(double cestLo, double cestHi, double weight,
                     double clo, double chi, double, double) {
    return (cestHi - cestLo)*weight/(cestHi*(chi - clo));
  }
};


/**
 * CentralityBinner contains a series of AnalysisObject of the same
 * quantity each in a different percentiles of another quantity.  For
 * example, a CentralityBinner may e.g. contain histograms of the
 * cross section differential in \f$ p_T \f$ in different centrality
 * regions for heavy ion collisions based on forward energy flow.
 **/
template <typename T = Histo1DPtr, typename MDist = MergeDistance>
class CentralityBinner: public ProjectionApplier {
  public:

  /// Create a new empty CentralityBinner. @a maxbins is the maximum
  /// number of bins used by the binner. Default is 1000, which is
  /// typically enough. @a wlim is the mximum allowed error allowed
  /// for the centrality limits before a warning is emitted.
  CentralityBinner(int maxbins = 200, double wlim = 0.02)
    : _currentCEst(-1.0), _maxBins(maxbins), _warnlimit(wlim), _weightsum(0.0) {
    _percentiles.insert(0.0);
    _percentiles.insert(1.0);
  }

  /// Set the centrality projection to be used. Note that this
  /// projection must have already been declared to Rivet.
  void setProjection(const CentralityEstimator & p, string pname) {
    declare(p, pname);
    _estimator = pname;
  }

  /// Return the class name.
  virtual std::string name() const {
    return "Rivet::CentralityBinner";
  }

  /// Add an AnalysisObject in the region between @a cmin and @a cmax to
  /// this set of CentralityBinners. The range represent
  /// percentiles and must be between 0 and 100. No overlaping bins
  /// are allowed.

  /// Note that (cmin=0, cmax=5), means the five percent MOST central
  /// events although the internal notation is reversed for
  /// convenience.

  /// Optionally supply corresponding limits @a cestmin and @a cestmax
  /// of the centrality extimator.

  void add(T t, double cmin, double cmax,
           double cestmin = -1.0, double cestmax = -1.0 ) {
    _percentiles.insert(max(1.0 - cmax/100.0, 0.0));
    _percentiles.insert(min(1.0 - cmin/100.0, 1.0));
    if ( _unfilled.empty() && _ready.empty() )
      _devnull = CentralityBinTraits<T>::clone(t);
    if ( cestmin < 0.0 )
      _unfilled.push_back(Bin(t, 1.0 - cmax/100.0, 1.0 - cmin/100.0));
    else
      _ready[t] = Bin(t, 1.0 - cmax/100.0, 1.0 - cmin/100.0, cestmin, cestmax);
  }

  /// Return one of the AnalysisObjects in the CentralityBinner for
  /// the given @a event. This version requires that a
  /// CentralityEstimator object has been assigned that can compute
  /// the value of the centrality estimator from the @a
  /// event. Optionally the @a weight of the event is given. This
  /// should be the weight that will be used to fill the
  /// AnalysisObject. If the centrality estimate is less than zero,
  /// the _devnull object will be returned.
  T select(const Event & event, double weight = 1.0) {
    return select(applyProjection<CentralityEstimator>
		  (event, _estimator).estimate(), weight);
  }

  /// Return one of the AnalysisObjecsts in the Setup the
  /// CentralityBinner depending on the value of the centrality
  /// estimator, @a cest. Optionally the @a weight of the event is
  /// given. This should be the weight that will be used to fill the
  /// AnalysisObject. If the centrality estimate is less than zero,
  /// the _devnull object will be returned.
  T select(double cest, double weight = 1.0);

  /// At the end of the run, calculate the percentiles and fill the
  /// AnalysisObjectss provided with the add() function. This is
  /// typically called from the finalize method in an Analysis, but
  /// can also be called earlier in which case the the select
  /// functions can be continued to run as before with the edges
  /// between the centrality regions now fixed.
  void finalize();

  /// Normalize each AnalysisObjects to the sum of event weights in the
  /// corresponding centrality bin.
  void normalizePerEvent() {
    for ( auto & b : _ready ) b.second.normalizePerEvent();
  }

  /// Return a map bin edges of the centrality extimator indexed by
  /// the corresponing percentile.
  map<double,double> edges() const {
    map<double,double> ret;
    for ( auto & b : _ready ) {
      ret[1.0 - b.second._centLo] = b.second._cestLo;
      ret[1.0 - b.second._centHi] = b.second._cestHi;
    }
    return ret;
  }

  /// Return the current AnalysisObject from the latest call to select().
  const T & current() const {
    return _currenT;
  }

  /// Return the value of the centrality estimator set in the latest
  /// call to select().
  double estimator() const {
    return _currentCEst;
  }

  vector<T> allObjects() {
    vector<T> ret;
    for ( auto & fb : _flexiBins ) ret.push_back(fb._t);
    if ( !ret.empty() ) return ret;
    for ( auto b : _ready ) ret.push_back(b.second._t);
    return ret;
  }

private:

  /// A flexible bin struct to be used to store temporary AnalysisObjects.
  struct FlexiBin {

    /// Construct with an initial centrality estimate and an event
    /// weight.
    FlexiBin(T & t, double cest = 0.0, double weight = 0.0)
      : _t(t), _cestLo(cest), _cestHi(cest), _weightsum(weight), _n(1), _m(0) {}

    /// Construct a temporary FlexiBin for finding a bin in a set.
    FlexiBin(double cest)
      : _cestLo(cest), _cestHi(cest), _weightsum(0.0), _n(0), _m(0) {}

    /// Merge in the contents of another FlexiBin into this.
    void merge(const FlexiBin & fb) {
      _cestLo = min(_cestLo, fb._cestLo);
      _cestHi = max(_cestHi, fb._cestHi);
      _weightsum += fb._weightsum;
      CentralityBinTraits<T>::add(_t, fb._t);
      _n += fb._n;
      _m += fb._m + 1;
    }

    /// Comparisons for containers.
    bool operator< (const FlexiBin & fb) const {
      return _cestLo < fb._cestLo;
    }

    /// Return true if the given centrality estimate is in the range
    /// of this bin.
    bool inRange(double cest) const {
      return cest == _cestLo || ( _cestLo < cest && cest < _cestHi );
    }

    /// The associated AnalysisObject.
    T _t;

    /// Current lower and upper edge of the centrality estimator for
    /// the fills in the associated AnalysiObject.
    double _cestLo, _cestHi;

    /// The sum of weights for all events entering the associated
    /// AnalysisObject.
    mutable double _weightsum;

    /// The number of times this bin has been selected.
    mutable int _n;

    /// The number of times this bin has been merged.
    mutable int _m;

  };

  struct Bin {

    /// Construct a completely empty bin.
    Bin()
      : _centLo(-1.0), _centHi(-1.0), _cestLo(-1.0), _cestHi(-1.0),
	_weightsum(0.0), _underflow(0.0), _overflow(0.0),
        _ambiguous(0), _ambweight(0.0) {}

    /// Constructor taking an AnalysisObject and centrality interval
    /// as argument. Optionally the interval in the estimator can be
    /// given, in which case this bin is considered to be
    /// "final".
    Bin(T t, double centLo, double centHi,
            double cestLo = -1.0, double cestHi = -1.0)
      : _t(t), _centLo(centLo), _centHi(centHi),
        _cestLo(cestLo), _cestHi(cestHi),
	_weightsum(0.0), _underflow(0.0), _overflow(0.0),
        _ambiguous(0.0), _ambweight(0.0) {}

    /// Return true if the given centrality estimate is in the range
    /// of this AnalysisObject.
    bool inRange(double cest) const {
      return _cestLo >= 0 && _cestLo <= cest &&
        ( _cestHi < 0.0 || cest <= _cestHi );
    }

    /// Normalise the AnalysisObject to the tital cross section.
    void normalizePerEvent() {
      CentralityBinTraits<T>::normalize(_t, _weightsum);
    }

    /// The AnalysisObject.
    T _t;

    /// The range in centrality.
    double _centLo, _centHi;

    /// The corresponding range in the centrality estimator.
    double _cestLo, _cestHi;

    /// The sum of event weights for this bin;
    double _weightsum;

    /// The weight in a final AnalysisObject that contains events
    /// below the centrality limit.
    double _underflow;

    /// The weight in a final AnalysisObject that contain events above
    /// the centrality limit.
    double _overflow;

    /// Number of ambiguous events in this bin.
    double _ambiguous;

    /// Sum of abmiguous weights.
    double _ambweight;

  };

protected:

  /// Convenient typedefs.
  typedef set<FlexiBin> FlexiBinSet;

  /// Find a bin corresponding to a given value of the centrality
  /// estimator.
  typename FlexiBinSet::iterator _findBin(double cest) {
    if ( _flexiBins.empty() ) return _flexiBins.end();
    auto it = _flexiBins.lower_bound(FlexiBin(cest));
    if ( it->_cestLo == cest ) return it;
    if ( it != _flexiBins.begin() ) --it;
    if ( it->_cestLo < cest && cest < it->_cestHi ) return it;
    return _flexiBins.end();
  }

  /// The name of the CentralityEstimator projection to be used.
  string _estimator;

  /// The current temporary AnalysisObject selected for the centrality
  /// estimator calculated from the event presented in setup().
  T _currenT;

  /// The current value of the centrality estimator.
  double _currentCEst;

  /// The oversampling of centrality bins. For each requested
  /// centrality bin this number of dynamic bins will be used.
  int _maxBins;

  /// If the fraction of events in a bin that comes from adjecent
  /// centrality bins exceeds this, emit a warning.
  double _warnlimit;

  /// The unfilled AnalysisObjectss where the esimator edges has not yet
  /// been determined.
  vector<Bin> _unfilled;

  /// The dynamic bins for ranges of centrality estimators.
  FlexiBinSet _flexiBins;

  /// The sum of all event weights so far.
  double _weightsum;

  /// The requested percentile limits.
  set<double> _percentiles;

  /// The filled AnalysisObjects where the estimator edges has been determined.
  map<T, Bin> _ready;

  /// A special AnalysisObject which will be filled if the centrality
  /// estimate is out of range (negative).
  T _devnull;

public:

  /// Print out the _flexiBins to cerr.
  void debug();
  void fulldebug();

};


/// Traits specialization for Profile histograms.
template <>
struct CentralityBinTraits<Profile1DPtr> {

  typedef Profile1DPtr T;

  /// Make a clone of the given object.
  static T clone(const T & t) {
    return Profile1DPtr(t->newclone());
  }

  /// Add the contents of @a o to @a t.
  static void add(T & t, const T & o) {
    *t += *o;
  }

  /// Scale the contents of a given object.
  static void scale(T & t, double f) {
    t->scaleW(f);
  }

  static void normalize(T & t, double sumw) {}

  /// Return the path of an AnalysisObject.
  static string path(T t) {
    return t->path();
  }

};


/// Traits specialization for Profile histograms.
template <>
struct CentralityBinTraits<Profile2DPtr> {

  typedef Profile2DPtr T;

  /// Make a clone of the given object.
  static T clone(const T & t) {
    return Profile2DPtr(t->newclone());
  }

  /// Add the contents of @a o to @a t.
  static void add(T & t, const T & o) {
    *t += *o;
  }

  /// Scale the contents of a given object.
  static void scale(T & t, double f) {
    t->scaleW(f);
  }

  static void normalize(T & t, double sumw) {}

  /// Return the name of an AnalysisObject.
  static string path(T t) {
    return t->path();
  }

};

template <typename T>
struct CentralityBinTraits< vector<T> > {

  /// Make a clone of the given object.
  static vector<T> clone(const vector<T> & tv) {
    vector<T> rtv;
    for ( auto t : tv ) rtv.push_back(CentralityBinTraits<T>::clone(t));
    return rtv;
  }

  /// Add the contents of @a o to @a t.
  static void add(vector<T> & tv, const vector<T> & ov) {
    for ( int i = 0, N = tv.size(); i < N; ++i )
      CentralityBinTraits::add(tv[i], ov[i]);
  }

  /// Scale the contents of a given object.
  static void scale(vector<T> & tv, double f) {
    for ( auto t : tv ) CentralityBinTraits<T>::scale(t, f);
  }

  static void normalize(vector<T> & tv, double sumw) {
    for ( auto t : tv ) CentralityBinTraits<T>::normalize(t, sumw);
  }

  /// Return the path of an AnalysisObject.
  static string path(const vector<T> & tv) {
    string ret = "(vector:";
    for ( auto t : tv ) {
      ret += " ";
      ret += CentralityBinTraits<T>::path(t);
    }
    ret += ")";
    return ret;
  }

};

template <size_t I, typename... Types>
struct TupleCentralityBinTraitsHelper {

  typedef tuple<Types...> Tuple;
  typedef typename tuple_element<I-1,Tuple>::type T;

  static void clone(Tuple & ret, const Tuple & tup) {
    get<I-1>(ret) = CentralityBinTraits<T>::clone(get<I-1>(tup));
    TupleCentralityBinTraitsHelper<I-1,Types...>::clone(ret, tup);
  }

  static void add(Tuple & tup, const Tuple & otup) {
    CentralityBinTraits<T>::add(get<I-1>(tup),get<I-1>(otup));
    TupleCentralityBinTraitsHelper<I-1,Types...>::add(tup, otup);
  }

  static void scale(Tuple & tup, double f) {
    CentralityBinTraits<T>::scale(get<I-1>(tup), f);
    TupleCentralityBinTraitsHelper<I-1,Types...>::scale(tup, f);
  }

  static void normalize(Tuple & tup, double sumw) {
    CentralityBinTraits<T>::normalize(get<I-1>(tup), sumw);
    TupleCentralityBinTraitsHelper<I-1,Types...>::normalize(tup, sumw);
  }

  static string path(const Tuple & tup) {
    return " " + CentralityBinTraits<T>::path(get<I-1>(tup))
      + TupleCentralityBinTraitsHelper<I-1,Types...>::path(tup);
  }
};

template <typename... Types>
struct TupleCentralityBinTraitsHelper<0,Types...> {

  typedef tuple<Types...> Tuple;

  static void clone(Tuple &, const Tuple &) {}
  static void add(Tuple & tup, const Tuple & otup) {}
  static void scale(Tuple & tup, double f) {}
  static void normalize(Tuple & tup, double sumw) {}
  static string path(const Tuple & tup) {return "";}

};

template <typename... Types>
struct CentralityBinTraits< tuple<Types...> > {

  typedef tuple<Types...> Tuple;
  static const size_t N = tuple_size<Tuple>::value;

  /// Make a clone of the given object.
  static Tuple clone(const Tuple & tup) {
    Tuple ret;
    TupleCentralityBinTraitsHelper<N,Types...>::clone(ret, tup);
    return ret;
  }

  /// Add the contents of @a o to @a t.
  static void add(Tuple & tup, const Tuple & otup) {
    TupleCentralityBinTraitsHelper<N,Types...>::add(tup, otup);
  }

  /// Scale the contents of a given object.
  static void scale(Tuple & tup, double f) {
    TupleCentralityBinTraitsHelper<N,Types...>::scale(tup, f);
  }

  static void normalize(Tuple & tup, double sumw) {
    TupleCentralityBinTraitsHelper<N,Types...>::normalize(tup, sumw);
  }

  /// Return the path of an AnalysisObject.
  static string path(const Tuple & tup) {
    string ret = "(tuple:";
    ret += TupleCentralityBinTraitsHelper<N,Types...>::path(tup);
    ret += ")";
    return ret;
  }

};

template <typename T, typename MDist>
T CentralityBinner<T,MDist>::select(double cest, double weight) {
  _currenT = _devnull;
  _currentCEst = cest;
  _weightsum += weight;

  // If estimator is negative, something has gone wrong.
  if ( _currentCEst < 0.0 ) return _currenT;

  // If we already have finalized the limits on the centrality
  // estimator, we just add the weights to their bins and return the
  // corresponding AnalysisObject.
  if ( _unfilled.empty() ) {
    for ( auto & b : _ready ) if ( b.second.inRange(_currentCEst) ) {
	b.second._weightsum += weight;
	return b.second._t;
      }
    return _currenT;
  }

  auto it = _findBin(cest);
  if ( it == _flexiBins.end() ) {
    _currenT = CentralityBinTraits<T>::clone(_unfilled.begin()->_t);
    it = _flexiBins.insert(FlexiBin(_currenT, _currentCEst, weight)).first;
  } else {
    it->_weightsum += weight;
    ++(it->_n);
    _currenT = it->_t;
  }

  if ( (int)_flexiBins.size() <= _maxBins ) return _currenT;


  set<double>::iterator citn = _percentiles.begin();
  set<double>::iterator cit0 = citn++;
  auto selectit = _flexiBins.end();
  double mindist = -1.0;
  double acc = 0.0;
  auto next = _flexiBins.begin();
  auto prev = next++;
  for ( ; next != _flexiBins.end(); prev = next++ ) {
    acc += prev->_weightsum/_weightsum;
    if ( acc > *citn ) {
      cit0 = citn++;
      continue;
    }
    if ( acc + next->_weightsum/_weightsum > *citn ) continue;
    double dist = MDist::dist(prev->_cestLo, next->_cestHi,
                              next->_weightsum + prev->_weightsum,
                              *cit0, *citn, next->_n + prev->_n,
                              next->_m + prev->_m);
    if ( mindist < 0.0 || dist < mindist ) {
      selectit = prev;
      mindist = dist;
    }
  }

  if ( selectit == _flexiBins.end() ) return _currenT;
  auto mergeit = selectit++;
  FlexiBin merged = *mergeit;
  merged.merge(*selectit);
  if ( merged.inRange(cest) || selectit->inRange(cest) )
    _currenT = merged._t;
  _flexiBins.erase(mergeit);
  _flexiBins.erase(selectit);
  _flexiBins.insert(merged);

  return _currenT;

}


template <typename T, typename MDist>
void CentralityBinner<T,MDist>::finalize() {

  // Take the contents of the dynamical binning and fill the original
  // AnalysisObjects.

  double clo = 0.0;
  for ( const FlexiBin & fb : _flexiBins ) {
    double chi = min(clo + fb._weightsum/_weightsum, 1.0);
    for ( Bin & bin : _unfilled ) {
      double olo = bin._centLo;
      double ohi = bin._centHi;
      if ( clo > ohi || chi <= olo ) continue;
      // If we only have partial overlap we need to scale
      double lo = max(olo, clo);
      double hi = min(ohi, chi);
      T t = CentralityBinTraits<T>::clone(fb._t);
      double frac = (hi - lo)/(chi - clo);
      CentralityBinTraits<T>::scale(t, frac);
      CentralityBinTraits<T>::add(bin._t, t);
      bin._weightsum += fb._weightsum*frac;
      if ( clo <= olo ) bin._cestLo = fb._cestLo +
                          (fb._cestHi - fb._cestLo)*(olo - clo)/(chi - clo);
      if ( clo < olo ) {
        bin._underflow = clo;
        bin._ambiguous += fb._n*frac;
        bin._ambweight += fb._weightsum*frac*(1.0 - frac);
      }
      if ( chi > ohi ) {
        bin._cestHi =
          fb._cestLo + (fb._cestHi - fb._cestLo)*(ohi - clo)/(chi - clo);
	bin._overflow = chi;
        bin._ambiguous += fb._n*frac;
        bin._ambweight += fb._weightsum*frac*(1.0 - frac);
      }
    }
    clo = chi;
  }
  _flexiBins.clear();
  for ( Bin & bin : _unfilled ) {
    if ( bin._overflow == 0.0 ) bin._overflow = 1.0;
    _ready[bin._t] = bin;
    if ( bin._ambweight/bin._weightsum >_warnlimit )
      MSG_WARNING("Analysis object \"" << CentralityBinTraits<T>::path(bin._t)
                  << "\", contains events with centralities between "
		  << bin._underflow*100.0
		  << " and " << bin._overflow*100.0 << "% ("
                  << int(bin._ambiguous + 0.5)
                  << " ambiguous events with effectively "
                  << 100.0*bin._ambweight/bin._weightsum
                  << "% of the weights)."
		  << "Consider increasing the number of bins.");

  }
  _unfilled.clear();

}

template <typename T, typename MDist>
void CentralityBinner<T,MDist>::fulldebug() {
  cerr <<  endl;
  double acc = 0.0;
  set<double>::iterator citn = _percentiles.begin();
  set<double>::iterator cit0 = citn++;
  int i = 0;
  for ( auto it = _flexiBins.begin(); it != _flexiBins.end(); ) {
    ++i;
    auto curr = it++;
    double w = curr->_weightsum/_weightsum;
    acc += w;
    if ( curr == _flexiBins.begin() || it == _flexiBins.end() || acc > *citn )
      cerr << "*";
    else
      cerr << " ";
    if ( acc > *citn ) cit0 = citn++;
    cerr << setw(6) << i
         << setw(12) << acc - w
	 << setw(12) << acc
	 << setw(8) << curr->_n
	 << setw(8) << curr->_m
	 << setw(12) << curr->_cestLo
	 << setw(12) << curr->_cestHi << endl;
  }
  cerr << "Number of sampler bins: " << _flexiBins.size() << endl;
}

template <typename T, typename MDist>
void CentralityBinner<T,MDist>::debug() {
  cerr <<  endl;
  double acc = 0.0;
  int i = 0;
  set<double>::iterator citn = _percentiles.begin();
  set<double>::iterator cit0 = citn++;
  for ( auto it = _flexiBins.begin(); it != _flexiBins.end(); ) {
    auto curr = it++;
    ++i;
    double w = curr->_weightsum/_weightsum;
    acc += w;
    if ( curr == _flexiBins.begin() || it == _flexiBins.end() || acc > *citn ) {
      if ( acc > *citn ) cit0 = citn++;
      cerr << setw(6) << i
           << setw(12) << acc - w
	   << setw(12) << acc
	   << setw(8) << curr->_n
	   << setw(8) << curr->_m
	   << setw(12) << curr->_cestLo
	   << setw(12) << curr->_cestHi << endl;

    }
  }
  cerr << "Number of sampler bins: " << _flexiBins.size() << endl;
}

/// Example of CentralityEstimator projection that the generated
/// centrality as given in the GenHeavyIon object in HepMC3.
class GeneratedCentrality: public CentralityEstimator {

public:

  /// Constructor.
  GeneratedCentrality() {}

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(GeneratedCentrality);

protected:

  /// Perform the projection on the Event
  void project(const Event& e) {
    _estimate = -1.0;
#if HEPMC_VERSION_CODE >= 3000000
    const HepMC::HeavyIon * hi = e.genEvent()->heavy_ion();
    if ( hi ) _estimate = 100.0 - hi->centrality; // @TODO We don't really know how to interpret this number!
#endif
  }

  /// Compare projections
  int compare(const Projection& p) const {
    return mkNamedPCmp(p, "GeneratedCentrality");
  }

};



}

#endif
