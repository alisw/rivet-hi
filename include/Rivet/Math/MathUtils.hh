// -*- C++ -*-
#ifndef RIVET_MathUtils_HH
#define RIVET_MathUtils_HH

#include "Rivet/Math/MathHeader.hh"
#include <type_traits>
#include <cassert>

namespace Rivet {


  /// @name Comparison functions for safe (floating point) equality tests
  //@{

  /// @brief Compare a number to zero
  ///
  /// This version for floating point types has a degree of fuzziness expressed
  /// by the absolute @a tolerance parameter, for floating point safety.
  template <typename NUM>
  inline typename std::enable_if<std::is_floating_point<NUM>::value, bool>::type
  isZero(NUM val, double tolerance=1e-8) {
    return fabs(val) < tolerance;
  }

  /// @brief Compare a number to zero
  ///
  /// SFINAE template specialisation for integers, since there is no FP
  /// precision issue.
  template <typename NUM>
  inline typename std::enable_if<std::is_integral<NUM>::value, bool>::type
  isZero(NUM val, double=1e-5) { //< NB. unused tolerance parameter for ints, still needs a default value!
    return val == 0;
  }


  /// @brief Compare two numbers for equality with a degree of fuzziness
  ///
  /// This version for floating point types (if any argument is FP) has a degree
  /// of fuzziness expressed by the fractional @a tolerance parameter, for
  /// floating point safety.
  template <typename N1, typename N2>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value &&
   (std::is_floating_point<N1>::value || std::is_floating_point<N2>::value), bool>::type
  fuzzyEquals(N1 a, N2 b, double tolerance=1e-5) {
    const double absavg = (std::abs(a) + std::abs(b))/2.0;
    const double absdiff = std::abs(a - b);
    const bool rtn = (isZero(a) && isZero(b)) || absdiff < tolerance*absavg;
    return rtn;
  }

  /// @brief Compare two numbers for equality with a degree of fuzziness
  ///
  /// Simpler SFINAE template specialisation for integers, since there is no FP
  /// precision issue.
  template <typename N1, typename N2>
  inline typename std::enable_if<
    std::is_integral<N1>::value && std::is_integral<N2>::value, bool>::type
    fuzzyEquals(N1 a, N2 b, double) { //< NB. unused tolerance parameter for ints, still needs a default value!
    return a == b;
  }


  /// @brief Compare two numbers for >= with a degree of fuzziness
  ///
  /// The @a tolerance parameter on the equality test is as for @c fuzzyEquals.
  template <typename N1, typename N2>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value, bool>::type
  fuzzyGtrEquals(N1 a, N2 b, double tolerance=1e-5) {
    return a > b || fuzzyEquals(a, b, tolerance);
  }


  /// @brief Compare two floating point numbers for <= with a degree of fuzziness
  ///
  /// The @a tolerance parameter on the equality test is as for @c fuzzyEquals.
  template <typename N1, typename N2>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value, bool>::type
  fuzzyLessEquals(N1 a, N2 b, double tolerance=1e-5) {
    return a < b || fuzzyEquals(a, b, tolerance);
  }

  //@}


  /// @name Ranges and intervals
  //@{

  /// Represents whether an interval is open (non-inclusive) or closed (inclusive).
  ///
  /// For example, the interval \f$ [0, \pi) \f$ is closed (an inclusive
  /// boundary) at 0, and open (a non-inclusive boundary) at \f$ \pi \f$.
  enum RangeBoundary { OPEN=0, SOFT=0, CLOSED=1, HARD=1 };

  /// @brief Determine if @a value is in the range @a low to @a high, for floating point numbers
  ///
  /// Interval boundary types are defined by @a lowbound and @a highbound.
  template <typename N1, typename N2, typename N3>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value && std::is_arithmetic<N3>::value, bool>::type
  inRange(N1 value, N2 low, N3 high,
          RangeBoundary lowbound=CLOSED, RangeBoundary highbound=OPEN) {
    if (lowbound == OPEN && highbound == OPEN) {
      return (value > low && value < high);
    } else if (lowbound == OPEN && highbound == CLOSED) {
      return (value > low && value <= high);
    } else if (lowbound == CLOSED && highbound == OPEN) {
      return (value >= low && value < high);
    } else { // if (lowbound == CLOSED && highbound == CLOSED) {
      return (value >= low && value <= high);
    }
  }

  /// @brief Determine if @a value is in the range @a low to @a high, for floating point numbers
  ///
  /// Interval boundary types are defined by @a lowbound and @a highbound.
  /// Closed intervals are compared fuzzily.
  template <typename N1, typename N2, typename N3>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value && std::is_arithmetic<N3>::value, bool>::type
  fuzzyInRange(N1 value, N2 low, N3 high,
               RangeBoundary lowbound=CLOSED, RangeBoundary highbound=OPEN) {
    if (lowbound == OPEN && highbound == OPEN) {
      return (value > low && value < high);
    } else if (lowbound == OPEN && highbound == CLOSED) {
      return (value > low && fuzzyLessEquals(value, high));
    } else if (lowbound == CLOSED && highbound == OPEN) {
      return (fuzzyGtrEquals(value, low) && value < high);
    } else { // if (lowbound == CLOSED && highbound == CLOSED) {
      return (fuzzyGtrEquals(value, low) && fuzzyLessEquals(value, high));
    }
  }

  /// Alternative version of inRange which accepts a pair for the range arguments.
  template <typename N1, typename N2, typename N3>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value && std::is_arithmetic<N3>::value, bool>::type
  inRange(N1 value, pair<N2, N3> lowhigh,
          RangeBoundary lowbound=CLOSED, RangeBoundary highbound=OPEN) {
    return inRange(value, lowhigh.first, lowhigh.second, lowbound, highbound);
  }


  // Alternative forms, with snake_case names and boundary types in names rather than as args -- from MCUtils

  /// @brief Boolean function to determine if @a value is within the given range
  ///
  /// @note The interval is closed (inclusive) at the low end, and open (exclusive) at the high end.
  template <typename N1, typename N2, typename N3>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value && std::is_arithmetic<N3>::value, bool>::type
  in_range(N1 val, N2 low, N3 high) {
    return inRange(val, low, high, CLOSED, OPEN);
  }

  /// @brief Boolean function to determine if @a value is within the given range
  ///
  /// @note The interval is closed at both ends.
  template <typename N1, typename N2, typename N3>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value && std::is_arithmetic<N3>::value, bool>::type
  in_closed_range(N1 val, N2 low, N3 high) {
    return inRange(val, low, high, CLOSED, CLOSED);
  }

  /// @brief Boolean function to determine if @a value is within the given range
  ///
  /// @note The interval is open at both ends.
  template <typename N1, typename N2, typename N3>
  inline typename std::enable_if<
    std::is_arithmetic<N1>::value && std::is_arithmetic<N2>::value && std::is_arithmetic<N3>::value, bool>::type
  in_open_range(N1 val, N2 low, N3 high) {
    return inRange(val, low, high, OPEN, OPEN);
  }

  /// @todo Add pair-based versions of the named range-boundary functions

  //@}


  /// @name Miscellaneous numerical helpers
  //@{

  /// Named number-type squaring operation.
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, NUM>::type
  sqr(NUM a) {
    return a*a;
  }

  /// @brief Named number-type addition in quadrature operation.
  ///
  /// @note Result has the sqrt operation applied.
  /// @todo When std::common_type can be used, generalise to multiple numeric types with appropriate return type.
  // template <typename N1, typename N2>
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, NUM>::type
  //std::common_type<N1, N2>::type
  add_quad(NUM a, NUM b) {
    return sqrt(a*a + b*b);
  }

  /// Named number-type addition in quadrature operation.
  ///
  /// @note Result has the sqrt operation applied.
  /// @todo When std::common_type can be used, generalise to multiple numeric types with appropriate return type.
  // template <typename N1, typename N2>
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, NUM>::type
  //std::common_type<N1, N2, N3>::type
  add_quad(NUM a, NUM b, NUM c) {
    return sqrt(a*a + b*b + c*c);
  }

  /// Return a/b, or @a fail if b = 0
  /// @todo When std::common_type can be used, generalise to multiple numeric types with appropriate return type.
  inline double safediv(double num, double den, double fail=0.0) {
    return (!isZero(den)) ? num/den : fail;
  }

  /// A more efficient version of pow for raising numbers to integer powers.
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, NUM>::type
  intpow(NUM val, unsigned int exp) {
    assert(exp >= 0);
    if (exp == 0) return (NUM) 1;
    else if (exp == 1) return val;
    return val * intpow(val, exp-1);
  }

  /// Find the sign of a number
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, int>::type
  sign(NUM val) {
    if (isZero(val)) return ZERO;
    const int valsign = (val > 0) ? PLUS : MINUS;
    return valsign;
  }

  //@}


  /// @name Physics statistical distributions
  //@{

  /// @brief CDF for the Breit-Wigner distribution
  inline double cdfBW(double x, double mu, double gamma) {
    // normalize to (0;1) distribution
    const double xn = (x - mu)/gamma;
    return std::atan(xn)/M_PI + 0.5;
  }

  /// @brief Inverse CDF for the Breit-Wigner distribution
  inline double invcdfBW(double p, double mu, double gamma) {
    const double xn = std::tan(M_PI*(p-0.5));
    return gamma*xn + mu;
  }

  //@}


  /// @name Binning helper functions
  //@{

  /// @brief Make a list of @a nbins + 1 values equally spaced between @a start and @a end inclusive.
  ///
  /// NB. The arg ordering and the meaning of the nbins variable is "histogram-like",
  /// as opposed to the Numpy/Matlab version.
  inline vector<double> linspace(size_t nbins, double start, double end, bool include_end=true) {
    assert(end >= start);
    assert(nbins > 0);
    vector<double> rtn;
    const double interval = (end-start)/static_cast<double>(nbins);
    for (size_t i = 0; i < nbins; ++i) {
      rtn.push_back(start + i*interval);
    }
    assert(rtn.size() == nbins);
    if (include_end) rtn.push_back(end); //< exact end, not result of n * interval
    return rtn;
  }


  /// @brief Make a list of @a nbins + 1 values exponentially spaced between @a start and @a end inclusive.
  ///
  /// NB. The arg ordering and the meaning of the nbins variable is "histogram-like",
  /// as opposed to the Numpy/Matlab version, and the start and end arguments are expressed
  /// in "normal" space, rather than as the logarithms of the start/end values as in Numpy/Matlab.
  inline vector<double> logspace(size_t nbins, double start, double end, bool include_end=true) {
    assert(end >= start);
    assert(start > 0);
    assert(nbins > 0);
    const double logstart = std::log(start);
    const double logend = std::log(end);
    const vector<double> logvals = linspace(nbins, logstart, logend, false);
    assert(logvals.size() == nbins);
    vector<double> rtn; rtn.reserve(nbins+1);
    rtn.push_back(start); //< exact start, not exp(log(start))
    for (size_t i = 1; i < logvals.size(); ++i) {
      rtn.push_back(std::exp(logvals[i]));
    }
    assert(rtn.size() == nbins);
    if (include_end) rtn.push_back(end); //< exact end, not exp(n * loginterval)
    return rtn;
  }


  /// @brief Make a list of @a nbins + 1 values spaced for equal area
  /// Breit-Wigner binning between @a start and @a end inclusive. @a
  /// mu and @a gamma are the Breit-Wigner parameters.
  ///
  /// @note The arg ordering and the meaning of the nbins variable is "histogram-like",
  /// as opposed to the Numpy/Matlab version, and the start and end arguments are expressed
  /// in "normal" space.
  inline vector<double> bwspace(size_t nbins, double start, double end, double mu, double gamma) {
    assert(end >= start);
    assert(nbins > 0);
    const double pmin = cdfBW(start, mu, gamma);
    const double pmax = cdfBW(end,   mu, gamma);
    const vector<double> edges = linspace(nbins, pmin, pmax);
    assert(edges.size() == nbins+1);
    vector<double> rtn;
    for (double edge : edges) {
      rtn.push_back(invcdfBW(edge, mu, gamma));
    }
    assert(rtn.size() == nbins+1);
    return rtn;
  }


  /// @brief Return the bin index of the given value, @a val, given a vector of bin edges
  ///
  /// An underflow always returns -1. If allow_overflow is false (default) an overflow
  /// also returns -1, otherwise it returns the Nedge-1, the index of an inclusive bin
  /// starting at the last edge.
  ///
  /// @note The @a binedges vector must be sorted
  /// @todo Use std::common_type<NUM1, NUM2>::type x = val; ?
  template <typename NUM, typename CONTAINER>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value && std::is_arithmetic<typename CONTAINER::value_type>::value, int>::type
  binIndex(NUM val, const CONTAINER& binedges, bool allow_overflow=false) {
    if (val < binedges.front()) return -1; ///< Below/out of histo range
    if (val >= binedges.back()) return allow_overflow ? int(binedges.size())-1 : -1; ///< Above/out of histo range
    auto it = std::upper_bound(binedges.begin(), binedges.end(), val);
    return std::distance(binedges.begin(), --it);
  }

  //@}


  /// @name Discrete statistics functions
  //@{

  /// Calculate the median of a sample
  /// @todo Support multiple container types via SFINAE
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, NUM>::type
  median(const vector<NUM>& sample) {
    if (sample.empty()) throw RangeError("Can't compute median of an empty set");
    vector<NUM> tmp = sample;
    std::sort(tmp.begin(), tmp.end());
    const size_t imid = tmp.size()/2; // len1->idx0, len2->idx1, len3->idx1, len4->idx2, ...
    if (sample.size() % 2 == 0) return (tmp.at(imid-1) + tmp.at(imid)) / 2.0;
    else return tmp.at(imid);
  }


  /// Calculate the mean of a sample
  /// @todo Support multiple container types via SFINAE
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, double>::type
  mean(const vector<NUM>& sample) {
    if (sample.empty()) throw RangeError("Can't compute mean of an empty set");
    double mean = 0.0;
    for (size_t i = 0; i < sample.size(); ++i) {
      mean += sample[i];
    }
    return mean/sample.size();
  }

  // Calculate the error on the mean, assuming Poissonian errors
  /// @todo Support multiple container types via SFINAE
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, double>::type
  mean_err(const vector<NUM>& sample) {
    if (sample.empty()) throw RangeError("Can't compute mean_err of an empty set");
    double mean_e = 0.0;
    for (size_t i = 0; i < sample.size(); ++i) {
      mean_e += sqrt(sample[i]);
    }
    return mean_e/sample.size();
  }


  /// Calculate the covariance (variance) between two samples
  /// @todo Support multiple container types via SFINAE
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, double>::type
  covariance(const vector<NUM>& sample1, const vector<NUM>& sample2) {
    if (sample1.empty() || sample2.empty()) throw RangeError("Can't compute covariance of an empty set");
    if (sample1.size() != sample2.size()) throw RangeError("Sizes of samples must be equal for covariance calculation");
    const double mean1 = mean(sample1);
    const double mean2 = mean(sample2);
    const size_t N = sample1.size();
    double cov = 0.0;
    for (size_t i = 0; i < N; i++) {
      const double cov_i = (sample1[i] - mean1)*(sample2[i] - mean2);
      cov += cov_i;
    }
    if (N > 1) return cov/(N-1);
    else return 0.0;
  }

  /// Calculate the error on the covariance (variance) of two samples, assuming poissonian errors
  /// @todo Support multiple container types via SFINAE
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, double>::type
  covariance_err(const vector<NUM>& sample1, const vector<NUM>& sample2) {
    if (sample1.empty() || sample2.empty()) throw RangeError("Can't compute covariance_err of an empty set");
    if (sample1.size() != sample2.size()) throw RangeError("Sizes of samples must be equal for covariance_err calculation");
    const double mean1 = mean(sample1);
    const double mean2 = mean(sample2);
    const double mean1_e = mean_err(sample1);
    const double mean2_e = mean_err(sample2);
    const size_t N = sample1.size();
    double cov_e = 0.0;
    for (size_t i = 0; i < N; i++) {
      const double cov_i = (sqrt(sample1[i]) - mean1_e)*(sample2[i] - mean2) +
        (sample1[i] - mean1)*(sqrt(sample2[i]) - mean2_e);
      cov_e += cov_i;
    }
    if (N > 1) return cov_e/(N-1);
    else return 0.0;
  }


  /// Calculate the correlation strength between two samples
  /// @todo Support multiple container types via SFINAE
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, double>::type
  correlation(const vector<NUM>& sample1, const vector<NUM>& sample2) {
    const double cov = covariance(sample1, sample2);
    const double var1 = covariance(sample1, sample1);
    const double var2 = covariance(sample2, sample2);
    const double correlation = cov/sqrt(var1*var2);
    const double corr_strength = correlation*sqrt(var2/var1);
    return corr_strength;
  }

  /// Calculate the error of the correlation strength between two samples assuming Poissonian errors
  /// @todo Support multiple container types via SFINAE
  template <typename NUM>
  inline typename std::enable_if<std::is_arithmetic<NUM>::value, double>::type
  correlation_err(const vector<NUM>& sample1, const vector<NUM>& sample2) {
    const double cov = covariance(sample1, sample2);
    const double var1 = covariance(sample1, sample1);
    const double var2 = covariance(sample2, sample2);
    const double cov_e = covariance_err(sample1, sample2);
    const double var1_e = covariance_err(sample1, sample1);
    const double var2_e = covariance_err(sample2, sample2);

    // Calculate the correlation
    const double correlation = cov/sqrt(var1*var2);
    // Calculate the error on the correlation
    const double correlation_err = cov_e/sqrt(var1*var2) -
      cov/(2*pow(3./2., var1*var2)) * (var1_e * var2 + var1 * var2_e);

    // Calculate the error on the correlation strength
    const double corr_strength_err = correlation_err*sqrt(var2/var1) +
      correlation/(2*sqrt(var2/var1)) * (var2_e/var1 - var2*var1_e/pow(2, var2));

    return corr_strength_err;
  }

  //@}


  /// @name Angle range mappings
  //@{

  /// @brief Reduce any number to the range [-2PI, 2PI]
  ///
  /// Achieved by repeated addition or subtraction of 2PI as required. Used to
  /// normalise angular measures.
  inline double _mapAngleM2PITo2Pi(double angle) {
    double rtn = fmod(angle, TWOPI);
    if (isZero(rtn)) return 0;
    assert(rtn >= -TWOPI && rtn <= TWOPI);
    return rtn;
  }

  /// Map an angle into the range (-PI, PI].
  inline double mapAngleMPiToPi(double angle) {
    double rtn = _mapAngleM2PITo2Pi(angle);
    if (isZero(rtn)) return 0;
    if (rtn > PI) rtn -= TWOPI;
    if (rtn <= -PI) rtn += TWOPI;
    assert(rtn > -PI && rtn <= PI);
    return rtn;
  }

  /// Map an angle into the range [0, 2PI).
  inline double mapAngle0To2Pi(double angle) {
    double rtn = _mapAngleM2PITo2Pi(angle);
    if (isZero(rtn)) return 0;
    if (rtn < 0) rtn += TWOPI;
    if (rtn == TWOPI) rtn = 0;
    assert(rtn >= 0 && rtn < TWOPI);
    return rtn;
  }

  /// Map an angle into the range [0, PI].
  inline double mapAngle0ToPi(double angle) {
    double rtn = fabs(mapAngleMPiToPi(angle));
    if (isZero(rtn)) return 0;
    assert(rtn > 0 && rtn <= PI);
    return rtn;
  }

  /// Map an angle into the enum-specified range.
  inline double mapAngle(double angle, PhiMapping mapping) {
    switch (mapping) {
    case MINUSPI_PLUSPI:
      return mapAngleMPiToPi(angle);
    case ZERO_2PI:
      return mapAngle0To2Pi(angle);
    case ZERO_PI:
      return mapAngle0To2Pi(angle);
    default:
      throw Rivet::UserError("The specified phi mapping scheme is not implemented");
    }
  }

  //@}


  /// @name Phase-space measure helpers
  //@{

  /// @brief Calculate the difference between two angles in radians
  ///
  /// Returns in the range [0, PI].
  inline double deltaPhi(double phi1, double phi2) {
    return mapAngle0ToPi(phi1 - phi2);
  }

  /// Calculate the abs difference between two pseudorapidities
  ///
  /// @note Just a cosmetic name for analysis code clarity.
  inline double deltaEta(double eta1, double eta2) {
    return fabs(eta1 - eta2);
  }

  /// Calculate the abs difference between two rapidities
  ///
  /// @note Just a cosmetic name for analysis code clarity.
  inline double deltaRap(double y1, double y2) {
    return fabs(y1 - y2);
  }

  /// Calculate the squared distance between two points in 2D rapidity-azimuthal
  /// ("\f$ \eta-\phi \f$") space. The phi values are given in radians.
  inline double deltaR2(double rap1, double phi1, double rap2, double phi2) {
    const double dphi = deltaPhi(phi1, phi2);
    return sqr(rap1-rap2) + sqr(dphi);
  }

  /// Calculate the distance between two points in 2D rapidity-azimuthal
  /// ("\f$ \eta-\phi \f$") space. The phi values are given in radians.
  inline double deltaR(double rap1, double phi1, double rap2, double phi2) {
    return sqrt(deltaR2(rap1, phi1, rap2, phi2));
  }

  /// Calculate a rapidity value from the supplied energy @a E and longitudinal momentum @a pz.
  inline double rapidity(double E, double pz) {
    if (isZero(E - pz)) {
      throw std::runtime_error("Divergent positive rapidity");
      return MAXDOUBLE;
    }
    if (isZero(E + pz)) {
      throw std::runtime_error("Divergent negative rapidity");
      return -MAXDOUBLE;
    }
    return 0.5*log((E+pz)/(E-pz));
  }

  //@}


  /// Calculate transverse mass of two vectors from provided pT and deltaPhi
  /// @note Several versions taking two vectors are found in Vector4.hh
  inline double mT(double pT1, double pT2, double dphi) {
    return sqrt(2*pT1*pT2 * (1 - cos(dphi)) );
  }


}


#endif
