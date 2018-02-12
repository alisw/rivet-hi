// -*- C++ -*-
#ifndef RIVET_JetSmearingFunctions_HH
#define RIVET_JetSmearingFunctions_HH

#include "Rivet/Jet.hh"
#include "Rivet/Tools/MomentumSmearingFunctions.hh"
#include "Rivet/Tools/ParticleSmearingFunctions.hh"
#include "Rivet/Tools/Random.hh"

namespace Rivet {


  /// @name Jet filtering, efficiency and smearing utils
  //@{

  /// @name Typedef for Jet smearing functions/functors
  typedef function<Jet(const Jet&)> JetSmearFn;

  /// @name Typedef for Jet efficiency functions/functors
  typedef function<double(const Jet&)> JetEffFn;



  /// Return a constant 0 given a Jet as argument
  inline double JET_EFF_ZERO(const Jet& p) { return 0; }
  /// Return a constant 1 given a Jet as argument
  inline double JET_EFF_ONE(const Jet& p) { return 1; }

  /// Take a Jet and return a constant efficiency
  struct JET_EFF_CONST {
    JET_EFF_CONST(double eff) : _eff(eff) {}
    double operator () (const Jet& )  const { return _eff; }
    double _eff;
  };


  /// Return 1 if the given Jet contains a b, otherwise 0
  inline double JET_BTAG_PERFECT(const Jet& j) { return j.bTagged() ? 1 : 0; }

  /// Return 1 if the given Jet contains a c, otherwise 0
  inline double JET_CTAG_PERFECT(const Jet& j) { return j.cTagged() ? 1 : 0; }


  /// @brief b-tagging efficiency functor, for more readable b-tag effs and mistag rates
  /// Note several constructors, allowing for optional specification of charm, tau, and light jet mistag rates
  struct JET_BTAG_EFFS {
    JET_BTAG_EFFS(double eff_b, double eff_light=0) : _eff_b(eff_b), _eff_c(-1), _eff_t(-1), _eff_l(eff_light) { }
    JET_BTAG_EFFS(double eff_b, double eff_c, double eff_light) : _eff_b(eff_b), _eff_c(eff_c), _eff_t(-1), _eff_l(eff_light) { }
    JET_BTAG_EFFS(double eff_b, double eff_c, double eff_tau, double eff_light) : _eff_b(eff_b), _eff_c(eff_c), _eff_t(eff_tau), _eff_l(eff_light) { }
    inline double operator () (const Jet& j) {
      if (j.bTagged()) return _eff_b;
      if (_eff_c >= 0 && j.cTagged()) return _eff_c;
      if (_eff_t >= 0 && j.tauTagged()) return _eff_t;
      return _eff_l;
    }
    double _eff_b, _eff_c, _eff_t, _eff_l;
  };


  /// Take a jet and return an unmodified copy
  /// @todo Modify constituent particle vectors for consistency
  /// @todo Set a null PseudoJet if the Jet is smeared?
  inline Jet JET_SMEAR_IDENTITY(const Jet& j) { return j; }


  /// @brief Functor for simultaneous efficiency-filtering and smearing of Jets
  ///
  /// A central element of the SmearedJets system
  ///
  /// @todo Include tagging efficiency functions?
  struct JetEffSmearFn {
    JetEffSmearFn(const JetSmearFn& s, const JetEffFn& e)
      : sfn(s), efn(e) {    }

    JetEffSmearFn(const JetEffFn& e, const JetSmearFn& s)
      : sfn(s), efn(e) {    }

    JetEffSmearFn(const JetSmearFn& s)
      : sfn(s), efn(JET_EFF_ONE) {    }

    JetEffSmearFn(const JetEffFn& e)
      : sfn(JET_SMEAR_IDENTITY), efn(e) {    }

    JetEffSmearFn(double eff)
      : JetEffSmearFn(JET_EFF_CONST(eff)) {    }

    /// Smear and calculate an efficiency for the given jet
    pair<Jet,double> operator() (const Jet& j) const {
      return make_pair(sfn(j), efn(j));
    }

    /// Compare to another, for use in the projection system
    int cmp(const JetEffSmearFn& other) const {
      // cout << "Eff hashes = " << get_address(efn) << "," << get_address(other.efn) << "; "
      //      << "smear hashes = " << get_address(sfn) << "," << get_address(other.sfn) << endl;
      if (get_address(sfn) == 0 || get_address(other.sfn) == 0) return UNDEFINED;
      if (get_address(efn) == 0 || get_address(other.efn) == 0) return UNDEFINED;
      return Rivet::cmp(get_address(sfn), get_address(other.sfn)) || Rivet::cmp(get_address(efn), get_address(other.efn));
    }

    /// Automatic conversion to a smearing function
    operator JetSmearFn () { return sfn; }
    /// Automatic conversion to an efficiency function
    /// @todo Ambiguity re. whether reco eff or a tagging efficiency...
    // operator JetEffFn () { return efn; }

    // Stored functions/functors
    JetSmearFn sfn;
    JetEffFn efn;
  };


  /// Return true if Jet @a j is chosen to survive a random efficiency selection
  template <typename FN>
  inline bool efffilt(const Jet& j, FN& feff) {
    return rand01() < feff(j);
  }

  /// A functor to return true if Jet @a j survives a random efficiency selection
  struct JetEffFilter {
    template <typename FN>
    JetEffFilter(const FN& feff) : _feff(feff) {}
    JetEffFilter(double eff) : JetEffFilter( [&](const Jet& j){return eff;} ) {}
    bool operator () (const Jet& j) const { return efffilt(j, _feff); }
  private:
    const JetEffFn _feff;
  };
  using jetEffFilter = JetEffFilter;

  //@}


}

#endif
