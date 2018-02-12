#include "Rivet/Tools/Cuts.hh"
#include "Rivet/Particle.hh"
#include "Rivet/Jet.hh"
#include "Rivet/Math/Vectors.hh"
#include "fastjet/PseudoJet.hh"
#include "HepMC/SimpleVector.h"

/// @todo Identify what can go into anonymous namespace

namespace Rivet {


  // Base for all wrapper classes that translate ClassToCheck to Cuttable
  class CuttableBase {
  public:
    virtual double getValue(Cuts::Quantity) const = 0;
    virtual ~CuttableBase() {}
  };


  // Cuttables can be directly passed to @ref _accept
  template <>
  bool CutBase::accept<CuttableBase>(const CuttableBase& t) const {
    return _accept(t);
  }


  // Open Cut singleton
  class Open_Cut : public CutBase {
  public:
    bool operator==(const Cut& c) const {
      std::shared_ptr<Open_Cut> cc = dynamic_pointer_cast<Open_Cut>(c);
      return bool(cc);
    }
  protected:
    // open cut accepts everything
    bool _accept(const CuttableBase&) const { return true; }
  };


  const Cut& Cuts::open() {
    // Only ever need one static open cut object
    static const Cut open = std::make_shared<Open_Cut>();
    return open;
  }

  // Constants for convenient access
  const Cut& Cuts::OPEN = Cuts::open();
  const Cut& Cuts::NOCUT = Cuts::open();


  // Cut constructor for ==
  class Cut_Eq : public CutBase {
  public:
    Cut_Eq(const Cuts::Quantity qty, double val) : _qty(qty), _val(val) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<Cut_Eq> cc = dynamic_pointer_cast<Cut_Eq>(c);
      return cc  &&  _qty == cc->_qty  &&  _val == cc->_val;
    }
  protected:
    bool _accept(const CuttableBase& o) const { return o.getValue(_qty) == _val; }
  private:
    Cuts::Quantity _qty;
    int _val;
  };


  // Cut constructor for !=
  class Cut_NEq : public CutBase {
  public:
    Cut_NEq(const Cuts::Quantity qty, double val) : _qty(qty), _val(val) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<Cut_NEq> cc = dynamic_pointer_cast<Cut_NEq>(c);
      return cc  &&  _qty == cc->_qty  &&  _val == cc->_val;
    }
  protected:
    bool _accept(const CuttableBase& o) const { return o.getValue(_qty) != _val; }
  private:
    Cuts::Quantity _qty;
    int _val;
  };


  // Cut constructor for >=
  class Cut_GtrEq : public CutBase {
  public:
    Cut_GtrEq(const Cuts::Quantity qty, double val) : _qty(qty), _val(val) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<Cut_GtrEq> cc = dynamic_pointer_cast<Cut_GtrEq>(c);
      return cc && _qty == cc->_qty  &&  _val == cc->_val;
    }
  protected:
    bool _accept(const CuttableBase& o) const { return o.getValue(_qty) >= _val; }
  private:
    Cuts::Quantity _qty;
    double _val;
  };


  // Cut constructor for <
  class Cut_Less : public CutBase {
  public:
    Cut_Less(const Cuts::Quantity qty, const double val) : _qty(qty), _val(val) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<Cut_Less> cc = dynamic_pointer_cast<Cut_Less>(c);
      return cc  && _qty == cc->_qty  &&  _val == cc->_val;
    }
  protected:
    bool _accept(const CuttableBase& o) const { return o.getValue(_qty) < _val; }
  private:
    Cuts::Quantity _qty;
    double _val;
  };


  // Cut constructor for >
  class Cut_Gtr : public CutBase {
  public:
    Cut_Gtr(const Cuts::Quantity qty, const double val) : _qty(qty), _val(val) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<Cut_Gtr> cc = dynamic_pointer_cast<Cut_Gtr>(c);
      return cc && _qty == cc->_qty  &&  _val == cc->_val;
    }
  protected:
    bool _accept(const CuttableBase& o) const { return o.getValue(_qty) > _val; }
  private:
    Cuts::Quantity _qty;
    double _val;
  };


  // Cut constructor for <=
  class Cut_LessEq : public CutBase {
  public:
    Cut_LessEq(const Cuts::Quantity qty, const double val) : _qty(qty), _val(val) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<Cut_LessEq> cc = dynamic_pointer_cast<Cut_LessEq>(c);
      return cc && _qty == cc->_qty  &&  _val == cc->_val;
    }
  protected:
    bool _accept(const CuttableBase& o) const { return o.getValue(_qty) <= _val; }
  private:
    Cuts::Quantity _qty;
    double _val;
  };


  template <typename T>
  inline Cut make_cut(T t) {
    return std::make_shared<T>(t);
  }

  Cut operator == (Cuts::Quantity qty, double n) {
    return make_cut(Cut_Eq(qty, n));
  }

  Cut operator != (Cuts::Quantity qty, double n) {
    return make_cut(Cut_NEq(qty, n));
  }

  Cut operator < (Cuts::Quantity qty, double n) {
    return make_cut(Cut_Less(qty, n));
  }

  Cut operator >= (Cuts::Quantity qty, double n) {
    return make_cut(Cut_GtrEq(qty, n));
  }

  Cut operator <= (Cuts::Quantity qty, double n) {
    return make_cut(Cut_LessEq(qty, n));
  }

  Cut operator > (Cuts::Quantity qty, double n) {
    return make_cut(Cut_Gtr(qty, n));
  }

  Cut Cuts::range(Cuts::Quantity qty, double m, double n) {
    if (m > n) std::swap(m,n);
    return (qty >= m) && (qty < n);
  }


  //////////////
  /// Combiners

  /// AND, OR, NOT, and XOR objects for combining cuts

  class CutsOr : public CutBase {
  public:
    CutsOr(const Cut& c1, const Cut& c2) : cut1(c1), cut2(c2) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<CutsOr> cc = dynamic_pointer_cast<CutsOr>(c);
      return cc && (   ( cut1 == cc->cut1  &&  cut2 == cc->cut2 )
                       || ( cut1 == cc->cut2  &&  cut2 == cc->cut1 ));
    }
  protected:
    bool _accept(const CuttableBase& o) const {
      return cut1->accept(o) || cut2->accept(o);
    }
  private:
    const Cut cut1;
    const Cut cut2;
  };


  class CutsAnd : public CutBase {
  public:
    CutsAnd(const Cut& c1, const Cut& c2) : cut1(c1), cut2(c2) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<CutsAnd> cc = dynamic_pointer_cast<CutsAnd>(c);
      return cc && (   ( cut1 == cc->cut1  &&  cut2 == cc->cut2 )
                       || ( cut1 == cc->cut2  &&  cut2 == cc->cut1 ));
    }
  protected:
    bool _accept(const CuttableBase& o) const {
      return cut1->accept(o) && cut2->accept(o);
    }
  private:
    const Cut cut1;
    const Cut cut2;
  };


  class CutInvert : public CutBase {
  public:
    CutInvert(const Cut& c1) : cut(c1) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<CutInvert> cc = dynamic_pointer_cast<CutInvert>(c);
      return cc && cut == cc->cut;
    }
  protected:
    bool _accept(const CuttableBase& o) const {
      return !cut->accept(o);
    }
  private:
    const Cut cut;
  };


  class CutsXor : public CutBase {
  public:
    CutsXor(const Cut& c1, const Cut& c2) : cut1(c1), cut2(c2) {}
    bool operator==(const Cut& c) const {
      std::shared_ptr<CutsXor> cc = dynamic_pointer_cast<CutsXor>(c);
      return cc && (   ( cut1 == cc->cut1  &&  cut2 == cc->cut2 )
                       || ( cut1 == cc->cut2  &&  cut2 == cc->cut1 ));
    }
  protected:
    bool _accept(const CuttableBase& o) const {
      bool A_and_B = cut1->accept(o) && cut2->accept(o);
      bool A_or_B  = cut1->accept(o) || cut2->accept(o);
      return A_or_B && (! A_and_B);
    }
  private:
    const Cut cut1;
    const Cut cut2;
  };


  ////////////
  ///Operators

  Cut operator && (const Cut& aptr, const Cut& bptr) {
    return make_cut(CutsAnd(aptr, bptr));
  }

  Cut operator || (const Cut& aptr, const Cut& bptr) {
    return make_cut(CutsOr(aptr, bptr));
  }

  Cut operator ! (const Cut& cptr) {
    return make_cut(CutInvert(cptr));
  }


  Cut operator & (const Cut& aptr, const Cut& bptr) {
    return make_cut(CutsAnd(aptr, bptr));
  }

  Cut operator | (const Cut& aptr, const Cut& bptr) {
    return make_cut(CutsOr(aptr, bptr));
  }

  Cut operator ~ (const Cut& cptr) {
    return make_cut(CutInvert(cptr));
  }

  Cut operator ^ (const Cut& aptr, const Cut& bptr) {
    return make_cut(CutsXor(aptr, bptr));
  }


  ///////////////////////
  /// Cuts

  // Non-functional Cuttable template class. Must be specialized below.
  template <typename T>
  class Cuttable;

  // Non-cuttables need to be wrapped into a Cuttable first
  #define SPECIALISE_ACCEPT(TYPENAME)                           \
    template <>                                                 \
    bool CutBase::accept<TYPENAME>(const TYPENAME& t) const {   \
      return _accept(Cuttable<TYPENAME>(t));                    \
    }                                                           \


  inline void qty_not_found() {
    throw Exception("Missing implementation for a Cuts::Quantity.");
  }


  template<>
  class Cuttable<Particle> : public CuttableBase {
  public:
    Cuttable(const Particle& p) : p_(p) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:         return p_.pT();
      case Cuts::Et:         return p_.Et();
      case Cuts::E:          return p_.E();
      case Cuts::mass:       return p_.mass();
      case Cuts::rap:        return p_.rap();
      case Cuts::absrap:     return p_.absrap();
      case Cuts::eta:        return p_.eta();
      case Cuts::abseta:     return p_.abseta();
      case Cuts::phi:        return p_.phi();
      case Cuts::pid:        return p_.pid();
      case Cuts::abspid:     return p_.abspid();
      case Cuts::charge:     return p_.charge();
      case Cuts::abscharge:  return p_.abscharge();
      case Cuts::charge3:    return p_.charge3();
      case Cuts::abscharge3: return p_.abscharge3();
      default: qty_not_found();
      }
      return -999.;
    }

  private:
    const Particle& p_;
  };

  SPECIALISE_ACCEPT(Particle)


  template<>
  class Cuttable<FourMomentum> : public CuttableBase {
  public:
    Cuttable(const FourMomentum& fm) : fm_(fm) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return fm_.pT();
      case Cuts::Et:     return fm_.Et();
      case Cuts::E:      return fm_.E();
      case Cuts::mass:   return fm_.mass();
      case Cuts::rap:    return fm_.rap();
      case Cuts::absrap: return fm_.absrap();
      case Cuts::eta:    return fm_.eta();
      case Cuts::abseta: return fm_.abseta();
      case Cuts::phi:    return fm_.phi();
      default: qty_not_found();
      }
      return -999.;
    }

  private:
    const FourMomentum& fm_;
  };

  SPECIALISE_ACCEPT(FourMomentum)


  template<>
  class Cuttable<Jet> : public CuttableBase {
  public:
    Cuttable(const Jet& jet) : jet_(jet) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return jet_.pT();
      case Cuts::Et:     return jet_.Et();
      case Cuts::E:      return jet_.E();
      case Cuts::mass:   return jet_.mass();
      case Cuts::rap:    return jet_.rapidity();
      case Cuts::absrap: return std::abs(jet_.rapidity());
      case Cuts::eta:    return jet_.pseudorapidity();
      case Cuts::abseta: return std::abs(jet_.pseudorapidity());
      case Cuts::phi:    return jet_.phi();
      default: qty_not_found();
      }
      return -999.;
    }

  private:
    const Jet& jet_;
  };

  SPECIALISE_ACCEPT(Jet)


  template<>
  class Cuttable<fastjet::PseudoJet> : public CuttableBase {
  public:
    Cuttable(const fastjet::PseudoJet& pjet) : pjet_(pjet) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return pjet_.perp();
      case Cuts::Et:     return pjet_.Et();
      case Cuts::mass:   return pjet_.m();
      case Cuts::rap:    return pjet_.rap();
      case Cuts::absrap: return std::abs(pjet_.rap());
      case Cuts::eta:    return pjet_.eta();
      case Cuts::abseta: return std::abs(pjet_.eta());
      case Cuts::phi:    return pjet_.phi();
      default: qty_not_found();
      }
      return -999.;
    }

  private:
    const fastjet::PseudoJet& pjet_;
  };

  SPECIALISE_ACCEPT(fastjet::PseudoJet)


  template<>
  class Cuttable<HepMC::FourVector> : public CuttableBase {
  public:
    Cuttable(const HepMC::FourVector& vec) : vec_(vec) {}
    double getValue(Cuts::Quantity qty) const {
      switch ( qty ) {
      case Cuts::pT:     return vec_.perp();
      case Cuts::E:      return vec_.e();
      case Cuts::Et:     return vec_.e() * sin(vec_.theta());
      case Cuts::mass:   return vec_.m();
      case Cuts::rap:    return 0.5*std::log((vec_.t()+vec_.z())/(vec_.t()-vec_.z()));
      case Cuts::absrap: return std::abs(getValue(Cuts::rap));
      case Cuts::eta:    return vec_.pseudoRapidity();
      case Cuts::abseta: return std::abs(vec_.pseudoRapidity());
      case Cuts::phi:    return vec_.phi();
      default: qty_not_found();
      }
      return -999.;
    }

  private:
    const HepMC::FourVector& vec_;
  };

  SPECIALISE_ACCEPT(HepMC::FourVector)


}
