#ifndef RIVET_Cuts_HH
#define RIVET_Cuts_HH

#include "Rivet/Tools/Cuts.fhh"
#include <memory>

namespace Rivet {


  class CutBase {
  public:

    /// Main work method, checking whether the cut is passed
    /// @internal Forwards the received object to @ref accept_, wrapped in the Cuttable converter
    template <typename ClassToCheck>
    bool accept(const ClassToCheck&) const;

    /// @brief Call operator alias for @a accept
    /// @note A bit subtle, because this gets wrapped in a shared_ptr so you need to dereference to get the functor
    template <typename ClassToCheck>
    bool operator () (const ClassToCheck& x) const { return accept(x); }

    /// Comparison to another Cut
    virtual bool operator == (const Cut&) const = 0;

    /// Default destructor
    virtual ~CutBase() {}

  protected:

    /// @internal Actual accept implementation, overloadable by various cut combiners
    virtual bool _accept(const CuttableBase&) const = 0;

  };


  /// Compare two cuts for equality, forwards to the cut-specific implementation
  inline bool operator == (const Cut& a, const Cut& b) { return *a == b; }


  /// Namespace used for ambiguous identifiers.
  namespace Cuts {

    /// Available categories of cut objects
    enum Quantity { pT=0, pt=0, Et=1, et=1, E=2, energy=2,
                    mass, rap, absrap, eta, abseta, phi,
                    pid, abspid, charge, abscharge, charge3, abscharge3 };

    /// Fully open cut singleton, accepts everything
    const Cut& open(); //< access by factory function

    extern const Cut& OPEN; //= open(); //< access by constant
    extern const Cut& NOCUT; //= open(); //< access by constant

    /// @name Shortcuts for common cuts, using the Quantity enums defined above
    //@{
    Cut range(Quantity, double m, double n);
    inline Cut ptIn(double m, double n) { return range(pT, m,n); }
    inline Cut etIn(double m, double n) { return range(Et, m,n); }
    inline Cut energyIn(double m, double n) { return range(energy, m,n); }
    inline Cut massIn(double m, double n) { return range(mass, m,n); }
    inline Cut rapIn(double m, double n) { return range(rap, m,n); }
    inline Cut absrapIn(double m, double n) { return range(absrap, m,n); }
    inline Cut etaIn(double m, double n) { return range(eta, m,n); }
    inline Cut absetaIn(double m, double n) { return range(abseta, m,n); }
    //@}

  }


  /// @name Cut constructors
  //@{
  Cut operator == (Cuts::Quantity, double);
  Cut operator != (Cuts::Quantity, double);
  Cut operator <  (Cuts::Quantity, double);
  Cut operator >  (Cuts::Quantity, double);
  Cut operator <= (Cuts::Quantity, double);
  Cut operator >= (Cuts::Quantity, double);

  /// @internal Overload helpers for integer arguments
  //@{
  inline Cut operator == (Cuts::Quantity qty, int i) { return qty ==  double(i); }
  inline Cut operator != (Cuts::Quantity qty, int i) { return qty !=  double(i); }
  // Cut operator == (Cuts::Quantity qty, int i);
  // Cut operator != (Cuts::Quantity qty, int i);
  inline Cut operator <  (Cuts::Quantity qty, int i) { return qty <  double(i); }
  inline Cut operator >  (Cuts::Quantity qty, int i) { return qty >  double(i); }
  inline Cut operator <= (Cuts::Quantity qty, int i) { return qty <= double(i); }
  inline Cut operator >= (Cuts::Quantity qty, int i) { return qty >= double(i); }
  //@}

  //@}


  /// @name Cut combiners
  //@{

  /// Logical AND operation on two cuts
  /// @note No comparison short-circuiting for overloaded &&!
  Cut operator && (const Cut & aptr, const Cut & bptr);
  /// Logical OR operation on two cuts
  /// @note No comparison short-circuiting for overloaded ||!
  Cut operator || (const Cut & aptr, const Cut & bptr);
  /// Logical NOT operation on a cut
  Cut operator ! (const Cut & cptr);

  /// Logical AND operation on two cuts
  Cut operator & (const Cut & aptr, const Cut & bptr);
  /// Logical OR operation on two cuts
  Cut operator | (const Cut & aptr, const Cut & bptr);
  /// Logical NOT operation on a cut
  Cut operator ~ (const Cut & cptr);
  /// Logical XOR operation on two cuts
  Cut operator ^ (const Cut & aptr, const Cut & bptr);

  //@}


}

#endif
