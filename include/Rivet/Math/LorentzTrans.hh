#ifndef RIVET_MATH_LORENTZTRANS
#define RIVET_MATH_LORENTZTRANS

#include "Rivet/Math/MathHeader.hh"
#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/MatrixN.hh"
#include "Rivet/Math/Matrix3.hh"
#include "Rivet/Math/Vector4.hh"
#include <iostream>

namespace Rivet {


  /// @brief Object implementing Lorentz transform calculations and boosts.
  ///
  /// @note These boosts are defined actively, i.e. as modifications of vectors
  /// rather than frame transformations. So the boost vector is the opposite of
  /// what you might expect.
  ///
  /// @todo Review the active/passive convention choice. Seems counterintuitive now...
  class LorentzTransform {
  public:

    /// @name Simple Lorentz factor conversions
    //@{

    /// Calculate the \f$ \gamma \f$ factor from \f$ \beta \f$
    static double beta2gamma(double beta) {
      return 1.0 / sqrt(1 - sqr(beta));
    }

    /// Calculate \f$ \beta \f$ from the \f$ \gamma \f$ factor
    static double gamma2beta(double gamma) {
      return sqrt(1 - sqr(1/gamma));
    }

    // /// Calculate the \f$ \gamma \f$ factor from \f$ \bar{\beta}^2 = 1-\beta^2 \f$
    // static double betabarsq2gamma(double betabarsq) {
    //   return 1.0 / sqrt(betabarsq);
    // }

    // /// Calculate \f$ \bar{\beta}^2 = 1-\beta^2 \f$ from the \f$ \gamma \f$ factor
    // static double gamma2betabarsq(double gamma) {
    //   return 1.0 / sqr(gamma);
    // }

    //@}


    /// @name Construction
    //@{

    /// Default (identity) constructor
    LorentzTransform() {
      _boostMatrix = Matrix<4>::mkIdentity();
    }


    /// Make an LT for an active boost (i.e. object velocity += in boost direction)
    static LorentzTransform mkObjTransformFromBeta(const Vector3& vbeta) {
      LorentzTransform rtn;
      return rtn.setBetaVec(vbeta);
    }

    /// Make an LT for a passive boost (i.e. object velocity -= in boost direction)
    static LorentzTransform mkFrameTransformFromBeta(const Vector3& vbeta) {
      LorentzTransform rtn;
      return rtn.setBetaVec(-vbeta);
    }

    /// Make an LT for an active boost (i.e. object velocity += in boost direction)
    static LorentzTransform mkObjTransformFromGamma(const Vector3& vgamma) {
      LorentzTransform rtn;
      if (!vgamma.isZero()) rtn.setGammaVec(vgamma);
      return rtn;
    }

    /// Make an LT for a passive boost (i.e. object velocity -= in boost direction)
    static LorentzTransform mkFrameTransformFromGamma(const Vector3& vgamma) {
      LorentzTransform rtn;
      if (!vgamma.isZero()) rtn.setGammaVec(-vgamma);
      return rtn;
    }

    //@}


    /// @name Boost vector and beta/gamma factors
    //@{

    /// Set up an active Lorentz boost from the \f$ \vec\beta \f$ vector
    LorentzTransform& setBetaVec(const Vector3& vbeta) {
      assert(vbeta.mod2() < 1);
      const double beta = vbeta.mod();
      const double gamma = beta2gamma(beta);
      _boostMatrix = Matrix<4>::mkIdentity();
      _boostMatrix.set(0, 0, gamma);
      _boostMatrix.set(1, 1, gamma);
      _boostMatrix.set(0, 1, +beta*gamma); //< +ve coeff since active boost
      _boostMatrix.set(1, 0, +beta*gamma); //< +ve coeff since active boost
      if (beta > 0) _boostMatrix = rotate(Vector3::mkX(), vbeta)._boostMatrix;
      return *this;
    }

    /// Get the \f$ \vec\beta \f$ vector for an active Lorentz boost
    Vector3 betaVec() const {
      FourMomentum boost(_boostMatrix.getColumn(0)); //< @todo WRONG?!
      //cout << "!!!" << boost << endl;
      if (boost.isZero()) return Vector3();
      assert(boost.E() > 0);
      const double beta = boost.p3().mod() / boost.E();
      return boost.p3().unit() * beta;
    }

    /// Get the \f$ \beta \f$ factor
    double beta() const {
      return betaVec().mod();
    }


    /// Set up an active Lorentz boost from the \f$ \vec\gamma \f$ vector
    LorentzTransform& setGammaVec(const Vector3& vgamma) {
      const double gamma = vgamma.mod();
      const double beta = gamma2beta(gamma);
      _boostMatrix = Matrix<4>::mkIdentity();
      _boostMatrix.set(0, 0, gamma);
      _boostMatrix.set(1, 1, gamma);
      _boostMatrix.set(0, 1, +beta*gamma); //< +ve coeff since active boost
      _boostMatrix.set(1, 0, +beta*gamma); //< +ve coeff since active boost
      if (beta > 0) _boostMatrix = rotate(Vector3::mkX(), vgamma)._boostMatrix;
      return *this;
    }

    /// Get the \f$ \vec\gamma \f$ vector for an active Lorentz boost
    Vector3 gammaVec() const {
      FourMomentum boost(_boostMatrix.getColumn(0)); //< @todo WRONG?!
      if (boost.isZero()) return Vector3();
      assert(boost.E() > 0);
      const double beta = boost.p3().mod() / boost.E();
      return boost.p3().unit() * beta;
    }

    /// Get the \f$ \gamma \f$ factor
    double gamma() const {
      return beta2gamma(beta());
    }

    //@}


    /// Apply this transformation to the given 4-vector
    FourVector transform(const FourVector& v4) const {
      return multiply(_boostMatrix, v4);
    }

    /// Apply this transformation to the given 4-mometum
    FourMomentum transform(const FourMomentum& v4) const {
      return multiply(_boostMatrix, v4);
    }


    /// @name Transform modifications
    //@{

    /// Rotate the transformation cf. the difference between vectors @a from and @a to
    LorentzTransform rotate(const Vector3& from, const Vector3& to) const {
      return rotate(Matrix3(from, to));
    }

    /// Rotate the transformation by @a angle radians about @a axis
    LorentzTransform rotate(const Vector3& axis, double angle) const {
      return rotate(Matrix3(axis, angle));
    }

    /// Rotate the transformation by the 3D rotation matrix @a rot
    LorentzTransform rotate(const Matrix3& rot) const {
      LorentzTransform lt = *this;
      const Matrix4 rot4 = _mkMatrix4(rot);
      const Matrix4 newlt = rot4 * _boostMatrix * rot4.inverse();
      lt._boostMatrix = newlt;
      return lt;
    }

    /// Calculate the inverse transform
    LorentzTransform inverse() const {
      LorentzTransform rtn;
      rtn._boostMatrix = _boostMatrix.inverse();
      return rtn;
    }

    /// Combine LTs, treating @a this as the LH matrix.
    LorentzTransform combine(const LorentzTransform& lt) const {
      LorentzTransform rtn;
      rtn._boostMatrix = _boostMatrix * lt._boostMatrix;
      return rtn;
    }

    /// Operator combination of two LTs
    LorentzTransform operator*(const LorentzTransform& lt) const {
      return combine(lt);
    }

    /// Pre-multiply m3 by this LT
    LorentzTransform preMult(const Matrix3& m3) {
      _boostMatrix = multiply(_mkMatrix4(m3),_boostMatrix);
      return *this;
    }

    /// Post-multiply m3 by this LT
    LorentzTransform postMult(const Matrix3& m3) {
      _boostMatrix *= _mkMatrix4(m3);
      return *this;
    }

    //@}


    /// Return the matrix form
    Matrix4 toMatrix() const {
      return _boostMatrix;
    }


  private:

    Matrix4 _mkMatrix4(const Matrix3& m3) const {
      Matrix4 m4 = Matrix4::mkIdentity();
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          m4.set(i+1, j+1, m3.get(i, j));
        }
      }
      return m4;
    }

    Matrix4 _boostMatrix;

  };



  inline LorentzTransform inverse(const LorentzTransform& lt) {
    return lt.inverse();
  }

  inline LorentzTransform combine(const LorentzTransform& a, const LorentzTransform& b) {
    return a.combine(b);
  }

  inline FourVector transform(const LorentzTransform& lt, const FourVector& v4) {
      return lt.transform(v4);
  }


  //////////////////////////


  inline string toString(const LorentzTransform& lt) {
    return toString(lt.toMatrix());
  }

  inline ostream& operator<<(std::ostream& out, const LorentzTransform& lt) {
    out << toString(lt);
    return out;
  }


}

#endif
