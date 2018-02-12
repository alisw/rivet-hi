#include "Rivet/Math/MathUtils.hh"
#include "Rivet/Math/Vectors.hh"
#include "Rivet/Math/Matrices.hh"
// #include "Rivet/Math/MatrixDiag.hh"
using namespace Rivet;

#include <iostream>
#include <limits>
#include <cassert>
using namespace std;

int main() {

  FourVector a(1,0,0,0);
  cout << a << ": interval = " << a.invariant() << endl;
  assert(fuzzyEquals(a.invariant(), 1));
  a.setZ(1);
  assert(isZero(a.invariant()));
  cout << a << ": interval = " << a.invariant() << endl;
  a.setY(2).setZ(3);
  cout << a << ": interval = " << a.invariant() << endl;
  assert(fuzzyEquals(a.invariant(), -12));
  cout << a << ": vector = " << a.vector3() << endl << endl;

  FourMomentum b(1,0,0,0);
  cout << b << ": mass = " << b.mass() << endl;
  assert(fuzzyEquals(b.mass2(), 1));
  b.setPz(1);
  cout << b << ": mass = " << b.mass() << endl;
  assert(isZero(b.mass2()));
  b.setPy(2).setPz(3).setE(6);
  cout << b << ": mass = " << b.mass2() << endl;
  assert(fuzzyEquals(b.mass2(), 23));
  cout << b << ": vector = " << b.vector3() << endl << endl;

  Matrix3 m;
  m.set(0, 0, 7/4.0);
  m.set(0, 1, 3 * sqrt(3)/4.0);
  m.set(1, 0, 3 * sqrt(3)/4.0);
  m.set(1, 1, 13/4.0);
  m.set(2, 2, 9);
  cout << m << endl << endl;
//  EigenSystem<3> es = diagonalize(m);

  cout << "Matrices:" << endl;
  cout << Matrix3() << endl;
  cout << Matrix3::mkIdentity() << endl;
  const Matrix3 I3 = Matrix3::mkIdentity();
  cout << Matrix3::mkIdentity() * m * I3 << endl;
  cout << "tr(0) & det(0): " << Matrix3().trace() << ", " << Matrix3().det() << endl;
  cout << "tr(I3) & det(I3): " << I3.trace() << ", " << I3.det() << endl;
  Matrix3 m1 = Matrix3::mkIdentity();
  Matrix3 m2 = m1;
  m1.setRow(1, Vector3(1,2,3));
  m2.setColumn(1, Vector3(3,2,1));
  Matrix3 m3 = Matrix3::mkZero();
  cout << m1 << " + " << m2 << " = " << m1 + m2 << endl;
  m3.setRow(0, Vector3(2,3,0)).setRow(1, Vector3(1,4,3)).setRow(2, Vector3(0,1,2));
  cout << m1+m2 << " == " << m3 << ": " << (m1+m2 == m3 ? "true" : "false") << endl;
  cout << endl;


  Vector3 v3(1,2,3);
  cout << "Vector: " << v3 << endl;
  cout << "Invert: " << v3 << " --> " << -v3 << endl;
  const Matrix3 rot90(Vector3(0,0,1), PI/2.0);
  const Matrix3 rot90m(Vector3(0,0,1), -PI/2.0);
  const Matrix3 rot180(Vector3(0,0,1), PI);
  const Matrix3 rot180m(Vector3(0,0,1), -PI);
  const Vector3 v3_90 = rot90*v3;
  cout << "Rot 90: " << v3 << " ---> " << v3_90 << endl;
  const Vector3 v3_90m = rot90m*v3;
  cout << "Rot -90: " << v3 << " ---> " << v3_90m << endl;
  const Vector3 v3_180 = rot180*v3;
  cout << "Rot 180: " << v3 << " ---> " << v3_180 << endl;
  const Vector3 v3_180m = rot180m*v3;
  cout << "Rot -180: " << v3 << " ---> " << v3_180m << endl;
  assert(fuzzyEquals(v3_180, v3_180m));

  const Vector3 v3_9090 = rot90*rot90*v3;
  cout << "Rot 2 x 90: " << v3 << " ---> " << v3_9090 << endl;
  assert(fuzzyEquals(v3_180, v3_9090));

  const Vector3 v3_90m90m = rot90m*rot90m*v3;
  cout << "Rot 2 x -90: " << v3 << " ---> " << v3_90m90m << endl;
  assert(fuzzyEquals(v3_180, v3_90m90m));

  const Vector3 v3_9090m = rot90*rot90m*v3;
  const Vector3 v3_90m90 = rot90m*rot90*v3;
  cout << "Rot 90*-90: "<< v3 << " ---> " << v3_9090m << endl;
  cout << "Rot -90*90: "<< v3 << " ---> " << v3_90m90 << endl;
  assert(fuzzyEquals(v3, v3_9090m));
  assert(fuzzyEquals(v3, v3_90m90));

  const Vector3 v3_90i = rot90.inverse()*v3;
  cout << "Rot (90)^-1: "<< v3 << " ---> " << v3_90i << endl;
  assert(fuzzyEquals(v3_90i, v3_90m));

  const Vector3 v3_9090i = rot90*rot90.inverse()*v3;
  const Vector3 v3_90i90 = rot90.inverse()*rot90*v3;
  cout << "Rot 90*(90)^-1: "<< v3 << " ---> " << v3_9090i << endl;
  cout << "Rot (90)^-1*90: "<< v3 << " ---> " << v3_90i90 << endl;
  assert(fuzzyEquals(v3, v3_9090i));
  assert(fuzzyEquals(v3, v3_90i90));

  const Matrix3 rot1(Vector3(0,1,0), PI/180.0);
  cout << "Rot 0 x 45 x 1: " << v3 << endl;
  for (size_t i = 0; i < 8; ++i) {
    for (size_t j = 0; j < 45; ++j) {
      v3 = rot1*v3;
    }
    cout << "Rot " << i+1 << " x 45 x 1: " << v3 << endl;
  }
  assert(fuzzyEquals(v3, Vector3(1,2,3)));
  cout << endl;

  cout << "Boosts:" << endl;
  LorentzTransform ltX = LorentzTransform::mkObjTransformFromBeta(Vector3(0.5, 0, 0));
  cout << "LTx: " << ltX << endl;
  cout << "I on LTx: " << ltX.rotate(Matrix3::mkIdentity()) << endl;
  cout << "Rot90 on LTx: " << ltX.rotate(rot90) << endl;
  cout << endl;

  cout << "X-boosts:" << endl;
  const FourMomentum p1 = FourMomentum(10, 0, 0, 1);
  const FourMomentum p2 = ltX.transform(p1);
  cout << p1 << " -> " << p2 << endl;
  cout << p2 << " -> " << ltX.inverse().transform(p2) << endl;
  //cout << p1.boostVector() << endl;
  const FourMomentum p3 = LorentzTransform::mkFrameTransformFromBeta(p1.boostVector()).transform(p1);
  cout << p1 << " -> " << p3 << endl;
  cout << endl;

  LorentzTransform ltY = LorentzTransform::mkObjTransformFromGamma(Vector3(0, 1.2, 0));
  cout << FourMomentum(1,0,0,1) << " -> " //<< "\n  "
       << (ltX * ltY).transform(FourMomentum(1,0,0,1)) << endl;
  cout << FourMomentum(1,0,0,1) << " -> " //<< "\n  "
       << (ltY * ltX).transform(FourMomentum(1,0,0,1)) << endl;
  cout << (ltX * ltY).betaVec() << endl;
  cout << (ltY * ltX).betaVec() << endl;
  cout << (ltX * ltX.inverse()).betaVec() << endl;

  // If we are already in the rest frame and there is no boost, then LT is trivial/identity
  LorentzTransform noBoost;
  cout << "Element  0,0 should be 1 and is " << noBoost.toMatrix().get(0,0) << endl;
  assert(noBoost.toMatrix().get(0,0)==1);
  cout << "Element  0,1 should be 0 and is " << noBoost.toMatrix().get(0,1) << endl;
  assert(noBoost.toMatrix().get(0,1)==0);
  cout << "Element  1,0 should be 0 and is " << noBoost.toMatrix().get(1,0) << endl;
  assert(noBoost.toMatrix().get(1,0)==0);
  cout << "Element  1,1 should be 1 and is " << noBoost.toMatrix().get(1,1) << endl;
  assert(noBoost.toMatrix().get(1,1)==1);

  return EXIT_SUCCESS;
}
