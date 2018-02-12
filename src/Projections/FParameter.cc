// -*- C++ -*-
#include "Rivet/Projections/FParameter.hh"

namespace Rivet {


  FParameter::FParameter(const FinalState& fsp) {
    setName("FParameter");
    addProjection(fsp, "FS");
    clear();
  }


  void FParameter::clear() {
    _lambdas = vector<double>(2, 0);
  }


  void FParameter::project(const Event& e) {
    const Particles prts = applyProjection<FinalState>(e, "FS").particles();
    calc(prts);
  }


  void FParameter::calc(const FinalState& fs) {
    calc(fs.particles());
  }

  void FParameter::calc(const vector<Particle>& fsparticles) {
    vector<Vector3> threeMomenta;
    threeMomenta.reserve(fsparticles.size());
    for (const Particle& p : fsparticles) {
      const Vector3 p3 = p.momentum().vector3();
      threeMomenta.push_back(p3);
    }
    _calcFParameter(threeMomenta);
  }

  void FParameter::calc(const vector<FourMomentum>& fsmomenta) {
    vector<Vector3> threeMomenta;
    threeMomenta.reserve(fsmomenta.size());
    for (const FourMomentum& v : fsmomenta) {
      threeMomenta.push_back(v.vector3());
    }
    _calcFParameter(threeMomenta);
  }

  void FParameter::calc(const vector<Vector3>& fsmomenta) {
    _calcFParameter(fsmomenta);
  }

  // Actually do the calculation
  void FParameter::_calcFParameter(const vector<Vector3>& fsmomenta) {

    // Return (with "safe nonsense" sphericity params) if there are no final state particles.
    if (fsmomenta.empty()) {
      MSG_DEBUG("No particles in final state...");
      clear();
      return;
    }

    // A small iteration over full momenta but set z-coord. to 0.0 to get transverse momenta
    vector <Vector3> fsperpmomenta;
    for (const Vector3& p : fsmomenta) {
      fsperpmomenta.push_back(Vector3(p.x(), p.y(), 0.0));
    }

    // Iterate over all the final state particles.
    Matrix<2> mMom;
    MSG_DEBUG("Number of particles = " << fsperpmomenta.size());
    for (const Vector3& p3 : fsperpmomenta) {

      double prefactor = 1.0/p3.mod();

      Matrix<2> mMomPart;
      for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 2; ++j) {
          mMomPart.set(i,j, p3[i]*p3[j]);
        }
      }
      mMom += prefactor * mMomPart;
    }

    MSG_DEBUG("Linearised transverse momentum tensor = " << mMom);

    // Check that the matrix is symmetric.
    const bool isSymm = mMom.isSymm();
    if (!isSymm) {
      MSG_ERROR("Error: momentum tensor not symmetric:");
      MSG_ERROR("[0,1] vs. [1,0]: " << mMom.get(0,1) << ", " << mMom.get(1,0));
    }
    // If not symmetric, something's wrong (we made sure the error msg appeared first).

    assert(isSymm);
    const double a = mMom.get(0,0);
    const double b = mMom.get(1,1);
    const double c = mMom.get(1,0);

    const double l1 = 0.5*(a+b+sqrt( (a-b)*(a-b) + 4 *c*c));
    const double l2 = 0.5*(a+b-sqrt( (a-b)*(a-b) + 4 *c*c));

    _lambdas = {l1, l2};

    // Debug output.
    MSG_DEBUG("Lambdas = ("
             << lambda1() << ", " << lambda2() << ")");
    MSG_DEBUG("Sum of lambdas = " << lambda1() + lambda2());
    MSG_DEBUG("F-Parameter = " << F());
  }
}
