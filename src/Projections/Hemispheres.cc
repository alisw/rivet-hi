// -*- C++ -*-
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Jet.hh"

namespace Rivet {


  void Hemispheres::project(const Event& e) {
    clear();

    // Get thrust axes.
    const AxesDefinition& ax = applyProjection<AxesDefinition>(e, "Axes");
    const Vector3 n = ax.axis1();
    const FinalState& fs = applyProjection<FinalState>(e, ax.getProjection("FS"));
    const Particles& particles = fs.particles();
    calc(n, particles);
  }


  void Hemispheres::calc(const Vector3& n, const Particles& particles) {
    vector<FourMomentum> p4s; p4s.reserve(particles.size());
    for (const Particle& p : particles) p4s.push_back(p.mom());
    calc(n, p4s);
  }


  void Hemispheres::calc(const Vector3& n, const Jets& jets) {
    vector<FourMomentum> p4s; p4s.reserve(jets.size());
    for (const Jet& j : jets) p4s.push_back(j.mom());
    calc(n, p4s);
  }


  void Hemispheres::calc(const Vector3& n, const vector<FourMomentum>& p4s) {
    MSG_DEBUG("Hemisphere axis = " << n);
    MSG_DEBUG("Number of constituents = " << p4s.size());

    FourMomentum p4With, p4Against;
    double Evis(0), broadWith(0), broadAgainst(0), broadDenom(0);
    for (const FourMomentum& p4 : p4s) {
      const Vector3 p3 = p4.vector3();
      const double p3Para = dot(p3, n);
      const double p3Trans = (p3 - p3Para * n).mod();

      // Update normalisations
      Evis += p4.E();
      broadDenom += 2.0 * p3.mod();

      // Update the mass and broadening variables
      if (p3Para > 0) {
        p4With += p4;
        broadWith += p3Trans;
      } else if (p3Para < 0) {
        p4Against += p4;
        broadAgainst += p3Trans;
      } else {
        // In the incredibly unlikely event that a particle goes exactly along the
        // thrust plane, add half to each hemisphere.
        MSG_WARNING("Particle split between hemispheres");
        p4With += 0.5 * p4;
        p4Against += 0.5 * p4;
        broadWith += 0.5 * p3Trans;
        broadAgainst += 0.5 * p3Trans;
      }
    }

    // Visible energy squared.
    _E2vis = sqr(Evis);

    // Calculate masses.
    const double mass2With = p4With.mass2();
    const double mass2Against = p4Against.mass2();
    _M2high = max(mass2With, mass2Against);
    _M2low = min(mass2With, mass2Against);

    // Calculate broadenings.
    broadWith /= broadDenom;
    broadAgainst /= broadDenom;
    _Bmax = max(broadWith, broadAgainst);
    _Bmin = min(broadWith, broadAgainst);

    // Calculate high-max correlation flag.
    const int maxMassID = (mass2With >= mass2Against);
    const int maxBroadID = (broadWith >= broadAgainst);
    _highMassEqMaxBroad = (maxMassID == maxBroadID);
  }


}
