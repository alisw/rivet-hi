// -*- C++ -*-
#include "Rivet/Projections/DISFinalState.hh"

namespace Rivet {


  void DISFinalState::project(const Event& e) {
    const DISKinematics& diskin = apply<DISKinematics>(e, "Kinematics");
    LorentzTransform hcmboost; //< Null boost = LAB frame by default
    if (_boosttype == HCM) hcmboost = diskin.boostHCM();
    else if (_boosttype == BREIT) hcmboost = diskin.boostBreit();

    const DISLepton& dislep = diskin.apply<DISLepton>(e, "Lepton");

    const FinalState& fs = apply<FinalState>(e, "FS");

    // Fill the particle list with all particles _other_ than the DIS scattered
    // lepton, with momenta boosted into the appropriate frame.
    _theParticles.clear();
    _theParticles.reserve(fs.particles().size()-1);
    const GenParticle* dislepGP = dislep.out().genParticle();
    // const GenParticle* dislepIN = dislep.in().genParticle();

    for (const Particle& p : fs.particles()) { ///< Ensure that we skip the DIS lepton
      Particle temp = p;
      if (_boosttype != LAB) temp.setMomentum(hcmboost.transform(temp.momentum()));
      if (p.genParticle() != dislepGP)  _theParticles.push_back(temp);
    }

  }


}
