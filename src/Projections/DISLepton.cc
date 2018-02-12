// -*- C++ -*-
#include "Rivet/Projections/DISLepton.hh"

namespace Rivet {


  int DISLepton::compare(const Projection& p) const {
    const DISLepton& other = pcast<DISLepton>(p);
    return mkNamedPCmp(other, "Beam") || mkNamedPCmp(other, "FS");
  }


  void DISLepton::project(const Event& e) {

    // Find incoming lepton beam
    const ParticlePair& inc = applyProjection<Beam>(e, "Beam").beams();
    bool firstIsLepton = PID::isLepton(inc.first.pid());
    bool secondIsLepton = PID::isLepton(inc.second.pid());
    if (firstIsLepton && !secondIsLepton) {
      _incoming = inc.first;
    } else if (!firstIsLepton && secondIsLepton) {
      _incoming = inc.second;
    } else {
      throw Error("DISLepton could not find the correct beam");
    }

    // // Find outgoing scattered lepton via HepMC graph
    // /// @todo Evidence that this doesn't work with Sherpa... FIX
    // const GenParticle* current_l = _incoming.genParticle();
    // bool found_next_vertex = true;
    // while (found_next_vertex) {
    //   found_next_vertex = false;
    //   if (!current_l->end_vertex()) break;
    //   // Get lists of outgoing particles consistent with a neutral (gamma/Z) or charged (W) DIS current
    //   /// @todo Avoid loops
    //   vector<const GenParticle*> out_n, out_c;
    //   for (const GenParticle* pp : particles_out(current_l, HepMC::children)) {
    //     if (current_l->pdg_id() == pp->pdg_id()) out_n.push_back(pp);
    //     if (std::abs(std::abs(current_l->pdg_id()) - std::abs(pp->pdg_id())) == 1) out_c.push_back(pp);
    //   }
    //   if (out_n.empty() && out_c.empty()) {
    //     MSG_WARNING("No lepton in the new vertex");
    //     break;
    //   }
    //   if (out_c.size() + out_n.size() > 1) {
    //     MSG_WARNING("More than one lepton in the new vertex");
    //     break;
    //   }
    //   current_l = out_c.empty() ? out_n.front() : out_c.front();
    //   found_next_vertex = true;
    // }
    // if (current_l != nullptr) {
    //   _outgoing = Particle(current_l);
    //   MSG_DEBUG("Found DIS lepton from event-record structure");
    //   return;
    // }

    // If no graph-connected scattered lepton, use the hardest (preferably same-flavour) prompt FS lepton in the event
    /// @todo Specify the charged or neutral current being searched for in the DISLepton constructor/API, and remove the guesswork
    const Particles fsleptons = applyProjection<FinalState>(e, "PromptFS").particles(isLepton, cmpMomByE);
    const Particles sfleptons = filter_select(fsleptons, Cuts::pid == _incoming.pid());
    MSG_DEBUG("SF leptons = " << sfleptons.size() << ", all leptons = " << fsleptons.size());
    if (!sfleptons.empty()) {
      _outgoing = sfleptons.front();
    } else if (!fsleptons.empty()) {
      _outgoing = fsleptons.front();
    } else {
      throw Error("Could not find the scattered lepton");
    }

  }


}
