#ifndef RIVET_RIVETFASTJET_HH
#define RIVET_RIVETFASTJET_HH

#include "Rivet/Config/RivetCommon.hh"

#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/tools/Filter.hh"

namespace Rivet {


  /// Unscoped awareness of FastJet's PseudoJet
  using fastjet::PseudoJet;

  /// Typedef for a collection of PseudoJet objects.
  typedef std::vector<PseudoJet> PseudoJets;


  /// Make a 3-momentum vector from a FastJet pseudojet
  inline Vector3 momentum3(const fastjet::PseudoJet& pj) {
    return Vector3(pj.px(), pj.py(), pj.pz());
  }

  /// Make a 4-momentum vector from a FastJet pseudojet
  inline FourMomentum momentum(const fastjet::PseudoJet& pj) {
    return FourMomentum(pj.E(), pj.px(), pj.py(), pj.pz());
  }


}

#endif
