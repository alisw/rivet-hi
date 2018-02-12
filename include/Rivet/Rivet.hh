#ifndef RIVET_Rivet_HH
#define RIVET_Rivet_HH

#include <string>

namespace Rivet {

  /// A function to get the Rivet version string
  std::string version();

}

// Import main user-interface Rivet classes
#include "Rivet/Particle.hh"
#include "Rivet/Event.hh"
#include "Rivet/Projection.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/AnalysisHandler.hh"

// Rivet path lookup stuff
#include "Rivet/Tools/RivetPaths.hh"

#endif
