#include "Rivet/AnalysisHandler.hh"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"

using namespace std;

int main() {
  Rivet::AnalysisHandler ah;
  Rivet::Log::setLevel("Rivet", Rivet::Log::TRACE);

  // Specify the analyses to be used
  ah.addAnalysis("EXAMPLE");
  ah.addAnalyses({{ "MC_JETS", "EXAMPLE_CUTS", "EXAMPLE_SMEAR" }});

  std::ifstream file("testApi.hepmc");
  HepMC::IO_GenEvent hepmcio(file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  double sum_of_weights = 0.0;
  while (evt) {
    // Analyse current event
    ah.analyze(*evt);
    sum_of_weights += evt->weights()[0];

    // Clean up and get next event
    delete evt; evt = nullptr;
    hepmcio >> evt;
  }
  file.close();

  ah.setCrossSection(1.0);
  ah.finalize();
  ah.writeData("out.yoda");

  return 0;
}
