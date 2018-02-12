#include "Rivet/AnalysisHandler.hh"
#include "Rivet/AnalysisLoader.hh"
#include "HepMC/IO_GenEvent.h"

using namespace std;

int main(int argc, char** argv) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <hepmcfile> <ana1> [<ana2> ...]" << endl;
    return 1;
  }

  foreach (const string& a, Rivet::AnalysisLoader::analysisNames())
    cout << a << endl;

  Rivet::AnalysisHandler ah;
  for (int i = 2; i < argc; ++i) {
    ah.addAnalysis(argv[i]);
  }

  std::ifstream file(argv[1]);
  HepMC::IO_GenEvent hepmcio(file);
  HepMC::GenEvent* evt = hepmcio.read_next_event();
  while (evt) {
    ah.analyze(*evt);
    delete evt; evt = 0;
    hepmcio >> evt;
  }
  file.close();

  ah.setCrossSection(1.0);
  ah.finalize();
  ah.writeData("Rivet.yoda");

  return 0;
}
