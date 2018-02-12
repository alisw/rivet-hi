#ifndef RIVET_RIVETYODA_HH
#define RIVET_RIVETYODA_HH

#include "Rivet/Config/RivetCommon.hh"
#include "YODA/AnalysisObject.h"
#include "YODA/Counter.h"
#include "YODA/Histo1D.h"
#include "YODA/Histo2D.h"
#include "YODA/Profile1D.h"
#include "YODA/Profile2D.h"
#include "YODA/Scatter1D.h"
#include "YODA/Scatter2D.h"
#include "YODA/Scatter3D.h"

namespace Rivet {

  typedef std::shared_ptr<YODA::AnalysisObject> AnalysisObjectPtr;
  typedef std::shared_ptr<YODA::Counter> CounterPtr;
  typedef std::shared_ptr<YODA::Histo1D> Histo1DPtr;
  typedef std::shared_ptr<YODA::Histo2D> Histo2DPtr;
  typedef std::shared_ptr<YODA::Profile1D> Profile1DPtr;
  typedef std::shared_ptr<YODA::Profile2D> Profile2DPtr;
  typedef std::shared_ptr<YODA::Scatter1D> Scatter1DPtr;
  typedef std::shared_ptr<YODA::Scatter2D> Scatter2DPtr;
  typedef std::shared_ptr<YODA::Scatter3D> Scatter3DPtr;

  using YODA::AnalysisObject;
  using YODA::Counter;
  using YODA::Histo1D;
  using YODA::HistoBin1D;
  using YODA::Histo2D;
  using YODA::HistoBin2D;
  using YODA::Profile1D;
  using YODA::ProfileBin1D;
  using YODA::Profile2D;
  using YODA::ProfileBin2D;
  using YODA::Scatter1D;
  using YODA::Point1D;
  using YODA::Scatter2D;
  using YODA::Point2D;
  using YODA::Scatter3D;
  using YODA::Point3D;


  /// Function to get a map of all the refdata in a paper with the given @a papername.
  map<string, AnalysisObjectPtr> getRefData(const string& papername);

  /// Get the file system path to the reference file for this paper.
  string getDatafilePath(const string& papername);


}

#endif
