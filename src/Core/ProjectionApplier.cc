// -*- C++ -*-
#include "Rivet/ProjectionApplier.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  // NB. Allow proj registration in constructor by default -- explicitly disable for Analysis
  ProjectionApplier::ProjectionApplier()
    : _allowProjReg(true), _owned(false),
      _projhandler(ProjectionHandler::getInstance())
  {  }


  ProjectionApplier::~ProjectionApplier() {
    if ( ! _owned )
      getProjHandler().removeProjectionApplier(*this);
  }


  const Projection& ProjectionApplier::_applyProjection(const Event& evt,
                                                        const string& name) const {
    return evt.applyProjection(getProjection(name));
  }


  const Projection& ProjectionApplier::_applyProjection(const Event& evt,
                                                        const Projection& proj) const {
    return evt.applyProjection(proj);
  }


  const Projection& ProjectionApplier::_declareProjection(const Projection& proj,
                                                          const string& name) {
    if (!_allowProjReg) {
      cerr << "Trying to register projection '"
           << proj.name() << "' outside init phase in '" << this->name() << "'." << endl;
      exit(2);
    }
    const Projection& reg = getProjHandler().registerProjection(*this, proj, name);
    return reg;
  }


}
