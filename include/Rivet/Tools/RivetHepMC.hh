// -*- C++ -*-
#ifndef RIVET_RivetHepMC_HH
#define RIVET_RivetHepMC_HH

#include "HepMC/GenEvent.h"
#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenRanges.h"
#include "HepMC/IO_GenEvent.h"
#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/Exceptions.hh"

namespace Rivet {


  using HepMC::GenEvent;
  using HepMC::GenParticle;
  using HepMC::GenVertex;

  #if HEPMC_VERSION_CODE >= 3000000
  using HepMC::GenEventPtr;
  using HepMC::GenParticlePtr;
  using HepMC::GenVertexPtr;
  #elif HEPMC_VERSION_CODE >= 2007000
  // HepMC 2.07 provides its own #defines
  #else
  #define GenEventPtr GenEvent*
  #define GenParticlePtr GenParticle*
  #define GenVertexPtr GenVertex*
  #endif


  /// @todo Use mcutils?


  inline std::vector<GenParticle const *> particles(const GenEvent* ge) {
    assert(ge);
    std::vector<const GenParticle*> rtn;
    for (GenEvent::particle_const_iterator pi = ge->particles_begin(); pi != ge->particles_end(); ++pi)
      rtn.push_back(*pi);
    return rtn;
  }

  inline std::vector<GenParticlePtr> particles(GenEvent* ge) {
    assert(ge);
    std::vector<GenParticle*> rtn;
    for (GenEvent::particle_iterator pi = ge->particles_begin(); pi != ge->particles_end(); ++pi)
      rtn.push_back(*pi);
    return rtn;
  }


  inline std::vector<const GenVertex*> vertices(const GenEvent* ge) {
    std::vector<GenVertex const *> rtn;
    for (GenEvent::vertex_const_iterator vi = ge->vertices_begin(); vi != ge->vertices_end(); ++vi)
      rtn.push_back(*vi);
    return rtn;
  }

  inline std::vector<GenVertex*> vertices(GenEvent* ge) {
    std::vector<GenVertex*> rtn;
    for (GenEvent::vertex_iterator vi = ge->vertices_begin(); vi != ge->vertices_end(); ++vi)
      rtn.push_back(*vi);
    return rtn;
  }


  //////////////////////////


  inline std::vector<const GenParticle*> particles(const GenVertex* gv, HepMC::IteratorRange range=HepMC::relatives) {
    std::vector<GenParticle const *> rtn;
    /// @todo A particle_const_iterator on GenVertex would be nice...
    // Before HepMC 2.7.0 there were no GV::particles_const_iterators and constness consistency was all screwed up :-/
    #if HEPMC_VERSION_CODE >= 2007000
    // for (GenVertex::particle_iterator pi = gv->particles_const_begin(range); pi != gv->particles_const_end(range); ++pi)
    for (GenVertex::particle_iterator pi = gv->particles_begin(range); pi != gv->particles_end(range); ++pi)
      rtn.push_back(*pi);
    #else
    GenVertex* gv2 = const_cast<GenVertex*>(gv);
    for (GenVertex::particle_iterator pi = gv2->particles_begin(range); pi != gv2->particles_end(range); ++pi)
      rtn.push_back(const_cast<const GenParticle*>(*pi));
    #endif
    return rtn;
  }

  inline std::vector<GenParticle*> particles(GenVertex* gv, HepMC::IteratorRange range=HepMC::relatives) {
    std::vector<GenParticle*> rtn;
    for (GenVertex::particle_iterator pi = gv->particles_begin(range); pi != gv->particles_end(range); ++pi)
      rtn.push_back(*pi);
    return rtn;
  }



  // Get iterator ranges as wrapped begin/end pairs
  /// @note GenVertex _in and _out iterators are actually, secretly the same types *sigh*
  struct GenVertexIterRangeC {
    typedef vector<GenParticle*>::const_iterator genvertex_particles_const_iterator;
    GenVertexIterRangeC(const genvertex_particles_const_iterator& begin, const genvertex_particles_const_iterator& end)
      : _begin(begin), _end(end) {  }
    const genvertex_particles_const_iterator& begin() { return _begin; }
    const genvertex_particles_const_iterator& end() { return _end; }
  private:
    const genvertex_particles_const_iterator _begin, _end;
  };

  inline GenVertexIterRangeC particles_in(const GenVertex* gv) {
    return GenVertexIterRangeC(gv->particles_in_const_begin(), gv->particles_in_const_end());
  }

  inline GenVertexIterRangeC particles_out(const GenVertex* gv) {
    return GenVertexIterRangeC(gv->particles_out_const_begin(), gv->particles_out_const_end());
  }



  #if HEPMC_VERSION_CODE >= 2007000

  // Get iterator ranges as wrapped begin/end pairs
  /// @note GenVertex _in and _out iterators are actually, secretly the same types *sigh*
  struct GenVertexIterRange {
    typedef vector<GenParticle*>::iterator genvertex_particles_iterator;
    GenVertexIterRange(const genvertex_particles_iterator& begin, const genvertex_particles_iterator& end)
      : _begin(begin), _end(end) {  }
    const genvertex_particles_iterator& begin() { return _begin; }
    const genvertex_particles_iterator& end() { return _end; }
  private:
    const genvertex_particles_iterator _begin, _end;
  };

  inline GenVertexIterRange particles_in(GenVertex* gv) {
    return GenVertexIterRange(gv->particles_in_begin(), gv->particles_in_end());
  }

  inline GenVertexIterRange particles_out(GenVertex* gv) {
    return GenVertexIterRange(gv->particles_out_begin(), gv->particles_out_end());
  }

  #endif


  //////////////////////////


  /// Get the direct parents or all-ancestors of GenParticle @a gp
  inline std::vector<const GenParticle*> particles_in(const GenParticle* gp, HepMC::IteratorRange range=HepMC::ancestors) {
    if (range != HepMC::parents && range != HepMC::ancestors)
      throw UserError("Requested particles_in(GenParticle*) with a non-'in' iterator range");
    if (!gp->production_vertex()) return std::vector<const GenParticle*>();
    #if HEPMC_VERSION_CODE >= 2007000
    return particles(gp->production_vertex(), range);
    #else
    // Before HepMC 2.7.0 the constness consistency of methods and their return types was all screwed up :-/
    std::vector<const GenParticle*> rtn;
    foreach (GenParticle* gp2, particles(gp->production_vertex(), range))
      rtn.push_back( const_cast<const GenParticle*>(gp2) );
    return rtn;
    #endif
  }

  /// Get the direct parents or all-ancestors of GenParticle @a gp
  inline std::vector<GenParticle*> particles_in(GenParticle* gp, HepMC::IteratorRange range=HepMC::ancestors) {
    if (range != HepMC::parents && range != HepMC::ancestors)
      throw UserError("Requested particles_in(GenParticle*) with a non-'in' iterator range");
    return (gp->production_vertex()) ? particles(gp->production_vertex(), range) : std::vector<GenParticle*>();
  }


  /// Get the direct children or all-descendents of GenParticle @a gp
  inline std::vector<const GenParticle*> particles_out(const GenParticle* gp, HepMC::IteratorRange range=HepMC::descendants) {
    if (range != HepMC::children && range != HepMC::descendants)
      throw UserError("Requested particles_out(GenParticle*) with a non-'out' iterator range");
    if (!gp->end_vertex()) return std::vector<const GenParticle*>();
    #if HEPMC_VERSION_CODE >= 2007000
    return particles(gp->end_vertex(), range);
    #else
    // Before HepMC 2.7.0 the constness consistency of methods and their return types was all screwed up :-/
    std::vector<const GenParticle*> rtn;
    foreach (GenParticle* gp2, particles(gp->end_vertex(), range))
      rtn.push_back( const_cast<const GenParticle*>(gp2) );
    return rtn;
    #endif
  }

  /// Get the direct children or all-descendents of GenParticle @a gp
  inline std::vector<GenParticle*> particles_out(GenParticle* gp, HepMC::IteratorRange range=HepMC::descendants) {
    if (range != HepMC::children && range != HepMC::descendants)
      throw UserError("Requested particles_out(GenParticle*) with a non-'out' iterator range");
    return (gp->end_vertex()) ? particles(gp->end_vertex(), range) : std::vector<GenParticle*>();
  }


  /// Get any relatives of GenParticle @a gp
  inline std::vector<const GenParticle*> particles(const GenParticle* gp, HepMC::IteratorRange range=HepMC::ancestors) {
    if (range == HepMC::parents || range == HepMC::ancestors)
      return particles_in(gp, range);
    if (range == HepMC::children || range == HepMC::descendants)
      return particles_in(gp, range);
    throw UserError("Requested particles(GenParticle*) with an unsupported iterator range");
  }


}

#endif
