// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/CentralityBinner.hh"

namespace Rivet {

/// Example of CentralityEstimator projection that uses summed Et in
/// the forward region.
class SumETFwdCentrality: public CentralityEstimator {

public:

  /// Constructor.
  SumETFwdCentrality() {
    declare(FinalState(Cuts::eta < -3.2 && Cuts::eta > -4.9 && Cuts::pT > 0.1*GeV),
	    "FSSumETFwdCentrality");
  }

    /// Clone on the heap.
    DEFAULT_RIVET_PROJ_CLONE(SumETFwdCentrality);

protected:

  /// Perform the projection on the Event
  void project(const Event& e) {
    const FinalState & fsfwd =
      apply<FinalState>(e, "FSSumETFwdCentrality");
    _estimate = 0.0;
    for ( const Particle & p : fsfwd.particles() ) {
      _estimate += p.Et();
    }
  }
  
  /// Compare projections
  int compare(const Projection& p) const {
    return mkNamedPCmp(p, "FSSumETFwdCentrality");
  }

};
    

  /// Generic analysis looking at various distributions of final state particles
  class MC_Centrality : public Analysis {
  public:

    /// Constructor
    MC_Centrality() : Analysis("MC_Centrality"), _cent(400, 0.02), _sumw(0.0) {}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      declare(ChargedFinalState(Cuts::eta > -2.7 && Cuts::eta < 2.7 &&
      				Cuts::pT > 0.1*GeV), "CFS");
      declare(FinalState(Cuts::eta < -3.2 && Cuts::eta > -4.9 &&
      			 Cuts::pT > 0.1*GeV), "CFSF");
      declare(FinalState(Cuts::eta > 2.09 && Cuts::eta < 3.84 &&
      			 Cuts::pT > 0.1*GeV), "MBB");
      declare(FinalState(Cuts::eta < -2.09 && Cuts::eta > -3.84 &&
      			 Cuts::pT > 0.1*GeV), "MBF");
      declare(GeneratedCentrality(), "GeneratedCentrality");
      _gencent.setProjection(GeneratedCentrality(), "GeneratedCentrality");

      
      // Histograms
      // The sum Et in the forward region.
      _hETfwd  = bookHisto1D("ETfwd", 400, 0.0, 200.0);

      vector<double> pclim =
	{100.0, 90.0, 60.0, 40.0, 30.0, 20.0, 10.0, 5.0, 1.0, 0.0 };
      vector<double> etlim =
	{ 0.0*GeV, 3.87*GeV, 12.71*GeV, 22.65*GeV, 31.49*GeV,
	  53.04*GeV, 64.64*GeV, 98.5*GeV, -1.0*GeV };

      for ( int i = 0; i < 8; ++i ) {
	_cent.add(bookHisto1D(2, 1, i+1),
		  pclim[i+1], pclim[i]);
	_fixcent.add(bookHisto1D(12, 1, i+1),
		     pclim[i+1], pclim[i], etlim[i], etlim[i+1]);
	_gencent.add(bookHisto1D(22, 1, i+1),
		     pclim[i+1], pclim[i], pclim[i], pclim[i+1]);
      }
            
      _centrue.clear();

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      if ( applyProjection<FinalState>(event,"MBF").particles().empty() ||
	   applyProjection<FinalState>(event,"MBB").particles().empty() )
	vetoEvent;

      _sumw += weight;

      const ChargedFinalState& parts =
	applyProjection<ChargedFinalState>(event,"CFS");
      const FinalState& forw =
	applyProjection<FinalState>(event,"CFSF");

      double sumEt = 0.0;
      for( const Particle & p : forw.particles() ) sumEt += p.Et();

      _hETfwd->fill(sumEt, weight);
      _centrue.insert(make_pair(sumEt, weight));

      Histo1DPtr ch = _cent.select(sumEt, weight);
      Histo1DPtr fch = _fixcent.select(sumEt, weight);
      Histo1DPtr gch = _gencent.select(event, weight);

      for ( const Particle & p : parts.particles() ) {
	if ( ch ) ch->fill(p.eta(), weight);
	if ( fch ) fch->fill(p.eta(), weight);
	if ( gch ) gch->fill(p.eta(), weight);
      }

    }


    /// Finalize
    void finalize() {

      if ( _sumw == 0.0 ) return;
        
      normalize(_hETfwd, _hETfwd->sumW()/_sumw);

      _cent.finalize();
      _cent.normalizePerEvent();
      _fixcent.finalize();
      _fixcent.normalizePerEvent();

      map<double,double> edges = _cent.edges();
      map<double,double> edgesf = _fixcent.edges();
      map<double,double> edges0 = edges;

      auto curr = edges0.rbegin();
      curr->second = 0.0;
      ++curr;
      pair<double,double> prev = *_centrue.begin();
      double acc = 0.0;
      for ( auto next: _centrue ) {
	double del = next.second/_sumw;
	if ( acc + del > 1.0 - curr->first ) {
	  curr->second = prev.first +
	    ((1.0 - curr->first) - acc)*(next.first - prev.first)/del;
	  ++curr;
	}
	prev = next;
	acc += del;
      }


      MSG_INFO("Cross check of centalty edges:\n");
      MSG_INFO("" << setw(10) << "%" << setw(10) << "true"
               << setw(10) << "estimate" << setw(10) << "fixed");
      for ( auto e : edges0 ) {
        MSG_INFO(""
		 << setw(10) << e.first << setw(10) << e.second
		 << setw(10) << edges[e.first] << setw(10) << edgesf[e.first]);
      }
      

    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hETfwd;
    CentralityBinner<Histo1DPtr> _cent;
    CentralityBinner<Histo1DPtr> _fixcent;
    CentralityBinner<Histo1DPtr> _gencent;
    //@}

    /// Keep track of the actually generated centralities.
    multimap<double, double> _centrue;
    double _sumw;


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_Centrality);

}
