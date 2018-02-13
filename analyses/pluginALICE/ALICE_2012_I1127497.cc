// -*- C++ -*-
#include "pluginALICE/HeavyIonAnalysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include <fstream>

#define NHISTOS 15

#define _USE_MATH_DEFINES
#include <math.h>

// maybe use boost

namespace Rivet {


  class ALICE_2012_I1127497 : public HeavyIonAnalysis {

  public:

    ALICE_2012_I1127497() :
      HeavyIonAnalysis("ALICE_2012_I1127497")
    {  }


    void init() {
      HeavyIonAnalysis::init();

      // Charged final states with |eta| < 0.5 and pT > 150 MeV
      const Cut& cut = Cuts::abseta < 0.5 && Cuts::pT > 150*MeV;
      const ChargedFinalState cfs(cut);
      addProjection(cfs,"CFS");

      // Set centrality method
      addCentralityMethod(HeavyIonAnalysis::ImpactParameter, 10000, "method1");

      // @debug Member "log" is declared below.
      //log.open("/data/mcplots/mcplots/scripts/mcprod/testarea/logs/log.log");

      // @note and @todo In principle, only ONE pp histogram needed since centrality does not make a difference here.
      // However, in some cases it had to be re-binned since some of the binnings of histograms in this analysis
      // differ from each other. This is therefore the "brute-force" but easy-to-implement way. Making this more
      // efficient this could be part of the optimization procedure in the end. Now, initialize AA and pp dummy
      // histograms.

        for( size_t d = 0; d < NHISTOS; ++d ) {
          m_sumOfWeights[1][d] = 0;
          _histNch[1][d] = bookHisto1D(d + 1, 1, 1);
          _histRAA[d] = bookScatter2D(d+16, 1, 1);

          m_sumOfWeights[0][d] = 0;
          // @note Give a proper name for the dummy histograms. They have to incorporate at least the name of the
          // analysis ALICE_2012_I1127497 and must of course not have the same path as another analysis object used
          // in this analysis.
          std::string name = _histNch[1][d]->name();
	  //cout << "PbPb name: " << _histNch[1][d]->name() << endl;
          std::string namePP = name + "-pp";
          _histNch[0][d] = bookHisto1D( namePP, refData(d + 1, 1, 1) ); // binning is also taken from ref data
	  //cout << "pp name: " << _histNch[0][d]->name() << endl;
        }

	_histWeights[0] = bookHisto1D( "weights_pp", 15, 0.0, 15.0, "weights_pp", "xlabel", "ylabel");
	_histWeights[1] = bookHisto1D( "weights_PbPb", 15, 0.0, 15.0, "weights_PbPb", "xlabel", "ylabel");

      // Centrality regions, 2 following numbers are boundaries for a certain region. Note, that some
      // regions overlap with other regions. The current implementation of centrality bins may not be
      // the most efficient and flexible one. Take this into account during the optimization procedure.
      m_centrRegions.clear();
      m_centrRegions += {{0., 0.05, 0.05, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5, 0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0., 0.1, 0., 0.2, 0.2, 0.4, 0.4, 0.6, 0.4, 0.8, 0.6, 0.8}};


    }

    void analyze(const Event& event) {

      const double weight = event.weight();

      // final state particles with at least pT = 50 MeV in eta range of |eta| < 0.8
      const ChargedFinalState& charged = applyProjection<ChargedFinalState>(event, "CFS");
      Particles chargedParticles = charged.particlesByPt();

      // @note Choose this convention here since JEWEL return a centrality value of -1 for a vacuum run.
      //float centr = -1;

      // @note It has to be checked explicitly, if the pointer is sane and is no null_ptr. The standard
      // constructor of HepMC does not initialize the pointer and returns a null_ptr in case no heavy-ion
      // was written to the HepMC file. This is e.g. the case for the current pp generators. Trying to
      // access this naivly may therefore lead to a segmentation fault.
      //if( event.genEvent()->heavy_ion() ) {
        // JEWEL saves the centrality as the impact parameter. This should be changed later in JEWEL directly.
        // Note, that changes of HepMC may be required as well. However, given the current implementation,
        // the centrality can be accessed via
        //centr = event.genEvent()->heavy_ion()->impact_parameter();
      //}

      // Check event type
      double centr = -1.;
      if (HeavyIonAnalysis::is_heavy_ion(event)) {
	centr = centrality(event, "method1") / 100.0;
      }
      std::cout << "Centrality: " << centr << std::endl;

      // Veto event for too large centralities since those are not used in the analysis at all.
      if( centr > 0.8 )
        vetoEvent;


      size_t pp_AA;
      float pT;
      vector< size_t > indices;
      indices.clear();

      if( centr < 0. ) {
        // Use this as a flag to decide later which histograms should be filled. pp_AA = 0 corresponds
        // to pp beam (and in case of JEWEL to a vacuum run).
	//std::cout << "Filling pp indices" << std::endl;
        pp_AA = 0;
        for(size_t i = 0; i < NHISTOS; ++i)
          // Push all indices in case of pp dummy histograms.
          indices.push_back(i);

      } else {
	//std::cout << "Filling PbPb indices" << std::endl;
        pp_AA = 1;
        for( size_t i = 0; i < NHISTOS; ++i ) {
          // Check centrality bins and push corresponding indices of AA histograms.
          if( inRange( centr, m_centrRegions[ 2 * i ], m_centrRegions[ 2 * i + 1] ) ) {
            indices.push_back(i);
          }
        }
      }


      // Fill the right histograms and add weights based on the flag pp_AA and indices pushed to "indices".
      for( size_t i = 0; i < indices.size(); ++i ) {
        m_sumOfWeights[pp_AA][indices.at(i)] += weight;
	_histWeights[pp_AA]->fill(indices.at(i), weight * 2. * M_PI * 1.6 );
        foreach (const Particle& p, chargedParticles) {
          pT = p.pT()/GeV;
          if(pT < 50.) {
	    //std::cout << "Filling histogram with pT=" << pT << " and weight=" << weight * ( 1 / pT ) << std::endl;
            _histNch[pp_AA][indices.at(i)]->fill( pT , weight * ( 1 / pT ) );
	    //cout << "Weight: " << weight << endl;
	    //cout << "pT: " << pT << endl;
          }
        }
      }

    }


    void finalize() {
      // Right scaling of the histograms with their individual weights.
      /*for( size_t pp_AA = 0; pp_AA < 2; ++pp_AA ) {
        for( size_t i = 0; i < NHISTOS; ++i ) {
          if( m_sumOfWeights[pp_AA][i] > 0 ) //@note is this necessary or can this be handeled by the scaling method?
            scale( _histNch[pp_AA][i], ( 1./m_sumOfWeights[pp_AA][i] / 2. /  M_PI / 1.6 ) ); // TODO PI
        }
	}*/
    }


    // @note First time a post() method is implemented and used by Rivet. This is only available on the
    // mcplots-alice-dev.cern.ch machine because the Rivet installation here is the only one which implements
    // Rivet::Analysis::post() and Rivet::AnalysisHandler::post()
    void post() {

      // Right scaling of the histograms with their individual weights.
      for( size_t pp_AA = 0; pp_AA < 2; ++pp_AA ) {
        for( size_t i = 0; i < NHISTOS; ++i ) {
	  std::cout << _histNch[pp_AA][i]->numEntries() << std::endl;
	  std::cout << _histWeights[pp_AA]->numEntries() << std::endl;
          if( m_sumOfWeights[pp_AA][i] > 0 ) //@note is this necessary or can this be handeled by the scaling method?
	    std::cout << "Scaling " << _histNch[pp_AA][i]->name() << " using weight " << _histWeights[pp_AA]->bin(i).sumW() << std::endl;
            scale( _histNch[pp_AA][i], ( 1. / _histWeights[pp_AA]->bin(i).sumW() ) ); // TODO PI
        }
      }

      // Divide AA and pp histograms to obtain R_AA histograms.
      for( size_t i = 0; i < NHISTOS; ++i) {
	//cout << "_histNch[1][" << i << "]->bin(0).numEntries(): " << _histNch[1][i]->bin(0).numEntries() << endl;
	//cout << "_histNch[0][" << i << "]->bin(0).numEntries(): " << _histNch[0][i]->bin(0).numEntries() << endl;
	//cout << "_histRAA[" << i << "]->bin(0): " << _histRAA[i]->bin(0) << endl << endl;
        divide( _histNch[1][i], _histNch[0][i], _histRAA[i] );
	//cout << "_histRAA[" << i << "]->bin(0): " << _histRAA[i]->bin(0) << endl << endl;
      }
      // @debug
      //log.close();
    }

  private:

    //std::ofstream validationFile;
    Histo1DPtr _histNch[2][NHISTOS];
    double m_sumOfWeights[2][NHISTOS];
    Histo1DPtr _histWeights[2];

    Scatter2DPtr _histRAA[NHISTOS];
    std::vector<float> m_centrRegions;
    // @debug
    //std::ofstream log;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALICE_2012_I1127497);


}
