// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/JetShape.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief D0 Run-1 jet shapes measurement
  class D0_1995_I398175 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(D0_1995_I398175);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs(-4.0, 4.0);
      declare(fs, "FS");
      //      FastJets jets(fs, FastJets::ANTIKT, 0.6);
      FastJets jets(fs, FastJets::D0ILCONE, 1.0);
      jets.useInvisibles();
      declare(jets, "Jets");


      // Specify jets pT bins
      _ptedges = {{ 45.0, 70.0, 105.0, 140.0, 1800.0}};

      // Book histograms
      for (size_t ptbin = 0; ptbin < 4; ++ptbin) {
	_jsnames_pT[ptbin] = "JetShape" + to_str(ptbin) ;
	const JetShape jsp(jets, 0.0, 1.0, 10, _ptedges[ptbin], _ptedges[ptbin+1], 0.0, 0.2, PSEUDORAPIDITY);
	declare(jsp, _jsnames_pT[ptbin]);
	_h_Rho_pT_central[ptbin] = bookProfile1D(ptbin+1, 1, 1);
      }

	const JetShape jspfwd0(jets, 0.0, 1.0, 10, 45, 70, 2.5, 3.5, PSEUDORAPIDITY);
	declare(jspfwd0, "JetShapeFwd0");
	const JetShape jspfwd1(jets, 0.0, 1.0, 10, 70, 105, 2.5, 3.5, PSEUDORAPIDITY);
	declare(jspfwd1, "JetShapeFwd1");
	_h_Rho_pT_forward[0] = bookProfile1D(5, 1, 1);
	_h_Rho_pT_forward[1] = bookProfile1D(6, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get jets and require at least one to pass pT and y cuts
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::ptIn(_ptedges.front()*GeV, _ptedges.back()*GeV) );
      MSG_DEBUG("Selecting jets with pT> "<<_ptedges.front());
      MSG_DEBUG("Jet multiplicity before cuts = " << jets.size());
      if (jets.size() == 0){
	MSG_DEBUG("No jets found in required pT and rapidity range");
	vetoEvent;
      }
      const double weight = event.weight();

      // Calculate and histogram jet shapes
      for (size_t ipt = 0; ipt < 4; ++ipt) {
	const JetShape& jsipt = apply<JetShape>(event, _jsnames_pT[ipt]);
	for (size_t ijet = 0; ijet < jsipt.numJets(); ++ijet) {
	  for (size_t rbin = 0; rbin < jsipt.numBins(); ++rbin) {
	    const double r_rho = jsipt.rBinMid(rbin);
	    MSG_DEBUG(ipt << " " << rbin << " (" << r_rho << ") " << jsipt.diffJetShape(ijet, rbin));
	    /// @note Bin width Jacobian factor of 0.7/0.1 = 7 in the differential shapes plot
	    //	    _profhistRho_pT[ipt]->fill(r_rho/0.7, (0.7/0.1)*jsipt.diffJetShape(ijet, rbin), weight);
	    const double r_Psi = jsipt.rBinMax(rbin);
	    MSG_DEBUG(ipt << " " << rbin << " (" << r_rho << ") " << jsipt.intJetShape(ijet, rbin));
	    _h_Rho_pT_central[ipt]->fill(r_Psi/1.0, jsipt.intJetShape(ijet, rbin), weight);
	  }
	}
      }


      const JetShape& jsiptfwd0 = apply<JetShape>(event, "JetShapeFwd0");
      for (size_t ijet = 0; ijet < jsiptfwd0.numJets(); ++ijet) {
	for (size_t rbin = 0; rbin < jsiptfwd0.numBins(); ++rbin) {
	  const double r_Psi = jsiptfwd0.rBinMax(rbin);
	  _h_Rho_pT_forward[0]->fill(r_Psi/1.0, jsiptfwd0.intJetShape(ijet, rbin), weight);
	}
      }

      const JetShape& jsiptfwd1 = apply<JetShape>(event, "JetShapeFwd1");
      for (size_t ijet = 0; ijet < jsiptfwd1.numJets(); ++ijet) {
        for (size_t rbin = 0; rbin < jsiptfwd1.numBins(); ++rbin) {
	  const double r_Psi = jsiptfwd1.rBinMax(rbin);
          _h_Rho_pT_forward[1]->fill(r_Psi/1.0, jsiptfwd1.intJetShape(ijet, rbin), weight);
	}
      }





    }

    /// Normalise histograms etc., after the run
    void finalize() {

      // scale(_h_YYYY, crossSection()/sumOfWeights()); // norm to cross section
      // normalize(_h_YYYY); // normalize to unity

    }

    //@}


  private:


    vector<double> _ptedges;
    string _jsnames_pT[4];
    /// @name Histograms
    //@{
    Profile1DPtr _h_Rho_pT_central[4];
    Profile1DPtr _h_Rho_pT_forward[2];

    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_1995_I398175);


}
