// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {

  /// @brief MPI sensitive di-jet balance variables for W->ejj or W->mujj events.
  class ATLAS_2013_I1216670 : public Analysis {
  public:

    /// @name Constructor
    ATLAS_2013_I1216670()
      : Analysis("ATLAS_2013_I1216670")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms, set up projections for W and jets
    void init() {

      _h_delta_jets_n = bookHisto1D(1, 1, 1);
      _h_delta_jets   = bookHisto1D(2, 1, 1);

      FinalState fs;

      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT >= 20*GeV;

      WFinder w_e_finder(fs, cuts, PID::ELECTRON, 40*GeV, MAXDOUBLE, 0.0*GeV, 0.0, WFinder::CLUSTERNODECAY, 
                         WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(w_e_finder, "W_E_FINDER");

      WFinder w_mu_finder(fs, cuts, PID::MUON, 40*GeV, MAXDOUBLE, 0.0*GeV, 0.0, WFinder::CLUSTERNODECAY, 
                          WFinder::NOTRACK, WFinder::TRANSMASS);
      declare(w_mu_finder, "W_MU_FINDER");

      VetoedFinalState jet_fs(fs);
      jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("W_E_FINDER"));
      jet_fs.addVetoOnThisFinalState(getProjection<WFinder>("W_MU_FINDER"));
      FastJets jets(jet_fs, FastJets::ANTIKT, 0.4);
      declare(jets, "JETS");

    }

    /// Do the analysis
    void analyze(const Event &e) {

      double weight = e.weight();

      const WFinder& w_e_finder  = apply<WFinder>(e, "W_E_FINDER" );
      const WFinder& w_mu_finder = apply<WFinder>(e, "W_MU_FINDER");
      Particle lepton, neutrino;
      Jets all_jets, jets;
  
      // Find exactly 1 W->e or W->mu boson
      if(w_e_finder.bosons().size() == 1 && w_mu_finder.bosons().size() == 0) {
        MSG_DEBUG(" Event identified as W->e nu."); 
        if( !(w_e_finder.mT() > 40*GeV && w_e_finder.constituentNeutrino().Et() > 25.0*GeV) )  vetoEvent;
        lepton = w_e_finder.constituentLepton();
      } else if(w_mu_finder.bosons().size() == 1 && w_e_finder.bosons().size() == 0) {
        MSG_DEBUG(" Event identified as W->mu nu.");  
        if( !(w_mu_finder.mT() > 40*GeV && w_mu_finder.constituentNeutrino().Et() > 25.0*GeV) )  vetoEvent;
        lepton = w_mu_finder.constituentLepton();
      } else {
        MSG_DEBUG(" No W found passing cuts.");    
        vetoEvent;
      }

      all_jets = apply<FastJets>(e, "JETS").jetsByPt(Cuts::pt>20.0*GeV && Cuts::absrap<2.8);

      // Remove jets DeltaR < 0.5 from W lepton
      for(Jets::iterator it = all_jets.begin(); it != all_jets.end(); ++it) {
        double distance = deltaR( lepton, (*it) );
        if(distance < 0.5) {
          MSG_DEBUG("   Veto jet DeltaR " << distance << " from W lepton"); 
        } else {
          jets.push_back(*it);
        }
      }
      
      // Exactly two jets required
      if( jets.size() != 2 )  vetoEvent;

      // Calculate analysis quantities from the two jets
      double delta_jets = (jets.front().momentum() + jets.back().momentum()).pT();
      double total_pt = jets.front().momentum().pT() + jets.back().momentum().pT();
      double delta_jets_n = delta_jets / total_pt;
      
      _h_delta_jets->fill(   delta_jets,   weight ); // Jet pT balance
      _h_delta_jets_n->fill( delta_jets_n, weight ); // Jet pT balance, normalised by scalar dijet pT
      
    }

    /// Finalize
    void finalize() {

      // Data is normalised to 0.03 and 3
      normalize(_h_delta_jets_n, 0.03);
      normalize(_h_delta_jets  , 3.0 );

    }

    //@}

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_delta_jets_n;
    Histo1DPtr _h_delta_jets;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1216670);

}
