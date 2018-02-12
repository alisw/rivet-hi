// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {

  // TODO this calculation needs checked!
  double impact(const FourMomentum& a, const FourMomentum& b) {
    const Vector3 a3 = a.vector3();
    const Vector3 b3 = b.vector3();

    double impact = 0;
    if (b3.polarRadius() !=0) {
      impact = (a3).cross((a3-b3)).polarRadius() / (b3).polarRadius();
    }
    return impact;
  } 
  
  /// @brief Add a short analysis description here
  class ALEPH_2016_I1492968 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ALEPH_2016_I1492968);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      const FinalState fs;
      addProjection(fs, "FS");

      FastJets jets(fs, FastJets::GENKTEE, 0.5, JetAlg::NO_MUONS, JetAlg::ALL_INVISIBLES);
      //FastJets jets(fs, FastJets::ANTIKT, 0.5, JetAlg::NO_MUONS, JetAlg::ALL_INVISIBLES);
      addProjection(jets, "Jets");

      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      addProjection(mu_id, "MUONS");
  
      addProjection(MissingMomentum(fs), "MissingMomenta");      
      // Book histograms
      //_h_costheta = bookHisto1D(2, 1, 1);
      _h_m_OS = bookHisto1D(3, 1, 1);
      _h_m_SS = bookHisto1D(5, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // B-jets
      const Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 5*GeV); // tODO jet eta?
      const Jets bjets = filter_select(jets,  [](const Jet& j) { return j.bTagged(); });
      if (bjets.size()<2) vetoEvent;

      // Muons
      const Particles all_muons = applyProjection<IdentifiedFinalState>(event, "MUONS").particles(Cuts::pT>2.5/GeV, cmpMomByE);
      const Particles b_muons = filter_select(all_muons, [](const Particle& m) {return cos(m.theta()) < 0.7; });
      if (b_muons.size()<2) vetoEvent;

      // Missing energy cut
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MissingMomenta");
      double Pmiss = met.missingMomentum().p();
      if (Pmiss/GeV>18) vetoEvent;
    
      // Impact paarameter considerations
      double b_muon_0_impactdistance = min(impact(b_muons[0].origin(), bjets[0].momentum()),impact(b_muons[0].origin(), bjets[1].momentum()));
      double b_muon_1_impactdistance = min(impact(b_muons[1].origin(), bjets[0].momentum()),impact(b_muons[1].origin(), bjets[1].momentum()));

      // Impact parameter cut
      if ((b_muon_0_impactdistance > 0.1) || (b_muon_1_impactdistance > 0.1)) vetoEvent;
      
      FourMomentum dimuon = b_muons[0].momentum() + b_muons[1].momentum();
      
      // Same sign
      if (b_muons[0].charge()*b_muons[1].charge()>0) {
        _h_m_SS->fill( dimuon.mass()/GeV, weight);
      }
      // Opposite sign
      else {
        _h_m_OS->fill( dimuon.mass()/GeV, weight);
        //
        //FourMomentum muonminus;
        //if (b_muons[0].charge() < 0)  muonminus = b_muons[0].momentum();
        //else muonminus = b_muons[1].momentum();
        
        //const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(-dimuon.betaVec());
        //FourMomentum boostedmuon = cms_boost.transform(muonminus);

        //double cosmuonboosted = boostedmuon.vector3().dot(cms_boost.betaVec())
          /// (boostedmuon.vector3().mod()*cms_boost.betaVec().mod());

        //_h_costheta->fill( cosmuonboosted, weight);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      //normalize(_h_costheta);

      // Normalize to data according to Arno.
      normalize(_h_m_OS, 1387);
      normalize(_h_m_SS, 1047);

    }

    //@}


    /// @name Histograms
    //@{
    //Histo1DPtr _h_costheta;
    Histo1DPtr _h_m_OS;
    Histo1DPtr _h_m_SS;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ALEPH_2016_I1492968);


}
