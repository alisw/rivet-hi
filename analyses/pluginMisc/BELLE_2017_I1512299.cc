// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BELLE_2017_I1512299 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(BELLE_2017_I1512299);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableFinalState(), "UFS");

      // Book histograms
      _h_w      = bookHisto1D(1, 1, 1);
      _h_costhv = bookHisto1D(2, 1, 1);
      _h_costhl = bookHisto1D(3, 1, 1);
      _h_chi    = bookHisto1D(4, 1, 1);

    }


    /// Perform the per-event analysis
    bool analyzeDecay(Particle mother, vector<int> ids) {
      // There is no point in looking for decays with less particles than to be analysed
      if (mother.children().size() == ids.size()) {
        bool decayfound = true;
        for (int id : ids) {
          if (!contains(mother, id)) decayfound = false;
        }
        return decayfound;
      }
      return false;
    }

    bool contains(Particle& mother, int id) {
      return any(mother.children(), HasPID(id));
    }


    double recoilW(const Particle& mother) {
      FourMomentum lepton, neutrino, meson, q;
      foreach(const Particle& c, mother.children()) {
        if (c.isNeutrino()) neutrino=c.mom();
        if (c.isLepton() &! c.isNeutrino()) lepton =c.mom();
        if (c.isHadron()) meson=c.mom();
      }
      q = lepton + neutrino; //no hadron before
      double mb2= mother.mom()*mother.mom();
      double mD2 = meson*meson;
      return (mb2 + mD2 - q*q )/ (2. * sqrt(mb2) * sqrt(mD2) );
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      FourMomentum pl, pnu, pB, pD, pDs, ppi;
      // Iterate of B0bar mesons
      for(const Particle& p : apply<UnstableFinalState>(event, "UFS").particles(Cuts::pid==-511)) {
        pB = p.momentum();
        // Find semileptonic decays
        if (analyzeDecay(p, {PID::DSTARPLUS,-12,11}) || analyzeDecay(p, {PID::DSTARPLUS,-14,13}) ) {
          _h_w->fill(recoilW(p), event.weight());
          // Get the necessary momenta for the angles
          bool foundDdecay=false;
          for (const Particle c : p.children()) {
            if ( (c.pid() == PID::DSTARPLUS)  && (analyzeDecay(c, {PID::PIPLUS, PID::D0}) || analyzeDecay(c, {PID::PI0, PID::DPLUS})) ) {
              foundDdecay=true;
              pDs = c.momentum();
              for (const Particle dc : c.children()) {
                if (dc.hasCharm()) pD = dc.momentum(); 
                else ppi = dc.momentum(); 
              }
            }
            if (c.pid() ==  11 || c.pid() ==  13) pl  = c.momentum();
            if (c.pid() == -12 || c.pid() == -14) pnu = c.momentum();
          }
          // This is the angle analysis
          if (foundDdecay) {

            // First boost all relevant momenta into the B-rest frame
            const LorentzTransform B_boost = LorentzTransform::mkFrameTransformFromBeta(pB.betaVec());
            // Momenta in B rest frame:
            FourMomentum lv_brest_Dstar = B_boost.transform(pDs);//lab2brest(gp_Dstar.particle.p());
            FourMomentum lv_brest_w     = B_boost.transform(pB - pDs); //lab2brest(p_lv_w);
            FourMomentum lv_brest_D     = B_boost.transform(pD); //lab2brest(gp_D.particle.p());
            FourMomentum lv_brest_lep   = B_boost.transform(pl); //lab2brest(gp_lep.p());
            
            const LorentzTransform Ds_boost = LorentzTransform::mkFrameTransformFromBeta(pDs.betaVec());
            FourMomentum lv_Dstarrest_D     = Ds_boost.transform(lv_brest_D);
            const LorentzTransform W_boost  = LorentzTransform::mkFrameTransformFromBeta((pB-pDs).betaVec());
            FourMomentum lv_wrest_lep       = W_boost.transform(lv_brest_lep);

            double cos_thetaV = cos(lv_brest_Dstar.p3().angle(lv_Dstarrest_D.p3()));
            _h_costhv->fill(cos_thetaV, event.weight());
            
            double cos_thetaL = cos(lv_brest_w.p3().angle(lv_wrest_lep.p3()));
            _h_costhl->fill(cos_thetaL, event.weight());

            Vector3 LTrans = lv_wrest_lep.p3()   - cos_thetaL*lv_wrest_lep.p3().perp()*lv_brest_w.p3().unit();
            Vector3 VTrans = lv_Dstarrest_D.p3() - cos_thetaV*lv_Dstarrest_D.p3().perp()*lv_brest_Dstar.p3().unit();
            float chi = atan2(LTrans.cross(VTrans).dot(lv_brest_w.p3().unit()), LTrans.dot(VTrans));
            if(chi < 0) chi += TWOPI; 

            _h_chi->fill(chi, event.weight());

            //const LorentzTransform W_boost = LorentzTransform::mkFrameTransformFromBeta((pl+pnu).betaVec());
            //const LorentzTransform D_boost = LorentzTransform::mkFrameTransformFromBeta((pD+ppi).betaVec());

            //FourMomentum pl_t = FourMomentum(W_boost.transform(pl));
            //FourMomentum pD_t = FourMomentum(D_boost.transform(pD));
            //double thetal = (pl+pnu).angle(pl_t);
            //double thetav = (pD+ppi).angle(pD_t);
            //_h_costhv->fill(cos(thetav), event.weight());
            //_h_costhl->fill(cos(thetal), event.weight());
          }
        }
      }
    }
        //else if (analyzeDecay(p, {413,-14,13}) ) {
          //_h_w->fill(recoilW(p), event.weight());
        //}

    /// Normalise histograms etc., after the run
    void finalize() {
    
      double GAMMA_B0 = 4.32e-13; // Total width in GeV, calculated from mean life time of 1.52 pico seconds
      double BR_B0_DSPLUS_ELL_NU = 0.0495; // Branching fraction from the same paper for B0bar to D*+ ell nu
      double NORM = GAMMA_B0 * BR_B0_DSPLUS_ELL_NU; // Normalise histos to partial width
      normalize(_h_w,      NORM); 
      normalize(_h_costhv, NORM); 
      normalize(_h_costhl, NORM); 
      normalize(_h_chi,    NORM);
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr _h_w;
    Histo1DPtr _h_costhv;
    Histo1DPtr _h_costhl;
    Histo1DPtr _h_chi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2017_I1512299);


}
