// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Track-based underlying event at 13 TeV in ATLAS
  class ATLAS_2017_I1509919 : public Analysis {
  public:

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2017_I1509919);


    // Pre-run histogram and projection booking
    void init() {

      declare(ChargedFinalState(Cuts::abseta < 2.5 && Cuts::pT > 500*MeV), "CFS500");

      for (size_t iR = 0; iR < NREGIONS; ++iR) {
        // Nch profiles vs. pT_lead
        _hist_nch[iR] = bookProfile1D(1, 1 + iR, 1);
        // pTsum profiles vs. pT_lead
        _hist_ptsum[iR] = bookProfile1D(2, 1 + iR, 1);
        // <pT> profiles vs pT_lead (not measured for trans diff)
        if (iR != kTransDiff)  _hist_ptavg[iR] = bookProfile1D(3, 1 + iR, 1);
        // <pT> profiles vs. Nch (not measured for trans diff)
        if (iR != kTransDiff)  _hist_dn_dpt[iR] = bookProfile1D(4, 1 + iR, 1);
        // Only measured for trans max/min
        if ( (iR == kTransMax) || (iR == kTransMin) )  _hist_dn_dpt2[iR] = bookProfile1D(5, 1 + iR, 1);
      }

      for (size_t iC = 0; iC < NCUTS; ++iC) {
        // Nch vs. Delta(phi) profiles
        _hist_N_vs_dPhi[iC] = bookProfile1D(6, 1 + iC, 1);
        // pT vs. Delta(phi) profiles
        _hist_pT_vs_dPhi[iC] = bookProfile1D(7, 1 + iC, 1);
        //ptLead histos only for 1 and 5 GeV cuts
        if ( (iC == 0) || (iC == 1) )  _hist_ptLead[iC] = bookHisto1D(8, 1 + iC, 1);
        // 
        _counters[iC] = bookCounter("Ctr_cut_" + toString(iC));
      }

    }


    void analyze(const Event& event) {

      // Get charged particles (tracks) with pT > 500 MeV
      const ChargedFinalState& charged500 = apply<ChargedFinalState>(event, "CFS500");
      const Particles& particlesAll = charged500.particlesByPt();
      MSG_DEBUG("Num tracks: " << particlesAll.size());

      const Cut& pcut = ( (Cuts::abspid != PID::SIGMAMINUS) && (Cuts::abspid != PID::SIGMAPLUS) &&
                          (Cuts::abspid != PID::XIMINUS)    && (Cuts::abspid != PID::OMEGAMINUS) );
      const Particles& particles = charged500.particlesByPt(pcut);
      MSG_DEBUG("Num tracks without strange baryons: " <<  particles.size());

      // Require at least one track in the event for pTlead histograms
      if (particles.empty()) vetoEvent;
      const double weight = event.weight();
      for (size_t iC = 0; iC < 2; ++iC) {
        if (particles[0].pT() < PTCUTS[iC]*GeV) continue;
        _counters[iC]->fill(weight);
        _hist_ptLead[iC]->fill( particles[0].pT()/GeV, weight);
      }

      // Require at least one track in the event with pT >= 1 GeV for the rest
      if (particles[0].pT() < 1*GeV) vetoEvent;

      // Identify leading track and its phi and pT
      const Particle& p_lead = particles[0];
      const double philead = p_lead.phi();
      const double etalead = p_lead.eta();
      const double pTlead  = p_lead.perp();
      MSG_DEBUG("Leading track: pT = " << pTlead << ", eta = " << etalead << ", phi = " << philead);

      // Iterate over all particles and count particles and scalar pTsum in three basic regions
      vector<double> num(NREGIONS, 0), ptSum(NREGIONS, 0.0), avgpt(NREGIONS, 0.0);

      // Temporary histos that bin Nch and pT in dPhi.
      Histo1D hist_num_dphi(*_hist_N_vs_dPhi[0].get(), "/hist_num_dphi");
      Histo1D hist_pt_dphi(*_hist_pT_vs_dPhi[0].get(), "/hist_pt_dphi");
      hist_num_dphi.reset();
      hist_pt_dphi .reset();

      int    tmpnch[2]   = {0,0};
      double tmpptsum[2] = {0,0};
      for (const Particle& p : particles) {
        const double pT   = p.pT()/GeV;
        const double dPhi = deltaPhi(philead, p.phi()); // in range (0,pi)
        const int    ir   = region_index(dPhi); // gives just toward/away/trans

        // Toward/away/trans region: just count
        num  [ir] += 1;
        ptSum[ir] += pT;

        // Determine which transverse side
        if (ir == kTrans) {
          const size_t iside = (mapAngleMPiToPi(p.phi() - philead) > 0) ? 0 : 1;
          tmpnch  [iside] += 1;
          tmpptsum[iside] += p.pT();
        }

        // Fill temp histos to bin Nch and pT in dPhi
        if (p.genParticle() != p_lead.genParticle()) { // We don't want to fill all those zeros from the leading track...
          hist_num_dphi.fill(dPhi/M_PI*180, 1);
          hist_pt_dphi .fill(dPhi/M_PI*180, pT/GeV);
        }
      }

      // Construct max/min/diff regions
      num[kTransMax ]   = std::max(tmpnch[0], tmpnch[1]);
      num[kTransMin ]   = std::min(tmpnch[0], tmpnch[1]);
      num[kTransDiff]   = num[kTransMax ] - num[kTransMin ];
      ptSum[kTransMax ] = std::max(tmpptsum[0], tmpptsum[1]);
      ptSum[kTransMin ] = std::min(tmpptsum[0], tmpptsum[1]);
      ptSum[kTransDiff] = ptSum[kTransMax ] - ptSum[kTransMin ];
      avgpt[kToward]    = (num[kToward] > 0 ) ? ptSum[kToward] / num[kToward] : 0. ;
      avgpt[kAway]      = (num[kAway  ] > 0 ) ? ptSum[kAway]   / num[kAway]   : 0. ;
      avgpt[kTrans]     = (num[kTrans ] > 0 ) ? ptSum[kTrans]  / num[kTrans]  : 0. ;
      // Avg pt max/min regions determined according sumpt max/min
      int sumptMaxRegID = (tmpptsum[0] >  tmpptsum[1]) ? 0 : 1 ;
      int sumptMinRegID = (sumptMaxRegID == 0) ? 1 : 0;
      avgpt[kTransMax ] = (tmpnch[sumptMaxRegID] > 0) ? tmpptsum[sumptMaxRegID] / tmpnch[sumptMaxRegID] : 0.;
      avgpt[kTransMin ] = (tmpnch[sumptMinRegID] > 0) ? tmpptsum[sumptMinRegID] / tmpnch[sumptMinRegID] : 0.;
      avgpt[kTransDiff] = ((tmpnch[sumptMaxRegID] > 0) && (tmpnch[sumptMinRegID] > 0)) ? avgpt[kTransMax ] - avgpt[kTransMin ] : 0.;


      // Now fill underlying event histograms

      // The densities are calculated by dividing the UE properties by dEta*dPhi
      // -- each basic region has a dPhi of 2*PI/3 and dEta is two times 2.5
      // min/max/diff regions are only half of that
      const double dEtadPhi[NREGIONS] = { 2*2.5 * 2*PI/3.0, 2*2.5 * 2*PI/3.0, 2*2.5 * 2*PI/3.0,
                                          2*2.5 *   PI/3.0, 2*2.5 *   PI/3.0, 2*2.5 *   PI/3.0 };
      for (size_t iR = 0; iR < NREGIONS; ++iR) {

        _hist_nch  [iR]->fill(pTlead/GeV, num[iR]   /dEtadPhi[iR]     , weight);
        _hist_ptsum[iR]->fill(pTlead/GeV, ptSum[iR] /GeV/dEtadPhi[iR] , weight);

        // <pT> profiles vs. pT_lead (first 3 are the same!)
        switch (iR) {
        case kToward    :
        case kAway      :
        case kTrans     :
          if (num[iR] > 0) _hist_ptavg[iR]->fill(pTlead/GeV, avgpt[iR]/GeV, weight);
          break;
        case kTransMax  :
          if (tmpnch[sumptMaxRegID] > 0) _hist_ptavg[iR]->fill(pTlead/GeV, avgpt[iR]/GeV, weight);
          break;
        case kTransMin  :
          if (tmpnch[sumptMinRegID] > 0) _hist_ptavg[iR]->fill(pTlead/GeV, avgpt[iR]/GeV, weight);
          break;
        case kTransDiff :
          break;
        default: //should not get here!!!
          MSG_WARNING("Unknown region in <pT> profiles vs.pt lead switch!!! : " << iR);
        }

        // <pT> profiles vs. Nch (first 3 are the same!)
        switch (iR) {
        case kToward    :
        case kAway      :
        case kTrans     :
          if (num[iR] > 0) _hist_dn_dpt[iR]->fill(num[iR] , avgpt[iR]/GeV, weight);
          break;
        case kTransMax  :
          if (tmpnch[sumptMaxRegID] > 0) {
            _hist_dn_dpt [iR]->fill(num[kTrans]          , avgpt[iR]/GeV, weight);
            _hist_dn_dpt2[iR]->fill(tmpnch[sumptMaxRegID], avgpt[iR]/GeV, weight);
          }
          break;
        case kTransMin  :
          if (tmpnch[sumptMinRegID] > 0) {
            _hist_dn_dpt [iR]->fill(num[kTrans]          , avgpt[iR]/GeV, weight);
            _hist_dn_dpt2[iR]->fill(tmpnch[sumptMinRegID], avgpt[iR]/GeV, weight);
          }
          break;
        case kTransDiff :
          break;
        default : //should not get here!!!
          MSG_INFO("unknown region in <pT> profiles vs. nch switch!!! : " <<  iR);
        }

      }


      // Update the "proper" dphi profile histograms
      // Note that we fill dN/dEtadPhi: dEta = 2*2.5, dPhi = 2*PI/nBins
      const double dEtadPhi2 = (2*2.5 * 2) * M_PI/180.;
      for (size_t i = 0; i < hist_num_dphi.numBins(); ++i) {

        // First Nch
        double mean = hist_num_dphi.bin(i).xMid() ;
        double value = 0.;
        if (hist_num_dphi.bin(i).numEntries() > 0) {
          mean  = hist_num_dphi.bin(i).xMean() ;
          value = hist_num_dphi.bin(i).area()/hist_num_dphi.bin(i).xWidth()/dEtadPhi2;
        }
        for (size_t iC = 0; iC < NCUTS; ++iC) {
          if (pTlead >= PTCUTS[iC]*GeV) _hist_N_vs_dPhi[iC] ->fill(mean, value, weight);
        }

        // Then pT
        mean = hist_pt_dphi.bin(i).xMid();
        value = 0.;
        if (hist_pt_dphi.bin(i).numEntries() > 0) {
          mean  = hist_pt_dphi.bin(i).xMean() ;
          value = hist_pt_dphi.bin(i).area()/hist_pt_dphi.bin(i).xWidth()/dEtadPhi2;
        }
        for (size_t iC = 0; iC < NCUTS; ++iC) {
          if (pTlead >= PTCUTS[iC]*GeV) _hist_pT_vs_dPhi[iC] ->fill(mean, value, weight);
        }
      }

    }


    void finalize() {
      for (size_t iC = 0; iC < NCUTS; ++iC) {
        if (iC == 0 || iC == 1) scale(_hist_ptLead[iC], 1.0/_counters[iC]->sumW());
      }
    }


  private:

    enum regionID {
      kToward = 0,
      kAway,
      kTrans,
      kTransMax,
      kTransMin,
      kTransDiff,
      NREGIONS
    };

    // Little helper function to identify basic Delta(phi) regions: toward/away/trans
    int region_index(double dphi) {
      assert(inRange(dphi, 0.0, PI, CLOSED, CLOSED));
      if (dphi < PI/3.0)    return kToward;
      if (dphi < 2*PI/3.0)  return kTrans;
      return kAway;
    }

    const static size_t NCUTS = 3;
    const vector<double> PTCUTS = { 1., 5., 10. };


    /// @name Histograms
    //@{

    // Nch, sumpT, avgpT profiles vs. pTlead
    Profile1DPtr _hist_nch   [NREGIONS]; //for regions: all 6 regions
    Profile1DPtr _hist_ptsum [NREGIONS]; //for regions: all 6 regions
    Profile1DPtr _hist_ptavg [NREGIONS]; //for regions: trans towards/away/all/min/max
    // Nch, sumpT, avgpT profiles vs. Nch
    Profile1DPtr _hist_dn_dpt [NREGIONS]; //regions: towards/away/ vs nch(region) & trans all/min/max vs nch(trans)
    Profile1DPtr _hist_dn_dpt2[NREGIONS]; //regions: trans min/max vs. nch(region)

    Profile1DPtr _hist_N_vs_dPhi [NCUTS];
    Profile1DPtr _hist_pT_vs_dPhi[NCUTS];
    Histo1DPtr   _hist_ptLead[NCUTS]; //for 1,5 GeV cuts only
    CounterPtr   _counters[NCUTS];

    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1509919);

}
