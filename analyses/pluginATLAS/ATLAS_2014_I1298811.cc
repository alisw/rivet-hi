// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2014_I1298811 : public Analysis {
  public:


    ATLAS_2014_I1298811()
      : Analysis("ATLAS_2014_I1298811") {    }


    void init() {
      // Configure projections
      const FinalState fs(-4.8, 4.8, 0*MeV);
      declare(fs, "FS");
      const FastJets jets(fs, FastJets::ANTIKT, 0.4);
      declare(jets, "Jets");

      // Book histograms
      for (size_t itopo = 0; itopo < 2; ++itopo) {
        // Profiles
        for (size_t iregion = 0; iregion < 3; ++iregion) {
          _p_ptsumch_vs_ptlead[itopo][iregion] = bookProfile1D(1+iregion, 1, itopo+1);
          _p_nch_vs_ptlead[itopo][iregion] = bookProfile1D(4+iregion, 1, itopo+1);
        }
        _p_etsum25_vs_ptlead_trans[itopo] = bookProfile1D(7, 1, itopo+1);
        _p_etsum48_vs_ptlead_trans[itopo] = bookProfile1D(8, 1, itopo+1);
        _p_chratio_vs_ptlead_trans[itopo] = bookProfile1D(9, 1, itopo+1);
        _p_ptmeanch_vs_ptlead_trans[itopo] = bookProfile1D(10, 1, itopo+1);
        // 1D histos
        for (size_t iregion = 0; iregion < 3; ++iregion) {
          for (size_t ipt = 0; ipt < 4; ++ipt) {
            _h_ptsumch[ipt][itopo][iregion] = bookHisto1D(13+3*ipt+iregion, 1, itopo+1);
            _h_nch[ipt][itopo][iregion] = bookHisto1D(25+3*ipt+iregion, 1, itopo+1);
          }
        }
      }
      _p_ptmeanch_vs_nch_trans[0] = bookProfile1D(11, 1, 1);
      _p_ptmeanch_vs_nch_trans[1] = bookProfile1D(12, 1, 1);

    }



    void analyze(const Event& event) {
      // Find the jets with pT > 20 GeV and *rapidity* within 2.8
      /// @todo Use Cuts instead rather than an eta cut in the proj and a y cut after
      const Jets alljets = apply<FastJets>(event, "Jets").jetsByPt(20*GeV);
      Jets jets;
      foreach (const Jet& j, alljets)
        if (j.absrap() < 2.8) jets.push_back(j);
      // Require at least one jet in the event
      if (jets.empty()) vetoEvent;

      // Get the event weight since we will be filling some histos
      const double weight = event.weight();

      // Identify the leading jet and its phi and pT
      const FourMomentum plead = jets[0].momentum();
      const double philead = plead.phi();
      const double etalead = plead.eta();
      const double ptlead  = plead.pT();
      MSG_DEBUG("Leading object: pT = " << ptlead << ", eta = " << etalead << ", phi = " << philead);

      // Sum particle properties in the transverse regions
      int tmpnch[2] = {0,0};
      double tmpptsum[2] = {0,0};
      double tmpetsum48[2] = {0,0};
      double tmpetsum25[2] = {0,0};
      const Particles particles = apply<FinalState>(event, "FS").particles();
      foreach (const Particle& p, particles) {
        // Only consider the transverse region(s), not toward or away
        if (!inRange(deltaPhi(p.phi(), philead), PI/3.0, TWOPI/3.0)) continue;
        // Work out which transverse side this particle is on
        const size_t iside = (mapAngleMPiToPi(p.phi() - philead) > 0) ? 0 : 1;
        MSG_TRACE(p.phi() << " vs. " << philead << ": " << iside);
        // Charged or neutral particle?
        const bool charged = PID::threeCharge(p.pdgId()) != 0;
        // Track observables
        if (charged && fabs(p.eta()) < 2.5 && p.pT() > 500*MeV) {
          tmpnch[iside] += 1;
          tmpptsum[iside] += p.pT();
        }
        // Cluster observables
        if ((charged && p.p3().mod() > 200*MeV) || (!charged && p.p3().mod() > 500*MeV)) {
          tmpetsum48[iside] += p.pT();
          if (fabs(p.eta()) < 2.5) tmpetsum25[iside] += p.pT();
        }
      }

      // Construct tot/max/min counts (for trans/max/min, indexed by iregion)
      const int nch[3] = { tmpnch[0] + tmpnch[1],
                           std::max(tmpnch[0], tmpnch[1]),
                           std::min(tmpnch[0], tmpnch[1]) };
      const double ptsum[3] = { tmpptsum[0] + tmpptsum[1],
                                std::max(tmpptsum[0], tmpptsum[1]),
                                std::min(tmpptsum[0], tmpptsum[1]) };
      const double etsum48[3] = { tmpetsum48[0] + tmpetsum48[1],
                                  std::max(tmpetsum48[0], tmpetsum48[1]),
                                  std::min(tmpetsum48[0], tmpetsum48[1]) };
      const double etsum25[3] = { tmpetsum25[0] + tmpetsum25[1],
                                  std::max(tmpetsum25[0], tmpetsum25[1]),
                                  std::min(tmpetsum25[0], tmpetsum25[1]) };


      //////////////////////////////////////////////////////////
      // Now fill the histograms with the computed quantities

      // phi sizes of each trans/max/min region (for indexing by iregion)
      const double dphi[3] = { 2*PI/3.0, PI/3.0, PI/3.0 };

      // Loop over inclusive jet and exclusive dijet configurations
      for (size_t itopo = 0; itopo < 2; ++itopo) {

        // Exit early if in the exclusive dijet iteration and the exclusive dijet cuts are not met
        if (itopo == 1) {
          if (jets.size() != 2) continue;
          const FourMomentum psublead = jets[1].momentum();
          // Delta(phi) cut
          const double phisublead = psublead.phi();
          if (deltaPhi(philead, phisublead) < 2.5) continue;
          // pT fraction cut
          const double ptsublead  = psublead.pT();
          if (ptsublead < 0.5*ptlead) continue;
          MSG_DEBUG("Exclusive dijet event");
        }

        // Plot profiles and distributions which have no max/min region definition
        _p_etsum25_vs_ptlead_trans[itopo]->fill(ptlead/GeV, etsum25[0]/5.0/dphi[0]/GeV, weight);
        _p_etsum48_vs_ptlead_trans[itopo]->fill(ptlead/GeV, etsum48[0]/9.6/dphi[0]/GeV, weight);
        if (etsum25[0] > 0) {
          _p_chratio_vs_ptlead_trans[itopo]->fill(ptlead/GeV, ptsum[0]/etsum25[0], weight);
        }
        const double ptmean = safediv(ptsum[0], nch[0], -1); ///< Return -1 if div by zero
        if (ptmean >= 0) {
          _p_ptmeanch_vs_ptlead_trans[itopo]->fill(ptlead/GeV, ptmean/GeV, weight);
          _p_ptmeanch_vs_nch_trans[itopo]->fill(nch[0], ptmean/GeV, weight);
        }

        // Plot remaining profile and 1D observables, which are defined in all 3 tot/max/min regions
        for (size_t iregion = 0; iregion < 3; ++iregion) {
          _p_ptsumch_vs_ptlead[itopo][iregion]->fill(ptlead/GeV, ptsum[iregion]/5.0/dphi[iregion]/GeV, weight);
          _p_nch_vs_ptlead[itopo][iregion]->fill(ptlead/GeV, nch[iregion]/5.0/dphi[iregion], weight);
          for (size_t ipt = 0; ipt < 4; ++ipt) {
            if (ipt == 1 && !inRange(ptlead/GeV, 20, 60)) continue;
            if (ipt == 2 && !inRange(ptlead/GeV, 60, 210)) continue;
            if (ipt == 3 && ptlead/GeV < 210) continue;
            _h_ptsumch[ipt][itopo][iregion]->fill(ptsum[iregion]/5.0/dphi[iregion]/GeV, weight);
            _h_nch[ipt][itopo][iregion]->fill(nch[iregion]/5.0/dphi[iregion], weight);
          }
        }
      }

    }



    void finalize() {
      for (size_t iregion = 0; iregion < 3; ++iregion) {
        for (size_t itopo = 0; itopo < 2; ++itopo) {
          for (size_t ipt = 0; ipt < 4; ++ipt) {
            normalize(_h_ptsumch[ipt][itopo][iregion], 1.0);
            normalize(_h_nch[ipt][itopo][iregion], 1.0);
          }
        }
      }
    }


  private:

    /// @name Histogram arrays
    //@{

    Profile1DPtr _p_ptsumch_vs_ptlead[2][3];
    Profile1DPtr _p_nch_vs_ptlead[2][3];
    Profile1DPtr _p_ptmeanch_vs_ptlead_trans[2];
    Profile1DPtr _p_etsum25_vs_ptlead_trans[2];
    Profile1DPtr _p_etsum48_vs_ptlead_trans[2];
    Profile1DPtr _p_chratio_vs_ptlead_trans[2];
    Profile1DPtr _p_ptmeanch_vs_nch_trans[2];
    Histo1DPtr _h_ptsumch[4][2][3];
    Histo1DPtr _h_nch[4][2][3];

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1298811);

}
