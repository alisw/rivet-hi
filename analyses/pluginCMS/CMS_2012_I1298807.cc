// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MergedFinalState.hh"

namespace Rivet {


  /// Inclusive ZZ production cross section and constraints on anomalous triple gauge couplings
  class CMS_2012_I1298807 : public Analysis {
  public:

    // Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2012_I1298807);


    /// Initialise projections and histograms
    void init() {

      // FinalState electrons(Cuts::abseta < 2.5 && Cuts::abspid == PID::ELECTRON);
      // FinalState muons(Cuts::abseta < 2.4 && Cuts::abspid == PID::MUON);
      // MergedFinalState leptons(electrons, muons);
      FinalState leptons((Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5) ||
                         (Cuts::abspid == PID::MUON && Cuts::abseta < 2.4));
      declare(leptons, "Leptons");

      Cut cut_el = Cuts::abseta < 2.5 && Cuts::pT > 7.0*GeV;
      Cut cut_mu = Cuts::abseta < 2.4 && Cuts::pT > 5.0*GeV;

      ZFinder zeefinder(FinalState(), cut_el, PID::ELECTRON, 60*GeV, 120*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
      declare(zeefinder, "ZeeFinder");
      ZFinder zmmfinder(FinalState(), cut_mu, PID::MUON, 60*GeV, 120*GeV, 0.1, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
      declare(zmmfinder, "ZmmFinder");

      VetoedFinalState fs_woZmm;
      fs_woZmm.addVetoOnThisFinalState(zmmfinder);
      VetoedFinalState fs_woZee;
      fs_woZee.addVetoOnThisFinalState(zeefinder);

      ZFinder zeefinder_woZee(fs_woZee, cut_el, PID::ELECTRON, 60*GeV, 120*GeV, 0.1, ZFinder::CLUSTERNODECAY);
      declare(zeefinder_woZee, "Zeefinder_WoZee");
      ZFinder zmmfinder_woZmm(fs_woZmm, cut_mu, PID::MUON, 60*GeV, 120*GeV, 0.1, ZFinder::CLUSTERNODECAY);
      declare(zmmfinder_woZmm, "Zmmfinder_WoZmm");

      // Book histograms
      _hist_pt_l1 = bookHisto1D(1, 1, 1);
      _hist_pt_z1 = bookHisto1D(1, 1, 2);
      _hist_pt_zz = bookHisto1D(1, 1, 3);
      _hist_m_zz = bookHisto1D(1, 1, 4);
      _hist_dphi_zz = bookHisto1D(1, 1, 5);
      _hist_dR_zz = bookHisto1D(1, 1, 6);

    }


    // Perform the per-event analysis
    void analyze(const Event& evt) {

      // Find leading leptons and apply cuts
      const Particles& leptons = apply<FinalState>(evt, "Leptons").particlesByPt();
      if (leptons.size() < 2) vetoEvent;
      const double leading_l_pt = leptons[0].pT();
      const double second_l_pt = leptons[1].pT();
      if (leading_l_pt < 20*GeV || second_l_pt < 10*GeV) vetoEvent;

      // Find acceptable ZZ combinations and build four-momenta, otherwise veto
      const ZFinder& zeefinder = applyProjection<ZFinder>(evt, "ZeeFinder");
      const ZFinder& zeefinder_woZee = applyProjection<ZFinder>(evt, "Zeefinder_WoZee");
      const ZFinder& zmmfinder = applyProjection<ZFinder>(evt, "ZmmFinder");
      const ZFinder& zmmfinder_woZmm = applyProjection<ZFinder>(evt, "Zmmfinder_WoZmm");

      FourMomentum pZ_a, pZ_b, pZ_1, pZ_2;
      FourMomentum pZZ, Z_a_l1, Z_a_l2, Z_b_l1, Z_b_l2;
      if (zeefinder.bosons().size() > 0 && zmmfinder.bosons().size() > 0) {
        pZ_a = zeefinder.bosons()[0];
        pZ_b = zmmfinder.bosons()[0];
        pZZ = pZ_a + pZ_b;
        pZ_1 = pZ_a;
        pZ_2 = pZ_b;
        Z_a_l1 = zeefinder.constituents()[0];
        Z_a_l2 = zeefinder.constituents()[1];
        Z_b_l1 = zmmfinder.constituents()[0];
        Z_b_l2 = zmmfinder.constituents()[1];
      } else if (zeefinder.bosons().size() > 0 && zeefinder_woZee.bosons().size() > 0) {
        pZ_a = zeefinder.bosons()[0];
        pZ_b = zeefinder_woZee.bosons()[0];
        pZZ = pZ_a + pZ_b;
        pZ_1 = pZ_a;
        pZ_2 = pZ_b;
        Z_a_l1 = zeefinder.constituents()[0];
        Z_a_l2 = zeefinder.constituents()[1];
        Z_b_l1 = zeefinder_woZee.constituents()[0];
        Z_b_l2 = zeefinder_woZee.constituents()[1];
      } else if (zmmfinder.bosons().size() > 0 && zmmfinder_woZmm.bosons().size() > 0) {
        pZ_a = zmmfinder.bosons()[0];
        pZ_b = zmmfinder_woZmm.bosons()[0];
        pZZ = pZ_a + pZ_b;
        pZ_1 = pZ_a;
        pZ_2 = pZ_b;
        Z_a_l1 = zmmfinder.constituents()[0];
        Z_a_l2 = zmmfinder.constituents()[1];
        Z_b_l1 = zmmfinder_woZmm.constituents()[0];
        Z_b_l2 = zmmfinder_woZmm.constituents()[1];
      } else {
        vetoEvent;
      }

      // Set ordered pT variables
      /// @todo Looks like there should be a nicer way than this
      double pt_l1 = Z_a_l1.pT();
      if (Z_a_l2.pT() > pt_l1) pt_l1 = Z_a_l2.pT();
      if (Z_b_l1.pT() > pt_l1) pt_l1 = Z_b_l1.pT();
      if (Z_b_l2.pT() > pt_l1) pt_l1 = Z_b_l2.pT();

      // Leading Z pT
      double pt_z1 = pZ_a.pT();
      if (pZ_b.pT() > pZ_a.pT()) {
        pt_z1 = pZ_b.pT();
        pZ_1 = pZ_b;
        pZ_2 = pZ_a;
      }

      // Fill histograms
      const double weight = evt.weight();
      _hist_pt_zz->fill(pZZ.pT()/GeV, weight);
      _hist_m_zz->fill(pZZ.mass()/GeV, weight);
      _hist_dphi_zz->fill(deltaPhi(pZ_a, pZ_b), weight);
      _hist_dR_zz->fill(deltaR(pZ_a, pZ_b, PSEUDORAPIDITY), weight);
      _hist_pt_z1->fill(pt_z1/GeV, weight);
      _hist_pt_l1->fill(pt_l1/GeV, weight);

    }


    /// Scale histograms
    /// @note This is all needed to undo bin width factor -- WHY DO PEOPLE USE UNPHYSICAL HISTOGRAMS?!?
    /// @todo If we introduce a "bar plot" or similar, it'd work better here
    void finalize() {

      double sum_height_pt_zz = 0;
      for (size_t i = 0; i < _hist_pt_zz->numBins(); i++) {
        _hist_pt_zz->bin(i).scaleW(1. / _hist_pt_zz->bin(i).width());
        sum_height_pt_zz += _hist_pt_zz->bin(i).height();
      }
      scale(_hist_pt_zz, 1. / sum_height_pt_zz);

      double sum_height_m_zz = 0;
      for (size_t i = 0; i < _hist_m_zz->numBins(); i++) {
        _hist_m_zz->bin(i).scaleW(1. / _hist_m_zz->bin(i).width());
        sum_height_m_zz += _hist_m_zz->bin(i).height();
      }
      scale(_hist_m_zz, 1. / sum_height_m_zz);

      double sum_height_dphi_zz = 0;
      for (size_t i = 0; i < _hist_dphi_zz->numBins(); i++) {
        _hist_dphi_zz->bin(i).scaleW(1. / _hist_dphi_zz->bin(i).width());
        sum_height_dphi_zz += _hist_dphi_zz->bin(i).height();
      }
      scale(_hist_dphi_zz, 1. / sum_height_dphi_zz);

      double sum_height_dR_zz = 0;
      for (size_t i = 0; i < _hist_dR_zz->numBins(); i++) {
        _hist_dR_zz->bin(i).scaleW(1. / _hist_dR_zz->bin(i).width());
        sum_height_dR_zz += _hist_dR_zz->bin(i).height();
      }
      scale(_hist_dR_zz, 1. / sum_height_dR_zz);

      double sum_height_pt_z1 = 0;
      for (size_t i = 0; i < _hist_pt_z1->numBins(); i++) {
        _hist_pt_z1->bin(i).scaleW(1. / _hist_pt_z1->bin(i).width());
        sum_height_pt_z1 += _hist_pt_z1->bin(i).height();
      }
      scale(_hist_pt_z1, 1. / sum_height_pt_z1);

      double sum_height_pt_l1 = 0;
      for (size_t i = 0; i < _hist_pt_l1->numBins(); i++) {
        _hist_pt_l1->bin(i).scaleW(1. / _hist_pt_l1->bin(i).width());
        sum_height_pt_l1 += _hist_pt_l1->bin(i).height();
      }
      scale(_hist_pt_l1, 1. / sum_height_pt_l1);
    }


    /// Histograms
    Histo1DPtr _hist_pt_zz, _hist_m_zz, _hist_dphi_zz, _hist_dR_zz, _hist_pt_z1, _hist_pt_l1;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2012_I1298807);

}
