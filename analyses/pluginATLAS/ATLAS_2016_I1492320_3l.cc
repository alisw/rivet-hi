// -*- C++ -*
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// @brief WWW cross-section at 8 TeV, 3L mode
  class ATLAS_2016_I1492320_3l : public Analysis {
  public:

    // Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1492320_3l);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Charged leptons within acceptance
      const PromptFinalState chLep_fid = PromptFinalState(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      const PromptFinalState photon_fs = PromptFinalState(Cuts::abspid == PID::PHOTON);
      const DressedLeptons dressed_leps(photon_fs, chLep_fid, 0.1, Cuts::pT > 20*GeV && Cuts::abseta < 2.5);
      declare(dressed_leps, "DressedLeptons");

      // Jets, anti-kt 0.4
      VetoedFinalState fsJets(FinalState(Cuts::abseta < 7.0)); //final state for jet finding: veto leptons and neutrinos
      fsJets.vetoNeutrinos();
      fsJets.addVetoOnThisFinalState(photon_fs);
      fsJets.addVetoOnThisFinalState(chLep_fid);
      declare(FastJets(fsJets, FastJets::ANTIKT, 0.4), "Jets");

      // Missing momentum
      declare(MissingMomentum(), "MET");

      // Histograms
      _h_fiducial_3l = bookCounter("d01-x01-y01");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the dressed leptons, sorted by pT of their constituent bare lepton (!!)
      vector<DressedLepton> _vbs_lep = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      if (_vbs_lep.size() != 3)  vetoEvent;
      std::sort(_vbs_lep.begin(), _vbs_lep.end(), [](const DressedLepton& l1, const DressedLepton& l2) {
        return (l1.constituentLepton().pT() > l2.constituentLepton().pT());
      });

      // Get the jets
      const Jets& _vbs_jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 4.5);
      if (_vbs_jets.size() > 1) vetoEvent;

      // Determine nsfos pairs for channel classification
      int nSFOS = 0;
      for (size_t i = 0; i < _vbs_lep.size(); ++i) {
        const double ch_l0 = _vbs_lep[i].charge();
        for (size_t j = i + 1; j < _vbs_lep.size(); ++j) {
          const double ch_l1 = _vbs_lep[j].charge();
          if (_vbs_lep[i].abspid() == _vbs_lep[j].abspid() && ch_l0*ch_l1 < 0) ++nSFOS;
        }
      }

      double minDRll = DBL_MAX, mSFOS_MinDiff = DBL_MAX, meeSS_MinDiff = DBL_MAX, mSF_min = DBL_MAX;
      for (size_t i = 0; i < _vbs_lep.size(); ++i) {
        const double ch_l0 = _vbs_lep[i].charge();
        for (size_t j = i + 1; j < _vbs_lep.size(); ++j) {
          const double ch_l1 = _vbs_lep[j].charge();
          const bool samesign = ch_l0*ch_l1 > 0;

          // Update min dR between leptons
          minDRll = min(minDRll, deltaR(_vbs_lep[i], _vbs_lep[j]));

          // Require same flavour
          if (_vbs_lep[i].abspid() != _vbs_lep[j].abspid()) continue;

          // SF dilepton mass (used several times)
          const double mSF = (_vbs_lep[i].momentum() + _vbs_lep[j].momentum()).mass();

          // Determine min for all same-flavor pairs
          mSF_min = min(mSF, mSF_min);

          // Determine min for all m_ee same-sign pairs
          if (_vbs_lep[i].abspid() == PID::ELECTRON && samesign) {
            if (fabs(mSF-ZMASS) < fabs(meeSS_MinDiff-ZMASS)) meeSS_MinDiff = mSF;
          }

          // Determine min for all mSFOS pairs
          if (!samesign && fabs(mSF-ZMASS) < abs(mSFOS_MinDiff-ZMASS)) mSFOS_MinDiff = mSF;
        }
      }

      if (minDRll < 0.1) vetoEvent;
      if (nSFOS == 0 && mSF_min < 20*GeV) vetoEvent;
      if (nSFOS == 0 && fabs(meeSS_MinDiff - ZMASS) < 15*GeV) vetoEvent;
      if (nSFOS == 1 && ((ZMASS - mSFOS_MinDiff) < 35*GeV && (mSFOS_MinDiff - ZMASS) < 20*GeV)) vetoEvent;
      if (nSFOS == 2 && fabs(mSFOS_MinDiff - ZMASS) < 20*GeV) vetoEvent;

      const Vector3& met = -1.0 * apply<MissingMomentum>(event, "MET").vectorEt();
      if (nSFOS == 1 && met.mod() < 45*GeV) vetoEvent;
      if (nSFOS == 2 && met.mod() < 55*GeV) vetoEvent;

      const double dPhi = deltaPhi((_vbs_lep[0].momentum() + _vbs_lep[1].momentum() + _vbs_lep[2].momentum()), met);
      if (dPhi < 2.5) vetoEvent;

      // Fill histo
      _h_fiducial_3l->fill(event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_fiducial_3l, crossSection()/sumOfWeights()/femtobarn);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    const double ZMASS = 91.1876*GeV;
    CounterPtr _h_fiducial_3l;

    //@}

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1492320_3l);

}
