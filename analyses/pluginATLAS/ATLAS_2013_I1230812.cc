// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"

namespace Rivet {


  /// Z + jets in pp at 7 TeV (combined channel / base class)
  /// @note This base class contains a "mode" variable for combined, e, and mu channel derived classes
  class ATLAS_2013_I1230812 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2013_I1230812(string name="ATLAS_2013_I1230812")
      : Analysis(name)
    {
      // This class uses the combined e+mu mode
      _mode = 1;
    }

    //@}


    /// Book histograms and initialise projections before the run
    void init() {
      // Determine the e/mu decay channels used (NB Prompt leptons only).
      /// @todo Note that Zs are accepted with any rapidity: the cuts are on the e/mu: is this correct?
      Cut pt20 = Cuts::pT >= 20.0*GeV;
      if (_mode == 1) {
        // Combined
        ZFinder zfinder(FinalState(Cuts::abseta < 2.5), pt20, PID::ELECTRON, 66*GeV, 116*GeV);
        declare(zfinder, "zfinder");
      } else if (_mode == 2) {
        // Electron
	    const Cut eta_e = Cuts::abseta < 1.37 || Cuts::absetaIn(1.52, 2.47);
        ZFinder zfinder(FinalState(eta_e), pt20, PID::ELECTRON, 66*GeV, 116*GeV);
        declare(zfinder, "zfinder");
      } else if (_mode == 3) {
        // Muon
        ZFinder zfinder(FinalState(Cuts::abseta < 2.4), pt20, PID::MUON, 66*GeV, 116*GeV);
        declare(zfinder, "zfinder");
      } else {
        MSG_ERROR("Unknown decay channel mode!!!");
      }

      // Define veto FS in order to prevent Z-decay products entering the jet algorithm
      VetoedFinalState had_fs;
      had_fs.addVetoOnThisFinalState(getProjection<ZFinder>("zfinder"));
      FastJets jets(had_fs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(true);
      declare(jets, "jets");

      _h_njet_incl              = bookHisto1D(  1, 1, _mode);
      _h_njet_incl_ratio        = bookScatter2D(2, 1, _mode, true);
      _h_njet_excl              = bookHisto1D(  3, 1, _mode);
      _h_njet_excl_ratio        = bookScatter2D(4, 1, _mode, true);
      _h_njet_excl_pt150        = bookHisto1D(  5, 1, _mode);
      _h_njet_excl_pt150_ratio  = bookScatter2D(6, 1, _mode, true);
      _h_njet_excl_vbf          = bookHisto1D ( 7, 1, _mode);
      _h_njet_excl_vbf_ratio    = bookScatter2D(8, 1, _mode, true);
      _h_ptlead                 = bookHisto1D(  9, 1, _mode);
      _h_ptseclead              = bookHisto1D( 10, 1, _mode);
      _h_ptthirdlead            = bookHisto1D( 11, 1, _mode);
      _h_ptfourthlead           = bookHisto1D( 12, 1, _mode);
      _h_ptlead_excl            = bookHisto1D( 13, 1, _mode);
      _h_pt_ratio               = bookHisto1D( 14, 1, _mode);
      _h_pt_z                   = bookHisto1D( 15, 1, _mode);
      _h_pt_z_excl              = bookHisto1D( 16, 1, _mode);
      _h_ylead                  = bookHisto1D( 17, 1, _mode);
      _h_yseclead               = bookHisto1D( 18, 1, _mode);
      _h_ythirdlead             = bookHisto1D( 19, 1, _mode);
      _h_yfourthlead            = bookHisto1D( 20, 1, _mode);
      _h_deltay                 = bookHisto1D( 21, 1, _mode);
      _h_mass                   = bookHisto1D( 22, 1, _mode);
      _h_deltaphi               = bookHisto1D( 23, 1, _mode);
      _h_deltaR                 = bookHisto1D( 24, 1, _mode);
      _h_ptthirdlead_vbf        = bookHisto1D( 25, 1, _mode);
      _h_ythirdlead_vbf         = bookHisto1D( 26, 1, _mode);
      _h_ht                     = bookHisto1D( 27, 1, _mode);
      _h_st                     = bookHisto1D( 28, 1, _mode);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const ZFinder& zfinder = apply<ZFinder>(event, "zfinder");
      MSG_DEBUG("# Z constituents = " << zfinder.constituents().size());
      if (zfinder.constituents().size() != 2) vetoEvent;

      FourMomentum z = zfinder.boson().momentum();
      FourMomentum lp = zfinder.constituents()[0].momentum();
      FourMomentum lm = zfinder.constituents()[1].momentum();
      if (deltaR(lp, lm) < 0.2) vetoEvent;

      Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::absrap < 4.4);
      ifilter_discard(jets, deltaRLess(lp, 0.5));
      ifilter_discard(jets, deltaRLess(lm, 0.5));

      const double weight = event.weight();

      // Fill jet multiplicities
      for (size_t ijet = 0; ijet <= jets.size(); ++ijet) {
        _h_njet_incl->fill(ijet, weight);
      }
      _h_njet_excl->fill(jets.size(), weight);

      // Require at least one jet
      if (jets.size() >= 1) {
        // Leading jet histos
        const double ptlead   = jets[0].pT()/GeV;
        const double yabslead = fabs(jets[0].rapidity());
        const double ptz   = z.pT()/GeV;
        _h_ptlead->fill(ptlead,   weight);
        _h_ylead ->fill(yabslead, weight);
        _h_pt_z  ->fill(ptz, weight);
        // Fill jet multiplicities
        if (ptlead > 150)  _h_njet_excl_pt150->fill(jets.size(), weight);

        // Loop over selected jets, fill inclusive distributions
        double st = 0;
        double ht = lp.pT()/GeV + lm.pT()/GeV;
        for (size_t ijet = 0; ijet < jets.size(); ++ijet) {
          ht += jets[ijet].pT()/GeV;
          st += jets[ijet].pT()/GeV;
        }
        _h_ht->fill(ht, weight);
        _h_st->fill(st, weight);

        // Require exactly one jet
        if (jets.size() == 1) {
          _h_ptlead_excl->fill(ptlead,   weight);
          _h_pt_z_excl  ->fill(ptz, weight);
        }
      }


      // Require at least two jets
      if (jets.size() >= 2) {
        // Second jet histos
        const double ptlead      = jets[0].pT()/GeV;
        const double pt2ndlead   = jets[1].pT()/GeV;
        const double ptratio     = pt2ndlead/ptlead;
        const double yabs2ndlead = fabs(jets[1].rapidity());
        _h_ptseclead->fill(pt2ndlead,   weight);
        _h_yseclead->fill( yabs2ndlead, weight);
        _h_pt_ratio->fill( ptratio,     weight);

        // Dijet histos
        const double deltaphi = fabs(deltaPhi(jets[1], jets[0]));
        const double deltarap = fabs(jets[0].rapidity() - jets[1].rapidity()) ;
        const double deltar   = fabs(deltaR(jets[0], jets[1], RAPIDITY));
        const double mass     = (jets[0].momentum() + jets[1].momentum()).mass()/GeV;
        _h_mass->fill(    mass,     weight);
        _h_deltay->fill(  deltarap, weight);
        _h_deltaphi->fill(deltaphi, weight);
        _h_deltaR->fill(  deltar,   weight);

        if (mass > 350 && deltarap > 3)  _h_njet_excl_vbf->fill(jets.size(), weight);
      }

      // Require at least three jets
      if (jets.size() >= 3) {
        // Third jet histos
        const double pt3rdlead   = jets[2].pT()/GeV;
        const double yabs3rdlead = fabs(jets[2].rapidity());
        _h_ptthirdlead->fill(pt3rdlead,   weight);
        _h_ythirdlead->fill( yabs3rdlead, weight);

        //Histos after VBF preselection
        const double deltarap = fabs(jets[0].rapidity() - jets[1].rapidity()) ;
        const double mass     = (jets[0].momentum() + jets[1].momentum()).mass();
        if (mass > 350 && deltarap > 3) {
          _h_ptthirdlead_vbf->fill(pt3rdlead,   weight);
          _h_ythirdlead_vbf->fill( yabs3rdlead, weight);
        }
      }

      // Require at least four jets
      if (jets.size() >= 4) {
        // Fourth jet histos
        const double pt4thlead   = jets[3].pT()/GeV;
        const double yabs4thlead = fabs(jets[3].rapidity());
        _h_ptfourthlead->fill(pt4thlead,   weight);
        _h_yfourthlead->fill( yabs4thlead, weight);
      }
    }

    /// @name Ratio calculator util functions
    //@{

    /// Calculate the efficiency error, being careful about div-by-zero
    double err_incl(const HistoBin1D &M, const HistoBin1D &N, bool hasWeights) {
      double r = safediv(M.sumW(), N.sumW());
      if (hasWeights) { // use F. James's approximation for weighted events
        return sqrt( safediv((1 - 2 * r) * M.sumW2() + r * r * N.sumW2(), N.sumW() * N.sumW()) );
      }
      return sqrt( safediv(r * (1 - r), N.sumW()) );
    }

    /// Calculate the ratio error, being careful about div-by-zero
    double err_excl(const HistoBin1D &A, const HistoBin1D &B) {
      double r = safediv(A.sumW(), B.sumW());
      double dAsquared = safediv(A.sumW2(), A.sumW() * A.sumW()); // squared relative error of A
      double dBsquared = safediv(B.sumW2(), B.sumW() * B.sumW()); // squared relative error of B
      return r * sqrt(dAsquared + dBsquared);
    }

    //@}


    void finalize() {
      bool hasWeights = _h_njet_incl->effNumEntries() != _h_njet_incl->numEntries();
      for (size_t i = 0; i < 6; ++i) {
        _h_njet_incl_ratio->point(i).setY(safediv(_h_njet_incl->bin(i + 1).sumW(), _h_njet_incl->bin(i).sumW()),
                                          err_incl(_h_njet_incl->bin(i + 1), _h_njet_incl->bin(i), hasWeights));
        _h_njet_excl_ratio->point(i).setY(safediv(_h_njet_excl->bin(i + 1).sumW(), _h_njet_excl->bin(i).sumW()),
                                          err_excl(_h_njet_excl->bin(i + 1), _h_njet_excl->bin(i)));
        if (i >= 1) {
          _h_njet_excl_pt150_ratio->point(i - 1).setY(safediv(_h_njet_excl_pt150->bin(i).sumW(), _h_njet_excl_pt150->bin(i - 1).sumW()),
                                                      err_excl(_h_njet_excl_pt150->bin(i), _h_njet_excl_pt150->bin(i - 1)));
          if (i >= 2) {
            _h_njet_excl_vbf_ratio->point(i - 2).setY(safediv(_h_njet_excl_vbf->bin(i).sumW(), _h_njet_excl_vbf->bin(i - 1).sumW()),
                                                      err_excl(_h_njet_excl_vbf->bin(i), _h_njet_excl_vbf->bin(i - 1)));
          }
        }
      }

      const double xs = crossSectionPerEvent()/picobarn;
      scale({_h_njet_incl, _h_njet_excl, _h_njet_excl_pt150, _h_njet_excl_vbf}, xs);
      scale({_h_ptlead, _h_ptseclead, _h_ptthirdlead, _h_ptfourthlead, _h_ptlead_excl}, xs);
      scale({_h_pt_ratio, _h_pt_z, _h_pt_z_excl}, xs);
      scale({_h_ylead, _h_yseclead, _h_ythirdlead, _h_yfourthlead}, xs);
      scale({_h_deltay, _h_mass, _h_deltaphi, _h_deltaR}, xs);
      scale({_h_ptthirdlead_vbf, _h_ythirdlead_vbf}, xs);
      scale({_h_ht, _h_st}, xs);
    }

    //@}


  protected:

    size_t _mode;


  private:

    Scatter2DPtr _h_njet_incl_ratio;
    Scatter2DPtr _h_njet_excl_ratio;
    Scatter2DPtr _h_njet_excl_pt150_ratio;
    Scatter2DPtr _h_njet_excl_vbf_ratio;
    Histo1DPtr _h_njet_incl;
    Histo1DPtr _h_njet_excl;
    Histo1DPtr _h_njet_excl_pt150;
    Histo1DPtr _h_njet_excl_vbf;
    Histo1DPtr _h_ptlead;
    Histo1DPtr _h_ptseclead;
    Histo1DPtr _h_ptthirdlead;
    Histo1DPtr _h_ptfourthlead;
    Histo1DPtr _h_ptlead_excl;
    Histo1DPtr _h_pt_ratio;
    Histo1DPtr _h_pt_z;
    Histo1DPtr _h_pt_z_excl;
    Histo1DPtr _h_ylead;
    Histo1DPtr _h_yseclead;
    Histo1DPtr _h_ythirdlead;
    Histo1DPtr _h_yfourthlead;
    Histo1DPtr _h_deltay;
    Histo1DPtr _h_mass;
    Histo1DPtr _h_deltaphi;
    Histo1DPtr _h_deltaR;
    Histo1DPtr _h_ptthirdlead_vbf;
    Histo1DPtr _h_ythirdlead_vbf;
    Histo1DPtr _h_ht;
    Histo1DPtr _h_st;
  };



  class ATLAS_2013_I1230812_EL : public ATLAS_2013_I1230812 {
  public:
    ATLAS_2013_I1230812_EL()
      : ATLAS_2013_I1230812("ATLAS_2013_I1230812_EL")
    {
      _mode = 2;
    }
  };



  class ATLAS_2013_I1230812_MU : public ATLAS_2013_I1230812 {
  public:
    ATLAS_2013_I1230812_MU()
      : ATLAS_2013_I1230812("ATLAS_2013_I1230812_MU")
    {
      _mode = 3;
    }
  };



  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1230812);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1230812_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1230812_MU);
}
