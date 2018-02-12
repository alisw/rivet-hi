// -*- C++ -*-
// ATLAS W+c analysis

//////////////////////////////////////////////////////////////////////////
/*
  Description of rivet analysis ATLAS_2014_I1282447 W+c production

  This rivet routine implements the ATLAS W+c analysis.
  Apart from those histograms, described and published on HEP Data, here
  are some helper histograms defined, these are:

  d02-x01-y01, d02-x01-y02 and d08-x01-y01 are ratios, the nominator ("_plus")
  and denominator ("_minus") histograms are also given, so that the ratios can
  be reconstructed if need be (e.g. when running on separate samples).

  d05 and d06 are ratios over inclusive W production.
  The routine has to be run on a sample for inclusive W production in order to
  make sure the denominator ("_winc") is correctly filled.

  The ratios can be constructed using the following sample code:
  python divideWCharm.py

  import yoda
  hists_wc   = yoda.read("Rivet_Wc.yoda")
  hists_winc = yoda.read("Rivet_Winc.yoda")

  ## division histograms --> ONLY for different plus minus runs
  # (merge before using yodamerge Rivet_plus.yoda Rivet_minus.yoda > Rivet_Wc.yoda)

  d02y01_plus = hists_wc["/ATLAS_2014_I1282447/d02-x01-y01_plus"]
  d02y01_minus = hists_wc["/ATLAS_2014_I1282447/d02-x01-y01_minus"]
  ratio_d02y01 =  d02y01_plus.divide(d02y01_minus)
  ratio_d02y01.path = "/ATLAS_2014_I1282447/d02-x01-y01"

  d02y02_plus = hists_wc["/ATLAS_2014_I1282447/d02-x01-y02_plus"]
  d02y02_minus = hists_wc["/ATLAS_2014_I1282447/d02-x01-y02_minus"]
  ratio_d02y02=  d02y02_plus.divide(d02y02_minus)
  ratio_d02y02.path = "/ATLAS_2014_I1282447/d02-x01-y02"

  d08y01_plus = hists_wc["/ATLAS_2014_I1282447/d08-x01-y01_plus"]
  d08y01_minus = hists_wc["/ATLAS_2014_I1282447/d08-x01-y01_minus"]
  ratio_d08y01=  d08y01_plus.divide(d08y01_minus)
  ratio_d08y01.path = "/ATLAS_2014_I1282447/d08-x01-y01"

  # inclusive cross section
  h_winc = hists_winc["/ATLAS_2014_I1282447/d05-x01-y01"]
  h_d    = hists_wc["/ATLAS_2014_I1282447/d01-x01-y02"]
  h_dstar= hists_wc["/ATLAS_2014_I1282447/d01-x01-y03"]

  ratio_wd      =  h_d.divide(h_winc)
  ratio_wd.path = "/ATLAS_2014_I1282447/d05-x01-y02"

  ratio_wdstar      =  h_d.divide(h_winc)
  ratio_wdstar.path = "/ATLAS_2014_I1282447/d05-x01-y03"

  # pT differential
  h_winc_plus  = hists_winc["/ATLAS_2014_I1282447/d06-x01-y01_winc"]
  h_winc_minus = hists_winc["/ATLAS_2014_I1282447/d06-x01-y02_winc"]

  h_wd_plus      = hists_wc["/ATLAS_2014_I1282447/d06-x01-y01_wplus"]
  h_wd_minus     = hists_wc["/ATLAS_2014_I1282447/d06-x01-y02_wminus"]
  h_wdstar_plus  = hists_wc["/ATLAS_2014_I1282447/d06-x01-y03_wplus"]
  h_wdstar_minus = hists_wc["/ATLAS_2014_I1282447/d06-x01-y04_wminus"]

  ratio_wd_plus       =  h_wd_plus.divide(h_winc_plus)
  ratio_wd_plus.path  = "/ATLAS_2014_I1282447/d06-x01-y01"
  ratio_wd_minus      =  h_wd_plus.divide(h_winc_minus)
  ratio_wd_minus.path = "/ATLAS_2014_I1282447/d06-x01-y02"

  ratio_wdstar_plus       =  h_wdstar_plus.divide(h_winc_plus)
  ratio_wdstar_plus.path  = "/ATLAS_2014_I1282447/d06-x01-y03"
  ratio_wdstar_minus      =  h_wdstar_plus.divide(h_winc_minus)
  ratio_wdstar_minus.path = "/ATLAS_2014_I1282447/d06-x01-y04"

  ratio_wd_plus =  h_wd_plus.divide(h_winc_plus)
  ratio_wd_plus.path = "/ATLAS_2014_I1282447/d06-x01-y01"
  ratio_wd_minus =  h_wd_plus.divide(h_winc_minus)
  ratio_wd_minus.path = "/ATLAS_2014_I1282447/d06-x01-y02"

  h_winc_plus= hists_winc["/ATLAS_2014_I1282447/d06-x01-y01_winc"]
  h_winc_minus= hists_winc["/ATLAS_2014_I1282447/d06-x01-y02_winc"]

  ## copy other histograms for plotting

  d01x01y01= hists_wc["/ATLAS_2014_I1282447/d01-x01-y01"]
  d01x01y01.path = "/ATLAS_2014_I1282447/d01-x01-y01"

  d01x01y02= hists_wc["/ATLAS_2014_I1282447/d01-x01-y02"]
  d01x01y02.path = "/ATLAS_2014_I1282447/d01-x01-y02"

  d01x01y03= hists_wc["/ATLAS_2014_I1282447/d01-x01-y03"]
  d01x01y03.path = "/ATLAS_2014_I1282447/d01-x01-y03"

  d03x01y01= hists_wc["/ATLAS_2014_I1282447/d03-x01-y01"]
  d03x01y01.path = "/ATLAS_2014_I1282447/d03-x01-y01"

  d03x01y02= hists_wc["/ATLAS_2014_I1282447/d03-x01-y02"]
  d03x01y02.path = "/ATLAS_2014_I1282447/d03-x01-y02"

  d04x01y01= hists_wc["/ATLAS_2014_I1282447/d04-x01-y01"]
  d04x01y01.path = "/ATLAS_2014_I1282447/d04-x01-y01"

  d04x01y02= hists_wc["/ATLAS_2014_I1282447/d04-x01-y02"]
  d04x01y02.path = "/ATLAS_2014_I1282447/d04-x01-y02"

  d04x01y03= hists_wc["/ATLAS_2014_I1282447/d04-x01-y03"]
  d04x01y03.path = "/ATLAS_2014_I1282447/d04-x01-y03"

  d04x01y04= hists_wc["/ATLAS_2014_I1282447/d04-x01-y04"]
  d04x01y04.path = "/ATLAS_2014_I1282447/d04-x01-y04"

  d07x01y01= hists_wc["/ATLAS_2014_I1282447/d07-x01-y01"]
  d07x01y01.path = "/ATLAS_2014_I1282447/d07-x01-y01"

  yoda.write([ratio_d02y01,ratio_d02y02,ratio_d08y01, ratio_wd ,ratio_wdstar,ratio_wd_plus,ratio_wd_minus ,ratio_wdstar_plus,ratio_wdstar_minus,d01x01y01,d01x01y02,d01x01y03,d03x01y01,d03x01y02,d04x01y01,d04x01y02,d04x01y03,d04x01y04,d07x01y01],"validation.yoda")

*/
//////////////////////////////////////////////////////////////////////////

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {



  class ATLAS_2014_I1282447 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1282447() : Analysis("ATLAS_2014_I1282447")
    {
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      /// @todo Initialise and register projections here
      UnstableFinalState fs;

      Cut cuts = Cuts::etaIn(-2.5, 2.5) & (Cuts::pT > 20*GeV);

      /// should use sample WITHOUT QED radiation off the electron
      WFinder wfinder_born_el(fs, cuts, PID::ELECTRON, 25*GeV, 8000*GeV, 15*GeV, 0.1, WFinder::CLUSTERALL, WFinder::TRACK);
      declare(wfinder_born_el, "WFinder_born_el");

      WFinder wfinder_born_mu(fs, cuts, PID::MUON    , 25*GeV, 8000*GeV, 15*GeV, 0.1, WFinder::CLUSTERALL, WFinder::TRACK);
      declare(wfinder_born_mu, "WFinder_born_mu");

      // all hadrons that could be coming from a charm decay --
      // -- for safety, use region -3.5 - 3.5
      declare(UnstableFinalState(Cuts::abseta <3.5), "hadrons");

      // Input for the jets: no neutrinos, no muons, and no electron which passed the electron cuts
      // also: NO electron, muon or tau (needed due to ATLAS jet truth reconstruction feature)
      VetoedFinalState veto;

      veto.addVetoOnThisFinalState(wfinder_born_el);
      veto.addVetoOnThisFinalState(wfinder_born_mu);
      veto.addVetoPairId(PID::ELECTRON);
      veto.addVetoPairId(PID::MUON);
      veto.addVetoPairId(PID::TAU);

      FastJets jets(veto, FastJets::ANTIKT, 0.4);
      declare(jets, "jets");

      // Book histograms

      // charge separated integrated cross sections
      _hist_wcjet_charge  = bookHisto1D("d01-x01-y01");
      _hist_wd_charge     = bookHisto1D("d01-x01-y02");
      _hist_wdstar_charge = bookHisto1D("d01-x01-y03");

      // charge integrated total cross sections
      _hist_wcjet_ratio = bookScatter2D("d02-x01-y01");
      _hist_wd_ratio    = bookScatter2D("d02-x01-y02");

      _hist_wcjet_minus = bookHisto1D("d02-x01-y01_minus");
      _hist_wd_minus    = bookHisto1D("d02-x01-y02_minus");

      _hist_wcjet_plus  = bookHisto1D("d02-x01-y01_plus");
      _hist_wd_plus     = bookHisto1D("d02-x01-y02_plus");

      // eta distributions
      _hist_wplus_wcjet_eta_lep   = bookHisto1D("d03-x01-y01");
      _hist_wminus_wcjet_eta_lep  = bookHisto1D("d03-x01-y02");

      _hist_wplus_wdminus_eta_lep = bookHisto1D("d04-x01-y01");
      _hist_wminus_wdplus_eta_lep = bookHisto1D("d04-x01-y02");
      _hist_wplus_wdstar_eta_lep  = bookHisto1D("d04-x01-y03");
      _hist_wminus_wdstar_eta_lep = bookHisto1D("d04-x01-y04");

      // ratio of cross section (WD over W inclusive) // postprocess!
      _hist_w_inc             = bookHisto1D("d05-x01-y01");
      _hist_wd_winc_ratio     = bookScatter2D("d05-x01-y02");
      _hist_wdstar_winc_ratio = bookScatter2D("d05-x01-y03");

      // ratio of cross section (WD over W inclusive -- function of pT of D meson)
      _hist_wplusd_wplusinc_pt_ratio       = bookScatter2D("d06-x01-y01");
      _hist_wminusd_wminusinc_pt_ratio     = bookScatter2D("d06-x01-y02");
      _hist_wplusdstar_wplusinc_pt_ratio   = bookScatter2D("d06-x01-y03");
      _hist_wminusdstar_wminusinc_pt_ratio = bookScatter2D("d06-x01-y04");

      // could use for postprocessing!
      _hist_wplusd_wplusinc_pt       = bookHisto1D("d06-x01-y01_wplus");
      _hist_wminusd_wminusinc_pt     = bookHisto1D("d06-x01-y02_wminus");
      _hist_wplusdstar_wplusinc_pt   = bookHisto1D("d06-x01-y03_wplus");
      _hist_wminusdstar_wminusinc_pt = bookHisto1D("d06-x01-y04_wminus");

      _hist_wplus_winc  = bookHisto1D("d06-x01-y01_winc");
      _hist_wminus_winc = bookHisto1D("d06-x01-y02_winc");

      // jet multiplicity of charge integrated W+cjet cross section (+0 or +1 jet in addition to the charm jet)
      _hist_wcjet_jets  = bookHisto1D("d07-x01-y01");

      // jet multiplicity of W+cjet cross section ratio (+0 or +1 jet in addition to the charm jet)
      _hist_wcjet_jets_ratio  = bookScatter2D("d08-x01-y01");
      _hist_wcjet_jets_plus   = bookHisto1D("d08-x01-y01_plus");
      _hist_wcjet_jets_minus  = bookHisto1D("d08-x01-y01_minus");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      double charge_weight = 0; // account for OS/SS events

      int    lepton_charge = 0;
      double lepton_eta    = 0.;

      /// Find leptons
      const WFinder& wfinder_born_el = apply<WFinder>(event, "WFinder_born_el");
      const WFinder& wfinder_born_mu = apply<WFinder>(event, "WFinder_born_mu");

      if (wfinder_born_el.empty() && wfinder_born_mu.empty()) {
        MSG_DEBUG("No W bosons found");
        vetoEvent;
      }

      bool keepevent = false;

      //check electrons
      if (!wfinder_born_el.empty()) {
        const FourMomentum nu = wfinder_born_el.constituentNeutrinos()[0];
        if (wfinder_born_el.mT() > 40*GeV && nu.pT() > 25*GeV) {
          keepevent = true;
          lepton_charge = wfinder_born_el.constituentLeptons()[0].charge();
          lepton_eta = fabs(wfinder_born_el.constituentLeptons()[0].pseudorapidity());
        }
      }

      //check muons
      if (!wfinder_born_mu.empty()) {
        const FourMomentum nu = wfinder_born_mu.constituentNeutrinos()[0];
        if (wfinder_born_mu.mT() > 40*GeV && nu.pT() > 25*GeV) {
          keepevent = true;
          lepton_charge = wfinder_born_mu.constituentLeptons()[0].charge();
          lepton_eta = fabs(wfinder_born_mu.constituentLeptons()[0].pseudorapidity());
        }
      }

      if (!keepevent) {
        MSG_DEBUG("Event does not pass mT and MET cuts");
        vetoEvent;
      }

      if (lepton_charge > 0) {
        _hist_wplus_winc->fill(10., weight);
        _hist_wplus_winc->fill(16., weight);
        _hist_wplus_winc->fill(30., weight);
        _hist_wplus_winc->fill(60., weight);
        _hist_w_inc->fill(+1, weight);
      }
      else if (lepton_charge < 0) {
        _hist_wminus_winc->fill(10., weight);
        _hist_wminus_winc->fill(16., weight);
        _hist_wminus_winc->fill(30., weight);
        _hist_wminus_winc->fill(60., weight);
        _hist_w_inc->fill(-1, weight);
      }

      // Find hadrons in the event
      const UnstableFinalState& fs = apply<UnstableFinalState>(event, "hadrons");

      /// FIND Different channels
      // 1: wcjet
      // get jets
      const Jets& jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT>25.0*GeV && Cuts::abseta<2.5);
      // loop over jets to select jets used to match to charm
      Jets js;
      int    matched_charmHadron = 0;
      double charm_charge = 0.;
      int    njets = 0;
      int    nj = 0;
      bool   mat_jet = false;

      double ptcharm = 0;
      if (matched_charmHadron > -1) {
        for (const Jet& j : jets) {
          mat_jet = false;
          njets += 1;
          for (const Particle& p : fs.particles()) {
            /// @todo Avoid touching HepMC!
            const GenParticle* part = p.genParticle();
            if (p.hasCharm()) {
              //if(isFromBDecay(p)) continue;
              if (p.fromBottom()) continue;
              if (p.pT() < 5*GeV ) continue;
              if (hasCharmedChildren(part)) continue;
              if (deltaR(p, j) < 0.3) {
                mat_jet = true;
                if (p.pT() > ptcharm) {
                  charm_charge = part->pdg_id();
                  ptcharm = p.pT();
                }
              }
            }
          }
          if (mat_jet) nj++;
        }

        if (charm_charge * lepton_charge > 0)  charge_weight = -1;
        else charge_weight = +1;

        if (nj == 1)  {
          if (lepton_charge > 0) {
            _hist_wcjet_charge        ->fill(         1, weight*charge_weight);
            _hist_wcjet_plus          ->fill(         0, weight*charge_weight);
            _hist_wplus_wcjet_eta_lep ->fill(lepton_eta, weight*charge_weight);
            _hist_wcjet_jets_plus     ->fill(njets-1   , weight*charge_weight);
          }
          else if (lepton_charge < 0) {
            _hist_wcjet_charge        ->fill(        -1, weight*charge_weight);
            _hist_wcjet_minus         ->fill(         0, weight*charge_weight);
            _hist_wminus_wcjet_eta_lep->fill(lepton_eta, weight*charge_weight);
            _hist_wcjet_jets_minus    ->fill(njets-1   , weight*charge_weight);
          }

          _hist_wcjet_jets->fill(njets-1, weight*charge_weight);
        }
      }

      // // 1/2: w+d(*) meson

      for (const Particle& p : fs.particles()) {

        /// @todo Avoid touching HepMC!
        const GenParticle* part = p.genParticle();
        if (p.pT() < 8*GeV)       continue;
        if (fabs(p.eta()) > 2.2)  continue;

        // W+D
        if (abs(part->pdg_id()) == 411) {
          if (lepton_charge * part->pdg_id() > 0)  charge_weight = -1;
          else                                     charge_weight = +1;

          // fill histos
          if (lepton_charge > 0) {
            _hist_wd_charge            ->fill(         1, weight*charge_weight);
            _hist_wd_plus              ->fill(         0, weight*charge_weight);
            _hist_wplus_wdminus_eta_lep->fill(lepton_eta, weight*charge_weight);
            _hist_wplusd_wplusinc_pt   ->fill(    p.pT(), weight*charge_weight);
          }
          else if (lepton_charge < 0) {
            _hist_wd_charge            ->fill(        -1, weight*charge_weight);
            _hist_wd_minus             ->fill(         0, weight*charge_weight);
            _hist_wminus_wdplus_eta_lep->fill(lepton_eta, weight*charge_weight);
            _hist_wminusd_wminusinc_pt ->fill(p.pT()    , weight*charge_weight);
          }
        }

        // W+Dstar
        if ( abs(part->pdg_id()) == 413 ) {
          if (lepton_charge*part->pdg_id() > 0) charge_weight = -1;
          else charge_weight = +1;

          if (lepton_charge > 0) {
            _hist_wdstar_charge->fill(+1, weight*charge_weight);
            _hist_wd_plus->fill( 0, weight*charge_weight);
            _hist_wplus_wdstar_eta_lep->fill( lepton_eta, weight*charge_weight);
            _hist_wplusdstar_wplusinc_pt->fill(  p.pT(), weight*charge_weight);
          }
          else if (lepton_charge < 0) {
            _hist_wdstar_charge->fill(-1, weight*charge_weight);
            _hist_wd_minus->fill(0, weight*charge_weight);
            _hist_wminus_wdstar_eta_lep->fill(lepton_eta, weight*charge_weight);
            _hist_wminusdstar_wminusinc_pt->fill(p.pT(), weight*charge_weight);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection() / sumOfWeights();

      // norm to cross section
      // d01
      scale(_hist_wcjet_charge,  sf);
      scale(_hist_wd_charge,     sf);
      scale(_hist_wdstar_charge, sf);

      //d02
      scale(_hist_wcjet_plus,  sf);
      scale(_hist_wcjet_minus, sf);
      scale(_hist_wd_plus,     sf);
      scale(_hist_wd_minus,    sf);

      divide(_hist_wcjet_plus, _hist_wcjet_minus, _hist_wcjet_ratio);
      divide(_hist_wd_plus,    _hist_wd_minus,    _hist_wd_ratio   );

      //d03
      scale(_hist_wplus_wcjet_eta_lep,  sf);
      scale(_hist_wminus_wcjet_eta_lep, sf);

      //d04
      scale(_hist_wplus_wdminus_eta_lep, crossSection()/sumOfWeights());
      scale(_hist_wminus_wdplus_eta_lep, crossSection()/sumOfWeights());

      scale(_hist_wplus_wdstar_eta_lep , crossSection()/sumOfWeights());
      scale(_hist_wminus_wdstar_eta_lep, crossSection()/sumOfWeights());

      //d05
      scale(_hist_w_inc, 0.01 * sf); // in percent --> /100
      divide(_hist_wd_charge,     _hist_w_inc, _hist_wd_winc_ratio    );
      divide(_hist_wdstar_charge, _hist_w_inc, _hist_wdstar_winc_ratio);

      //d06, in percentage!
      scale(_hist_wplusd_wplusinc_pt,       sf);
      scale(_hist_wminusd_wminusinc_pt,     sf);
      scale(_hist_wplusdstar_wplusinc_pt,   sf);
      scale(_hist_wminusdstar_wminusinc_pt, sf);

      scale(_hist_wplus_winc,  0.01 * sf); // in percent --> /100
      scale(_hist_wminus_winc, 0.01 * sf); // in percent --> /100

      divide(_hist_wplusd_wplusinc_pt,       _hist_wplus_winc , _hist_wplusd_wplusinc_pt_ratio      );
      divide(_hist_wminusd_wminusinc_pt,     _hist_wminus_winc, _hist_wminusd_wminusinc_pt_ratio    );
      divide(_hist_wplusdstar_wplusinc_pt,   _hist_wplus_winc , _hist_wplusdstar_wplusinc_pt_ratio  );
      divide(_hist_wminusdstar_wminusinc_pt, _hist_wminus_winc, _hist_wminusdstar_wminusinc_pt_ratio);


      //d07
      scale(_hist_wcjet_jets, sf);

      //d08
      scale(_hist_wcjet_jets_minus, sf);
      scale(_hist_wcjet_jets_plus,  sf);
      divide(_hist_wcjet_jets_plus, _hist_wcjet_jets_minus , _hist_wcjet_jets_ratio);
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here

    // Check whether particle comes from b-decay
    /// @todo Use built-in method and avoid HepMC
    bool isFromBDecay(const Particle& p) {

      bool isfromB = false;

      if (p.genParticle() == nullptr)  return false;

      const GenParticle* part = p.genParticle();
      const GenVertex* ivtx = const_cast<const GenVertex*>(part->production_vertex());
      while (ivtx) {
        if (ivtx->particles_in_size() < 1) {
          isfromB = false;
          break;
        }
        const HepMC::GenVertex::particles_in_const_iterator iPart_invtx = ivtx->particles_in_const_begin();
        part = (*iPart_invtx);
        if (!part) {
          isfromB = false;
          break;
        }
        isfromB = PID::hasBottom(part->pdg_id());
        if (isfromB == true)  break;
        ivtx = const_cast<const GenVertex*>(part->production_vertex());
        if ( part->pdg_id() == 2212 || !ivtx )  break; // reached beam
      }
      return isfromB;
    }


    // Check whether particle has charmed children
    /// @todo Use built-in method and avoid HepMC!
    bool hasCharmedChildren(const GenParticle *part) {

      bool hasCharmedChild = false;
      if (part == nullptr)  return false;

      const GenVertex* ivtx = const_cast<const GenVertex*>(part->end_vertex());
      if (ivtx == nullptr)  return false;

      // if (ivtx->particles_out_size() < 2) return false;
      HepMC::GenVertex::particles_out_const_iterator iPart_invtx = ivtx->particles_out_const_begin();
      HepMC::GenVertex::particles_out_const_iterator end_invtx = ivtx->particles_out_const_end();

      for ( ; iPart_invtx != end_invtx; iPart_invtx++ ) {
        const GenParticle* p2 = (*iPart_invtx);
        if (p2 == part)  continue;
        hasCharmedChild = PID::hasCharm(p2->pdg_id());
        if (hasCharmedChild == true)  break;
        hasCharmedChild = hasCharmedChildren(p2);
        if (hasCharmedChild == true)  break;
      }
      return hasCharmedChild;
    }


  private:

    /// @name Histograms
    //@{

    //d01-x01-
    Histo1DPtr   _hist_wcjet_charge;
    Histo1DPtr   _hist_wd_charge;
    Histo1DPtr   _hist_wdstar_charge;

    //d02-x01-
    Scatter2DPtr _hist_wcjet_ratio;
    Scatter2DPtr _hist_wd_ratio;
    Histo1DPtr _hist_wcjet_plus;
    Histo1DPtr _hist_wd_plus;

    Histo1DPtr _hist_wcjet_minus;
    Histo1DPtr _hist_wd_minus;

    //d03-x01-
    Histo1DPtr _hist_wplus_wcjet_eta_lep;
    Histo1DPtr _hist_wminus_wcjet_eta_lep;

    //d04-x01-
    Histo1DPtr _hist_wplus_wdminus_eta_lep;
    Histo1DPtr _hist_wminus_wdplus_eta_lep;

    //d05-x01-
    Histo1DPtr _hist_wplus_wdstar_eta_lep;
    Histo1DPtr _hist_wminus_wdstar_eta_lep;


    // postprocessing histos
    //d05-x01
    Histo1DPtr _hist_w_inc;
    Scatter2DPtr _hist_wd_winc_ratio;
    Scatter2DPtr _hist_wdstar_winc_ratio;

    //d06-x01
    Histo1DPtr _hist_wplus_winc;
    Histo1DPtr _hist_wminus_winc;

    Scatter2DPtr _hist_wplusd_wplusinc_pt_ratio;
    Scatter2DPtr _hist_wminusd_wminusinc_pt_ratio;
    Scatter2DPtr _hist_wplusdstar_wplusinc_pt_ratio;
    Scatter2DPtr _hist_wminusdstar_wminusinc_pt_ratio;

    Histo1DPtr _hist_wplusd_wplusinc_pt ;
    Histo1DPtr _hist_wminusd_wminusinc_pt;
    Histo1DPtr _hist_wplusdstar_wplusinc_pt;
    Histo1DPtr _hist_wminusdstar_wminusinc_pt;

    // d07-x01
    Histo1DPtr _hist_wcjet_jets ;

    //d08-x01
    Scatter2DPtr  _hist_wcjet_jets_ratio ;
    Histo1DPtr    _hist_wcjet_jets_plus ;
    Histo1DPtr    _hist_wcjet_jets_minus;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1282447);

}
