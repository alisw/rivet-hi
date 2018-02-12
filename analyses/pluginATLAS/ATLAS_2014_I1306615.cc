// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief ATLAS H->yy differential cross-sections measurement
  ///
  /// @author Michaela Queitsch-Maitland <michaela.queitsch-maitland@cern.ch>
  //
  // arXiv: http://arxiv.org/abs/ARXIV:1407.4222
  // HepData: http://hepdata.cedar.ac.uk/view/ins1306615
  class ATLAS_2014_I1306615 : public Analysis {
  public:

    // Constructor
    ATLAS_2014_I1306615()
      : Analysis("ATLAS_2014_I1306615")
    {    }


    // Book histograms and initialise projections before the run
    void init() {

      // Final state
      // All particles within |eta| < 5.0
      const FinalState FS(Cuts::abseta<5.0);
      declare(FS,"FS");

      // Project photons with pT > 25 GeV and |eta| < 2.37
      IdentifiedFinalState ph_FS(Cuts::abseta<2.37 && Cuts::pT>25.0*GeV);
      ph_FS.acceptIdPair(PID::PHOTON);
      declare(ph_FS, "PH_FS");

      // Project photons for dressing
      IdentifiedFinalState ph_dressing_FS(FS);
      ph_dressing_FS.acceptIdPair(PID::PHOTON);

      // Project bare electrons
      IdentifiedFinalState el_bare_FS(FS);
      el_bare_FS.acceptIdPair(PID::ELECTRON);
      declare(el_bare_FS,"el_bare_FS");

      // Project dressed electrons with pT > 15 GeV and |eta| < 2.47
      DressedLeptons el_dressed_FS(ph_dressing_FS, el_bare_FS, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 15*GeV);
      declare(el_dressed_FS,"EL_DRESSED_FS");

      // Project bare muons
      IdentifiedFinalState mu_bare_FS(FS);
      mu_bare_FS.acceptIdPair(PID::MUON);

      // Project dressed muons with pT > 15 GeV and |eta| < 2.47
      //DressedLeptons mu_dressed_FS(ph_dressing_FS, mu_bare_FS, 0.1, true, -2.47, 2.47, 15.0*GeV, false);
      DressedLeptons mu_dressed_FS(ph_dressing_FS, mu_bare_FS, 0.1, Cuts::abseta < 2.47 && Cuts::pT > 15*GeV);
      declare(mu_dressed_FS,"MU_DRESSED_FS");

      // Final state excluding muons and neutrinos (for jet building and photon isolation)
      VetoedFinalState veto_mu_nu_FS(FS);
      veto_mu_nu_FS.vetoNeutrinos();
      veto_mu_nu_FS.addVetoPairId(PID::MUON);
      declare(veto_mu_nu_FS, "VETO_MU_NU_FS");

      // Build the anti-kT R=0.4 jets, using FinalState particles (vetoing muons and neutrinos)
      FastJets jets(veto_mu_nu_FS, FastJets::ANTIKT, 0.4);
      declare(jets, "JETS");

      // Book histograms
      // 1D distributions
      _h_pT_yy         = bookHisto1D(1,1,1);
      _h_y_yy          = bookHisto1D(2,1,1);
      _h_Njets30       = bookHisto1D(3,1,1);
      _h_Njets50       = bookHisto1D(4,1,1);
      _h_pT_j1         = bookHisto1D(5,1,1);
      _h_y_j1          = bookHisto1D(6,1,1);
      _h_HT            = bookHisto1D(7,1,1);
      _h_pT_j2         = bookHisto1D(8,1,1);
      _h_Dy_jj         = bookHisto1D(9,1,1);
      _h_Dphi_yy_jj    = bookHisto1D(10,1,1);
      _h_cosTS_CS      = bookHisto1D(11,1,1);
      _h_cosTS_CS_5bin = bookHisto1D(12,1,1);
      _h_Dphi_jj       = bookHisto1D(13,1,1);
      _h_pTt_yy        = bookHisto1D(14,1,1);
      _h_Dy_yy         = bookHisto1D(15,1,1);
      _h_tau_jet       = bookHisto1D(16,1,1);
      _h_sum_tau_jet   = bookHisto1D(17,1,1);
      _h_y_j2          = bookHisto1D(18,1,1);
      _h_pT_j3         = bookHisto1D(19,1,1);
      _h_m_jj          = bookHisto1D(20,1,1);
      _h_pT_yy_jj      = bookHisto1D(21,1,1);

      // 2D distributions of cosTS_CS x pT_yy
      _h_cosTS_pTyy_low  = bookHisto1D(22,1,1);
      _h_cosTS_pTyy_high = bookHisto1D(22,1,2);
      _h_cosTS_pTyy_rest = bookHisto1D(22,1,3);

      // 2D distributions of Njets x pT_yy
      _h_pTyy_Njets0 = bookHisto1D(23,1,1);
      _h_pTyy_Njets1 = bookHisto1D(23,1,2);
      _h_pTyy_Njets2 = bookHisto1D(23,1,3);

      _h_pTj1_excl = bookHisto1D(24,1,1);

      // Fiducial regions
      _h_fidXSecs = bookHisto1D(29,1,1);
    }

    // Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();
      _weight = weight;

      // Get final state particles
      const ParticleVector& FS_ptcls          = apply<FinalState>(event, "FS").particles();
      const ParticleVector& ptcls_veto_mu_nu  = apply<VetoedFinalState>(event, "VETO_MU_NU_FS").particles();
      const ParticleVector& photons           = apply<IdentifiedFinalState>(event, "PH_FS").particlesByPt();
      const vector<DressedLepton>& el_dressed = apply<DressedLeptons>(event, "EL_DRESSED_FS").dressedLeptons();
      const vector<DressedLepton>& mu_dressed = apply<DressedLeptons>(event, "MU_DRESSED_FS").dressedLeptons();

      // For isolation calculation
      float dR_iso    = 0.4;
      float ETcut_iso = 14.0;
      FourMomentum ET_iso;

      // Fiducial selection: pT > 25 GeV, |eta| < 2.37 and isolation (in cone deltaR = 0.4) is < 14 GeV
      vector<const Particle*> fid_photons;
      foreach (const Particle& ph, photons) {

	// Veto photons from hadron or tau decay
	if ( fromHadronDecay(ph) ) continue;

	// Calculate isolation
	ET_iso = - ph.momentum();
	// Loop over fs truth particles (excluding muons and neutrinos)
	foreach (const Particle& p, ptcls_veto_mu_nu) {
	  // Check if the truth particle is in a cone of 0.4
	  if ( deltaR(ph.momentum(), p.momentum()) < dR_iso )
	    ET_iso += p.momentum();
	}

	// Check isolation
	if ( ET_iso.Et() > ETcut_iso ) continue;

	// Fill vector of photons passing fiducial selection
	fid_photons.push_back(&ph);
      }

      if(fid_photons.size() < 2) vetoEvent;

      const FourMomentum& y1 = fid_photons[0]->momentum();
      const FourMomentum& y2 = fid_photons[1]->momentum();

      double m_yy = (y1 + y2).mass();

      // Relative pT cuts
      if ( y1.pT() < 0.35 * m_yy || y2.pT() < 0.25 * m_yy ) vetoEvent;

      // Mass window cut
      if ( m_yy < 105 || m_yy > 160 ) vetoEvent;

      // -------------------------------------------- //
      // Passed diphoton baseline fiducial selection! //
      // -------------------------------------------- //

      // Electron selection
      vector<const Particle*> good_el;
      foreach(const DressedLepton& els, el_dressed)  {

	const Particle& el = els.constituentLepton();
	if ( el.momentum().pT()        < 15   ) continue;
	if ( fabs(el.momentum().eta()) > 2.47 ) continue;
	if ( deltaR(el.momentum(), y1) < 0.4  ) continue;
	if ( deltaR(el.momentum(), y2) < 0.4  ) continue;
	if ( fromHadronDecay(el)              ) continue; // Veto electrons from hadron or tau decay
	good_el.push_back(&el);
      }

      // Muon selection
      vector<const Particle*> good_mu;
      foreach(const DressedLepton& mus, mu_dressed)  {

	const Particle& mu = mus.constituentLepton();
	if ( mu.momentum().pT()        < 15   ) continue;
	if ( fabs(mu.momentum().eta()) > 2.47 ) continue;
	if ( deltaR(mu.momentum(), y1) < 0.4  ) continue;
	if ( deltaR(mu.momentum(), y2) < 0.4  ) continue;
	if ( fromHadronDecay(mu)              ) continue; // Veto muons from hadron or tau decay
	good_mu.push_back(&mu);
      }

      // Find prompt, invisible particles for missing ET calculation
      // Based on VisibleFinalState projection
      FourMomentum invisible(0,0,0,0);
      foreach (const Particle& p, FS_ptcls) {

	// Veto non-prompt particles (from hadron or tau decay)
        if ( fromHadronDecay(p) ) continue;
	// Charged particles are visible
	if ( PID::threeCharge( p.pid() ) != 0 ) continue;
	// Neutral hadrons are visible
	if ( PID::isHadron( p.pid() ) ) continue;
        // Photons are visible
	if ( p.pid() == PID::PHOTON ) continue;
        // Gluons are visible (for parton level analyses)
	if ( p.pid() == PID::GLUON ) continue;
	// Everything else is invisible
	invisible += p.momentum();
      }
      double MET = invisible.Et();

      // Jet selection
      // Get jets with pT > 25 GeV and |rapidity| < 4.4
      //const Jets& jets = apply<FastJets>(event, "JETS").jetsByPt(25.0*GeV, MAXDOUBLE, -4.4, 4.4, RAPIDITY);
      const Jets& jets = apply<FastJets>(event, "JETS").jetsByPt(Cuts::pT>25.0*GeV && Cuts::absrap <4.4);

      vector<const Jet*> jets_25;
      vector<const Jet*> jets_30;
      vector<const Jet*> jets_50;

      foreach (const Jet& jet, jets) {

	bool passOverlap = true;
	// Overlap with leading photons
	if ( deltaR(y1, jet.momentum()) < 0.4 ) passOverlap = false;
	if ( deltaR(y2, jet.momentum()) < 0.4 ) passOverlap = false;

	// Overlap with good electrons
	foreach (const Particle* el, good_el)
	  if ( deltaR(el->momentum(), jet.momentum()) < 0.2 ) passOverlap = false;

	if ( ! passOverlap ) continue;

	if ( fabs(jet.momentum().eta()) < 2.4 || ( fabs(jet.momentum().eta()) > 2.4 && jet.momentum().pT() > 30 ) ) jets_25.push_back(&jet);
	if ( jet.momentum().pT() > 30 ) jets_30.push_back(&jet);
	if ( jet.momentum().pT() > 50 ) jets_50.push_back(&jet);
      }

      // Fiducial regions
      _h_fidXSecs->fill(1,_weight);
      if ( jets_30.size() >= 1 ) _h_fidXSecs->fill(2, _weight);
      if ( jets_30.size() >= 2 ) _h_fidXSecs->fill(3, _weight);
      if ( jets_30.size() >= 3 ) _h_fidXSecs->fill(4, _weight);
      if ( jets_30.size() >= 2 && passVBFCuts(y1 + y2, jets_30.at(0)->momentum(), jets_30.at(1)->momentum()) ) _h_fidXSecs->fill(5, _weight);
      if ( (good_el.size() + good_mu.size()) > 0 ) _h_fidXSecs->fill(6, _weight);
      if ( MET > 80 ) _h_fidXSecs->fill(7, _weight);

      // Fill histograms
      // Inclusive variables
      _pT_yy    = (y1 + y2).pT();
      _y_yy     = fabs( (y1 + y2).rapidity() );
      _cosTS_CS = cosTS_CS(y1, y2);
      _pTt_yy   = pTt(y1, y2);
      _Dy_yy    = fabs( deltaRap(y1, y2) );

      _Njets30 = jets_30.size() > 3 ? 3 : jets_30.size();
      _Njets50 = jets_50.size() > 3 ? 3 : jets_50.size();
      _h_Njets30->fill(_Njets30, _weight);
      _h_Njets50->fill(_Njets50, _weight);

      _pT_j1 = jets_30.size() > 0 ? jets_30.at(0)->momentum().pT() : 0.;
      _pT_j2 = jets_30.size() > 1 ? jets_30.at(1)->momentum().pT() : 0.;
      _pT_j3 = jets_30.size() > 2 ? jets_30.at(2)->momentum().pT() : 0.;

      _HT = 0.0;
      foreach (const Jet* jet, jets_30)
	_HT += jet->momentum().pT();

      _tau_jet     = tau_jet_max(y1 + y2, jets_25);
      _sum_tau_jet = sum_tau_jet(y1 + y2, jets_25);

      _h_pT_yy        ->fill(_pT_yy    ,_weight);
      _h_y_yy         ->fill(_y_yy     ,_weight);
      _h_pT_j1        ->fill(_pT_j1    ,_weight);
      _h_cosTS_CS     ->fill(_cosTS_CS ,_weight);
      _h_cosTS_CS_5bin->fill(_cosTS_CS ,_weight);
      _h_HT           ->fill(_HT       ,_weight);
      _h_pTt_yy       ->fill(_pTt_yy   ,_weight);
      _h_Dy_yy        ->fill(_Dy_yy    ,_weight);
      _h_tau_jet      ->fill(_tau_jet  ,_weight);
      _h_sum_tau_jet  ->fill(_sum_tau_jet,_weight);

      // >=1 jet variables
      if ( jets_30.size() >= 1 ) {
	FourMomentum j1 = jets_30[0]->momentum();
	_y_j1 = fabs( j1.rapidity() );

	_h_pT_j2->fill(_pT_j2 ,_weight);
	_h_y_j1 ->fill(_y_j1  ,_weight);
      }

      // >=2 jet variables
      if ( jets_30.size() >= 2 ) {
	FourMomentum j1 = jets_30[0]->momentum();
	FourMomentum j2 = jets_30[1]->momentum();

	_Dy_jj      = fabs( deltaRap(j1, j2) );
	_Dphi_jj    = fabs( deltaPhi(j1, j2) );
	_Dphi_yy_jj = fabs( deltaPhi(y1 + y2, j1 + j2) );
	_m_jj       = (j1 + j2).mass();
	_pT_yy_jj   = (y1 + y2 + j1 + j2).pT();
	_y_j2       = fabs( j2.rapidity() );

	_h_Dy_jj      ->fill(_Dy_jj     ,_weight);
	_h_Dphi_jj    ->fill(_Dphi_jj   ,_weight);
	_h_Dphi_yy_jj ->fill(_Dphi_yy_jj,_weight);
	_h_m_jj       ->fill(_m_jj      ,_weight);
	_h_pT_yy_jj   ->fill(_pT_yy_jj  ,_weight);
	_h_pT_j3      ->fill(_pT_j3     ,_weight);
	_h_y_j2       ->fill(_y_j2      ,_weight);
      }

      // 2D distributions of cosTS_CS x pT_yy
      if ( _pT_yy < 80 )
	_h_cosTS_pTyy_low->fill(_cosTS_CS, _weight);
      else if ( _pT_yy > 80 && _pT_yy < 200 )
	_h_cosTS_pTyy_high->fill(_cosTS_CS,_weight);
      else if ( _pT_yy > 200 )
	_h_cosTS_pTyy_rest->fill(_cosTS_CS,_weight);

      // 2D distributions of pT_yy x Njets
      if ( _Njets30 == 0 )
	_h_pTyy_Njets0->fill(_pT_yy, _weight);
      else if ( _Njets30 == 1 )
	_h_pTyy_Njets1->fill(_pT_yy, _weight);
      else if ( _Njets30 >= 2 )
	_h_pTyy_Njets2->fill(_pT_yy, _weight);

      if ( _Njets30 == 1 ) _h_pTj1_excl->fill(_pT_j1, _weight);

    }

    // Normalise histograms after the run
    void finalize() {

      const double xs = crossSectionPerEvent()/femtobarn;

      scale(_h_pT_yy, xs);
      scale(_h_y_yy, xs);
      scale(_h_pT_j1, xs);
      scale(_h_y_j1, xs);
      scale(_h_HT, xs);
      scale(_h_pT_j2, xs);
      scale(_h_Dy_jj, xs);
      scale(_h_Dphi_yy_jj, xs);
      scale(_h_cosTS_CS, xs);
      scale(_h_cosTS_CS_5bin, xs);
      scale(_h_Dphi_jj, xs);
      scale(_h_pTt_yy, xs);
      scale(_h_Dy_yy, xs);
      scale(_h_tau_jet, xs);
      scale(_h_sum_tau_jet, xs);
      scale(_h_y_j2, xs);
      scale(_h_pT_j3, xs);
      scale(_h_m_jj, xs);
      scale(_h_pT_yy_jj, xs);
      scale(_h_cosTS_pTyy_low, xs);
      scale(_h_cosTS_pTyy_high, xs);
      scale(_h_cosTS_pTyy_rest, xs);
      scale(_h_pTyy_Njets0, xs);
      scale(_h_pTyy_Njets1, xs);
      scale(_h_pTyy_Njets2, xs);
      scale(_h_pTj1_excl, xs);
      scale(_h_Njets30, xs);
      scale(_h_Njets50, xs);
      scale(_h_fidXSecs, xs);
    }

    // Trace event record to see if particle came from a hadron (or a tau from a hadron decay)
    // Based on fromDecay() function
    bool fromHadronDecay(const Particle& p ) {
      if (p.genParticle() == NULL) return false;
      const GenVertex* prodVtx = p.genParticle()->production_vertex();
      if (prodVtx == NULL) return false;
      foreach (const GenParticle* ancestor, particles(prodVtx, HepMC::ancestors)) {
        const PdgId pid = ancestor->pdg_id();
        if (ancestor->status() == 2 && PID::isHadron(pid)) return true;
        if (ancestor->status() == 2 && (abs(pid) == PID::TAU && fromHadronDecay(ancestor))) return true;
      }
      return false;
    }

    // VBF-enhanced dijet topology selection cuts
    bool passVBFCuts(const FourMomentum &H, const FourMomentum &j1, const FourMomentum &j2) {
      return ( fabs(deltaRap(j1, j2)) > 2.8 && (j1 + j2).mass() > 400 && fabs(deltaPhi(H, j1 + j2)) > 2.6 );
    }

    // Cosine of the decay angle in the Collins-Soper frame
    double cosTS_CS(const FourMomentum &y1, const FourMomentum &y2) {
      return fabs( ( (y1.E() + y1.pz())* (y2.E() - y2.pz()) - (y1.E() - y1.pz()) * (y2.E() + y2.pz()) )
		   / ((y1 + y2).mass() * sqrt(pow((y1 + y2).mass(), 2) + pow((y1 + y2).pt(), 2)) ) );
    }

    // Diphoton pT along thrust axis
    double pTt(const FourMomentum &y1, const FourMomentum &y2) {
      return fabs(y1.px() * y2.py() - y2.px() * y1.py()) / (y1 - y2).pT()*2;
    }

    // Tau of jet  (see paper for description)
    // tau_jet = mT/(2*cosh(y*)), where mT = pT (+) m, and y* = rapidty in Higgs rest frame
    double tau_jet( const FourMomentum &H, const FourMomentum &jet ) {
      return sqrt( pow(jet.pT(),2) + pow(jet.mass(),2) ) / (2.0 * cosh( jet.rapidity() - H.rapidity() ) );
    }

    // Maximal (leading) tau_jet (see paper for description)
    double tau_jet_max( const FourMomentum &H, const vector<const Jet*> jets, double tau_jet_cut = 8. ) {
      double max_tj = 0;
      for (size_t i=0; i < jets.size(); ++i) {
	FourMomentum jet = jets[i]->momentum();
	if ( tau_jet(H, jet) > tau_jet_cut )
	  max_tj = max( tau_jet(H, jet), max_tj );
      }
      return max_tj;
    }

    // Scalar sum of tau for all jets (see paper for description)
    double sum_tau_jet( const FourMomentum &H, const vector<const Jet*> jets, double tau_jet_cut = 8. ) {
      double sum_tj = 0;
      for (size_t i=0; i < jets.size(); ++i) {
	FourMomentum jet = jets[i]->momentum();
	if ( tau_jet(H, jet) > tau_jet_cut )
	  sum_tj += tau_jet(H, jet);
      }
      return sum_tj;
    }

  private:

    Histo1DPtr _h_pT_yy;
    Histo1DPtr _h_y_yy;
    Histo1DPtr _h_Njets30;
    Histo1DPtr _h_Njets50;
    Histo1DPtr _h_pT_j1;
    Histo1DPtr _h_y_j1;
    Histo1DPtr _h_HT;
    Histo1DPtr _h_pT_j2;
    Histo1DPtr _h_Dy_jj;
    Histo1DPtr _h_Dphi_yy_jj;
    Histo1DPtr _h_cosTS_CS;
    Histo1DPtr _h_cosTS_CS_5bin;
    Histo1DPtr _h_Dphi_jj;
    Histo1DPtr _h_pTt_yy;
    Histo1DPtr _h_Dy_yy;
    Histo1DPtr _h_tau_jet;
    Histo1DPtr _h_sum_tau_jet;
    Histo1DPtr _h_y_j2;
    Histo1DPtr _h_pT_j3;
    Histo1DPtr _h_m_jj;
    Histo1DPtr _h_pT_yy_jj;
    Histo1DPtr _h_cosTS_pTyy_low;
    Histo1DPtr _h_cosTS_pTyy_high;
    Histo1DPtr _h_cosTS_pTyy_rest;
    Histo1DPtr _h_pTyy_Njets0;
    Histo1DPtr _h_pTyy_Njets1;
    Histo1DPtr _h_pTyy_Njets2;
    Histo1DPtr _h_pTj1_excl;
    Histo1DPtr _h_fidXSecs;

    double _weight;
    int _Njets30;
    int _Njets50;
    double _pT_yy;
    double _y_yy;
    double _cosTS_CS;
    double _pT_j1;
    double _m_jj;
    double _y_j1;
    double _HT;
    double _pT_j2;
    double _y_j2;
    double _Dphi_yy_jj;
    double _pT_yy_jj;
    double _Dphi_jj;
    double _Dy_jj;
    double _pT_j3;
    double _pTt_yy;
    double _Dy_yy;
    double _tau_jet;
    double _sum_tau_jet;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1306615);

}
