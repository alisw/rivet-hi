#include "Rivet/Analysis.hh"
#include "Rivet/Math/LorentzTrans.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  namespace { //< only visible in this compilation unit

    /// @brief Pseudo top finder
    ///
    /// Find top quark in the particle level.
    /// The definition is based on the agreement at the LHC working group.
    class PseudoTop : public FinalState {
    public:
      /// @name Standard constructors and destructors.
      //@{

      /// The default constructor. May specify the minimum and maximum
      /// pseudorapidity \f$ \eta \f$ and the min \f$ p_T \f$ (in GeV).
      PseudoTop(double lepR = 0.1, double lepMinPt = 20, double lepMaxEta = 2.4,
                double jetR = 0.4, double jetMinPt = 30, double jetMaxEta = 4.7)
        : FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV),
          _lepR(lepR), _lepMinPt(lepMinPt), _lepMaxEta(lepMaxEta),
          _jetR(jetR), _jetMinPt(jetMinPt), _jetMaxEta(jetMaxEta)
      {
        setName("PseudoTop");
      }

      enum TTbarMode {CH_NONE=-1, CH_FULLHADRON = 0, CH_SEMILEPTON, CH_FULLLEPTON};
      enum DecayMode {CH_HADRON = 0, CH_MUON, CH_ELECTRON};

      TTbarMode mode() const {
        if (!_isValid) return CH_NONE;
        if (_mode1 == CH_HADRON && _mode2 == CH_HADRON) return CH_FULLHADRON;
        else if ( _mode1 != CH_HADRON && _mode2 != CH_HADRON) return CH_FULLLEPTON;
        else return CH_SEMILEPTON;
      }
      DecayMode mode1() const {return _mode1;}
      DecayMode mode2() const {return _mode2;}

      /// Clone on the heap.
      virtual unique_ptr<Projection> clone() const {
        return unique_ptr<Projection>(new PseudoTop(*this));
      }

      //@}

    public:
      Particle t1() const {return _t1;}
      Particle t2() const {return _t2;}
      Particle b1() const {return _b1;}
      Particle b2() const {return _b2;}
      ParticleVector wDecays1() const {return _wDecays1;}
      ParticleVector wDecays2() const {return _wDecays2;}
      Jets jets() const {return _jets;}
      Jets bjets() const {return _bjets;}
      Jets ljets() const {return _ljets;}

    protected:
      // Apply the projection to the event
      void project(const Event& e); // override; ///< @todo Re-enable when C++11 allowed
      void cleanup(std::map<double, std::pair<size_t, size_t> >& v, const bool doCrossCleanup=false) const;

    private:
      const double _lepR, _lepMinPt, _lepMaxEta;
      const double _jetR, _jetMinPt, _jetMaxEta;

      //constexpr ///< @todo Re-enable when C++11 allowed
      static double _tMass; // = 172.5*GeV; ///< @todo Re-enable when C++11 allowed
      //constexpr ///< @todo Re-enable when C++11 allowed
      static double _wMass; // = 80.4*GeV; ///< @todo Re-enable when C++11 allowed

    private:
      bool _isValid;
      DecayMode _mode1, _mode2;

      Particle _t1, _t2;
      Particle _b1, _b2;
      ParticleVector _wDecays1, _wDecays2;
      Jets _jets, _bjets, _ljets;

    };

    // More implementation below the analysis code

  }



  /// Pseudo-top analysis from CMS
  class CMS_2015_I1370682 : public Analysis {
  public:

    CMS_2015_I1370682()
      : Analysis("CMS_2015_I1370682"),
        _applyCorrection(true),
        _doShapeOnly(true)
    {    }


    void init() {
      declare(PseudoTop(0.1, 20, 2.4, 0.5, 30, 2.4), "ttbar");

      // Lepton + Jet channel
      _hSL_topPt         = bookHisto1D("d15-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hSL_topPtTtbarSys = bookHisto1D("d16-x01-y01"); // 1/sigma dsigma/dpt*(top)
      _hSL_topY          = bookHisto1D("d17-x01-y01"); // 1/sigma dsigma/dy(top)
      _hSL_ttbarDelPhi   = bookHisto1D("d18-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
      _hSL_topPtLead     = bookHisto1D("d19-x01-y01"); // 1/sigma dsigma/dpt(t1)
      _hSL_topPtSubLead  = bookHisto1D("d20-x01-y01"); // 1/sigma dsigma/dpt(t2)
      _hSL_ttbarPt       = bookHisto1D("d21-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
      _hSL_ttbarY        = bookHisto1D("d22-x01-y01"); // 1/sigma dsigma/dy(ttbar)
      _hSL_ttbarMass     = bookHisto1D("d23-x01-y01"); // 1/sigma dsigma/dm(ttbar)

      // Dilepton channel
      _hDL_topPt         = bookHisto1D("d24-x01-y01"); // 1/sigma dsigma/dpt(top)
      _hDL_topPtTtbarSys = bookHisto1D("d25-x01-y01"); // 1/sigma dsigma/dpt*(top)
      _hDL_topY          = bookHisto1D("d26-x01-y01"); // 1/sigma dsigma/dy(top)
      _hDL_ttbarDelPhi   = bookHisto1D("d27-x01-y01"); // 1/sigma dsigma/ddeltaphi(t,tbar)
      _hDL_topPtLead     = bookHisto1D("d28-x01-y01"); // 1/sigma dsigma/dpt(t1)
      _hDL_topPtSubLead  = bookHisto1D("d29-x01-y01"); // 1/sigma dsigma/dpt(t2)
      _hDL_ttbarPt       = bookHisto1D("d30-x01-y01"); // 1/sigma dsigma/dpt(ttbar)
      _hDL_ttbarY        = bookHisto1D("d31-x01-y01"); // 1/sigma dsigma/dy(ttbar)
      _hDL_ttbarMass     = bookHisto1D("d32-x01-y01"); // 1/sigma dsigma/dm(ttbar)

    }


    void analyze(const Event& event) {

      // Get the ttbar candidate
      const PseudoTop& ttbar = apply<PseudoTop>(event, "ttbar");
      if ( ttbar.mode() == PseudoTop::CH_NONE ) vetoEvent;

      const FourMomentum& t1P4 = ttbar.t1().momentum();
      const FourMomentum& t2P4 = ttbar.t2().momentum();
      const double pt1 = std::max(t1P4.pT(), t2P4.pT());
      const double pt2 = std::min(t1P4.pT(), t2P4.pT());
      const double dPhi = deltaPhi(t1P4, t2P4);
      const FourMomentum ttP4 = t1P4 + t2P4;
      const FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttP4.betaVec()).transform(t1P4);

      const double weight = event.weight();

      if ( ttbar.mode() == PseudoTop::CH_SEMILEPTON ) {
        const Particle lCand1 = ttbar.wDecays1()[0]; // w1 dau0 is the lepton in the PseudoTop
        if (lCand1.pT() < 33*GeV || lCand1.abseta() > 2.1) vetoEvent;
        _hSL_topPt->fill(t1P4.pT(), weight);
        _hSL_topPt->fill(t2P4.pT(), weight);
        _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
        _hSL_topY->fill(t1P4.rapidity(), weight);
        _hSL_topY->fill(t2P4.rapidity(), weight);
        _hSL_ttbarDelPhi->fill(dPhi, weight);
        _hSL_topPtLead->fill(pt1, weight);
        _hSL_topPtSubLead->fill(pt2, weight);
        _hSL_ttbarPt->fill(ttP4.pT(), weight);
        _hSL_ttbarY->fill(ttP4.rapidity(), weight);
        _hSL_ttbarMass->fill(ttP4.mass(), weight);
      }
      else if ( ttbar.mode() == PseudoTop::CH_FULLLEPTON ) {
        const Particle lCand1 = ttbar.wDecays1()[0]; // dau0 are the lepton in the PseudoTop
        const Particle lCand2 = ttbar.wDecays2()[0]; // dau0 are the lepton in the PseudoTop
        if (lCand1.pT() < 20*GeV || lCand1.abseta() > 2.4) vetoEvent;
        if (lCand2.pT() < 20*GeV || lCand2.abseta() > 2.4) vetoEvent;
        _hDL_topPt->fill(t1P4.pT(), weight);
        _hDL_topPt->fill(t2P4.pT(), weight);
        _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
        _hDL_topY->fill(t1P4.rapidity(), weight);
        _hDL_topY->fill(t2P4.rapidity(), weight);
        _hDL_ttbarDelPhi->fill(dPhi, weight);
        _hDL_topPtLead->fill(pt1, weight);
        _hDL_topPtSubLead->fill(pt2, weight);
        _hDL_ttbarPt->fill(ttP4.pT(), weight);
        _hDL_ttbarY->fill(ttP4.rapidity(), weight);
        _hDL_ttbarMass->fill(ttP4.mass(), weight);
      }

    }


    void finalize() {
      if ( _applyCorrection ) {
        // Correction functions for TOP-12-028 paper, (parton bin height)/(pseudotop bin height)
        const double ch15[] = { 5.473609, 4.941048, 4.173346, 3.391191, 2.785644, 2.371346, 2.194161, 2.197167, };
        const double ch16[] = { 5.470905, 4.948201, 4.081982, 3.225532, 2.617519, 2.239217, 2.127878, 2.185918, };
        const double ch17[] = { 10.003667, 4.546519, 3.828115, 3.601018, 3.522194, 3.524694, 3.600951, 3.808553, 4.531891, 9.995370, };
        const double ch18[] = { 4.406683, 4.054041, 3.885393, 4.213646, };
        const double ch19[] = { 6.182537, 5.257703, 4.422280, 3.568402, 2.889408, 2.415878, 2.189974, 2.173210, };
        const double ch20[] = { 5.199874, 4.693318, 3.902882, 3.143785, 2.607877, 2.280189, 2.204124, 2.260829, };
        const double ch21[] = { 6.053523, 3.777506, 3.562251, 3.601356, 3.569347, 3.410472, };
        const double ch22[] = { 11.932351, 4.803773, 3.782709, 3.390775, 3.226806, 3.218982, 3.382678, 3.773653, 4.788191, 11.905338, };
        const double ch23[] = { 7.145255, 5.637595, 4.049882, 3.025917, 2.326430, 1.773824, 1.235329, };

        const double ch24[] = { 2.268193, 2.372063, 2.323975, 2.034655, 1.736793, };
        const double ch25[] = { 2.231852, 2.383086, 2.341894, 2.031318, 1.729672, 1.486993, };
        const double ch26[] = { 3.993526, 2.308249, 2.075136, 2.038297, 2.036302, 2.078270, 2.295817, 4.017713, };
        const double ch27[] = { 2.205978, 2.175010, 2.215376, 2.473144, };
        const double ch28[] = { 2.321077, 2.371895, 2.338871, 2.057821, 1.755382, };
        const double ch29[] = { 2.222707, 2.372591, 2.301688, 1.991162, 1.695343, };
        const double ch30[] = { 2.599677, 2.026855, 2.138620, 2.229553, };
        const double ch31[] = { 5.791779, 2.636219, 2.103642, 1.967198, 1.962168, 2.096514, 2.641189, 5.780828, };
        const double ch32[] = { 2.006685, 2.545525, 2.477745, 2.335747, 2.194226, 2.076500, };

        applyCorrection(_hSL_topPt, ch15);
        applyCorrection(_hSL_topPtTtbarSys, ch16);
        applyCorrection(_hSL_topY, ch17);
        applyCorrection(_hSL_ttbarDelPhi, ch18);
        applyCorrection(_hSL_topPtLead, ch19);
        applyCorrection(_hSL_topPtSubLead, ch20);
        applyCorrection(_hSL_ttbarPt, ch21);
        applyCorrection(_hSL_ttbarY, ch22);
        applyCorrection(_hSL_ttbarMass, ch23);

        applyCorrection(_hDL_topPt, ch24);
        applyCorrection(_hDL_topPtTtbarSys, ch25);
        applyCorrection(_hDL_topY, ch26);
        applyCorrection(_hDL_ttbarDelPhi, ch27);
        applyCorrection(_hDL_topPtLead, ch28);
        applyCorrection(_hDL_topPtSubLead, ch29);
        applyCorrection(_hDL_ttbarPt, ch30);
        applyCorrection(_hDL_ttbarY, ch31);
        applyCorrection(_hDL_ttbarMass, ch32);
      }

      if ( _doShapeOnly ) {
        normalize(_hSL_topPt        );
        normalize(_hSL_topPtTtbarSys);
        normalize(_hSL_topY         );
        normalize(_hSL_ttbarDelPhi  );
        normalize(_hSL_topPtLead    );
        normalize(_hSL_topPtSubLead );
        normalize(_hSL_ttbarPt      );
        normalize(_hSL_ttbarY       );
        normalize(_hSL_ttbarMass    );

        normalize(_hDL_topPt        );
        normalize(_hDL_topPtTtbarSys);
        normalize(_hDL_topY         );
        normalize(_hDL_ttbarDelPhi  );
        normalize(_hDL_topPtLead    );
        normalize(_hDL_topPtSubLead );
        normalize(_hDL_ttbarPt      );
        normalize(_hDL_ttbarY       );
        normalize(_hDL_ttbarMass    );
      }
      else {
        const double s = 1./sumOfWeights();
        scale(_hSL_topPt        , s);
        scale(_hSL_topPtTtbarSys, s);
        scale(_hSL_topY         , s);
        scale(_hSL_ttbarDelPhi  , s);
        scale(_hSL_topPtLead    , s);
        scale(_hSL_topPtSubLead , s);
        scale(_hSL_ttbarPt      , s);
        scale(_hSL_ttbarY       , s);
        scale(_hSL_ttbarMass    , s);
        scale(_hDL_topPt        , s);
        scale(_hDL_topPtTtbarSys, s);
        scale(_hDL_topY         , s);
        scale(_hDL_ttbarDelPhi  , s);
        scale(_hDL_topPtLead    , s);
        scale(_hDL_topPtSubLead , s);
        scale(_hDL_ttbarPt      , s);
        scale(_hDL_ttbarY       , s);
        scale(_hDL_ttbarMass    , s);
      }

    }


    void applyCorrection(Histo1DPtr h, const double* cf) {
      vector<YODA::HistoBin1D>& bins = h->bins();
      for (size_t i=0, n=bins.size(); i<n; ++i ) {
        const double s = cf[i];
        YODA::HistoBin1D& bin = bins[i];
        bin.scaleW(s);
      }
    }


  private:

    const bool _applyCorrection, _doShapeOnly;
    Histo1DPtr _hSL_topPt, _hSL_topPtTtbarSys, _hSL_topY, _hSL_ttbarDelPhi, _hSL_topPtLead,
      _hSL_topPtSubLead, _hSL_ttbarPt, _hSL_ttbarY, _hSL_ttbarMass;
    Histo1DPtr _hDL_topPt, _hDL_topPtTtbarSys, _hDL_topY, _hDL_ttbarDelPhi, _hDL_topPtLead,
      _hDL_topPtSubLead, _hDL_ttbarPt, _hDL_ttbarY, _hDL_ttbarMass;

  };



  DECLARE_RIVET_PLUGIN(CMS_2015_I1370682);


  ///////////////

  // More PseudoTop implementation
  namespace {


    double PseudoTop::_tMass = 172.5*GeV;
    double PseudoTop::_wMass = 80.4*GeV;


    void PseudoTop::cleanup(map<double, pair<size_t, size_t> >& v, const bool doCrossCleanup) const {
      vector<map<double, pair<size_t, size_t> >::iterator> toErase;
      set<size_t> usedLeg1, usedLeg2;
      if ( !doCrossCleanup ) {
        /// @todo Reinstate when C++11 allowed: for (auto key = v.begin(); key != v.end(); ++key) {
        for (map<double, pair<size_t, size_t> >::iterator key = v.begin(); key != v.end(); ++key) {
          const size_t leg1 = key->second.first;
          const size_t leg2 = key->second.second;
          if (usedLeg1.find(leg1) == usedLeg1.end() and
              usedLeg2.find(leg2) == usedLeg2.end()) {
            usedLeg1.insert(leg1);
            usedLeg2.insert(leg2);
          } else {
            toErase.push_back(key);
          }
        }
      }
      else {
        /// @todo Reinstate when C++11 allowed: for (auto key = v.begin(); key != v.end(); ++key) {
        for (map<double, pair<size_t, size_t> >::iterator key = v.begin(); key != v.end(); ++key) {
          const size_t leg1 = key->second.first;
          const size_t leg2 = key->second.second;
          if (usedLeg1.find(leg1) == usedLeg1.end() and
              usedLeg1.find(leg2) == usedLeg1.end()) {
            usedLeg1.insert(leg1);
            usedLeg1.insert(leg2);
          } else {
            toErase.push_back(key);
          }
        }
      }
      /// @todo Reinstate when C++11 allowed:  for (auto& key : toErase) v.erase(key);
      for (size_t i = 0; i < toErase.size(); ++i) v.erase(toErase[i]);
    }


    void PseudoTop::project(const Event& e) {
      // Leptons : do the lepton clustering anti-kt R=0.1 using stable photons and leptons not from hadron decay
      // Neutrinos : neutrinos not from hadron decay
      // MET : vector sum of all invisible particles in x-y plane
      // Jets : anti-kt R=0.4 using all particles excluding neutrinos and particles used in lepton clustering
      //        add ghost B hadrons during the jet clustering to identify B jets.

      // W->lv : dressed lepton and neutrino pairs
      // W->jj : light flavored dijet
      // W candidate : select lv or jj pairs which minimise |mW1-80.4|+|mW2-80.4|
      //               lepton-neutrino pair will be selected with higher priority

      // t->Wb : W candidate + b jet
      // t candidate : select Wb pairs which minimise |mtop1-172.5|+|mtop2-172.5|

      _isValid = false;
      _theParticles.clear();
      _wDecays1.clear();
      _wDecays2.clear();
      _jets.clear();
      _bjets.clear();
      _ljets.clear();
      _mode1 = _mode2 = CH_HADRON;

      // Collect final state particles
      Particles pForLep, pForJet;
      Particles neutrinos; // Prompt neutrinos
      /// @todo Avoid this unsafe jump into HepMC -- all this can be done properly via VisibleFS and HeavyHadrons projections
      for (const GenParticle* p : Rivet::particles(e.genEvent())) {
        const int status = p->status();
        const int pdgId = p->pdg_id();
        if (status == 1) {
          Particle rp = *p;
          if (!PID::isHadron(pdgId) && !rp.fromHadron()) {
            // Collect particles not from hadron decay
            if (rp.isNeutrino()) {
              // Prompt neutrinos are kept in separate collection
              neutrinos.push_back(rp);
            } else if (pdgId == 22 || rp.isLepton()) {
              // Leptons and photons for the dressing
              pForLep.push_back(rp);
            }
          } else if (!rp.isNeutrino()) {
            // Use all particles from hadron decay
            pForJet.push_back(rp);
          }
        } else if (PID::isHadron(pdgId) && PID::hasBottom(pdgId)) {
          // NOTE: Consider B hadrons with pT > 5GeV - not in CMS proposal
          //if ( p->momentum().perp() < 5 ) continue;

          // Do unstable particles, to be used in the ghost B clustering
          // Use last B hadrons only
          bool isLast = true;
          for (const GenParticlePtr pp : Rivet::particles(p->end_vertex(), HepMC::children)) {
            if (PID::hasBottom(pp->pdg_id())) {
              isLast = false;
              break;
            }
          }
          if (!isLast) continue;

          // Rescale momentum by 10^-20
          Particle ghost(pdgId, FourMomentum(p->momentum())*1e-20/p->momentum().rho());
          pForJet.push_back(ghost);
        }
      }

      // Start object building from trivial thing - prompt neutrinos
      sortByPt(neutrinos);

      // Proceed to lepton dressing
      FastJets fjLep(FinalState(), FastJets::ANTIKT, _lepR);
      fjLep.calc(pForLep);

      Jets leptons;
      vector<int> leptonsId;
      set<int> dressedIdxs;
      for (const Jet& lep : fjLep.jetsByPt(_lepMinPt)) {
        if (lep.abseta() > _lepMaxEta) continue;
        double leadingPt = -1;
        int leptonId = 0;
        for (const Particle& p : lep.particles()) {
          /// @warning Barcodes aren't future-proof in HepMC
          dressedIdxs.insert(p.genParticle()->barcode());
          if (p.isLepton() && p.pT() > leadingPt) {
            leadingPt = p.pT();
            leptonId = p.pid();
          }
        }
        if (leptonId == 0) continue;
        leptons.push_back(lep);
        leptonsId.push_back(leptonId);
      }

      // Re-use particles not used in lepton dressing
      for (const Particle& rp : pForLep) {
        /// @warning Barcodes aren't future-proof in HepMC
        const int barcode = rp.genParticle()->barcode();
        // Skip if the particle is used in dressing
        if (dressedIdxs.find(barcode) != dressedIdxs.end()) continue;
        // Put back to be used in jet clustering
        pForJet.push_back(rp);
      }

      // Then do the jet clustering
      FastJets fjJet(FinalState(), FastJets::ANTIKT, _jetR);
      //fjJet.useInvisibles(); // NOTE: CMS proposal to remove neutrinos (AB: wouldn't work anyway, since they were excluded from clustering inputs)
      fjJet.calc(pForJet);
      for (const Jet& jet : fjJet.jetsByPt(_jetMinPt)) {
        if (jet.abseta() > _jetMaxEta) continue;
        _jets.push_back(jet);
        bool isBJet = false;
        for (const Particle& rp : jet.particles()) {
          if (PID::hasBottom(rp.pdgId())) {
            isBJet = true;
            break;
          }
        }
        if ( isBJet ) _bjets.push_back(jet);
        else _ljets.push_back(jet);
      }

      // Every building blocks are ready. Continue to pseudo-W and pseudo-top combination

      if (_bjets.size() < 2) return; // Ignore single top for now
      map<double, pair<size_t, size_t> > wLepCandIdxs;
      map<double, pair<size_t, size_t> > wHadCandIdxs;

      // Collect leptonic-decaying W's
      for (size_t iLep = 0, nLep = leptons.size(); iLep < nLep; ++iLep) {
        const Jet& lep = leptons.at(iLep);
        for (size_t iNu = 0, nNu = neutrinos.size(); iNu < nNu; ++iNu) {
          const Particle& nu = neutrinos.at(iNu);
          const double m = (lep.momentum()+nu.momentum()).mass();
          const double dm = std::abs(m-_wMass);
          wLepCandIdxs[dm] = make_pair(iLep, iNu);
        }
      }

      // Continue to hadronic decaying W's
      for (size_t i = 0, nLjet = _ljets.size(); i < nLjet; ++i) {
        const Jet& ljet1 = _ljets[i];
        for (size_t j = i+1; j < nLjet; ++j) {
          const Jet& ljet2 = _ljets[j];
          const double m = (ljet1.momentum()+ljet2.momentum()).mass();
          const double dm = std::abs(m-_wMass);
          wHadCandIdxs[dm] = make_pair(i, j);
        }
      }

      // Cleanup W candidate, choose pairs with minimum dm if they share decay products
      cleanup(wLepCandIdxs);
      cleanup(wHadCandIdxs, true);
      const size_t nWLepCand = wLepCandIdxs.size();
      const size_t nWHadCand = wHadCandIdxs.size();

      if (nWLepCand + nWHadCand < 2) return; // We skip single top

      int w1Q = 1, w2Q = -1;
      int w1dau1Id = 1, w2dau1Id = -1;
      FourMomentum w1dau1LVec, w1dau2LVec;
      FourMomentum w2dau1LVec, w2dau2LVec;
      if (nWLepCand == 0) { // Full hadronic case
        const pair<size_t, size_t>& idPair1 = wHadCandIdxs.begin()->second;
        const pair<size_t, size_t>& idPair2 = (++wHadCandIdxs.begin())->second;  ///< @todo Reinstate std::next
        const Jet& w1dau1 = _ljets[idPair1.first];
        const Jet& w1dau2 = _ljets[idPair1.second];
        const Jet& w2dau1 = _ljets[idPair2.first];
        const Jet& w2dau2 = _ljets[idPair2.second];
        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
      } else if (nWLepCand == 1) { // Semi-leptonic case
        const pair<size_t, size_t>& idPair1 = wLepCandIdxs.begin()->second;
        const pair<size_t, size_t>& idPair2 = wHadCandIdxs.begin()->second;
        const Jet& w1dau1 = leptons[idPair1.first];
        const Particle& w1dau2 = neutrinos[idPair1.second];
        const Jet& w2dau1 = _ljets[idPair2.first];
        const Jet& w2dau2 = _ljets[idPair2.second];
        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
        w1dau1Id = leptonsId[idPair1.first];
        w1Q = w1dau1Id > 0 ? -1 : 1;
        w2Q = -w1Q;
        switch (w1dau1Id) {
        case 13: case -13: _mode1 = CH_MUON; break;
        case 11: case -11: _mode1 = CH_ELECTRON; break;
        }
      } else { // Full leptonic case
        const pair<size_t, size_t>& idPair1 = wLepCandIdxs.begin()->second;
        const pair<size_t, size_t>& idPair2 = (++wLepCandIdxs.begin())->second;  ///< @todo Reinstate std::next
        const Jet& w1dau1 = leptons[idPair1.first];
        const Particle& w1dau2 = neutrinos[idPair1.second];
        const Jet& w2dau1 = leptons[idPair2.first];
        const Particle& w2dau2 = neutrinos[idPair2.second];
        w1dau1LVec = w1dau1.momentum();
        w1dau2LVec = w1dau2.momentum();
        w2dau1LVec = w2dau1.momentum();
        w2dau2LVec = w2dau2.momentum();
        w1dau1Id = leptonsId[idPair1.first];
        w2dau1Id = leptonsId[idPair2.first];
        w1Q = w1dau1Id > 0 ? -1 : 1;
        w2Q = w2dau1Id > 0 ? -1 : 1;
        switch (w1dau1Id) {
        case 13: case -13: _mode1 = CH_MUON; break;
        case 11: case -11: _mode1 = CH_ELECTRON; break;
        }
        switch (w2dau1Id) {
        case 13: case -13: _mode2 = CH_MUON; break;
        case 11: case -11: _mode2 = CH_ELECTRON; break;
        }
      }
      const FourMomentum w1LVec = w1dau1LVec+w1dau2LVec;
      const FourMomentum w2LVec = w2dau1LVec+w2dau2LVec;

      // Combine b jets
      double sumDm = 1e9;
      FourMomentum b1LVec, b2LVec;
      for (size_t i = 0, n = _bjets.size(); i < n; ++i) {
        const Jet& bjet1 = _bjets[i];
        const double mtop1 = (w1LVec+bjet1.momentum()).mass();
        const double dmtop1 = std::abs(mtop1-_tMass);
        for (size_t j=0; j<n; ++j) {
          if (i == j) continue;
          const Jet& bjet2 = _bjets[j];
          const double mtop2 = (w2LVec+bjet2.momentum()).mass();
          const double dmtop2 = std::abs(mtop2-_tMass);

          if (sumDm <= dmtop1+dmtop2) continue;

          sumDm = dmtop1+dmtop2;
          b1LVec = bjet1.momentum();
          b2LVec = bjet2.momentum();
        }
      }
      if (sumDm >= 1e9) return; // Failed to make top, but this should not happen.

      const FourMomentum t1LVec = w1LVec + b1LVec;
      const FourMomentum t2LVec = w2LVec + b2LVec;

      // Put all of them into candidate collection
      _t1 = Particle(w1Q*6, t1LVec);
      _b1 = Particle(w1Q*5, b1LVec);
      _wDecays1.push_back(Particle(w1dau1Id, w1dau1LVec));
      _wDecays1.push_back(Particle(-w1dau1Id+w1Q, w1dau2LVec));

      _t2 = Particle(w2Q*6, t2LVec);
      _b2 = Particle(w2Q*5, b2LVec);
      _wDecays2.push_back(Particle(w2dau1Id, w2dau1LVec));
      _wDecays2.push_back(Particle(-w2dau1Id+w2Q, w2dau2LVec));

      _isValid = true;
    }

  }


}
