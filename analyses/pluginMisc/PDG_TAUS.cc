// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/TauFinder.hh"

namespace Rivet {


  class PDG_TAUS : public Analysis {
  public:

    /// Constructor
    PDG_TAUS()
      : Analysis("PDG_TAUS"),
        _weights_had(0),
        _weights_mu(0),
        _weights_el(0)
    {   }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      TauFinder tauleptonic(TauFinder::LEPTONIC); // open cuts, leptonic decays
      declare(tauleptonic, "TauLeptonic");

      TauFinder tauhadronic(TauFinder::HADRONIC); // open cuts, hadronic decays
      declare(tauhadronic, "TauHadronic");

      populateDecayMap();

      _h_ratio_mu        = bookHisto1D(1, 1, 1);
      _h_ratio_el        = bookHisto1D(1, 1, 2);
      _h_1prong_pinu     = bookHisto1D(2, 1, 1);
      _h_1prong_Kpnu     = bookHisto1D(2, 1, 2);
      _h_1prong_pipinu   = bookHisto1D(2, 1, 3);
      _h_1prong_Kppinu   = bookHisto1D(2, 1, 4);
      _h_1prong_pipipinu = bookHisto1D(2, 1, 5);
      _h_1prong_Knpinu   = bookHisto1D(2, 1, 6);
      _h_3prong_pipipinu = bookHisto1D(2, 2, 1);
      _h_5prong          = bookHisto1D(2, 3, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& e) {
      const double weight = e.weight();

      const TauFinder& taulep = apply<TauFinder>(e, "TauLeptonic");
      const TauFinder& tauhad = apply<TauFinder>(e, "TauHadronic");

      // Hadronic tau decays --- prong decays
      foreach(const Particle& tau, tauhad.taus()) {
        _weights_had += weight;
        int prongs = countProngs(tau); // number of charged particles among decay products
        // Only do 1 prong decays here
        if (prongs == 1) {
          ////// Exclusive decay modes "1-prong"
          if (analyzeDecay(tau,   decay_pids["pinu"], true))     _h_1prong_pinu->fill(1, weight);
          if (analyzeDecay(tau,   decay_pids["Kpnu"], true))     _h_1prong_Kpnu->fill(1, weight);
          if (analyzeDecay(tau, decay_pids["pipinu"], true))     _h_1prong_pipinu->fill(1, weight);
          if (analyzeDecay(tau, decay_pids["Kppinu"]  , true))   _h_1prong_Kppinu->fill(1, weight);
          if (analyzeDecay(tau, decay_pids["pipipinu"], true))   _h_1prong_pipipinu->fill(1, weight);
          // Kshort, Klong --- (twice) filling the K0 labelled PDG histo
          if (analyzeDecay(tau, decay_pids["KSpinu"]  , true))   _h_1prong_Knpinu->fill(1, weight);
          if (analyzeDecay(tau, decay_pids["KLpinu"]  , true))   _h_1prong_Knpinu->fill(1, weight);
        }
        else if (prongs == 3) {
          if (analyzeDecay(tau, decay_pids["3pipipinu"], true))  _h_3prong_pipipinu->fill(1, weight);
        }
        else if (prongs == 5 && !any(tau.children(), HasAbsPID(310))) _h_5prong->fill(1, weight);
      }

      // Leptonic tau decays --- look for radiative and non-radiative 1 prong decays
      foreach(const Particle& tau, taulep.taus()) {
        int prongs = countProngs(tau); // number of charged particles among decay products
        // Only do 1 prong decays here
        if (prongs == 1) {
          analyzeRadiativeDecay(tau, decay_pids["muids"], _weights_mu,  weight, true, _h_ratio_mu);
          analyzeRadiativeDecay(tau, decay_pids["elids"], _weights_el,  weight, true, _h_ratio_el);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_ratio_mu, 1./_weights_mu);
      scale(_h_ratio_el, 1./_weights_el);

      const double norm = _weights_had + _weights_mu + _weights_el;
      scale(_h_1prong_pinu,     1./norm);
      scale(_h_1prong_Kpnu,     1./norm);
      scale(_h_1prong_pipinu,   1./norm);
      scale(_h_1prong_Kppinu,   1./norm);
      scale(_h_1prong_pipipinu, 1./norm);
      scale(_h_1prong_Knpinu,   1./norm);
      scale(_h_3prong_pipipinu, 1./norm);
      scale(_h_5prong,          1./norm);
    }


    // Short hand
    bool contains(Particle& mother, int id, bool abs=false) {
      if (abs) return any(mother.children(), HasAbsPID(id));
      return any(mother.children(), HasPID(id));
    }


    // Count charged decay products
    int countProngs(Particle mother) {
      int n_prongs = 0;
      foreach(Particle p, mother.children())
        if (p.threeCharge()!=0) ++n_prongs;
      return n_prongs;
    }


    // Set up a lookup table for decays
    void populateDecayMap() {
      decay_pids["muids"]     = {{ 13,14,16 }};
      decay_pids["elids"]     = {{ 11,12,16 }};
      decay_pids["pinu"]      = {{ 211,16 }};
      decay_pids["Kpnu"]      = {{ 321,16 }};
      decay_pids["pipinu"]    = {{ 111,211,16 }};
      decay_pids["Kppinu"]    = {{ 111,321,16 }};
      decay_pids["pipipinu"]  = {{ 111,111,211,16 }};
      decay_pids["KSpinu"]    = {{ 211,310,16 }};
      decay_pids["KLpinu"]    = {{ 211,130,16 }};
      decay_pids["3pipipinu"] = {{ 211,211,211,16 }};
    }


    bool analyzeDecay(Particle mother, vector<int> ids, bool absolute) {
      // There is no point in looking for decays with less particles than to be analysed
      if (mother.children().size() == ids.size()) {
        bool decayfound = true;
        foreach (int id, ids) {
          if (!contains(mother, id, absolute)) decayfound = false;
        }
        return decayfound;
      } // end of first if
      return false;
    }


    // Look for radiative (and non-radiative) tau decays to fill a ratio histo
    void analyzeRadiativeDecay(Particle mother, vector<int> ids, double &w_incl, double e_weight, bool absolute, Histo1DPtr h_ratio) {
      // w_incl   ... reference to a global weight counter for all leptonic tau decays
      // e_weight ... the current event weight
      // h_ratio  ... pointer to ratio histo --- filled with e_weight in case of radiative events only

      // There is no point in looking for decays with less particles than to be analysed
      if (mother.children().size() >= ids.size()) {
        bool decayfound = true;
        foreach (int id, ids) {
          if (!contains(mother, id, absolute)) decayfound = false;
        }
        // Do not increment counters if the specified decay products were not found
        if (decayfound) {
          w_incl += e_weight; // the (global) weight counter for leptonic decays
          bool radiative = any(mother.children(), HasPID(PID::PHOTON));

          // Only fill the histo if there is a radiative decay
          if (radiative) {
            // Iterate over decay products to find photon with 5 MeV energy
            foreach (const Particle& son, mother.children()) {
              if (son.pid() == PID::PHOTON) {
                // Require photons to have at least 5 MeV energy in the rest frame of the tau
                // boosted taus
                if (!mother.momentum().betaVec().isZero()) {
                  LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(mother.momentum().betaVec());
                  if (cms_boost.transform(son.momentum())[0]/MeV > 5.) {
                    h_ratio->fill(1, e_weight);
                    break;
                  }
                }
                // not boosted taus
                else {
                  if (son.momentum()[0]/MeV > 5.) {
                    h_ratio->fill(1, e_weight);
                    break;
                  }
                }
              }
            } // end loop over decay products
          } // end of radiative
        } // end of decayfound
      } // end of first if
    }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_ratio_mu, _h_ratio_el;
    Histo1DPtr _h_1prong_pinu, _h_1prong_Kpnu, _h_1prong_Kppinu, _h_1prong_pipinu, _h_1prong_pipipinu, _h_1prong_Knpinu;
    Histo1DPtr _h_3prong_pipipinu;
    Histo1DPtr _h_5prong;
    //@}

    double _weights_had, _weights_mu, _weights_el;
    map<string, vector<int> > decay_pids;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(PDG_TAUS);

}
