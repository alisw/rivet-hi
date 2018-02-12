// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Event.hh"

namespace Rivet {


  class EHS_1988_I265504 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(EHS_1988_I265504);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(ChargedFinalState(), "CFS");
      declare(Beam(),"Beam");

      switch ( beamIds().first ) {
      case PID::PIPLUS:
        _h_cpos_xF = bookHisto1D(1, 1, 1);
        _h_cpos_eta = bookHisto1D(3, 1, 1);
        _h_cpos_pT2 = bookHisto1D(5, 1, 1);
        _h_cneg_xF = bookHisto1D(2, 1, 1);
        _h_cneg_eta = bookHisto1D(4, 1, 1);
        _h_cneg_pT2 = bookHisto1D(6, 1, 1);
        break;

      case PID::KPLUS:
        _h_cpos_xF = bookHisto1D(1, 1, 2);
        _h_cpos_eta = bookHisto1D(3, 1, 2);
        _h_cpos_pT2 = bookHisto1D(5, 1, 2);
        _h_cneg_xF = bookHisto1D(2, 1, 2);
        _h_cneg_eta = bookHisto1D(4, 1, 2);
        _h_cneg_pT2 = bookHisto1D(6, 1, 2);
        break;

      case PID::PROTON:
        _h_cpos_xF = bookHisto1D(1, 1, 3);
        _h_cpos_eta = bookHisto1D(3, 1, 3);
        _h_cpos_pT2 = bookHisto1D(5, 1, 3);
        _h_cneg_xF = bookHisto1D(2, 1, 3);
        _h_cneg_eta = bookHisto1D(4, 1, 3);
        _h_cneg_pT2 = bookHisto1D(6, 1, 3);
        break;
      }

      // Calculate boost from lab to CM frame
      _beamboost = cmsTransform( beams() );
      MSG_DEBUG("Boost vector: " << _beamboost );

      // Transform beam into CMS frame
      Particle _beam_cm = beams().first;
      _beam_cm.transformBy(_beamboost);
      // Beam momentum in CM frame defines Feynman-x
      _pz_max = _beam_cm.pz();

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const FinalState& fs = apply<FinalState>(event, "CFS");
      for (const Particle& p: fs.particles()) {
        // Only interested in pi- or positively charged
        if (p.charge() < 0 && p.pid() != PID::PIMINUS) continue;
        // Slow proton cut: reject lab momenta < 1.2GeV
        if (p.pid() == PID::PROTON && p.p() < 1.2*GeV) continue;
        // Transform to cm frame
        const FourMomentum pcm = _beamboost.transform(p);
        const double xF = pcm.pz()/_pz_max;

        if (p.charge() > 0) {
          _h_cpos_xF->fill( xF, weight);
          _h_cpos_pT2->fill( p.pT2(), weight);
          _h_cpos_eta->fill( p.eta(), weight);
        } else if (p.pid() == PID::PIMINUS) {
          _h_cneg_xF->fill( xF, weight);
          _h_cneg_pT2->fill( p.pT2(), weight);
          _h_cneg_eta->fill( p.eta(), weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale( {_h_cpos_xF, _h_cpos_pT2,_h_cpos_eta, _h_cneg_xF, _h_cneg_eta, _h_cneg_pT2},
             crossSection()/millibarn/sumOfWeights() );
    }

    //@}


    /// @name Histograms
    //@{
    LorentzTransform _beamboost;
    double _pz_max;
    Histo1DPtr _h_cpos_xF, _h_cpos_eta, _h_cpos_pT2;
    Histo1DPtr _h_cneg_xF, _h_cneg_eta, _h_cneg_pT2;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(EHS_1988_I265504);

}
