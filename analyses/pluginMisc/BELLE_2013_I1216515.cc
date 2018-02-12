// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief BELLE pion and kaon continuum production
  /// @author Peter Richardson
  class BELLE_2013_I1216515 : public Analysis {
  public:

    BELLE_2013_I1216515()
      : Analysis("BELLE_2013_I1216515")
    { }


    void analyze(const Event& e) {
      const double weight = e.weight();

      // Loop through charged FS particles and look for charmed mesons/baryons
      const ChargedFinalState& fs = apply<ChargedFinalState>(e, "FS");

      const Beam beamproj = apply<Beam>(e, "Beams");
      const ParticlePair& beams = beamproj.beams();
      const FourMomentum mom_tot = beams.first.momentum() + beams.second.momentum();
      const LorentzTransform cms_boost = LorentzTransform::mkFrameTransformFromBeta(mom_tot.betaVec());
      MSG_DEBUG("CMS energy sqrt s = " << beamproj.sqrtS());

      foreach (const Particle& p, fs.particles()) {
        // energy in CMS frame
        const double en = cms_boost.transform(p.momentum()).t();
	const double z = 2.*en/beamproj.sqrtS();
        const int PdgId = p.abspid();
        MSG_DEBUG("pdgID = " << PdgId << "  Energy = " << en);
        switch (PdgId) {
	case PID::PIPLUS:
	  _histPion->fill(z,weight);
	  break;
	case PID::KPLUS:
	  _histKaon->fill(z,weight);
	  break;
	default :
	  break;
        }
      }
    } // analyze


    void finalize() {

      scale(_histPion,crossSection()/femtobarn/sumOfWeights());
      scale(_histKaon,crossSection()/femtobarn/sumOfWeights());
    } // finalize


    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      _histPion = bookHisto1D(1,1,1);
      _histKaon = bookHisto1D(1,1,2);

    } // init

  private:

    //@{
    // Histograms for continuum data (sqrt(s) = 10.52 GeV)
    Histo1DPtr _histPion;
    Histo1DPtr _histKaon;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BELLE_2013_I1216515);

}
