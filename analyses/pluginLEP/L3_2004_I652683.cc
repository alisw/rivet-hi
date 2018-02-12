// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InitialQuarks.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  /// Jet rates and event shapes at LEP I+II
  class L3_2004_I652683 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(L3_2004_I652683);
    // L3_2004_I652683() : Analysis("L3_2004_I652683")
    // {    }


    /// Book histograms and initialise projections before the run
    void init() {

      // Projections to use
      const FinalState FS;
      declare(FS, "FS");
      declare(Beam(), "beams");
      const ChargedFinalState CFS;
      declare(CFS, "CFS");
      const Thrust thrust(FS);
      declare(thrust, "thrust");
      declare(ParisiTensor(FS), "Parisi");
      declare(Hemispheres(thrust), "Hemispheres");
      declare(InitialQuarks(), "initialquarks");

      // Book the histograms
      _h_Thrust_udsc             = bookHisto1D(47, 1, 1);
      _h_Thrust_bottom           = bookHisto1D(47, 1, 2);
      _h_heavyJetmass_udsc       = bookHisto1D(48, 1, 1);
      _h_heavyJetmass_bottom     = bookHisto1D(48, 1, 2);
      _h_totalJetbroad_udsc      = bookHisto1D(49, 1, 1);
      _h_totalJetbroad_bottom    = bookHisto1D(49, 1, 2);
      _h_wideJetbroad_udsc       = bookHisto1D(50, 1, 1);
      _h_wideJetbroad_bottom     = bookHisto1D(50, 1, 2);
      _h_Cparameter_udsc         = bookHisto1D(51, 1, 1);
      _h_Cparameter_bottom       = bookHisto1D(51, 1, 2);
      _h_Dparameter_udsc         = bookHisto1D(52, 1, 1);
      _h_Dparameter_bottom       = bookHisto1D(52, 1, 2);
      _h_Ncharged                = bookHisto1D(59, 1, 1);
      _h_Ncharged_udsc           = bookHisto1D(59, 1, 2);
      _h_Ncharged_bottom         = bookHisto1D(59, 1, 3);
      _h_scaledMomentum          = bookHisto1D(65, 1, 1);
      _h_scaledMomentum_udsc     = bookHisto1D(65, 1, 2);
      _h_scaledMomentum_bottom   = bookHisto1D(65, 1, 3);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beam average momentum
      const ParticlePair& beams = apply<Beam>(event, "beams").beams();
      const double beamMomentum = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;

      // InitialQuarks projection to have udsc events separated from b events
      /// @todo Yuck!!! Eliminate when possible...
      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "initialquarks");
      Particles quarks;
      if ( iqf.particles().size() == 2 ) {
        flavour = iqf.particles().front().abspid();
        quarks  = iqf.particles();
      } else {
        map<int, Particle> quarkmap;
        for (const Particle& p : iqf.particles()) {
          if (quarkmap.find(p.pid()) == quarkmap.end()) quarkmap[p.pid()] = p;
          else if (quarkmap[p.pid()].E() < p.E()) quarkmap[p.pid()] = p;
        }
        double max_energy = 0.;
        for (int i = 1; i <= 5; ++i) {
          double energy = 0.;
          if (quarkmap.find(i) != quarkmap.end())
            energy += quarkmap[ i].E();
          if (quarkmap.find(-i) != quarkmap.end())
            energy += quarkmap[-i].E();
          if (energy > max_energy)
            flavour = i;
        }
        if (quarkmap.find(flavour) != quarkmap.end())
          quarks.push_back(quarkmap[flavour]);
        if (quarkmap.find(-flavour) != quarkmap.end())
          quarks.push_back(quarkmap[-flavour]);
      }

      // Flavour label
      /// @todo Change to a bool?
      const int iflav = (flavour == PID::DQUARK || flavour == PID::UQUARK || flavour == PID::SQUARK || flavour == PID::CQUARK) ? 1 : (flavour == PID::BQUARK) ? 5 : 0;

      // Update weight sums
      const double weight = event.weight();
      if (iflav == 1) {
        _sumW_udsc += weight;
      } else if (iflav == 5) {
        _sumW_b += weight;
      }
      _sumW_ch += weight;

      // Charged multiplicity
      const FinalState& cfs = applyProjection<FinalState>(event, "CFS");
      _h_Ncharged->fill(cfs.size(), weight);
      if (iflav == 1) {
        _sumW_ch_udsc += weight;
        _h_Ncharged_udsc->fill(cfs.size(), weight);
      } else if (iflav == 5) {
        _sumW_ch_b += weight;
        _h_Ncharged_bottom->fill(cfs.size(), weight);
      }

      // Scaled momentum
      const Particles& chparticles = cfs.particlesByPt();
      for (const Particle& p : chparticles) {
        const Vector3 momentum3 = p.p3();
        const double mom = momentum3.mod();
        const double scaledMom = mom/beamMomentum;
        const double logScaledMom = std::log(scaledMom);
        _h_scaledMomentum->fill(-logScaledMom, weight);
        if (iflav == 1) {
          _h_scaledMomentum_udsc->fill(-logScaledMom, weight);
        } else if (iflav == 5) {
          _h_scaledMomentum_bottom->fill(-logScaledMom, weight);
        }
      }

      // Thrust
      const Thrust& thrust = applyProjection<Thrust>(event, "thrust");
      if (iflav == 1) {
        _h_Thrust_udsc->fill(thrust.thrust(), weight);
      } else if (iflav == 5) {
        _h_Thrust_bottom->fill(thrust.thrust(), weight);
      }

      // C and D Parisi parameters
      const ParisiTensor& parisi = applyProjection<ParisiTensor>(event, "Parisi");
      if (iflav == 1) {
        _h_Cparameter_udsc->fill(parisi.C(), weight);
        _h_Dparameter_udsc->fill(parisi.D(), weight);
      } else if (iflav == 5) {
        _h_Cparameter_bottom->fill(parisi.C(), weight);
        _h_Dparameter_bottom->fill(parisi.D(), weight);
      }

      // The hemisphere variables
      const Hemispheres& hemisphere = applyProjection<Hemispheres>(event, "Hemispheres");
      if (iflav == 1) {
        _h_heavyJetmass_udsc->fill(hemisphere.scaledM2high(), weight);
        _h_totalJetbroad_udsc->fill(hemisphere.Bsum(), weight);
        _h_wideJetbroad_udsc->fill(hemisphere.Bmax(), weight);
      } else if (iflav == 5) {
        _h_heavyJetmass_bottom->fill(hemisphere.scaledM2high(), weight);
        _h_totalJetbroad_bottom->fill(hemisphere.Bsum(), weight);
        _h_wideJetbroad_bottom->fill(hemisphere.Bmax(), weight);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale({_h_Thrust_udsc, _h_heavyJetmass_udsc, _h_totalJetbroad_udsc, _h_wideJetbroad_udsc, _h_Cparameter_udsc, _h_Dparameter_udsc}, 1/_sumW_udsc);
      scale({_h_Thrust_bottom, _h_heavyJetmass_bottom, _h_totalJetbroad_bottom, _h_wideJetbroad_bottom, _h_Cparameter_bottom, _h_Dparameter_bottom}, 1./_sumW_b);
      scale(_h_Ncharged, 2/_sumW_ch);
      scale(_h_Ncharged_udsc, 2/_sumW_ch_udsc);
      scale(_h_Ncharged_bottom, 2/_sumW_ch_b);
      scale(_h_scaledMomentum, 1/_sumW_ch);
      scale(_h_scaledMomentum_udsc, 1/_sumW_ch_udsc);
      scale(_h_scaledMomentum_bottom, 1/_sumW_ch_b);
    }


    /// Weight counters
    double _sumW_udsc = 0, _sumW_b = 0, _sumW_ch = 0, _sumW_ch_udsc = 0, _sumW_ch_b = 0;

    /// @name Histograms
    //@{
    Histo1DPtr _h_Thrust_udsc, _h_Thrust_bottom;
    Histo1DPtr _h_heavyJetmass_udsc, _h_heavyJetmass_bottom;
    Histo1DPtr _h_totalJetbroad_udsc, _h_totalJetbroad_bottom;
    Histo1DPtr _h_wideJetbroad_udsc, _h_wideJetbroad_bottom;
    Histo1DPtr _h_Cparameter_udsc, _h_Cparameter_bottom;
    Histo1DPtr _h_Dparameter_udsc, _h_Dparameter_bottom;
    Histo1DPtr _h_Ncharged, _h_Ncharged_udsc, _h_Ncharged_bottom;
    Histo1DPtr _h_scaledMomentum, _h_scaledMomentum_udsc, _h_scaledMomentum_bottom;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(L3_2004_I652683);

}
