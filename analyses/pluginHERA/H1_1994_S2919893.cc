// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Math/Constants.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"

namespace Rivet {


  /// @brief H1 energy flow and charged particle spectra
  /// @author Peter Richardson
  /// Based on the equivalent HZTool analysis
  class H1_1994_S2919893 : public Analysis {
  public:

    /// Constructor
    H1_1994_S2919893()
      : Analysis("H1_1994_S2919893")
    {

      // Initialise member variables
      _w77  = make_pair(0.0, 0.0);
      _w122 = make_pair(0.0, 0.0);
      _w169 = make_pair(0.0, 0.0);
      _w117 = make_pair(0.0, 0.0);
      _wEnergy = make_pair(0.0, 0.0);
    }


    /// @name Analysis methods
    //@{

    /// Initialise projections and histograms
    void init() {
      // Projections
      declare(DISLepton(), "Lepton");
      declare(DISKinematics(), "Kinematics");
      declare(FinalState(), "FS");

      // Histos
      _histEnergyFlowLowX =  bookHisto1D(1, 1, 1);
      _histEnergyFlowHighX = bookHisto1D(1, 1, 2);

      _histEECLowX = bookHisto1D(2, 1, 1);
      _histEECHighX = bookHisto1D(2, 1, 2);

      _histSpectraW77 = bookHisto1D(3, 1, 1);
      _histSpectraW122 = bookHisto1D(3, 1, 2);
      _histSpectraW169 = bookHisto1D(3, 1, 3);
      _histSpectraW117 = bookHisto1D(3, 1, 4);

      _histPT2 = bookProfile1D(4, 1, 1);
    }


    /// Analyse each event
    void analyze(const Event& event) {

      // Get the DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      const double x  = dk.x();
      const double w2 = dk.W2();
      const double w = sqrt(w2);

      // Momentum of the scattered lepton
      const DISLepton& dl = apply<DISLepton>(event,"Lepton");
      const FourMomentum leptonMom = dl.out();
      const double ptel = leptonMom.pT();
      const double enel = leptonMom.E();
      const double thel = leptonMom.angle(dk.beamHadron().mom())/degree;

      // Extract the particles other than the lepton
      const FinalState& fs = apply<FinalState>(event, "FS");
      Particles particles;
      particles.reserve(fs.particles().size());
      const GenParticle* dislepGP = dl.out().genParticle();
      foreach (const Particle& p, fs.particles()) {
        const GenParticle* loopGP = p.genParticle();
        if (loopGP == dislepGP) continue;
        particles.push_back(p);
      }

      // Cut on the forward energy
      double efwd = 0.0;
      foreach (const Particle& p, particles) {
        const double th = p.angle(dk.beamHadron())/degree;
        if (inRange(th, 4.4, 15)) efwd += p.E();
      }

      // Apply the cuts
      // Lepton energy and angle, w2 and forward energy
      MSG_DEBUG("enel/GeV = " << enel/GeV << ", thel = " << thel
                << ", w2 = " << w2 << ", efwd/GeV = " << efwd/GeV);
      bool cut = enel/GeV > 14. && thel > 157. && thel < 172.5 && w2 >= 3000. && efwd/GeV > 0.5;
      if (!cut) vetoEvent;

      // Weight of the event
      const double weight = event.weight();
      (x < 1e-3 ? _wEnergy.first : _wEnergy.second) += weight;

      // Boost to hadronic CM
      const LorentzTransform hcmboost = dk.boostHCM();
      // Loop over the particles
      long ncharged(0);
      for (size_t ip1 = 0; ip1 < particles.size(); ++ip1) {
        const Particle& p = particles[ip1];

        const double th = p.angle(dk.beamHadron().momentum()) / degree;
        // Boost momentum to lab
        const FourMomentum hcmMom = hcmboost.transform(p.momentum());
        // Angular cut
        if (th <= 4.4) continue;

        // Energy flow histogram
        const double et = fabs(hcmMom.Et());
        const double eta = hcmMom.eta();
        (x < 1e-3 ? _histEnergyFlowLowX : _histEnergyFlowHighX)->fill(eta, et*weight);
        if (PID::threeCharge(p.pid()) != 0) {
          /// @todo Use units in w comparisons... what are the units?
          if (w > 50. && w <= 200.) {
            double xf= 2 * hcmMom.z() / w;
            double pt2 = hcmMom.pT2();
            if (w > 50. && w <= 100.) {
              _histSpectraW77 ->fill(xf, weight);
            } else if (w > 100. && w <= 150.) {
              _histSpectraW122->fill(xf, weight);
            } else if (w > 150. && w <= 200.) {
              _histSpectraW169->fill(xf, weight);
            }
            _histSpectraW117->fill(xf, weight);
            /// @todo Is this profile meant to be filled with 2 weight factors?
            _histPT2->fill(xf, pt2*weight/GeV2, weight);
            ++ncharged;
          }
        }


        // Energy-energy correlation
        if (th <= 8.) continue;
        double phi1 = p.phi(ZERO_2PI);
        double eta1 = p.eta();
        double et1 = fabs(p.momentum().Et());
        for (size_t ip2 = ip1+1; ip2 < particles.size(); ++ip2) {
          const Particle& p2 = particles[ip2];

          //double th2 = beamAngle(p2.momentum(), order);
          double th2 = p2.angle(dk.beamHadron().momentum()) / degree;
          if (th2 <= 8.) continue;
          double phi2 = p2.phi(ZERO_2PI);

          /// @todo Use angle function
          double deltaphi = phi1 - phi2;
          if (fabs(deltaphi) > PI) deltaphi = fabs(fabs(deltaphi) - TWOPI);
          double eta2 = p2.eta();
          double omega = sqrt(sqr(eta1-eta2) + sqr(deltaphi));
          double et2 = fabs(p2.momentum().Et());
          double wt = et1*et2 / sqr(ptel) * weight;
          (x < 1e-3 ? _histEECLowX : _histEECHighX)->fill(omega, wt);
        }
      }

      // Factors for normalization
      if (w > 50. && w <= 200.) {
        if (w <= 100.) {
          _w77.first  += ncharged*weight;
          _w77.second += weight;
        } else if (w <= 150.) {
          _w122.first  += ncharged*weight;
          _w122.second += weight;
        } else {
          _w169.first  += ncharged*weight;
          _w169.second += weight;
        }
        _w117.first  += ncharged*weight;
        _w117.second += weight;
      }
    }


    // Normalize inclusive single particle distributions to the average number of charged particles per event.
    void finalize() {
      normalize(_histSpectraW77, _w77.first/_w77.second);
      normalize(_histSpectraW122, _w122.first/_w122.second);
      normalize(_histSpectraW169, _w169.first/_w169.second);
      normalize(_histSpectraW117, _w117.first/_w117.second);

      scale(_histEnergyFlowLowX , 1./_wEnergy.first );
      scale(_histEnergyFlowHighX, 1./_wEnergy.second);

      scale(_histEECLowX , 1./_wEnergy.first );
      scale(_histEECHighX, 1./_wEnergy.second);
    }

    //@}


  private:

    /// Polar angle with right direction of the beam
    inline double beamAngle(const FourVector& v, bool order) {
      double thel = v.polarAngle()/degree;
      if (thel < 0) thel += 180.;
      if (!order) thel = 180 - thel;
      return thel;
    }

    /// @name Histograms
    //@{
    Histo1DPtr _histEnergyFlowLowX, _histEnergyFlowHighX;
    Histo1DPtr _histEECLowX, _histEECHighX;
    Histo1DPtr _histSpectraW77, _histSpectraW122, _histSpectraW169, _histSpectraW117;
    Profile1DPtr _histPT2;
    //@}

    /// @name Storage of weights to calculate averages for normalisation
    //@{
    pair<double,double> _w77, _w122, _w169, _w117, _wEnergy;
    //@}

  };



  DECLARE_RIVET_PLUGIN(H1_1994_S2919893);

}
