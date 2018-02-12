// -*- C++ -*-
#ifndef RIVET_SmearingFunctions_HH
#define RIVET_SmearingFunctions_HH

#include "Rivet/Tools/MomentumSmearingFunctions.hh"
#include "Rivet/Tools/ParticleSmearingFunctions.hh"
#include "Rivet/Tools/JetSmearingFunctions.hh"

namespace Rivet {


  /// @name Electron efficiency and smearing functions
  //@{

  /// ATLAS Run 1 electron reconstruction efficiency
  /// @todo Include reco eff (but no e/y discrimination) in forward region
  /// @todo How to use this in combination with tracking eff?
  inline double ELECTRON_EFF_ATLAS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 10*GeV) return 0;
    return (e.abseta() < 1.5) ? 0.95 : 0.85;
  }

  /// ATLAS Run 2 electron reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_EFF_ATLAS_RUN2(const Particle& e) {
    return ELECTRON_EFF_ATLAS_RUN1(e);
  }


  /// @brief ATLAS Run 2 'loose' electron identification/selection efficiency
  ///
  /// Values read from Fig 3 of ATL-PHYS-PUB-2015-041
  /// @todo What about faking by jets or non-electrons?
  inline double ELECTRON_IDEFF_ATLAS_RUN2_LOOSE(const Particle& e) {

    // Manually symmetrised eta eff histogram
    const static vector<double> edges_eta = { 0.0,   0.1,   0.8,   1.37,  1.52,  2.01,  2.37,  2.47 };
    const static vector<double> effs_eta  = { 0.950, 0.965, 0.955, 0.885, 0.950, 0.935, 0.90 };
    // Et eff histogram (10-20 is a guess)
    const static vector<double> edges_et = { 0,   10,   20,   25,   30,   35,   40,    45,    50,   60,  80 };
    const static vector<double> effs_et  = { 0.0, 0.90, 0.91, 0.92, 0.94, 0.95, 0.955, 0.965, 0.97, 0.98 };

    if (e.abseta() > 2.47) return 0.0; // no ID outside the tracker

    const int i_eta = binIndex(e.abseta(), edges_eta);
    const int i_et = binIndex(e.Et()/GeV, edges_et, true);
    const double eff = effs_et[i_et] * effs_eta[i_eta] / 0.95; //< norm factor as approximate double differential
    return min(eff, 1.0);
  }


  /// @brief ATLAS Run 1 'medium' electron identification/selection efficiency
  inline double ELECTRON_IDEFF_ATLAS_RUN1_MEDIUM(const Particle& e) {

    const static vector<double> eta_edges_10 = {0.000, 0.049, 0.454, 1.107, 1.46, 1.790, 2.277, 2.500};
    const static vector<double> eta_vals_10  = {0.730, 0.757, 0.780, 0.771, 0.77, 0.777, 0.778};

    const static vector<double> eta_edges_15 = {0.000, 0.053, 0.456, 1.102, 1.463, 1.783, 2.263, 2.500};
    const static vector<double> eta_vals_15  = {0.780, 0.800, 0.819, 0.759, 0.749, 0.813, 0.829};

    const static vector<double> eta_edges_20 = {0.000, 0.065, 0.362, 0.719, 0.980, 1.289, 1.455, 1.681, 1.942, 2.239, 2.452, 2.500};
    const static vector<double> eta_vals_20  = {0.794, 0.806, 0.816, 0.806, 0.797, 0.774, 0.764, 0.788, 0.793, 0.806, 0.825};

    const static vector<double> eta_edges_25 = {0.000, 0.077, 0.338, 0.742, 1.004, 1.265, 1.467, 1.692, 1.940, 2.227, 2.452, 2.500};
    const static vector<double> eta_vals_25  = {0.833, 0.843, 0.853, 0.845, 0.839, 0.804, 0.790, 0.825, 0.830, 0.833, 0.839};

    const static vector<double> eta_edges_30 = {0.000, 0.077, 0.350, 0.707, 0.980, 1.289, 1.479, 1.681, 1.942, 2.239, 2.441, 2.500};
    const static vector<double> eta_vals_30  = {0.863, 0.872, 0.881, 0.874, 0.870, 0.824, 0.808, 0.847, 0.845, 0.840, 0.842};

    const static vector<double> eta_edges_35 = {0.000, 0.058, 0.344, 0.700, 1.009, 1.270, 1.458, 1.685, 1.935, 2.231, 2.468, 2.500};
    const static vector<double> eta_vals_35  = {0.878, 0.889, 0.901, 0.895, 0.893, 0.849, 0.835, 0.868, 0.863, 0.845, 0.832};

    const static vector<double> eta_edges_40 = {0.000, 0.047, 0.355, 0.699, 0.983, 1.280, 1.446, 1.694, 1.943, 2.227, 2.441, 2.500};
    const static vector<double> eta_vals_40  = {0.894, 0.901, 0.909, 0.905, 0.904, 0.875, 0.868, 0.889, 0.876, 0.848, 0.827};

    const static vector<double> eta_edges_45 = {0.000, 0.058, 0.356, 0.712, 0.997, 1.282, 1.459, 1.686, 1.935, 2.220, 2.444, 2.500};
    const static vector<double> eta_vals_45  = {0.900, 0.911, 0.923, 0.918, 0.917, 0.897, 0.891, 0.904, 0.894, 0.843, 0.796};

    const static vector<double> eta_edges_50 = {0.000, 0.059, 0.355, 0.711, 0.983, 1.280, 1.469, 1.682, 1.919, 2.227, 2.441, 2.500};
    const static vector<double> eta_vals_50  = {0.903, 0.913, 0.923, 0.922, 0.923, 0.903, 0.898, 0.908, 0.895, 0.831, 0.774};

    const static vector<double> eta_edges_60 = {0.000, 0.053, 0.351, 0.720, 1.006, 1.291, 1.469, 1.696, 1.946, 2.243, 2.455, 2.500};
    const static vector<double> eta_vals_60  = {0.903, 0.917, 0.928, 0.924, 0.927, 0.915, 0.911, 0.915, 0.899, 0.827, 0.760};

    const static vector<double> eta_edges_80 = {0.000, 0.053, 0.351, 0.720, 0.994, 1.292, 1.482, 1.708, 1.934, 2.220, 2.458, 2.500};
    const static vector<double> eta_vals_80  = {0.936, 0.942, 0.952, 0.956, 0.956, 0.934, 0.931, 0.944, 0.933, 0.940, 0.948};

    const static vector<double> et_edges = { 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80 };
    const static vector< vector<double> > et_eta_edges = { eta_edges_10, eta_edges_15, eta_edges_20, eta_edges_25, eta_edges_30, eta_edges_35, eta_edges_40, eta_edges_45, eta_edges_50, eta_edges_60, eta_edges_80 };
    const static vector< vector<double> > et_eta_vals  = { eta_vals_10, eta_vals_15, eta_vals_20, eta_vals_25, eta_vals_30, eta_vals_35, eta_vals_40, eta_vals_45, eta_vals_50, eta_vals_60, eta_vals_80 };

    if (e.abseta() > 2.5 || e.Et() < 10*GeV) return 0.0;
    const int i_et = binIndex(e.Et()/GeV, et_edges, true);
    const int i_eta = binIndex(e.abseta(), et_eta_edges[i_et]);
    return et_eta_vals[i_et][i_eta];
  }

  /// @brief ATLAS Run 2 'medium' electron identification/selection efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_IDEFF_ATLAS_RUN2_MEDIUM(const Particle& e) {
    return ELECTRON_IDEFF_ATLAS_RUN1_MEDIUM(e);
  }


  /// @brief ATLAS Run 1 'tight' electron identification/selection efficiency
  inline double ELECTRON_IDEFF_ATLAS_RUN1_TIGHT(const Particle& e) {

    const static vector<double> eta_edges_10 = {0.000, 0.049, 0.459, 1.100, 1.461, 1.789, 2.270, 2.500};
    const static vector<double> eta_vals_10  = {0.581, 0.632, 0.668, 0.558, 0.548, 0.662, 0.690};

    const static vector<double> eta_edges_15 = {0.000, 0.053, 0.450, 1.096, 1.463, 1.783, 2.269, 2.500};
    const static vector<double> eta_vals_15 =  {0.630, 0.678, 0.714, 0.633, 0.616, 0.700, 0.733};

    const static vector<double> eta_edges_20 = {0.000, 0.065, 0.362, 0.719, 0.992, 1.277, 1.479, 1.692, 1.930, 2.227, 2.464, 2.500};
    const static vector<double> eta_vals_20 =  {0.653, 0.695, 0.735, 0.714, 0.688, 0.635, 0.625, 0.655, 0.680, 0.691, 0.674};

    const static vector<double> eta_edges_25 = {0.000, 0.077, 0.362, 0.719, 0.992, 1.300, 1.479, 1.692, 1.942, 2.227, 2.464, 2.500};
    const static vector<double> eta_vals_25 =  {0.692, 0.732, 0.768, 0.750, 0.726, 0.677, 0.667, 0.692, 0.710, 0.706, 0.679};

    const static vector<double> eta_edges_30 = {0.000, 0.053, 0.362, 0.719, 1.004, 1.277, 1.467, 1.681, 1.954, 2.239, 2.452, 2.500};
    const static vector<double> eta_vals_30 =  {0.724, 0.763, 0.804, 0.789, 0.762, 0.702, 0.690, 0.720, 0.731, 0.714, 0.681};

    const static vector<double> eta_edges_35 = {0.000, 0.044, 0.342, 0.711, 0.971, 1.280, 1.456, 1.683, 1.944, 2.218, 2.442, 2.500};
    const static vector<double> eta_vals_35 =  {0.736, 0.778, 0.824, 0.811, 0.784, 0.730, 0.718, 0.739, 0.743, 0.718, 0.678};

    const static vector<double> eta_edges_40 = {0.000, 0.047, 0.355, 0.699, 0.983, 1.268, 1.457, 1.671, 1.931, 2.204, 2.453, 2.500};
    const static vector<double> eta_vals_40 =  {0.741, 0.774, 0.823, 0.823, 0.802, 0.764, 0.756, 0.771, 0.771, 0.734, 0.684};

    const static vector<double> eta_edges_45 = {0.000, 0.056, 0.354, 0.711, 0.984, 1.280, 1.458, 1.684, 1.945, 2.207, 2.442, 2.500};
    const static vector<double> eta_vals_45 =  {0.758, 0.792, 0.841, 0.841, 0.823, 0.792, 0.786, 0.796, 0.794, 0.734, 0.663};

    const static vector<double> eta_edges_50 = {0.000, 0.059, 0.355, 0.699, 0.983, 1.268, 1.446, 1.682, 1.943, 2.216, 2.453, 2.500};
    const static vector<double> eta_vals_50 =  {0.771, 0.806, 0.855, 0.858, 0.843, 0.810, 0.800, 0.808, 0.802, 0.730, 0.653};

    const static vector<double> eta_edges_60 = {0.000, 0.050, 0.350, 0.707, 0.981, 1.278, 1.468, 1.694, 1.944, 2.242, 2.453, 2.500};
    const static vector<double> eta_vals_60 =  {0.773, 0.816, 0.866, 0.865, 0.853, 0.820, 0.812, 0.817, 0.804, 0.726, 0.645};

    const static vector<double> eta_edges_80 = {0.000, 0.051, 0.374, 0.720, 0.981, 1.279, 1.468, 1.707, 1.945, 2.207, 2.457, 2.500};
    const static vector<double> eta_vals_80 =  {0.819, 0.855, 0.899, 0.906, 0.900, 0.869, 0.865, 0.873, 0.869, 0.868, 0.859};

    const static vector<double> et_edges = { 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 80 };
    const static vector< vector<double> > et_eta_edges = { eta_edges_10, eta_edges_15, eta_edges_20, eta_edges_25, eta_edges_30, eta_edges_35, eta_edges_40, eta_edges_45, eta_edges_50, eta_edges_60, eta_edges_80 };
    const static vector< vector<double> > et_eta_vals  = { eta_vals_10, eta_vals_15, eta_vals_20, eta_vals_25, eta_vals_30, eta_vals_35, eta_vals_40, eta_vals_45, eta_vals_50, eta_vals_60, eta_vals_80 };

    if (e.abseta() > 2.5 || e.Et() < 10*GeV) return 0.0;
    const int i_et = binIndex(e.Et()/GeV, et_edges, true);
    const int i_eta = binIndex(e.abseta(), et_eta_edges[i_et]);
    return et_eta_vals[i_et][i_eta];
  }

  /// @brief ATLAS Run 2 'tight' electron identification/selection efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_IDEFF_ATLAS_RUN2_TIGHT(const Particle& e) {
    return ELECTRON_IDEFF_ATLAS_RUN1_TIGHT(e);
  }



  /// ATLAS Run 1 electron reco smearing
  inline Particle ELECTRON_SMEAR_ATLAS_RUN1(const Particle& e) {
    static const vector<double> edges_eta = {0., 2.5, 3.};
    static const vector<double> edges_pt = {0., 0.1, 25.};
    static const vector<double> e2s = {0.000, 0.015, 0.005,
                                       0.005, 0.005, 0.005,
                                       0.107, 0.107, 0.107};
    static const vector<double> es = {0.00, 0.00, 0.05,
                                      0.05, 0.05, 0.05,
                                      2.08, 2.08, 2.08};
    static const vector<double> cs = {0.00, 0.00, 0.25,
                                      0.25, 0.25, 0.25,
                                      0.00, 0.00, 0.00};

    const int i_eta = binIndex(e.abseta(), edges_eta, true);
    const int i_pt = binIndex(e.pT()/GeV, edges_pt, true);
    const int i = i_eta*edges_pt.size() + i_pt;

    // Calculate absolute resolution in GeV
    const double c1 = sqr(e2s[i]), c2 = sqr(es[i]), c3 = sqr(cs[i]);
    const double resolution = sqrt(c1*e.E2() + c2*e.E() + c3) * GeV;

    // normal_distribution<> d(e.E(), resolution);
    // const double mass = e.mass2() > 0 ? e.mass() : 0; //< numerical carefulness...
    // const double smeared_E = max(d(gen), mass); //< can't let the energy go below the mass!
    // return Particle(e.pid(), FourMomentum::mkEtaPhiME(e.eta(), e.phi(), mass, smeared_E));
    return Particle(e.pid(), P4_SMEAR_E_GAUSS(e, resolution));
  }


  /// ATLAS Run 2 electron reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle ELECTRON_SMEAR_ATLAS_RUN2(const Particle& e) {
    return ELECTRON_SMEAR_ATLAS_RUN1(e);
  }


  /// @todo Add charge flip efficiency?



  /// CMS Run 1 electron reconstruction efficiency
  /// @todo How to use this in combination with tracking eff?
  inline double ELECTRON_EFF_CMS_RUN1(const Particle& e) {
    if (e.abseta() > 2.5) return 0;
    if (e.pT() < 10*GeV) return 0;
    return (e.abseta() < 1.5) ? 0.95 : 0.85;
  }


  /// CMS Run 2 electron reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double ELECTRON_EFF_CMS_RUN2(const Particle& e) {
    return ELECTRON_EFF_CMS_RUN1(e);
  }


  /// @brief CMS electron energy smearing, preserving direction
  ///
  /// Calculate resolution
  /// for pT > 0.1 GeV, E resolution = |eta| < 0.5 -> sqrt(0.06^2 + pt^2 * 1.3e-3^2)
  ///                                  |eta| < 1.5 -> sqrt(0.10^2 + pt^2 * 1.7e-3^2)
  ///                                  |eta| < 2.5 -> sqrt(0.25^2 + pt^2 * 3.1e-3^2)
  inline Particle ELECTRON_SMEAR_CMS_RUN1(const Particle& e) {
    // Calculate absolute resolution in GeV from functional form
    double resolution = 0;
    const double abseta = e.abseta();
    if (e.pT() > 0.1*GeV && abseta < 2.5) { //< should be a given from efficiencies
      if (abseta < 0.5) {
        resolution = add_quad(0.06, 1.3e-3 * e.pT()/GeV) * GeV;
      } else if (abseta < 1.5) {
        resolution = add_quad(0.10, 1.7e-3 * e.pT()/GeV) * GeV;
      } else { // still |eta| < 2.5
        resolution = add_quad(0.25, 3.1e-3 * e.pT()/GeV) * GeV;
      }
    }

    // normal_distribution<> d(e.E(), resolution);
    // const double mass = e.mass2() > 0 ? e.mass() : 0; //< numerical carefulness...
    // const double smeared_E = max(d(gen), mass); //< can't let the energy go below the mass!
    // return Particle(e.pid(), FourMomentum::mkEtaPhiME(e.eta(), e.phi(), mass, smeared_E));
    return Particle(e.pid(), P4_SMEAR_E_GAUSS(e, resolution));
  }


  /// CMS Run 2 electron reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle ELECTRON_SMEAR_CMS_RUN2(const Particle& e) {
    return ELECTRON_SMEAR_CMS_RUN1(e);
  }

  //@}



  /// @name Photon efficiency and smearing functions
  //@{

  /// ATLAS Run 1 photon reco efficiency
  /// @todo Currently identical to CMS, cf. Delphes
  inline double PHOTON_EFF_ATLAS_RUN1(const Particle& y) {
    if (y.pT() < 10*GeV || y.abseta() > 2.5) return 0;
    return (y.abseta() < 1.5) ? 0.95 : 0.85;
  }

  /// ATLAS Run 2 photon reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double PHOTON_EFF_ATLAS_RUN2(const Particle& y) {
    return PHOTON_EFF_ATLAS_RUN1(y);
  }

  /// CMS Run 1 photon reco efficiency
  /// @todo Currently identical to ATLAS, cf. Delphes
  inline double PHOTON_EFF_CMS_RUN1(const Particle& y) {
    if (y.pT() < 10*GeV || y.abseta() > 2.5) return 0;
    return (y.abseta() < 1.5) ? 0.95 : 0.85;
  }

  /// CMS Run 2 photon reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double PHOTON_EFF_CMS_RUN2(const Particle& y) {
    return PHOTON_EFF_CMS_RUN1(y);
  }

  //@}



  /// @name Muon efficiency and smearing functions
  //@{

  /// ATLAS Run 1 muon reco efficiency
  inline double MUON_EFF_ATLAS_RUN1(const Particle& m) {
    if (m.abseta() > 2.7) return 0;
    if (m.pT() < 10*GeV) return 0;
    return (m.abseta() < 1.5) ? 0.95 : 0.85;
  }

  /// ATLAS Run 2 muon reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double MUON_EFF_ATLAS_RUN2(const Particle& m) {
    return MUON_EFF_ATLAS_RUN1(m);
  }


  /// ATLAS Run 1 muon reco smearing
  inline Particle MUON_SMEAR_ATLAS_RUN1(const Particle& m) {
    static const vector<double> edges_eta = {0, 1.5, 2.5};
    static const vector<double> edges_pt = {0, 0.1, 1.0, 10., 200.};
    static const vector<double> res = {0., 0.03, 0.02, 0.03, 0.05,
                                       0., 0.04, 0.03, 0.04, 0.05};

    const int i_eta = binIndex(m.abseta(), edges_eta, true);
    const int i_pt = binIndex(m.pT()/GeV, edges_pt, true);
    const int i = i_eta*edges_pt.size() + i_pt;

    const double resolution = res[i];

    // Smear by a Gaussian centered on the current pT, with width given by the resolution
    // normal_distribution<> d(m.pT(), resolution*m.pT());
    // const double smeared_pt = max(d(gen), 0.);
    // const double mass = m.mass2() > 0 ? m.mass() : 0; //< numerical carefulness...
    // return Particle(m.pid(), FourMomentum::mkEtaPhiMPt(m.eta(), m.phi(), mass, smeared_pt));
    return Particle(m.pid(), P4_SMEAR_PT_GAUSS(m, resolution*m.pT()));
  }

  /// ATLAS Run 2 muon reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle MUON_SMEAR_ATLAS_RUN2(const Particle& m) {
    return MUON_SMEAR_ATLAS_RUN1(m);
  }




  /// CMS Run 1 muon reco efficiency
  inline double MUON_EFF_CMS_RUN1(const Particle& m) {
    if (m.abseta() > 2.4) return 0;
    if (m.pT() < 10*GeV) return 0;
    return 0.95 * (m.abseta() < 1.5 ? 1 : exp(0.5 - 5e-4*m.pT()/GeV));
  }

  /// CMS Run 2 muon reco efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double MUON_EFF_CMS_RUN2(const Particle& m) {
    return MUON_EFF_CMS_RUN1(m);
  }


  /// CMS Run 1 muon reco smearing
  inline Particle MUON_SMEAR_CMS_RUN1(const Particle& m) {
    // Calculate fractional resolution
    // for pT > 0.1 GeV, mom resolution = |eta| < 0.5 -> sqrt(0.01^2 + pt^2 * 2.0e-4^2)
    //                                    |eta| < 1.5 -> sqrt(0.02^2 + pt^2 * 3.0e-4^2)
    //                                    |eta| < 2.5 -> sqrt(0.05^2 + pt^2 * 2.6e-4^2)
    double resolution = 0;
    const double abseta = m.abseta();
    if (m.pT() > 0.1*GeV && abseta < 2.5) {
      if (abseta < 0.5) {
        resolution = add_quad(0.01, 2.0e-4 * m.pT()/GeV);
      } else if (abseta < 1.5) {
        resolution = add_quad(0.02, 3.0e-4 * m.pT()/GeV);
      } else { // still |eta| < 2.5... but isn't CMS' mu acceptance < 2.4?
        resolution = add_quad(0.05, 2.6e-4 * m.pT()/GeV);
      }
    }

    // Smear by a Gaussian centered on the current pT, with width given by the resolution
    // normal_distribution<> d(m.pT(), resolution*m.pT());
    // const double smeared_pt = max(d(gen), 0.);
    // const double mass = m.mass2() > 0 ? m.mass() : 0; //< numerical carefulness...
    // return Particle(m.pid(), FourMomentum::mkEtaPhiMPt(m.eta(), m.phi(), mass, smeared_pt));
    return Particle(m.pid(), P4_SMEAR_PT_GAUSS(m, resolution*m.pT()));
  }

  /// CMS Run 2 muon reco smearing
  /// @todo Currently just a copy of the Run 1 version: fix!
  inline Particle MUON_SMEAR_CMS_RUN2(const Particle& m) {
    return MUON_SMEAR_CMS_RUN1(m);
  }

  //@}



  /// @name Tau efficiency and smearing functions
  //@{

  /// @brief ATLAS Run 1 8 TeV tau efficiencies (medium working point)
  ///
  /// Taken from http://arxiv.org/pdf/1412.7086.pdf
  ///   20-40 GeV 1-prong LMT eff|mis = 0.66|1/10, 0.56|1/20, 0.36|1/80
  ///   20-40 GeV 3-prong LMT eff|mis = 0.45|1/60, 0.38|1/100, 0.27|1/300
  ///   > 40 GeV 1-prong LMT eff|mis = 0.66|1/15, 0.56|1/25, 0.36|1/80
  ///   > 40 GeV 3-prong LMT eff|mis = 0.45|1/250, 0.38|1/400, 0.27|1/1300
  inline double TAU_EFF_ATLAS_RUN1(const Particle& t) {
    if (t.abseta() > 2.5) return 0; //< hmm... mostly
    double pThadvis = 0;
    Particles chargedhadrons;
    for (const Particle& p : t.children()) {
      if (p.isHadron()) {
        pThadvis += p.pT(); //< right definition? Paper is unclear
        if (p.charge3() != 0 && p.abseta() < 2.5 && p.pT() > 1*GeV) chargedhadrons += p;
      }
    }
    if (chargedhadrons.empty()) return 0; //< leptonic tau
    if (pThadvis < 20*GeV) return 0; //< below threshold
    if (pThadvis < 40*GeV) {
      if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.56 : 1/20.;
      if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.38 : 1/100.;
    } else {
      if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.56 : 1/25.;
      if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.38 : 1/400.;
    }
    return 0;
  }


  /// @brief ATLAS Run 2 13 TeV tau efficiencies (medium working point)
  ///
  /// From https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/PUBNOTES/ATL-PHYS-PUB-2015-045/ATL-PHYS-PUB-2015-045.pdf
  ///   LMT 1 prong efficiency/mistag = 0.6|1/30, 0.55|1/50, 0.45|1/120
  ///   LMT 3 prong efficiency/mistag = 0.5|1/30, 0.4|1/110, 0.3|1/300
  inline double TAU_EFF_ATLAS_RUN2(const Particle& t) {
    if (t.abseta() > 2.5) return 0; //< hmm... mostly
    double pThadvis = 0;
    Particles chargedhadrons;
    for (const Particle& p : t.children()) {
      if (p.isHadron()) {
        pThadvis += p.pT(); //< right definition? Paper is unclear
        if (p.charge3() != 0 && p.abseta() < 2.5 && p.pT() > 1*GeV) chargedhadrons += p;
      }
    }
    if (chargedhadrons.empty()) return 0; //< leptonic tau
    if (pThadvis < 20*GeV) return 0; //< below threshold
    if (chargedhadrons.size() == 1) return (t.abspid() == PID::TAU) ? 0.55 : 1/50.;
    if (chargedhadrons.size() == 3) return (t.abspid() == PID::TAU) ? 0.40 : 1/110.;
    return 0;
  }


  /// ATLAS Run 1 tau smearing
  /// @todo Currently a copy of the crappy jet smearing that is probably wrong...
  inline Particle TAU_SMEAR_ATLAS_RUN1(const Particle& t) {
    // Const fractional resolution for now
    static const double resolution = 0.03;

    // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
    /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
    const double fsmear = max(randnorm(1., resolution), 0.);
    const double mass = t.mass2() > 0 ? t.mass() : 0; //< numerical carefulness...
    return Particle(t.pid(), FourMomentum::mkXYZM(t.px()*fsmear, t.py()*fsmear, t.pz()*fsmear, mass));
  }


  /// ATLAS Run 2 tau smearing
  /// @todo Currently a copy of the Run 1 version
  inline Particle TAU_SMEAR_ATLAS_RUN2(const Particle& t) {
    return TAU_SMEAR_ATLAS_RUN1(t);
  }


  /// CMS Run 2 tau efficiency
  ///
  /// @todo Needs work; this is the dumb version from Delphes 3.3.2
  inline double TAU_EFF_CMS_RUN2(const Particle& t) {
    return (t.abspid() == PID::TAU) ? 0.6 : 0;
  }

  /// CMS Run 1 tau efficiency
  ///
  /// @todo Needs work; this is just a copy of the Run 2 version in Delphes 3.3.2
  inline double TAU_EFF_CMS_RUN1(const Particle& t) {
    return TAU_EFF_CMS_RUN2(t);
  }


  /// CMS Run 1 tau smearing
  /// @todo Currently a copy of the crappy ATLAS one
  inline Particle TAU_SMEAR_CMS_RUN1(const Particle& t) {
    return TAU_SMEAR_ATLAS_RUN1(t);
  }


  /// CMS Run 2 tau smearing
  /// @todo Currently a copy of the Run 1 version
  inline Particle TAU_SMEAR_CMS_RUN2(const Particle& t) {
    return TAU_SMEAR_CMS_RUN1(t);
  }

  //@}



  /// @name Jet efficiency and smearing functions
  //@{

  /// Return the ATLAS Run 1 jet flavour tagging efficiency for the given Jet
  inline double JET_BTAG_ATLAS_RUN1(const Jet& j) {
    /// @todo This form drops past ~100 GeV, asymptotically to zero efficiency... really?!
    if (j.abseta() > 2.5) return 0;
    const auto ftagsel = [&](const Particle& p){ return p.pT() > 5*GeV && deltaR(p,j) < 0.3; };
    if (j.bTagged(ftagsel)) return 0.80*tanh(0.003*j.pT()/GeV)*(30/(1+0.0860*j.pT()/GeV));
    if (j.cTagged(ftagsel)) return 0.20*tanh(0.020*j.pT()/GeV)*( 1/(1+0.0034*j.pT()/GeV));
    return 0.002 + 7.3e-6*j.pT()/GeV;
  }
  /// Return the ATLAS Run 2 MC2c20 jet flavour tagging efficiency for the given Jet
  inline double JET_BTAG_ATLAS_RUN2_MV2C20(const Jet& j) {
    if (j.abseta() > 2.5) return 0;
    if (j.bTagged(Cuts::pT > 5*GeV)) return 0.77;
    if (j.cTagged(Cuts::pT > 5*GeV)) return 1/4.5;
    return 1/140.;
  }
  /// Return the ATLAS Run 2 MC2c10 jet flavour tagging efficiency for the given Jet
  inline double JET_BTAG_ATLAS_RUN2_MV2C10(const Jet& j) {
    if (j.abseta() > 2.5) return 0;
    if (j.bTagged(Cuts::pT > 5*GeV)) return 0.77;
    if (j.cTagged(Cuts::pT > 5*GeV)) return 1/6.0;
    return 1/134.;
  }


  /// ATLAS Run 1 jet smearing
  inline Jet JET_SMEAR_ATLAS_RUN1(const Jet& j) {
    // Jet energy resolution lookup
    //   Implemented by Matthias Danninger for GAMBIT, based roughly on
    //   https://atlas.web.cern.ch/Atlas/GROUPS/PHYSICS/CONFNOTES/ATLAS-CONF-2015-017/
    //   Parameterisation can be still improved, but eta dependence is minimal
    /// @todo Also need a JES uncertainty component?
    static const vector<double> binedges_pt = {0., 50., 70., 100., 150., 200., 1000., 10000.};
    static const vector<double> jer = {0.145, 0.115, 0.095, 0.075, 0.07, 0.05, 0.04, 0.04}; //< note overflow value
    const int ipt = binIndex(j.pT()/GeV, binedges_pt, true);
    if (ipt < 0) return j;
    const double resolution = jer.at(ipt);

    // Smear by a Gaussian centered on 1 with width given by the (fractional) resolution
    /// @todo Is this the best way to smear? Should we preserve the energy, or pT, or direction?
    const double fsmear = max(randnorm(1., resolution), 0.);
    const double mass = j.mass2() > 0 ? j.mass() : 0; //< numerical carefulness...
    return Jet(FourMomentum::mkXYZM(j.px()*fsmear, j.py()*fsmear, j.pz()*fsmear, mass));
  }

  /// ATLAS Run 2 jet smearing
  /// @todo Just a copy of the Run 1 one: improve!!
  inline Jet JET_SMEAR_ATLAS_RUN2(const Jet& j) {
    return JET_SMEAR_ATLAS_RUN1(j);
  }

  /// CMS Run 2 jet smearing
  /// @todo Just a copy of the suboptimal ATLAS one: improve!!
  inline Jet JET_SMEAR_CMS_RUN2(const Jet& j) {
    return JET_SMEAR_ATLAS_RUN1(j);
  }

  //@}


  /// @name ETmiss smearing functions
  //@{

  inline Vector3 MET_SMEAR_IDENTITY(const Vector3& met, double) { return met; }

  /// @brief ATLAS Run 1 ETmiss smearing
  ///
  /// Based on https://arxiv.org/pdf/1108.5602v2.pdf, Figs 14 and 15
  inline Vector3 MET_SMEAR_ATLAS_RUN1(const Vector3& met, double set) {
    // Linearity offset (Fig 14)
    Vector3 smeared_met = met;
    if (met.mod()/GeV < 25*GeV) smeared_met *= 1.05;
    else if (met.mod()/GeV < 40*GeV) smeared_met *= (1.05 - (0.04/15)*(met.mod()/GeV - 25)); //< linear decrease
    else smeared_met *= 1.01;

    // Smear by a Gaussian with width given by the resolution(sumEt) ~ 0.45 sqrt(sumEt) GeV
    const double resolution = 0.45 * sqrt(set/GeV) * GeV;
    const double metsmear = max(randnorm(smeared_met.mod(), resolution), 0.);
    smeared_met = metsmear * smeared_met.unit();

    return smeared_met;
  }

  /// ATLAS Run 2 ETmiss smearing
  /// @todo Just a copy of the Run 1 one: improve!!
  inline Vector3 MET_SMEAR_ATLAS_RUN2(const Vector3& met, double set) {
    return MET_SMEAR_ATLAS_RUN1(met, set);
  }

  /// CMS Run 1 ETmiss smearing
  /// @todo Just a copy of the ATLAS one: improve!!
  inline Vector3 MET_SMEAR_CMS_RUN1(const Vector3& met, double set) {
    return MET_SMEAR_ATLAS_RUN1(met, set);
  }

  /// CMS Run 2 ETmiss smearing
  /// @todo Just a copy of the ATLAS one: improve!!
  inline Vector3 MET_SMEAR_CMS_RUN2(const Vector3& met, double set) {
    return MET_SMEAR_ATLAS_RUN2(met, set);
  }

  //@}


  /// @name Tracking efficiency and smearing functions
  //@{

  /// ATLAS Run 1 tracking efficiency
  inline double TRK_EFF_ATLAS_RUN1(const Particle& p) {
    if (p.charge3() == 0) return 0;
    if (p.abseta() > 2.5) return 0;
    if (p.pT() < 0.1*GeV) return 0;

    if (p.abspid() == PID::ELECTRON) {
      if (p.abseta() < 1.5) {
        if (p.pT() < 1*GeV) return 0.73;
        if (p.pT() < 100*GeV) return 0.95;
        return 0.99;
      } else {
        if (p.pT() < 1*GeV) return 0.50;
        if (p.pT() < 100*GeV) return 0.83;
        else return 0.90;
      }
    } else if (p.abspid() == PID::MUON) {
      if (p.abseta() < 1.5) {
        return (p.pT() < 1*GeV) ? 0.75 : 0.99;
      } else {
        return (p.pT() < 1*GeV) ? 0.70 : 0.98;
      }
    } else { // charged hadrons
      if (p.abseta() < 1.5) {
        return (p.pT() < 1*GeV) ? 0.70 : 0.95;
      } else {
        return (p.pT() < 1*GeV) ? 0.60 : 0.85;
      }
    }
  }

  /// ATLAS Run 2 tracking efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double TRK_EFF_ATLAS_RUN2(const Particle& p) {
    return TRK_EFF_ATLAS_RUN1(p);
  }


  /// CMS Run 1 tracking efficiency
  inline double TRK_EFF_CMS_RUN1(const Particle& p) {
    if (p.charge3() == 0) return 0;
    if (p.abseta() > 2.5) return 0;
    if (p.pT() < 0.1*GeV) return 0;

    if (p.abspid() == PID::ELECTRON) {
      if (p.abseta() < 1.5) {
        if (p.pT() < 1*GeV) return 0.73;
        if (p.pT() < 100*GeV) return 0.95;
        return 0.99;
      } else {
        if (p.pT() < 1*GeV) return 0.50;
        if (p.pT() < 100*GeV) return 0.83;
        else return 0.90;
      }
    } else if (p.abspid() == PID::MUON) {
      if (p.abseta() < 1.5) {
        return (p.pT() < 1*GeV) ? 0.75 : 0.99;
      } else {
        return (p.pT() < 1*GeV) ? 0.70 : 0.98;
      }
    } else { // charged hadrons
      if (p.abseta() < 1.5) {
        return (p.pT() < 1*GeV) ? 0.70 : 0.95;
      } else {
        return (p.pT() < 1*GeV) ? 0.60 : 0.85;
      }
    }
  }

  /// CMS Run 2 tracking efficiency
  /// @todo Currently just a copy of Run 1: fix!
  inline double TRK_EFF_CMS_RUN2(const Particle& p) {
    return TRK_EFF_CMS_RUN1(p);
  }

  //@}


}

#endif
