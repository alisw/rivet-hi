// -*- C++ -*-
#ifndef RIVET_MC_ParticleAnalysis_HH
#define RIVET_MC_ParticleAnalysis_HH

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Base class providing common functionality for MC particle species validation analyses
  class MC_ParticleAnalysis : public Analysis {
  public:

    /// Default constructor.
    MC_ParticleAnalysis(const string& name,
                        size_t nparticles,
                        const string& particle_name);


    /// @name Analysis methods
    //@{

    /// Bookings as usual, using the particle name specified by the derived classe
    virtual void init();

    /// To be implemented by derived classes, making particle selection then calling _analyze
    virtual void analyze(const Event& event) = 0;

    /// Normalization, division, etc.
    virtual void finalize();

    /// For derived classes to call, passing the sorted particle collection that they wish to analyse
    virtual void _analyze(const Event& event, const Particles& particles);
    //@}


  protected:

    /// The number of particles for which histograms are to be initialised
    size_t _nparts;

    /// The name of the particle species/group being analysed
    std::string _pname;

    /// @name Histograms
    //@{
    std::vector<Histo1DPtr> _h_pt;
    std::vector<Histo1DPtr> _h_eta;
    std::vector<Histo1DPtr> _h_eta_plus, _h_eta_minus;
    std::vector<Histo1DPtr> _h_rap;
    std::vector<Histo1DPtr> _h_rap_plus, _h_rap_minus;
    std::map<std::pair<size_t, size_t>, Histo1DPtr> _h_deta;
    std::map<std::pair<size_t, size_t>, Histo1DPtr> _h_dphi;
    std::map<std::pair<size_t, size_t>, Histo1DPtr> _h_dR;
    Histo1DPtr _h_multi_exclusive, _h_multi_inclusive;
    Histo1DPtr _h_multi_exclusive_prompt, _h_multi_inclusive_prompt;
    Scatter2DPtr _h_multi_ratio, _h_multi_ratio_prompt;
    //@}

  };

}

#endif
