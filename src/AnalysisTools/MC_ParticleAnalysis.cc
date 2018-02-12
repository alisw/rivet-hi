// -*- C++ -*-
#include "Rivet/Analyses/MC_ParticleAnalysis.hh"

namespace Rivet {

  


  MC_ParticleAnalysis::MC_ParticleAnalysis(const string& name,
                                           size_t nparticles,
                                           const string& particle_name)
    : Analysis(name),
      _nparts(nparticles), _pname(particle_name),
      _h_pt(nparticles),
      _h_eta(nparticles), _h_eta_plus(nparticles), _h_eta_minus(nparticles),
      _h_rap(nparticles), _h_rap_plus(nparticles), _h_rap_minus(nparticles)
  {
    setNeedsCrossSection(true); // legitimate use, since a base class has no .info file!
  }



  // Book histograms
  void MC_ParticleAnalysis::init() {

    for (size_t i = 0; i < _nparts; ++i) {
      const string ptname = _pname + "_pt_" + to_str(i+1);
      const double ptmax = 1.0/(double(i)+2.0) * (sqrtS()>0.?sqrtS():14000.)/GeV/2.0;
      const int nbins_pt = 100/(i+1);
      _h_pt[i] = bookHisto1D(ptname, logspace(nbins_pt, 1.0, ptmax));

      const string etaname = _pname + "_eta_" + to_str(i+1);
      _h_eta[i] = bookHisto1D(etaname, i > 1 ? 25 : 50, -5.0, 5.0);
      _h_eta_plus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));
      _h_eta_minus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));

      const string rapname = _pname + "_y_" + to_str(i+1);
      _h_rap[i] = bookHisto1D(rapname, i > 1 ? 25 : 50, -5.0, 5.0);
      _h_rap_plus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));
      _h_rap_minus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));

      for (size_t j = i+1; j < min(size_t(3), _nparts); ++j) {
        const pair<size_t, size_t> ij = std::make_pair(i, j);

        string detaname = _pname + "s_deta_" + to_str(i+1) + to_str(j+1);
        _h_deta.insert(make_pair(ij, bookHisto1D(detaname, 25, -5.0, 5.0)));

        string dphiname = _pname + "s_dphi_" + to_str(i+1) + to_str(j+1);
        _h_dphi.insert(make_pair(ij, bookHisto1D(dphiname, 25, 0.0, M_PI)));

        string dRname = _pname + "s_dR_" + to_str(i+1) + to_str(j+1);
        _h_dR.insert(make_pair(ij, bookHisto1D(dRname, 25, 0.0, 5.0)));
      }
    }

    _h_multi_exclusive = bookHisto1D(_pname + "_multi_exclusive", _nparts+3, -0.5, _nparts+3-0.5);
    _h_multi_inclusive = bookHisto1D(_pname + "_multi_inclusive", _nparts+3, -0.5, _nparts+3-0.5);
    _h_multi_ratio = bookScatter2D(_pname + "_multi_ratio");

    _h_multi_exclusive_prompt = bookHisto1D(_pname + "_multi_exclusive_prompt", _nparts+3, -0.5, _nparts+3-0.5);
    _h_multi_inclusive_prompt = bookHisto1D(_pname + "_multi_inclusive_prompt", _nparts+3, -0.5, _nparts+3-0.5);
    _h_multi_ratio_prompt = bookScatter2D(_pname + "_multi_ratio_prompt");
  }


  // Do the analysis
  void MC_ParticleAnalysis::_analyze(const Event& event, const Particles& particles) {
    const double weight = event.weight();
    Particles promptparticles;
    foreach (const Particle& p, particles)
      if (!p.fromDecay()) promptparticles += p;

    for (size_t i = 0; i < _nparts; ++i) {
      if (particles.size() < i+1) continue;
      _h_pt[i]->fill(particles[i].pt()/GeV, weight);

      // Eta
      const double eta_i = particles[i].eta();
      _h_eta[i]->fill(eta_i, weight);
      (eta_i > 0.0 ? _h_eta_plus : _h_eta_minus)[i]->fill(fabs(eta_i), weight);

      // Rapidity
      const double rap_i = particles[i].rapidity();
      _h_rap[i]->fill(rap_i, weight);
      (rap_i > 0.0 ? _h_rap_plus : _h_rap_minus)[i]->fill(fabs(rap_i), weight);

      // Inter-particle properties
      for (size_t j = i+1; j < min(size_t(3),_nparts); ++j) {
        if (particles.size() < j+1) continue;
        std::pair<size_t, size_t> ij = std::make_pair(i, j);
        double deta = particles[i].eta() - particles[j].eta();
        double dphi = deltaPhi(particles[i].momentum(), particles[j].momentum());
        double dR = deltaR(particles[i].momentum(), particles[j].momentum());
        _h_deta[ij]->fill(deta, weight);
        _h_dphi[ij]->fill(dphi, weight);
        _h_dR[ij]->fill(dR, weight);
      }
    }

    // Multiplicities
    _h_multi_exclusive->fill(particles.size(), weight);
    _h_multi_exclusive_prompt->fill(promptparticles.size(), weight);
    for (size_t i = 0; i < _nparts+2; ++i) {
      if (particles.size() >= i) _h_multi_inclusive->fill(i, weight);
      if (promptparticles.size() >= i) _h_multi_inclusive_prompt->fill(i, weight);
    }

  }


  // Finalize
  void MC_ParticleAnalysis::finalize() {
    for (size_t i = 0; i < _nparts; ++i) {
      scale(_h_pt[i], crossSection()/sumOfWeights());
      scale(_h_eta[i], crossSection()/sumOfWeights());
      scale(_h_rap[i], crossSection()/sumOfWeights());

      // Create eta/rapidity ratio plots
      divide(*_h_eta_plus[i], *_h_eta_minus[i], bookScatter2D(_pname + "_eta_pmratio_" + to_str(i+1)));
      divide(*_h_rap_plus[i], *_h_rap_minus[i], bookScatter2D(_pname + "_y_pmratio_" + to_str(i+1)));
    }

    // Scale the d{eta,phi,R} histograms
    typedef map<pair<size_t, size_t>, Histo1DPtr> HistMap;
    foreach (HistMap::value_type& it, _h_deta) scale(it.second, crossSection()/sumOfWeights());
    foreach (HistMap::value_type& it, _h_dphi) scale(it.second, crossSection()/sumOfWeights());
    foreach (HistMap::value_type& it, _h_dR) scale(it.second, crossSection()/sumOfWeights());

    // Fill inclusive multi ratios
    for (size_t i = 0; i < _h_multi_inclusive->numBins()-1; ++i) {
      _h_multi_ratio->addPoint(i+1, 0, 0.5, 0);
      if (_h_multi_inclusive->bin(i).sumW() > 0.0) {
        const double ratio = _h_multi_inclusive->bin(i+1).sumW() / _h_multi_inclusive->bin(i).sumW();
        const double relerr_i = _h_multi_inclusive->bin(i).relErr();
        const double relerr_j = _h_multi_inclusive->bin(i+1).relErr();
        const double err = ratio * (relerr_i + relerr_j);
        _h_multi_ratio->point(i).setY(ratio, err);
      }
    }
    for (size_t i = 0; i < _h_multi_inclusive_prompt->numBins()-1; ++i) {
      _h_multi_ratio_prompt->addPoint(i+1, 0, 0.5, 0);
      if (_h_multi_inclusive_prompt->bin(i).sumW() > 0.0) {
        const double ratio = _h_multi_inclusive_prompt->bin(i+1).sumW() / _h_multi_inclusive_prompt->bin(i).sumW();
        const double relerr_i = _h_multi_inclusive_prompt->bin(i).relErr();
        const double relerr_j = _h_multi_inclusive_prompt->bin(i+1).relErr();
        const double err = ratio * (relerr_i + relerr_j);
        _h_multi_ratio_prompt->point(i).setY(ratio, err);
      }
    }

    scale(_h_multi_exclusive, crossSection()/sumOfWeights());
    scale(_h_multi_exclusive_prompt, crossSection()/sumOfWeights());
    scale(_h_multi_inclusive, crossSection()/sumOfWeights());
    scale(_h_multi_inclusive_prompt, crossSection()/sumOfWeights());
  }


}
