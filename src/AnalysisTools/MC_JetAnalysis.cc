// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  


  MC_JetAnalysis::MC_JetAnalysis(const string& name,
                                 size_t njet,
                                 const string& jetpro_name,
                                 double jetptcut)
    : Analysis(name), _njet(njet), _jetpro_name(jetpro_name), _jetptcut(jetptcut),
      _h_pT_jet(njet),
      _h_eta_jet(njet), _h_eta_jet_plus(njet), _h_eta_jet_minus(njet),
      _h_rap_jet(njet), _h_rap_jet_plus(njet), _h_rap_jet_minus(njet),
      _h_mass_jet(njet)
  {
    setNeedsCrossSection(true); // legitimate use, since a base class has no .info file!
  }



  // Book histograms
  void MC_JetAnalysis::init() {
    const double sqrts = sqrtS() ? sqrtS() : 14000.*GeV;

    for (size_t i = 0; i < _njet; ++i) {
      const string pTname = "jet_pT_" + to_str(i+1);
      const double pTmax = 1.0/(double(i)+2.0) * sqrts/GeV/2.0;
      const int nbins_pT = 100/(i+1);
      if (pTmax > 10) { // Protection aginst logspace exception, needed for LEP
        _h_pT_jet[i] = bookHisto1D(pTname, logspace(nbins_pT, 10.0, pTmax));
      }

      const string massname = "jet_mass_" + to_str(i+1);
      const double mmax = 100.0;
      const int nbins_m = 100/(i+1);
      _h_mass_jet[i] = bookHisto1D(massname, logspace(nbins_m, 1.0, mmax));

      const string etaname = "jet_eta_" + to_str(i+1);
      _h_eta_jet[i] = bookHisto1D(etaname, i > 1 ? 25 : 50, -5.0, 5.0);
      _h_eta_jet_plus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));
      _h_eta_jet_minus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));

      const string rapname = "jet_y_" + to_str(i+1);
      _h_rap_jet[i] = bookHisto1D(rapname, i>1 ? 25 : 50, -5.0, 5.0);
      _h_rap_jet_plus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));
      _h_rap_jet_minus[i].reset(new Histo1D(i > 1 ? 15 : 25, 0, 5));

      for (size_t j = i+1; j < min(size_t(3), _njet); ++j) {
        const std::pair<size_t, size_t> ij = std::make_pair(i, j);

        string detaname = "jets_deta_" + to_str(i+1) + to_str(j+1);
        _h_deta_jets.insert(make_pair(ij, bookHisto1D(detaname, 25, -5.0, 5.0)));

        string dphiname = "jets_dphi_" + to_str(i+1) + to_str(j+1);
        _h_dphi_jets.insert(make_pair(ij, bookHisto1D(dphiname, 25, 0.0, M_PI)));

        string dRname = "jets_dR_" + to_str(i+1) + to_str(j+1);
        _h_dR_jets.insert(make_pair(ij, bookHisto1D(dRname, 25, 0.0, 5.0)));
      }
    }

    _h_jet_multi_exclusive = bookHisto1D("jet_multi_exclusive", _njet+3, -0.5, _njet+3-0.5);
    _h_jet_multi_inclusive = bookHisto1D("jet_multi_inclusive", _njet+3, -0.5, _njet+3-0.5);
    _h_jet_multi_ratio = bookScatter2D("jet_multi_ratio");
    _h_jet_HT = bookHisto1D("jet_HT", logspace(50, _jetptcut, sqrts/GeV/2.0));
    _h_mjj_jets = bookHisto1D("jets_mjj", 40, 0.0, sqrts/GeV/2.0);
  }


  // Do the analysis
  void MC_JetAnalysis::analyze(const Event & e) {
    const double weight = e.weight();

    const Jets& jets = apply<FastJets>(e, _jetpro_name).jetsByPt(_jetptcut);

    for (size_t i = 0; i < _njet; ++i) {
      if (jets.size() < i+1) continue;
      _h_pT_jet[i]->fill(jets[i].pT()/GeV, weight);
      // Check for numerical precision issues with jet masses
      double m2_i = jets[i].mass2();
      if (m2_i < 0) {
        if (m2_i < -1e-4) {
          MSG_WARNING("Jet mass2 is negative: " << m2_i << " GeV^2. "
                      << "Truncating to 0.0, assuming numerical precision is to blame.");
        }
        m2_i = 0.0;
      }

      // Jet mass
      _h_mass_jet[i]->fill(sqrt(m2_i)/GeV, weight);

      // Jet eta
      const double eta_i = jets[i].eta();
      _h_eta_jet[i]->fill(eta_i, weight);
      (eta_i > 0.0 ? _h_eta_jet_plus : _h_eta_jet_minus)[i]->fill(fabs(eta_i), weight);

      // Jet rapidity
      const double rap_i = jets[i].rapidity();
      _h_rap_jet[i]->fill(rap_i, weight);
      (rap_i > 0.0 ? _h_rap_jet_plus : _h_rap_jet_minus)[i]->fill(fabs(rap_i), weight);

      // Inter-jet properties
      for (size_t j = i+1; j < min(size_t(3),_njet); ++j) {
        if (jets.size() < j+1) continue;
        std::pair<size_t, size_t> ij = std::make_pair(i, j);
        double deta = jets[i].eta()-jets[j].eta();
        double dphi = deltaPhi(jets[i].momentum(),jets[j].momentum());
        double dR = deltaR(jets[i].momentum(), jets[j].momentum());
        _h_deta_jets[ij]->fill(deta, weight);
        _h_dphi_jets[ij]->fill(dphi, weight);
        _h_dR_jets[ij]->fill(dR, weight);
      }
    }

    // Multiplicities
    _h_jet_multi_exclusive->fill(jets.size(), weight);
    for (size_t i = 0; i < _njet+2; ++i) {
      if (jets.size() >= i) {
        _h_jet_multi_inclusive->fill(i, weight);
      }
    }

    // HT
    double HT = 0.0;
    foreach (const Jet& jet, jets) {
      HT += jet.pT();
    }
    _h_jet_HT->fill(HT, weight);

    // mjj
    if (jets.size() > 1) {
      double mjj = (jets[0].momentum() + jets[1].momentum()).mass();
      _h_mjj_jets->fill(mjj, weight);
    }
  }


  // Finalize
  void MC_JetAnalysis::finalize() {
    for (size_t i = 0; i < _njet; ++i) {
      scale(_h_pT_jet[i], crossSection()/sumOfWeights());
      scale(_h_mass_jet[i], crossSection()/sumOfWeights());
      scale(_h_eta_jet[i], crossSection()/sumOfWeights());
      scale(_h_rap_jet[i], crossSection()/sumOfWeights());

      // Create eta/rapidity ratio plots
      divide(*_h_eta_jet_plus[i], *_h_eta_jet_minus[i], bookScatter2D("jet_eta_pmratio_" + to_str(i+1)));
      divide(*_h_rap_jet_plus[i], *_h_rap_jet_minus[i], bookScatter2D("jet_y_pmratio_" + to_str(i+1)));
    }

    // Scale the d{eta,phi,R} histograms
    typedef map<pair<size_t, size_t>, Histo1DPtr> HistMap;
    foreach (HistMap::value_type& it, _h_deta_jets) scale(it.second, crossSection()/sumOfWeights());
    foreach (HistMap::value_type& it, _h_dphi_jets) scale(it.second, crossSection()/sumOfWeights());
    foreach (HistMap::value_type& it, _h_dR_jets) scale(it.second, crossSection()/sumOfWeights());

    // Fill inclusive jet multi ratio
    int Nbins = _h_jet_multi_inclusive->numBins();
    for (int i = 0; i < Nbins-1; ++i) {
      _h_jet_multi_ratio->addPoint(i+1, 0, 0.5, 0);
      if (_h_jet_multi_inclusive->bin(i).sumW() > 0.0) {
        const double ratio = _h_jet_multi_inclusive->bin(i+1).sumW()/_h_jet_multi_inclusive->bin(i).sumW();
        const double relerr_i = _h_jet_multi_inclusive->bin(i).relErr();
        const double relerr_j = _h_jet_multi_inclusive->bin(i+1).relErr();
        const double err = ratio * (relerr_i + relerr_j);
        _h_jet_multi_ratio->point(i).setY(ratio, err);
      }
    }

    scale(_h_jet_multi_exclusive, crossSection()/sumOfWeights());
    scale(_h_jet_multi_inclusive, crossSection()/sumOfWeights());
    scale(_h_jet_HT, crossSection()/sumOfWeights());
    scale(_h_mjj_jets, crossSection()/sumOfWeights());
  }


}
