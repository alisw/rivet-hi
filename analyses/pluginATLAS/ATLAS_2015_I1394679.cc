// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Differential multijet cross-section measurement in different variables
  class ATLAS_2015_I1394679 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1394679);


    /// @name Analysis methods
    //@

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections here
      const FinalState fs;
      declare(fs, "FinalState");
      FastJets fj04(fs, FastJets::ANTIKT, 0.4, JetAlg::ALL_MUONS, JetAlg::DECAY_INVISIBLES);
      declare(fj04, "AntiKt4jets");

      // Histograms
      _h["pt1"] = bookHisto1D(1, 1, 1);
      _h["pt2"] = bookHisto1D(2, 1, 1);
      _h["pt3"] = bookHisto1D(3, 1, 1);
      _h["pt4"] = bookHisto1D(4, 1, 1);
      _h["HT"]  = bookHisto1D(5, 1, 1);
      _h["M4j"] = bookHisto1D(6, 1, 1);

      // Histograms with different pt/m4j cuts
      for (size_t i_hist = 0; i_hist < 4; ++i_hist) {
        _h["M2jratio_"+to_str(i_hist)] = bookHisto1D( 7 + i_hist, 1, 1);
        _h["dPhiMin2j_"+to_str(i_hist)] = bookHisto1D(11 + i_hist, 1, 1);
        _h["dPhiMin3j_"+to_str(i_hist)] = bookHisto1D(15 + i_hist, 1, 1);
        _h["dYMin2j_"+to_str(i_hist)] = bookHisto1D(19 + i_hist, 1, 1);
        _h["dYMin3j_"+to_str(i_hist)] = bookHisto1D(23 + i_hist, 1, 1);
        _h["dYMax2j_"+to_str(i_hist)] = bookHisto1D(27 + i_hist, 1, 1);
        for (size_t ygap = 0; ygap < 4; ++ygap) {
          _h["sumPtCent_"+to_str(ygap)+to_str(i_hist)] = bookHisto1D(31 + i_hist + ygap * 4, 1, 1);
        }
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& alljetskt4 = apply<FastJets>(event, "AntiKt4jets").jetsByPt(Cuts::pT > 60*GeV && Cuts::absrap < 2.8);

      // Require 4 jets with rap < 2.8 and passing pT cuts
      int nJets = alljetskt4.size();
      if (nJets < 4) vetoEvent;
      Jets leadingJetskt4 = alljetskt4; leadingJetskt4.resize(4);
      Jet jet1(leadingJetskt4[0]), jet2(leadingJetskt4[1]), jet3(leadingJetskt4[2]), jet4(leadingJetskt4[3]);
      if (jet1.pT() < 100*GeV) vetoEvent;
      if (jet4.pT() <  64*GeV) vetoEvent;

      // dR cut
      const double dRcut = 0.65;
      double drmin = 9999;
      for (int ijet = 0; ijet < 4; ++ijet) {
        for (int jjet = ijet + 1; jjet < 4; ++jjet) {
          double myDR = deltaR(alljetskt4[ijet], alljetskt4[jjet], RAPIDITY);
          if (myDR < drmin) drmin = myDR;
        }
      }
      if (drmin < dRcut) vetoEvent;

      // Variables for calculation in loops over jets
      FourMomentum sum_alljets;
      double HT = 0; // scalar sum of pt of 4 leading jets
      double Mjj = 99999; // minimum combined mass of 2 jets
      double minDphi_ij = 999, minDphi_ijk = 999; // minimum azimuthal distance btw 2 & 3 jets
      double maxDrap_ij = -999;  // maximum rapidity distance btw 2 jets
      double minDrap_ij = 999, minDrap_ijk = 999;  // minimum rapidity distance btw 2 & 3 jets
      size_t maxY_i = -1, maxY_j = -1;

      // Loop over 4 leading jets
      for (size_t ij = 0; ij< 4; ++ij) {
        Jet& jeti = leadingJetskt4.at(ij);
        sum_alljets += jeti.mom();
        HT += jeti.pT();

        for (size_t jj = 0; jj< 4; ++jj) {
          if ( ij == jj )  continue;
          Jet& jetj = leadingJetskt4.at(jj);

          const double auxDphi = fabs(deltaPhi(jeti, jetj));
          minDphi_ij = std::min(auxDphi, minDphi_ij);

          const double auxDrap = fabs(deltaRap(jeti, jetj));
          minDrap_ij = std::min(auxDrap, minDrap_ij);
          if (auxDrap > maxDrap_ij) {
            maxDrap_ij = auxDrap;
            maxY_i = ij;
            maxY_j = jj;
          }

          FourMomentum sum_twojets = jeti.mom() + jetj.mom();
          Mjj = std::min(Mjj, sum_twojets.mass());

          for (size_t kj = 0; kj < 4; ++kj) {
            if (kj == ij || kj == jj) continue;
            Jet& jetk = leadingJetskt4.at(kj);

            const double auxDphi2 = auxDphi + fabs(deltaPhi(jeti, jetk));
            minDphi_ijk = std::min(auxDphi2, minDphi_ijk);

            const double auxDrap2 = auxDrap + fabs(deltaRap(jeti, jetk));
            minDrap_ijk = std::min(auxDrap2, minDrap_ijk);
          }
        }
      } //end loop over 4 leading jets


      // Combined mass of 4 leading jets
      const double Mjjjj = sum_alljets.mass();

      // Sum of central jets pT
      double sumpt_twojets_cent = 0; // Scalar sum pt of central jets, with different rapidity gaps
      for (size_t ijet = 0; ijet < 4; ++ijet) {
        if (ijet == maxY_i || ijet == maxY_j) continue; // these are the forward jets
        sumpt_twojets_cent += leadingJetskt4.at(ijet).pT();
      }


      // Fill histos
      // Mass and pT cuts in which the analysis tables are split; values are in GeV and cuts are inclusive
      const double m4jcuts[4]   = {500, 1000, 1500, 2000};
      const double pt1cutA[4]   = {100,  400,  700, 1000};
      const double pt1cutB[4]   = {100,  250,  400,  550};
      const double rapGapCut[4] = {1, 2, 3, 4};
      const double weight = event.weight();

      _h["pt1"]->fill(jet1.pt(), weight);
      _h["pt2"]->fill(jet2.pt(), weight);
      _h["pt3"]->fill(jet3.pt(), weight);
      _h["pt4"]->fill(jet4.pt(), weight);
      _h["HT"] ->fill(HT,        weight);
      _h["M4j"]->fill(Mjjjj,     weight);

      for (size_t i_cut = 0; i_cut < 4; ++i_cut) {
        const string icutstr = to_str(i_cut);

        if (Mjjjj > m4jcuts[i_cut])
          _h["M2jratio_"+icutstr]->fill( Mjj/Mjjjj , weight);

        if (jet1.pT() > pt1cutA[i_cut]) {
          _h["dPhiMin2j_"+icutstr]->fill(minDphi_ij , weight);
          _h["dPhiMin3j_"+icutstr]->fill(minDphi_ijk, weight);
          _h["dYMin2j_"+icutstr]->fill(minDrap_ij , weight);
          _h["dYMin3j_"+icutstr]->fill(minDrap_ijk, weight);
        }

        if (jet1.pt() > pt1cutB[i_cut]) {
          _h["dYMax2j_"+icutstr]->fill( maxDrap_ij , weight);
          for (size_t yy = 0; yy < 4; ++yy) {
            if (maxDrap_ij > rapGapCut[yy])
              _h["sumPtCent_"+to_str(yy)+icutstr]->fill(sumpt_twojets_cent, weight);
          }
        }

      } //end loop over pt/m4j cuts

    }



    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = (crossSection()/femtobarn) / sumOfWeights();
      /// @todo Migrate to C++11 range-for loop
      for (map<string, Histo1DPtr>::iterator hit = _h.begin(); hit != _h.end(); ++hit) {
        scale(hit->second, sf);
      }
    }

    //@


  private:

    /// @name Histograms
    //@{
    map<string, Histo1DPtr> _h;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1394679);

}
