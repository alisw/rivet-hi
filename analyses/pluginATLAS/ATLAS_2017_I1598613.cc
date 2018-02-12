// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"

namespace Rivet {


  class ATLAS_2017_I1598613 : public Analysis {
  public:

    /// Constructor
    ATLAS_2017_I1598613(string name="ATLAS_2017_I1598613")
      : Analysis(name)
    {
      //default to the 3-muon mode:
      _mode = 0;
    }


    struct HistoHandler {
      Histo1DPtr histo;
      Scatter2DPtr scatter;
      unsigned int d, x, y;

      HistoHandler() {}

      void fill(double value, double weight) {
        histo->fill(value, weight);
      }
    };


    /// Book histograms and initialise projections before the run
    void init() {

      // Get the particles needed for each running mode:
      if (_mode == 0) {
        // Get photons to dress leptons
        FinalState photons(Cuts::abspid == PID::PHOTON);
        FinalState muons(Cuts::abspid == PID::MUON);
        Cut eta_lep = Cuts::abseta < 2.5;
        DressedLeptons dressedmuons(photons, muons, 0.1, eta_lep && Cuts::pT >= 6*GeV, true);
        declare(dressedmuons, "dressedmuons");
      } else {
        declare(HeavyHadrons(Cuts::absrap < 2.4 && Cuts::pT > 15.5*GeV), "BHadrons");
      }

      //Book the histograms:
      _h["dPhi"]      = bookHandler(1,  1, 1);
      _h["dy"]        = bookHandler(1,  2, 1);
      _h["yboost"]    = bookHandler(1,  3, 1);
      _h["dR"]        = bookHandler(1,  4, 1);
      _h["lowpT_dR"]  = bookHandler(1,  5, 1);
      _h["highpT_dR"] = bookHandler(1,  6, 1);
      _h["M"]         = bookHandler(1,  7, 1);
      _h["pT"]        = bookHandler(1,  8, 1);
      _h["MopT"]      = bookHandler(1,  9, 1);
      _h["pToM"]      = bookHandler(1, 10, 1);
    }


    HistoHandler bookHandler(unsigned int id_d, unsigned int id_x, unsigned int id_y) {
      HistoHandler dummy;
      dummy.histo = bookHisto1D(id_d, id_x, id_y);
      if (_mode) {
        dummy.scatter = bookScatter2D(id_d, id_x, id_y, true);
        dummy.d = id_d;
        dummy.x = id_x;
        dummy.y = id_y + 1;
      }
      return dummy;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      if (_mode == 1) { // make the 2-b-hadron-level plots
        const Particles& bHadrons = apply<HeavyHadrons>(event, "BHadrons").bHadrons();
        if (bHadrons.size() > 1) {
          sortBy(bHadrons, cmpMomByPt);

          float dphiBB = deltaPhi(bHadrons[0], bHadrons[1]);
          float dRBB = deltaR(bHadrons[0], bHadrons[1], RAPIDITY);
          float dyBB = fabs(bHadrons[0].rapidity() - bHadrons[1].rapidity());
          float yboostBB = 0.5*fabs(bHadrons[0].rapidity() + bHadrons[1].rapidity());
          FourMomentum systemBB = bHadrons[0].momentum() +  bHadrons[1].momentum();
          // Due to the additional particles produced in the decay,
          // the 3 muons carry only a fraction of the momentum of the b-hadrons,
          // scale down mass and pT to match 3-muon-level more closely
          float MBB = systemBB.mass()/1.75;
          float pTBB = systemBB.pT()/1.75;

          _h["dPhi"].fill(dphiBB, weight);
          _h["dy"].fill(dyBB, weight);
          _h["yboost"].fill(yboostBB, weight);
          _h["dR"].fill(dRBB, weight);
          _h["M"].fill(MBB/GeV, weight);
          _h["pT"].fill(pTBB/GeV, weight);
          _h["MopT"].fill(MBB/pTBB, weight);
          _h["pToM"].fill(pTBB/MBB, weight);

          if (pTBB >= 20*GeV)  _h["highpT_dR"].fill(dRBB, weight);
          else                 _h["lowpT_dR"].fill(dRBB,  weight);
        }
      }


      if (_mode == 0) { // the 3-muon-level analysis

        // First, simply check that we have enough muons
        const vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
        if (muons.size() < 3)  vetoEvent;

        // Not sure if this is going to work, but ..
        vector<DressedLepton> Jpsi_muons, third_muons;
        for (DressedLepton mu : muons) {
          if (mu.constituentLepton().fromBottom() && mu.constituentLepton().hasAncestor(PID::JPSI)) {
            Jpsi_muons.push_back(mu);
          }
          else if (mu.constituentLepton().fromBottom()) {
            third_muons.push_back(mu);
          }
        }

        // Veto events without enough muons:
        if (Jpsi_muons.size() < 2)  vetoEvent;

        // At this point, we must always have a Jpsi. So get its 4-vector:
        FourMomentum Jpsi = Jpsi_muons[0].momentum() + Jpsi_muons[1].momentum();

        // If there is more than one J/psi, take the one closest to PDG mass,
        // and push all the other muons back to the 3rd muon list
        size_t mu1 = 0, mu2 = 1;
        if (Jpsi_muons.size() > 2) {
          for (size_t i = 0; i < Jpsi_muons.size(); ++i) {
            for (size_t j = i; j < Jpsi_muons.size(); ++j) {
              FourMomentum temp = Jpsi_muons[i].momentum() + Jpsi_muons[j].momentum();
              if (fabs(temp.mass() - 3.096) < fabs(Jpsi.mass() - 3.096)) {
                Jpsi = temp;
                mu1 = i;
                mu2 = j;
              }
            }
          }

          for (size_t i = 0; i < Jpsi_muons.size(); ++i) {
            if (i == mu1 || i == mu2)  continue;
            third_muons.push_back(Jpsi_muons[i]);
          }
        }

        // There *is* more than one Jpsi:
        if (Jpsi_muons[mu1].abseta() >= 2.3) vetoEvent;
        if (Jpsi_muons[mu2].abseta() >= 2.3) vetoEvent;

        // We should now have the full list of 3rd muons to consider. Make sure we have one:
        if (third_muons.empty())  vetoEvent;

        // Sort the third muons by pT and pick highest one
        std::sort(third_muons.begin(), third_muons.end(), [](const DressedLepton &l1, const DressedLepton &l2) {
          return (l1.pT() > l2.pT());
        });
        FourMomentum third_mu = third_muons[0].momentum();

        // Finally, make the plots!
        float dphi = deltaPhi(Jpsi, third_mu);
        float dR = deltaR(Jpsi, third_mu, RAPIDITY);
        float dy = fabs(Jpsi.rapidity() - third_mu.rapidity());
        float yboost = 0.5*fabs(Jpsi.rapidity() + third_mu.rapidity());
        FourMomentum system = Jpsi +  third_mu;
        float M = system.mass();
        float pT = system.pT();

        _h["dPhi"].fill(dphi, weight);
        _h["dy"].fill(dy, weight);
        _h["yboost"].fill(yboost, weight);
        _h["dR"].fill(dR, weight);
        if (pT >= 20*GeV)  _h["highpT_dR"].fill(dR, weight);
        else  _h["lowpT_dR"].fill(dR, weight);

        _h["M"].fill(M, weight);
        _h["pT"].fill(pT, weight);
        _h["MopT"].fill(M/pT, weight);
        _h["pToM"].fill(pT/M, weight);

      } //< end of muon analysis.
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (map<string, HistoHandler>::iterator hit = _h.begin(); hit != _h.end(); ++hit) {
        normalize(hit->second.histo);
        if (_mode == 1) {
          applyTransferFnAndNorm(hit->second);
          removeAnalysisObject(hit->second.histo);
        }
      }
    }


    void applyTransferFnAndNorm(HistoHandler handler) { ///< @todo Pass as const reference?
      // Load transfer function from reference data file
      const YODA::Scatter2D& myTransferFn = refData(handler.d, handler.x, handler.y);
      double area = 0.0;
      for (size_t i = 0; i < handler.scatter->numPoints(); ++i) {
        const Point2D& f = myTransferFn.point(i);
        Point2D& p = handler.scatter->point(i);
        const HistoBin1D&  b = handler.histo->bin(i);
        double newy;
        try {
          newy = b.height();
        } catch (const Exception&) { // LowStatsError or WeightError
          newy = 0;
        }
        double newey;
        try {
          newey = b.heightErr();
        } catch (const Exception&) { // LowStatsError or WeightError
          newey = 0;
        }
        // apply transfer function here
        newy *= f.y(); newey *= f.y();
        double rp = safediv(newey, newy);
        double rf = safediv(f.yErrAvg(), f.y());
        newey = newy * sqrt(rp*rp + rf*rf);
        // set new values
        p.setY(newy);
        p.setYErrMinus(newey);
        p.setYErrPlus(newey);
        area += newy * (p.xMax() - p.xMin());
      }
      if (area > 0.) { // normalise to unity
        for (size_t i = 0; i < handler.scatter->numPoints(); ++i)
          handler.scatter->point(i).scaleY(1.0 / area);
      }
    }


  protected:

    /// Analysis-mode switch
    int _mode;

    /// Histograms
    map<string, HistoHandler> _h;

  };



  /// Specialised subclass for the BB analysis
  struct ATLAS_2017_I1598613_BB : public ATLAS_2017_I1598613 {
    ATLAS_2017_I1598613_BB() : ATLAS_2017_I1598613("ATLAS_2017_I1598613_BB") { _mode = 1; }
  };


  // Hooks for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1598613);
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1598613_BB);


}
