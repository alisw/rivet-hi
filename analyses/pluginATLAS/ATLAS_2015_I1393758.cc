#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  class ATLAS_2015_I1393758 : public Analysis {

  public:
  
    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1393758);
  
  public:

    void init() {

      addProjection(FastJets(FinalState(), FastJets::ANTIKT, 0.4), "Jets");

      forward_kappa3 = bookProfile1D(1, 1, 1);
      forward_kappa5 = bookProfile1D(2, 1, 1);
      forward_kappa7 = bookProfile1D(3, 1, 1);

      central_kappa3 = bookProfile1D(4, 1, 1);
      central_kappa5 = bookProfile1D(5, 1, 1);
      central_kappa7 = bookProfile1D(6, 1, 1);

      forwardRMS_kappa3 = bookScatter2D("d07-x01-y01", true);
      forwardRMS_kappa5 = bookScatter2D("d08-x01-y01", true);
      forwardRMS_kappa7 = bookScatter2D("d09-x01-y01", true);

      centralRMS_kappa3 = bookScatter2D("d10-x01-y01", true);
      centralRMS_kappa5 = bookScatter2D("d11-x01-y01", true);
      centralRMS_kappa7 = bookScatter2D("d12-x01-y01", true);

    }
  
    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();
      
      Jets m_goodJets = applyProjection<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.1);

      if (m_goodJets.size() < 2)        vetoEvent;
      if (m_goodJets[0].pT() < 50*GeV)  vetoEvent;
      if (m_goodJets[1].pT() < 50*GeV)  vetoEvent;
      if (fabs(1.0 - m_goodJets[0].pT()/m_goodJets[1].pT()) > 0.5)  vetoEvent;

      bool check = m_goodJets[0].abseta() < m_goodJets[1].abseta();
      int pos_f = int(check);
      int pos_c = int(!check);

      double kappa3_f   = CalculateJetCharge(m_goodJets[pos_f], 0.3, 0.5, 1.8);
      double kappa5_f   = CalculateJetCharge(m_goodJets[pos_f], 0.5, 0.5, 1.2);
      double kappa7_f   = CalculateJetCharge(m_goodJets[pos_f], 0.7, 0.5, 0.9);
      double pT_f = m_goodJets[pos_f].pT();
      
      double kappa3_c = CalculateJetCharge(m_goodJets[pos_c], 0.3, 0.5, 1.8);
      double kappa5_c   = CalculateJetCharge(m_goodJets[pos_c], 0.5, 0.5, 1.2);
      double kappa7_c   = CalculateJetCharge(m_goodJets[pos_c], 0.7, 0.5, 0.9);
      double pT_c = m_goodJets[pos_c].pT();

      forward_kappa3->fill(pT_f, kappa3_f, weight);
      forward_kappa5->fill(pT_f, kappa5_f, weight);
      forward_kappa7->fill(pT_f, kappa7_f, weight);

      central_kappa3->fill(pT_c, kappa3_c, weight);
      central_kappa5->fill(pT_c, kappa5_c, weight);
      central_kappa7->fill(pT_c, kappa7_c, weight);
    }

    double CalculateJetCharge(Jet& jet, double kappa=0.5, double pTcut=0.5, double Qmax=1.2) {
      double PTkap = pow(jet.momentum().pT(),kappa);
      double jetcharge = 0.;
      foreach (const Particle& p, jet.particles()) {
        if (p.pT() < pTcut)  continue;
        if (p.threeCharge()) jetcharge += (p.threeCharge()/3.)*pow(p.pT(),kappa)/PTkap;
      }
      //Overflow and underflow
      if (jetcharge > Qmax) jetcharge = Qmax*0.9999;
      if (jetcharge < -Qmax) jetcharge = -Qmax*0.9999;
      return jetcharge;
    }
  
    /// Normalise histograms etc., after the run
    void finalize() {
    
      if (numEvents() > 2) {
        for (unsigned int i = 0; i < forward_kappa3->numBins(); ++i) {
	        double stdv_fkappa3 = forward_kappa3->bin(i).numEntries() > 1? forward_kappa3->bin(i).stdDev() : 0.0;
	        //See Eq. 3 for the factor of two: https://web.eecs.umich.edu/~fessler/papers/files/tr/stderr.pdf
          double yerr_fkappa3  = safediv(sqrt(forward_kappa3->bin(i).sumW2()), 2.*forward_kappa3->bin(i).sumW());
	        forwardRMS_kappa3->point(i).setY(stdv_fkappa3, yerr_fkappa3);

          double stdv_fkappa5 = forward_kappa5->bin(i).numEntries() > 1? forward_kappa5->bin(i).stdDev() : 0.0;
          double yerr_fkappa5  = safediv(sqrt(forward_kappa5->bin(i).sumW2()), 2.*forward_kappa5->bin(i).sumW());
          forwardRMS_kappa5->point(i).setY(stdv_fkappa5, yerr_fkappa5);

          double stdv_fkappa7 = forward_kappa7->bin(i).numEntries() > 1? forward_kappa7->bin(i).stdDev() : 0.0;
          double yerr_fkappa7  = safediv(sqrt(forward_kappa7->bin(i).sumW2()), 2.*forward_kappa7->bin(i).sumW());
          forwardRMS_kappa7->point(i).setY(stdv_fkappa7, yerr_fkappa7);

          double stdv_ckappa3 = central_kappa3->bin(i).numEntries() > 1? central_kappa3->bin(i).stdDev() : 0.0;
          double yerr_ckappa3  = safediv(sqrt(central_kappa3->bin(i).sumW2()), 2.*central_kappa3->bin(i).sumW());
          centralRMS_kappa3->point(i).setY(stdv_ckappa3, yerr_ckappa3);

          double stdv_ckappa5 = central_kappa5->bin(i).numEntries() > 1? central_kappa5->bin(i).stdDev() : 0.0;
          double yerr_ckappa5  = safediv(sqrt(central_kappa5->bin(i).sumW2()), 2.*central_kappa5->bin(i).sumW());
          centralRMS_kappa5->point(i).setY(stdv_ckappa5, yerr_ckappa5);

          double stdv_ckappa7 = central_kappa7->bin(i).numEntries() > 1? central_kappa7->bin(i).stdDev() : 0.0;
          double yerr_ckappa7  = safediv(sqrt(central_kappa7->bin(i).sumW2()), 2.*central_kappa7->bin(i).sumW());
          centralRMS_kappa7->point(i).setY(stdv_ckappa7, yerr_ckappa7);

        }
      }

    }


  private:
  
    Profile1DPtr forward_kappa3;
    Profile1DPtr forward_kappa5;
    Profile1DPtr forward_kappa7;

    Profile1DPtr central_kappa3;
    Profile1DPtr central_kappa5;
    Profile1DPtr central_kappa7;

    Scatter2DPtr forwardRMS_kappa3;
    Scatter2DPtr forwardRMS_kappa5;
    Scatter2DPtr forwardRMS_kappa7;

    Scatter2DPtr centralRMS_kappa3;
    Scatter2DPtr centralRMS_kappa5;
    Scatter2DPtr centralRMS_kappa7;
  
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1393758);
}
