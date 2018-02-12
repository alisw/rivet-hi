// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief CMS underlying event in leading track events at 7 TeV
  /// @author Paolo Gunnellini (DESY)
  ///
  /// CMS measurement of the underlying event in "leading track" events.
  class CMS_PAS_FSQ_12_020 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_PAS_FSQ_12_020);


    /// Book histograms and initialise projections before the run
    void init() {

      const ChargedFinalState cfs(Cuts::abseta < 0.8 && Cuts::pT > 0.5*GeV);
      addProjection(cfs, "Tracks");

      _NchgPDFden1 = bookProfile1D(7,1,1);
      _NchgPMNden1 = bookProfile1D(6,1,1);
      _NchgPMXden1 = bookProfile1D(5,1,1);

      _PTsumPDFden1 = bookProfile1D(10,1,1);
      _PTsumPMNden1 = bookProfile1D(9,1,1);
      _PTsumPMXden1 = bookProfile1D(8,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {


      // Require at least one track in the event with pT >= 0.5 GeV
      const FinalState& cfs = applyProjection<ChargedFinalState>(event, "Tracks");
      if (cfs.empty()) vetoEvent;
      const Particles trks = cfs.particlesByPt();

      // Identify leading track and its phi and pT
      const Particle p_lead = trks[0];
      const double philead = p_lead.momentum().phi();
      const double ptlead  = p_lead.momentum().pT();

      // Loop over particles and build transverse side variables
      double NchgP1 = 0, NchgP2 = 0, PTsumP1 = 0, PTsumP2 = 0;
      for (const Particle& p : trks) {

        // Region definition -- if not in transverse region, ignore
        const double dphi = mapAngle0To2Pi(p.phi() - philead);
        if (!inRange(dphi, PI/3, 2*PI/3) && !inRange(dphi, 4*PI/3, 5*PI/3)) continue;

        // Transverse region 1
        if (inRange(dphi, PI/3, 2*PI/3)) {
          NchgP1 += 1;
          PTsumP1 += p.pT();
        }
        // Transverse region 2
        else if (inRange(dphi, 4*PI/3, 5*PI/3)) {
          NchgP2 += 1;
          PTsumP2 += p.pT();
        }
      }

      // Calculate total variables
      // const double NchgPtot = (NchgP1 + NchgP2)/2;
      const double NchgPmax = max(NchgP1,NchgP2);
      const double NchgPmin = min(NchgP1,NchgP2);
      // const double PTsumPtot = (PTsumP1 + PTsumP2)/2;
      const double PTsumPmax = max(PTsumP1,PTsumP2);
      const double PTsumPmin = min(PTsumP1,PTsumP2);
      //
      const double PTsumPMXden = PTsumPmax/AREA;
      const double PTsumPMNden = PTsumPmin/AREA;
      const double NchgPMXden = NchgPmax/AREA;
      const double NchgPMNden = NchgPmin/AREA;
      //
      const double NchgPDFden = NchgPMXden - NchgPMNden;
      const double PTsumPDFden = PTsumPMXden - PTsumPMNden;

      // Fill histograms
      const double weight = event.weight();
      _NchgPMXden1->fill(ptlead/GeV, NchgPmax/AREA, weight);
      _NchgPMNden1->fill(ptlead/GeV, NchgPmin/AREA, weight);
      _NchgPDFden1->fill(ptlead/GeV, NchgPDFden, weight);
      _PTsumPMXden1->fill(ptlead/GeV, PTsumPmax/AREA, weight);
      _PTsumPMNden1->fill(ptlead/GeV, PTsumPmin/AREA, weight);
      _PTsumPDFden1->fill(ptlead/GeV, PTsumPDFden, weight);

    }


    /// eta-phi area of the transverse region
    constexpr static double AREA = 2*0.8 * M_PI/3;

    /// Histograms
    Profile1DPtr _NchgPden1, _NchgPMXden1, _NchgPMNden1, _NchgPDFden1, _PTsumPden1, _PTsumPMXden1, _PTsumPMNden1, _PTsumPDFden1;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_PAS_FSQ_12_020);

}
