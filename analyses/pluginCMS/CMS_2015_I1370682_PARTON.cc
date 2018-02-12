#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  class CMS_2015_I1370682_PARTON : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1370682_PARTON);


    /// Book projections and histograms
    void init() {
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicPartonTops");
      declare(PartonicTops(PartonicTops::HADRONIC),    "HadronicPartonTops");

      _hSL_topPt         = bookHisto1D(15, 1, 1);
      _hSL_topPtTtbarSys = bookHisto1D(16, 1, 1);
      _hSL_topY          = bookHisto1D(17, 1, 1);
      _hSL_ttbarDelPhi   = bookHisto1D(18, 1, 1);
      _hSL_topPtLead     = bookHisto1D(19, 1, 1);
      _hSL_topPtSubLead  = bookHisto1D(20, 1, 1);
      _hSL_ttbarPt       = bookHisto1D(21, 1, 1);
      _hSL_ttbarY        = bookHisto1D(22, 1, 1);
      _hSL_ttbarMass     = bookHisto1D(23, 1, 1);

      _hDL_topPt           = bookHisto1D(24, 1, 1);
      _hDL_topPtTtbarSys   = bookHisto1D(25, 1, 1);
      _hDL_topY            = bookHisto1D(26, 1, 1);
      _hDL_ttbarDelPhi     = bookHisto1D(27, 1, 1);
      _hDL_topPtLead       = bookHisto1D(28, 1, 1);
      _hDL_topPtSubLead    = bookHisto1D(29, 1, 1);
      _hDL_ttbarPt         = bookHisto1D(30, 1, 1);
      _hDL_ttbarY          = bookHisto1D(31, 1, 1);
      _hDL_ttbarMass       = bookHisto1D(32, 1, 1);
      }


      void analyze(const Event& event) {

        // Do the analysis only for the ttbar full leptonic or semileptonic channel, without tau decay
        const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
        const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
        const bool isSemilepton = (leptonicpartontops.size() == 1 && hadronicpartontops.size() == 1);
        const bool isDilepton = (leptonicpartontops.size() == 2 && hadronicpartontops.size() == 0);
        if (!isSemilepton && !isDilepton) vetoEvent;

        // Parton level at full phase space
        // Fill top quarks defined in the parton level, full phase space
        const FourMomentum t1P4 = leptonicpartontops[0];
        const FourMomentum t2P4 = isSemilepton ? hadronicpartontops[0] : leptonicpartontops[1];
        const double t1Pt = t1P4.pT(), t2Pt = t2P4.pT();
        const FourMomentum ttbarP4 = t1P4 + t2P4;
        const FourMomentum t1P4AtCM = LorentzTransform::mkFrameTransformFromBeta(ttbarP4.betaVec()).transform(t1P4);
        const double dPhi = deltaPhi(t1P4.phi(), t2P4.phi());

        const double weight = event.weight();
        if (isSemilepton) {
          _hSL_topPt->fill(t1Pt, weight);
          _hSL_topPt->fill(t2Pt, weight);
          _hSL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hSL_topY->fill(t1P4.rapidity(), weight);
          _hSL_topY->fill(t2P4.rapidity(), weight);
          _hSL_ttbarDelPhi->fill(dPhi, weight);
          _hSL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hSL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hSL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hSL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hSL_ttbarMass->fill(ttbarP4.mass(), weight);
        } else { // if (isDilepton) {
          _hDL_topPt->fill(t1Pt, weight);
          _hDL_topPt->fill(t2Pt, weight);
          _hDL_topPtTtbarSys->fill(t1P4AtCM.pT(), weight);
          _hDL_topY->fill(t1P4.rapidity(), weight);
          _hDL_topY->fill(t2P4.rapidity(), weight);
          _hDL_ttbarDelPhi->fill(dPhi, weight);
          _hDL_topPtLead->fill(std::max(t1Pt, t2Pt), weight);
          _hDL_topPtSubLead->fill(std::min(t1Pt, t2Pt), weight);
          _hDL_ttbarPt->fill(ttbarP4.pT(), weight);
          _hDL_ttbarY->fill(ttbarP4.rapidity(), weight);
          _hDL_ttbarMass->fill(ttbarP4.mass(), weight);
        }
      }


      void finalize() {
        normalize({_hSL_topPt, _hSL_topPtTtbarSys, _hSL_topY, _hSL_ttbarDelPhi, _hSL_topPtLead,
              _hSL_topPtSubLead, _hSL_ttbarPt, _hSL_ttbarY, _hSL_ttbarMass});
        normalize({_hDL_topPt, _hDL_topPtTtbarSys, _hDL_topY, _hDL_ttbarDelPhi, _hDL_topPtLead,
              _hDL_topPtSubLead, _hDL_ttbarPt, _hDL_ttbarY, _hDL_ttbarMass});
      }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _hSL_topPt, _hSL_topPtTtbarSys, _hSL_topY, _hSL_ttbarDelPhi, _hSL_topPtLead,
      _hSL_topPtSubLead, _hSL_ttbarPt, _hSL_ttbarY, _hSL_ttbarMass;
    Histo1DPtr _hDL_topPt, _hDL_topPtTtbarSys, _hDL_topY, _hDL_ttbarDelPhi, _hDL_topPtLead,
      _hDL_topPtSubLead, _hDL_ttbarPt, _hDL_ttbarY, _hDL_ttbarMass;
    //@}

  };


  DECLARE_RIVET_PLUGIN(CMS_2015_I1370682_PARTON);

}
