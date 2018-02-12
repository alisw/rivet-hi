// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/PrimaryHadrons.hh"
#include "Rivet/Projections/HeavyHadrons.hh"

namespace Rivet {




  class MC_HFJETS : public Analysis {
  public:

    /// Constructor
    MC_HFJETS()
      : Analysis("MC_HFJETS")
    {    }


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FastJets fj(FinalState(-5, 5), FastJets::ANTIKT, 0.6);
      fj.useInvisibles();
      declare(fj, "Jets");
      declare(HeavyHadrons(Cuts::abseta < 5 && Cuts::pT > 500*MeV), "BCHadrons");

      _h_ptCJetLead = bookHisto1D("ptCJetLead", linspace(5, 0, 20, false) + logspace(25, 20, 200));
      _h_ptCHadrLead = bookHisto1D("ptCHadrLead", linspace(5, 0, 10, false) + logspace(25, 10, 200));
      _h_ptFracC = bookHisto1D("ptfracC", 50, 0, 1.5);
      _h_eFracC = bookHisto1D("efracC", 50, 0, 1.5);

      _h_ptBJetLead = bookHisto1D("ptBJetLead", linspace(5, 0, 20, false) + logspace(25, 20, 200));
      _h_ptBHadrLead = bookHisto1D("ptBHadrLead", linspace(5, 0, 10, false) + logspace(25, 10, 200));
      _h_ptFracB = bookHisto1D("ptfracB", 50, 0, 1.5);
      _h_eFracB = bookHisto1D("efracB", 50, 0, 1.5);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Get jets and heavy hadrons
      const Jets& jets = apply<JetAlg>(event, "Jets").jetsByPt();
      const Particles bhadrons = sortByPt(apply<HeavyHadrons>(event, "BCHadrons").bHadrons());
      const Particles chadrons = sortByPt(apply<HeavyHadrons>(event, "BCHadrons").cHadrons());
      MSG_DEBUG("# b hadrons = " << bhadrons.size() << ", # c hadrons = " << chadrons.size());

      // Max HF hadron--jet axis dR to be regarded as a jet tag
      const double MAX_DR = 0.3;

      // Tag the leading b and c jets with a deltaR < 0.3 match
      // b-tagged jet are excluded from also being considered as c-tagged
      /// @todo Do this again with the ghost match?
      MSG_DEBUG("Getting b/c-tags");
      bool gotLeadingB = false, gotLeadingC = false;;
      foreach (const Jet& j, jets) {
        if (!gotLeadingB) {
          FourMomentum leadBJet, leadBHadr;
          double dRmin = MAX_DR;
          foreach (const Particle& b, bhadrons) {
            const double dRcand = min(dRmin, deltaR(j, b));
            if (dRcand < dRmin) {
              dRmin = dRcand;
              leadBJet = j.momentum();
              leadBHadr = b.momentum();
              MSG_DEBUG("New closest b-hadron jet tag candidate: dR = " << dRmin
                        << " for jet pT = " << j.pT()/GeV << " GeV, "
                        << " b hadron pT = " << b.pT()/GeV << " GeV, PID = " << b.pid());
            }
          }
          if (dRmin < MAX_DR) {
            // A jet has been tagged, so fill the histos and break the loop
            _h_ptBJetLead->fill(leadBJet.pT()/GeV, weight);
            _h_ptBHadrLead->fill(leadBHadr.pT()/GeV, weight);
            _h_ptFracB->fill(leadBHadr.pT() / leadBJet.pT(), weight);
            _h_eFracB->fill(leadBHadr.E() / leadBJet.E(), weight);
            gotLeadingB = true;
            continue; // escape this loop iteration so the same jet isn't c-tagged
          }
        }
        if (!gotLeadingC) {
          FourMomentum leadCJet, leadCHadr;
          double dRmin = MAX_DR;
          foreach (const Particle& c, chadrons) {
            const double dRcand = min(dRmin, deltaR(j, c));
            if (dRcand < dRmin) {
              dRmin = dRcand;
              leadCJet = j.momentum();
              leadCHadr = c.momentum();
              MSG_DEBUG("New closest c-hadron jet tag candidate: dR = " << dRmin
                        << " for jet pT = " << j.pT()/GeV << " GeV, "
                        << " c hadron pT = " << c.pT()/GeV << " GeV, PID = " << c.pid());
            }
          }
          if (dRmin < MAX_DR) {
            // A jet has been tagged, so fill the histos and break the loop
            _h_ptCJetLead->fill(leadCJet.pT()/GeV, weight);
            _h_ptCHadrLead->fill(leadCHadr.pT()/GeV, weight);
            _h_ptFracC->fill(leadCHadr.pT() / leadCJet.pT(), weight);
            _h_eFracC->fill(leadCHadr.E() / leadCJet.E(), weight);
            gotLeadingB = true;
          }
        }
        // If we've found both a leading b and a leading c jet, break the loop over jets
        if (gotLeadingB && gotLeadingC) break;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_ptCJetLead);
      normalize(_h_ptCHadrLead);
      normalize(_h_ptFracC);
      normalize(_h_eFracC);
      normalize(_h_ptBJetLead);
      normalize(_h_ptBHadrLead);
      normalize(_h_ptFracB);
      normalize(_h_eFracB);
    }

    //@}

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_ptCJetLead, _h_ptCHadrLead, _h_ptFracC, _h_eFracC;
    Histo1DPtr _h_ptBJetLead, _h_ptBHadrLead, _h_ptFracB, _h_eFracB;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HFJETS);

}
