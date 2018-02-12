// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/HadronicFinalState.hh"
#include "Rivet/Jet.hh"
#include "Rivet/Tools/BinnedHistogram.hh"
#include "fastjet/JetDefinition.hh"

namespace fastjet {

class P_scheme : public JetDefinition::Recombiner {
 public:
  std::string description() const {return "";}
  void recombine(const PseudoJet & pa, const PseudoJet & pb,
                         PseudoJet & pab) const {
    PseudoJet tmp = pa + pb;
    double E = sqrt(tmp.px()*tmp.px() + tmp.py()*tmp.py() + tmp.pz()*tmp.pz());
    pab.reset_momentum(tmp.px(), tmp.py(), tmp.pz(), E);
  }
  void preprocess(PseudoJet & p) const {
    double E = sqrt(p.px()*p.px() + p.py()*p.py() + p.pz()*p.pz());
    p.reset_momentum(p.px(), p.py(), p.pz(), E);
  }
  ~P_scheme() { }
};

}

namespace Rivet {

  class OPAL_2004_I631361_qq : public Analysis {
  public:

    /// Constructor
    OPAL_2004_I631361_qq()
      : Analysis("OPAL_2004_I631361_qq"), _sumWEbin(7,0.)
    {     }

    /// @name Analysis methods
    //@{
    int getEbin(double E_glue) {
      int ih = -1;
      if (inRange(E_glue/GeV, 5.0, 5.5)) {
        ih = 0;
      } else if (inRange(E_glue/GeV, 5.5, 6.5)) {
        ih = 1;
      } else if (inRange(E_glue/GeV, 6.5, 7.5)) {
        ih = 2;
      } else if (inRange(E_glue/GeV, 7.5, 9.5)) {
        ih = 3;
      } else if (inRange(E_glue/GeV, 9.5, 13.0)) {
        ih = 4;
      } else if (inRange(E_glue/GeV, 13.0, 16.0)) {
        ih = 5;
      } else if (inRange(E_glue/GeV, 16.0, 20.0)) {
        ih = 6;
      }
      assert(ih >= 0);
      return ih;
    }

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      addProjection(fs, "FS");
      addProjection(HadronicFinalState(fs), "HFS");

      const ChargedFinalState cfs;
      addProjection(cfs, "CFS");
      addProjection(HadronicFinalState(cfs), "HCFS");

      _h_chMult.addHistogram(5.0, 5.5, bookHisto1D(1,1,1));
      _h_chMult.addHistogram(5.5, 6.5, bookHisto1D(1,1,2));
      _h_chMult.addHistogram(6.5, 7.5, bookHisto1D(1,1,3));
      _h_chMult.addHistogram(7.5, 9.5, bookHisto1D(2,1,1));
      _h_chMult.addHistogram(9.5, 13.0, bookHisto1D(2,1,2));
      _h_chMult.addHistogram(13.0, 16.0, bookHisto1D(3,1,1));
      _h_chMult.addHistogram(16.0, 20.0, bookHisto1D(3,1,2));

      _h_chFragFunc.addHistogram(13.0, 16.0, bookHisto1D(5,1,1));
      _h_chFragFunc.addHistogram(16.0, 20.0, bookHisto1D(5,1,2));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // cut on the number of charged particles
      const Particles& chParticles = applyProjection<FinalState>(event, "CFS").particles();
      if(chParticles.size() < 5) vetoEvent;  
      // cluster the jets
      const Particles& particles = applyProjection<FinalState>(event, "FS").particles();
      fastjet::JetDefinition ee_kt_def(fastjet::ee_kt_algorithm, &p_scheme);
      PseudoJets pParticles;
      foreach(Particle p, particles) {
        PseudoJet temp = p.pseudojet();
        if(p.fromBottom()) {
          temp.set_user_index(5);
        }
        pParticles.push_back(temp);
      }
      fastjet::ClusterSequence cluster(pParticles, ee_kt_def);
      // rescale energys to just keep the directions of the jets
      // and keep track of b tags
      PseudoJets pJets = sorted_by_E(cluster.exclusive_jets_up_to(3));
      if(pJets.size() < 3) vetoEvent;
      array<Vector3, 3> dirs;
      for(int i=0; i<3; i++) {
        dirs[i] = Vector3(pJets[i].px(),pJets[i].py(),pJets[i].pz()).unit();
      }
      array<bool, 3> bTagged;
      Jets jets;
      for(int i=0; i<3; i++) {
        double Ejet = sqrtS()*sin(angle(dirs[(i+1)%3],dirs[(i+2)%3])) /
          (sin(angle(dirs[i],dirs[(i+1)%3])) + sin(angle(dirs[i],dirs[(i+2)%3])) + sin(angle(dirs[(i+1)%3],dirs[(i+2)%3])));
        jets.push_back(FourMomentum(Ejet,Ejet*dirs[i].x(),Ejet*dirs[i].y(),Ejet*dirs[i].z()));
        bTagged[i] = false;
        foreach(PseudoJet particle, pJets[i].constituents()) {
          if(particle.user_index() > 1 and !bTagged[i]) {
            bTagged[i] = true;
          }
        }
      }

      int QUARK1 = 0, QUARK2 = 1, GLUON = 2; 
      
      if(jets[QUARK2].E() > jets[QUARK1].E()) swap(QUARK1, QUARK2);
      if(jets[GLUON].E() > jets[QUARK1].E())  swap(QUARK1,  GLUON);
      if(!bTagged[QUARK2]) {
        if(!bTagged[GLUON]) vetoEvent;
        else swap(QUARK2, GLUON);
      }
      if(bTagged[GLUON]) vetoEvent;

      // exclude collinear or soft jets
      double k1 = jets[QUARK1].E()*min(angle(jets[QUARK1].momentum(),jets[QUARK2].momentum()),
				       angle(jets[QUARK1].momentum(),jets[GLUON].momentum())); 
      double k2 = jets[QUARK2].E()*min(angle(jets[QUARK2].momentum(),jets[QUARK1].momentum()),
				       angle(jets[QUARK2].momentum(),jets[GLUON].momentum()));
      if(k1<8.0*GeV || k2<8.0*GeV) vetoEvent;

      double sqg = (jets[QUARK1].momentum()+jets[GLUON].momentum()).mass2();
      double sgq = (jets[QUARK2].momentum()+jets[GLUON].momentum()).mass2();
      double s = (jets[QUARK1].momentum()+jets[QUARK2].momentum()+jets[GLUON].momentum()).mass2();

      double Eg = 0.5*sqrt(sqg*sgq/s);

      if(Eg < 5.0 || Eg > 20.0) { vetoEvent; }
      else if(Eg > 9.5) {
        //requirements for experimental reconstructability raise as energy raises
        if(!bTagged[QUARK1]) {
          vetoEvent;
        }
      }

      // all cuts applied, increment sum of weights
      const double weight = event.weight();
      _sumWEbin[getEbin(Eg)] += weight;


      // transform to frame with event in y-z and glue jet in z direction
      Matrix3 glueTOz(jets[GLUON].momentum().vector3(), Vector3(0,0,1));
      Vector3 transQuark = glueTOz*jets[QUARK2].momentum().vector3();
      Matrix3 quarksTOyz(Vector3(transQuark.x(), transQuark.y(), 0), Vector3(0,1,0));

      // work out transformation to symmetric frame
      array<double, 3> x_cm;
      array<double, 3> x_cm_y;
      array<double, 3> x_cm_z;
      array<double, 3> x_pr;
      for(int i=0; i<3; i++) {
        x_cm[i] = 2*jets[i].E()/sqrt(s);
        Vector3 p_transf = quarksTOyz*glueTOz*jets[i].p3();
        x_cm_y[i] = 2*p_transf.y()/sqrt(s);
        x_cm_z[i] = 2*p_transf.z()/sqrt(s);
      }
      x_pr[GLUON] = sqrt(4*(1-x_cm[QUARK1])*(1-x_cm[QUARK2])/(3+x_cm[GLUON]));
      x_pr[QUARK1] = x_pr[GLUON]/(1-x_cm[QUARK1]);
      x_pr[QUARK2] = x_pr[GLUON]/(1-x_cm[QUARK2]);
      double gamma = (x_pr[QUARK1] + x_pr[GLUON] + x_pr[QUARK2])/2;
      double beta_z = x_pr[GLUON]/(gamma*x_cm[GLUON]) - 1;
      double beta_y = (x_pr[QUARK2]/gamma - x_cm[QUARK2] - beta_z*x_cm_z[QUARK2])/x_cm_y[QUARK2];

      LorentzTransform toSymmetric = LorentzTransform::mkObjTransformFromBeta(Vector3(0.,beta_y,beta_z)).
	postMult(quarksTOyz*glueTOz);
      
      FourMomentum transGlue = toSymmetric.transform(jets[GLUON].momentum());
      double cutAngle = angle(toSymmetric.transform(jets[QUARK2].momentum()), transGlue)/2;

      int nCh = 0;
      foreach(const Particle& chP, chParticles ) {
        FourMomentum pSymmFrame = toSymmetric.transform(FourMomentum(chP.p3().mod(), chP.px(), chP.py(), chP.pz()));
	if(angle(pSymmFrame, transGlue) < cutAngle) {
          _h_chFragFunc.fill(Eg, pSymmFrame.E()*sin(cutAngle)/Eg, weight);
          nCh++;
        }
      }
      _h_chMult.fill(Eg, nCh, weight);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      foreach(Histo1DPtr hist, _h_chMult.getHistograms()) {
        normalize(hist);
      }
      for(int i=0; i<2; i++) {
        if(!isZero(_sumWEbin[i+5])) {
          scale(_h_chFragFunc.getHistograms()[i], 1./_sumWEbin[i+5]);
        }
      }
    }

    //@}


  private:

    // Data members like post-cuts event weight counters go here
    vector<double> _sumWEbin;

    // p scheme jet definition
    fastjet::P_scheme p_scheme;

    /// @name Histograms
    //@{
    BinnedHistogram<double> _h_chMult;
    BinnedHistogram<double> _h_chFragFunc;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(OPAL_2004_I631361_qq);
}
