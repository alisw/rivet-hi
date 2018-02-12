// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  namespace {

    /// Number of event shape variables
    /// @todo Move into the EventShape class
    const int NEVTVAR = 5;

    /// Number of leading jet pT thresholds
    /// @todo Move into the analysis class
    const int NJETPTMN = 5;
    /// Leading jet pT thresholds
    /// @todo Move into the analysis class
    const double LEADINGPTTHRESHOLD[NJETPTMN] = { 110.0, 170.0, 250.0, 320.0, 390.0 };


    // Helpers for event shape calculations in hidden namespace; implementation at bottom of file
    /// @todo Why a class? Improve/remove this junk
    class EventShape {
    public:

      /// Constructor from vectors of four-vectors as input objects in the event to calculate the event shapes
      EventShape(const vector<double>& px_vector, const vector<double>& py_vector, const vector<double>& pz_vector,
                 const vector<double>& e_vector, double eta_central, int irap, int nmn)
        : _object_px(px_vector), _object_py(py_vector), _object_pz(pz_vector),
          _object_e(e_vector), _eta_c(eta_central), _irap(irap), _nmnjet(nmn)
      {   }

      /// @brief Returns the values of the five event shapes
      ///
      /// Event shape indices:
      /// 0. central transverse thrust
      /// 1. central total jet broadening
      /// 2. central total jet mass
      /// 3. central total transverse jet mass
      /// 4. central three-jet resolution threshold
      vector<double> getEventShapes() {
        _calculate(); ///< @todo There should be some test for success/failure!!
        return _event_shapes;
      }

      /// Returns the global thrust axis Nx, Ny, Nz=0
      vector<double> getThrustAxis() {
        _calculate(); ///< @todo There should be some test for success/failure!!
        return _thrust_axis;
      }

      /// Returns the central thrust axis Nx, Ny, Nz=0
      vector<double> getThrustAxisC() {
        _calculate(); ///< @todo There should be some test for success/failure!!
        return _thrust_axis_c;
      }

      // /// @brief Choice of the central region
      // void setEtaC(double eta_central) { _eta_c = eta_central; }

      // // Whether to use the rapidity y (rap==1)  or the pseudorapidity eta (rap==0)
      // void setRapType(int irap) { _irap = irap; }


    private:

      /// Calculate everything
      int _calculate();

      /// Returns the difference in phi between two vectors
      double _delta_phi(double, double);

      /// The Lorentz scalar product
      double _lorentz_sp(const vector<double>&, const vector<double>&);

      // Calculates the three-jet resolutions
      double _three_jet_res(const vector<double>&, const vector<double>&, const vector<double>&, const vector<double>&, int);

      // Calculates the thrust axis and the tau values
      vector<double> _thrust(const vector<double>&, const vector<double>&);


      vector<double> _object_px, _object_py, _object_pz, _object_p;
      vector<double> _object_pt, _object_e, _object_phi, _object_eta;
      vector<double> _event_shapes;
      vector<double> _thrust_axis, _thrust_axis_c;

      double _eta_c;
      int _irap;
      size_t _nmnjet;

    };

  }




  class CMS_2014_I1305624 : public Analysis {
  public:

    /// Constructor
    CMS_2014_I1305624()
      : Analysis("CMS_2014_I1305624")
    {    }


    /// @name Analysis methods

    /// Book histograms and initialise projections before the run
    void init() {
      const FastJets jets(FinalState(Cuts::abseta < 2.6), FastJets::ANTIKT, 0.5);
      declare(jets, "Jets");

      for (int ij=0; ij < NJETPTMN; ij++) {
        _h_thrustc[ij] = bookHisto1D(1, 1, ij+1);
        _h_broadt[ij] = bookHisto1D(1, 2, ij+1);
        _h_tot3dmass[ij] = bookHisto1D(1, 3, ij+1);
        _h_tottrnsmass[ij] = bookHisto1D(1, 4, ij+1);
        _h_y23c[ij] = bookHisto1D(1, 5, ij+1);
        //
        _alow1[ij] = _h_thrustc[ij]->xMin();
        _alow2[ij] = _h_broadt[ij]->xMin();
        _alow3[ij] = _h_tot3dmass[ij]->xMin();
        _alow4[ij] = _h_tottrnsmass[ij]->xMin();
        _alow5[ij] = _h_y23c[ij]->xMin();
        //
        _ahgh1[ij] = _h_thrustc[ij]->xMax();
        _ahgh2[ij] = _h_broadt[ij]->xMax();
        _ahgh3[ij] = _h_tot3dmass[ij]->xMax();
        _ahgh4[ij] = _h_tottrnsmass[ij]->xMax();
        _ahgh5[ij] = _h_y23c[ij]->xMax();
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(30.0*GeV);
      if (jets.size() < 2) vetoEvent;
      if (jets[0].abseta() > 2.4 || jets[1].abseta() > 2.4) vetoEvent;

      const double leadingpt = jets[0].pT();
      if (leadingpt < 110*GeV) vetoEvent;

      vector<double> jtpx, jtpy, jtpz, jten;
      foreach (const Jet& j, jets) {
        if (j.abseta() < 2.4) {
          jtpx.push_back(j.px());
          jtpy.push_back(j.py());
          jtpz.push_back(j.pz());
          jten.push_back(j.E());
        }
      }

      EventShape eventshape(jtpx, jtpy, jtpz, jten, 2.4, 0, 2);
      const vector<double> eventvar = eventshape.getEventShapes();
      if (eventvar[NEVTVAR] < 0) vetoEvent; // Jets are not only one hemisphere

      const double weight = event.weight();
      for (int ij = NJETPTMN-1; ij >= 0; --ij) {
        if (leadingpt/GeV > LEADINGPTTHRESHOLD[ij]) {
          if (inRange(eventvar[0], _alow1[ij], _ahgh1[ij])) _h_thrustc[ij]->fill(eventvar[0], weight);
          if (inRange(eventvar[2], _alow3[ij], _ahgh3[ij])) _h_tot3dmass[ij]->fill(eventvar[2], weight);
          if (inRange(eventvar[3], _alow4[ij], _ahgh4[ij])) _h_tottrnsmass[ij]->fill(eventvar[3], weight);
          if (eventvar[NEVTVAR] >= 3) {
            if (inRange(eventvar[1], _alow2[ij], _ahgh2[ij])) _h_broadt[ij]->fill(eventvar[1], weight);
            if (inRange(eventvar[4], _alow5[ij], _ahgh5[ij])) _h_y23c[ij]->fill(eventvar[4], weight);
          }
          break;
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (int ij = 0; ij < NJETPTMN; ij++) {
        normalize(_h_thrustc[ij]);
        normalize(_h_broadt[ij]);
        normalize(_h_tot3dmass[ij]);
        normalize(_h_tottrnsmass[ij]);
        normalize(_h_y23c[ij]);
      }
    }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_thrustc[NJETPTMN];
    Histo1DPtr _h_broadt[NJETPTMN];
    Histo1DPtr _h_tot3dmass[NJETPTMN];
    Histo1DPtr _h_tottrnsmass[NJETPTMN];
    Histo1DPtr _h_y23c[NJETPTMN];
    //@}

    // Data members
    double _alow1[NJETPTMN], _alow2[NJETPTMN], _alow3[NJETPTMN], _alow4[NJETPTMN], _alow5[NJETPTMN];
    double _ahgh1[NJETPTMN], _ahgh2[NJETPTMN], _ahgh3[NJETPTMN], _ahgh4[NJETPTMN], _ahgh5[NJETPTMN];

  };


  DECLARE_RIVET_PLUGIN(CMS_2014_I1305624);



  /////////////////////


  namespace {

    // EventShape helper class method implementations:

    int EventShape::_calculate() {
      if (!_event_shapes.empty() && !_thrust_axis.empty() && !_thrust_axis_c.empty())
        return 1; //< return success if this appears to already have been run

      const size_t length = (size_t) _object_px.size();

      if (((size_t) _object_py.size() != length) ||
          ((size_t) _object_pz.size() != length) ||
          ((size_t) _object_e.size() != length)) {
        /// @todo Change to exception or assert
        // cout << "ERROR!!!! Input vectors differ in size! Change that please!" << endl;
        // cout<<"py_size: "<<_object_py.size()<<" ,pz_size: "<<_object_pz.size()
        //     <<" ,px_size: "<<_object_px.size()<<" ,E_size: "<<_object_e.size()<<endl;
        return 0;
      }

      if (!_object_p.empty()) {
        _object_p.clear();
        _object_pt.clear();
        _object_eta.clear();
        _object_phi.clear();
        _event_shapes.clear();
        _thrust_axis.clear();
        _thrust_axis_c.clear();
      }

      for (size_t j = 0; j < length; j++) {
        _object_p.push_back(0.);
        _object_pt.push_back(0.);
        _object_eta.push_back(0.);
        _object_phi.push_back(0.);
      }

      for (int j = 0; j < NEVTVAR; j++) {
        _event_shapes.push_back(-50.);
      }

      _event_shapes.push_back(double(_object_px.size())); //< WTF?

      for (size_t j = 0; j < 3; j++) {
        _thrust_axis.push_back(0.);
        _thrust_axis_c.push_back(0.);
      }

      double theta = 0;

      for (size_t k = 0; k < length; k++) {
        _object_p[k] = sqrt(pow(_object_px[k],2) + pow(_object_py[k],2) + pow(_object_pz[k],2));
        _object_pt[k] = sqrt(pow(_object_px[k],2) + pow(_object_py[k],2));
        if (_object_p[k] > _object_e[k] + 1e-4) {
          /// @todo Change to exception or assert
          // cout << "ERROR!!! object " << k <<" has P = " << _object_p[k]
          //      << " which is bigger than E = " << _object_e[k] <<" "
          //      << _object_px[k] <<" "<< _object_py[k] <<" "
          //      << _object_pz[k] <<" of total length "<< length
          //      << endl;
          return 0;
        }

        //to prevent a division by zero
        if (_irap == 0) {
          if (fabs(_object_pz[k]) > 1e-5) {
            theta = atan(_object_pt[k]/(_object_pz[k]));
          } else {
            theta = M_PI/2;
          }
          if (theta < 0.) theta = theta + M_PI;
          _object_eta[k] = -log(tan(0.5*theta));
        }
        if (_irap == 1) {
          if (_object_pz[k] == _object_e[k]) {
            /// @todo Change to exception or assert
            // cout << "ERROR!!! object "<<k<<" has Pz "<< _object_pz[k] <<" which is equal to E = "<< _object_e[k] <<endl;
            return 0;
          }
          _object_eta[k]=0.5*log((_object_e[k]+_object_pz[k])/(_object_e[k]-_object_pz[k]));
        }
        if (_irap != 0 && _irap != 1) {
          /// @todo Change to exception or assert
          // cout << "ERROR!!!, The choice to use the rapidity y or the pseudorapidity eta is not set correctly! Change that please!" << endl;
          return 0;
        }
        _object_phi[k] = atan2(_object_py[k], _object_px[k]);
      }

      vector<double> object_px_in, object_py_in, object_pz_in, object_pt_in, object_e_in, object_et_in, object_eta_in;
      vector<double> object_px_out, object_py_out, object_pz_out, object_e_out, object_pt_out, object_eta_out;
      if (!object_px_in.empty()) { //< FFS, this is impossible: it's only just been created!
        object_px_in.clear();
        object_py_in.clear();
        object_pz_in.clear();
        object_pt_in.clear();
        object_e_in.clear();
        object_et_in.clear();
        object_eta_in.clear();
        object_px_out.clear();
        object_py_out.clear();
        object_pz_out.clear();
        object_pt_out.clear();
        object_e_out.clear();
        object_eta_out.clear();
      }

      size_t nin = 0;

      for (size_t j = 0; j < length; j++) {
        if (fabs(_object_eta[j]) < _eta_c) {
          object_px_in.push_back(_object_px[j]);
          object_py_in.push_back(_object_py[j]);
          object_pz_in.push_back(_object_pz[j]);
          object_e_in.push_back(_object_e[j]);
          object_pt_in.push_back(sqrt(pow(_object_px[j],2)+pow(_object_py[j],2)));
          object_et_in.push_back(sqrt((pow(_object_e[j],2)*pow(_object_pt[j],2))/(pow(_object_pt[j],2)+pow(_object_pz[j],2))));
          object_eta_in.push_back(_object_eta[j]);
          nin += 1;
      } else {
          object_px_out.push_back(_object_px[j]);
          object_py_out.push_back(_object_py[j]);
          object_pz_out.push_back(_object_pz[j]);
          object_e_out.push_back(_object_e[j]);
          object_pt_out.push_back(sqrt(pow(_object_px[j],2)+pow(_object_py[j],2)));
          object_eta_out.push_back(_object_eta[j]);
        }
      }

      if (object_px_in.size() != nin) {
        /// @todo Change to exception or assert
        cout<<"ERROR!!! wrong dimension of 'in' momenta"<<endl;
        //return 0; ///< @todo Why not do this?
      }
      const size_t nout = length - nin;

      if (nin < _nmnjet) {
        for (int i = 0; i < NEVTVAR; i++) {
          _event_shapes[i] = -50.0;
        }
      }

      _event_shapes[NEVTVAR] = nin;

      if (nin >= _nmnjet) {
        double p_sum_c = 0; //GMA
        double pt_sum_c = 0;
        double eta_cw=0;
        double px_sum_in = 0;
        double py_sum_in = 0;
        for (size_t j = 0; j < nin; j++) {
          pt_sum_c += object_pt_in[j];
          p_sum_c += sqrt(pow(object_pt_in[j],2.) + pow(object_pz_in[j], 2.0)); //GMA
          eta_cw += object_pt_in[j]*object_eta_in[j];
          px_sum_in += object_px_in[j];
          py_sum_in += object_py_in[j];
        }
        eta_cw /= pt_sum_c;

        double expTerm = 0;
        for (size_t j = 0; j < nout; j++) {
          expTerm += object_pt_out[j] * exp(-fabs(object_eta_out[j]-eta_cw));
        }
        expTerm /= pt_sum_c;

        //the central global transverse thrust centrthr is calculated
        double centrthr = 0;
        vector<double> thrust_central = _thrust(object_px_in, object_py_in);

        for (size_t l=0; l<3; l++) _thrust_axis_c[l] = thrust_central[l];
        //the variable which gets resummed is not thrust
        //but tau = 1 - thrust - see calculation
        centrthr = thrust_central[3];
        _event_shapes[0] = centrthr;

        double alpha_c = atan2(_thrust_axis_c[1], _thrust_axis_c[0]);
        //central jet masses
        //define two jet masses in region U and D
        double cenjm_up = 0;
        double cenjm_down= 0;
        double dot_product = 0;

        vector<double> up_sum;
        vector<double> down_sum;
        for (size_t j=0; j<4;j++) {
          up_sum.push_back(0.);
          down_sum.push_back(0.);
        }
        for (size_t i=0;i<nin;i++) {
          dot_product = object_px_in[i] * _thrust_axis_c[0] + object_py_in[i] * _thrust_axis_c[1];
          if (dot_product >= 0) {
            up_sum[0]+=object_px_in[i];
            up_sum[1]+=object_py_in[i];
            up_sum[2]+=object_pz_in[i];
            up_sum[3]+=object_e_in[i];
          } else {
            down_sum[0]+=object_px_in[i];
            down_sum[1]+=object_py_in[i];
            down_sum[2]+=object_pz_in[i];
            down_sum[3]+=object_e_in[i];
          }
        }
        cenjm_up = _lorentz_sp(up_sum, up_sum) / pow(p_sum_c, 2.); //GMA pow(pt_sum_c,2);
        cenjm_down = _lorentz_sp(down_sum, down_sum) / pow(p_sum_c, 2.); //GMA pow(pt_sum_c,2);

        //central total jet mass centotjm
        double centotjm=0;
        centotjm = cenjm_up + cenjm_down;

        _event_shapes[2]=centotjm;

        double centrjm_up=0, centrjm_down=0;
        vector<double> upsum;
        vector<double> downsum;
        for (size_t j = 0; j < 3; j++) {
          upsum.push_back(0.);
          downsum.push_back(0.);
        }
        for (size_t i = 0; i < nin; i++) {
          dot_product = object_px_in[i]*_thrust_axis_c[0]+object_py_in[i]*_thrust_axis_c[1];
          if (dot_product >= 0) {
            upsum[0] += object_px_in[i];
            upsum[1] += object_py_in[i];
            upsum[2] += object_et_in[i];
          } else {
            downsum[0] += object_px_in[i];
            downsum[1] += object_py_in[i];
            downsum[2] += object_et_in[i];
          }
        }
        centrjm_up = _lorentz_sp(upsum, upsum) / pow(pt_sum_c, 2);
        centrjm_down = _lorentz_sp(downsum, downsum) / pow(pt_sum_c, 2);
        double centottrjm = centrjm_up + centrjm_down;

        _event_shapes[3] = centottrjm;

        //central three-jet resolution threshold
        double ceny3=0;
        if (nin < 3) {
          ceny3 = -1.0;
        } else {
          ceny3 = _three_jet_res(object_px_in, object_py_in, object_pz_in, object_e_in, _irap);
        }

        _event_shapes[4] = ceny3;

        //the central jet broadenings in the up and down region
        double cenbroad_up=0;
        double cenbroad_down=0;

        double eta_up=0;
        size_t num_up=0;
        double eta_down =0;
        size_t num_down =0;
        double phi_temp =0;
        double phi_up_aver =0;
        double phi_down_aver =0;
        double pt_sum_up =0;
        double pt_sum_down =0;
        double dot_product_b =0;
        vector<double> phi_up;
        vector<double> phi_down;
        double py_rot =0;
        double px_rot =0;

        for (size_t j = 0; j < 4; j++) {
          up_sum.push_back(0.);
          down_sum.push_back(0.);
        }

        for (size_t i=0;i<nin;i++) {
          dot_product_b =sqrt(object_px_in[i]*_thrust_axis_c[0] + object_py_in[i]*_thrust_axis_c[1]);
          if (dot_product_b>=0){
            pt_sum_up += object_pt_in[i];
            //rotate the coordinate system so that
            //the central thrust axis is e_x
            px_rot = cos(alpha_c)*object_px_in[i]+sin(alpha_c)*object_py_in[i];
            py_rot = - sin(alpha_c)*object_px_in[i]+cos(alpha_c)*object_py_in[i];
            //calculate the eta and phi in the rotated system
            eta_up += object_pt_in[i]*object_eta_in[i];
            phi_temp = atan2(py_rot,px_rot);

            if(phi_temp > M_PI/2){
              phi_temp = phi_temp - M_PI/2;
            }
            if (phi_temp < -M_PI/2){
              phi_temp = phi_temp + M_PI/2;
            }
            phi_up.push_back(phi_temp);
            phi_up_aver += object_pt_in[i]*phi_temp;
            num_up += 1;
          } else {
            eta_down += object_pt_in[i]*object_eta_in[i];
            pt_sum_down += object_pt_in[i];
            px_rot = cos(alpha_c)*object_px_in[i]+sin(alpha_c)*object_py_in[i];
            py_rot = - sin(alpha_c)*object_px_in[i]+cos(alpha_c)*object_py_in[i];
            phi_temp = atan2(py_rot,px_rot);
            if (phi_temp > M_PI/2) {
              //if phi is bigger than pi/2 in the new system calculate
              //the difference to the thrust axis
              phi_temp = M_PI -phi_temp;
            }
            if (phi_temp<-M_PI/2) {
              //if phi is smaller than
              phi_temp = -M_PI-phi_temp;
            }
            phi_down.push_back(phi_temp);
            //calculate the pt-weighted phi
            phi_down_aver += object_pt_in[i]*phi_temp;
            num_down += 1;
          }
        }
        if (num_up!=0){
          eta_up = eta_up/pt_sum_up;
          phi_up_aver = phi_up_aver/pt_sum_up;
        }
        if (num_down!=0) {
          eta_down = eta_down/pt_sum_down;
          phi_down_aver = phi_down_aver/pt_sum_down;
        }

        size_t index_up=0, index_down=0;
        for (size_t i = 0; i < nin; i++) {
          dot_product_b = object_px_in[i]*_thrust_axis_c[0] + object_py_in[i]*_thrust_axis_c[1];
          if (dot_product_b >= 0) {
            //calculate the broadenings of the regions with the rotated system
            //and the pt-weighted average of phi in the rotated system
            cenbroad_up += object_pt_in[i]*sqrt(pow(object_eta_in[i]-eta_up, 2) +
                                                pow(_delta_phi(phi_up[index_up], phi_up_aver), 2));
            index_up += 1;
          } else {
            cenbroad_down += object_pt_in[i]*sqrt(pow(object_eta_in[i]-eta_down, 2)+
                                                  pow(_delta_phi(phi_down[index_down], phi_down_aver), 2));
            index_down += 1;
          }
        }

        if (index_up == 0 || index_down ==0) _event_shapes[NEVTVAR] *= -1;

        cenbroad_up=cenbroad_up/(2*pt_sum_c);
        cenbroad_down=cenbroad_down/(2*pt_sum_c);

        //central total jet broadening
        double centotbroad = 0;
        centotbroad = cenbroad_up + cenbroad_down;

        _event_shapes[1] = centotbroad;

        for (int ij = 0; ij < 5; ij++) {
          if (_event_shapes[ij] < 1.e-20) _event_shapes[ij] = 1.e-20;
          _event_shapes[ij] = log(_event_shapes[ij]);
        }
      }

      return 1;
    }


    double EventShape::_three_jet_res(const vector<double>& in_object_px, const vector<double>& in_object_py, const vector<double>& in_object_pz, const vector<double>& in_object_e, int irap) {

      size_t y3_length = (size_t)in_object_px.size();
      if (((size_t) in_object_py.size()!=y3_length) ||
          ((size_t) in_object_pz.size()!=y3_length) ||
          (in_object_e.size()!=y3_length)) {
        // cout << "ERROR!!!! Input vectors differ in size! Change that please!" << endl;
        // cout<<"py_size: "<<in_object_py.size()<<" ,pz_size: "<<in_object_pz.size()
        //     <<" ,px_size: "<<in_object_px.size()<<" , E_size: "<<in_object_e.size() <<endl;
        return 0.0;
      }

      vector<double> in_object_p, in_object_pt, in_object_eta, in_object_phi;
      if (!in_object_p.empty()) {
        in_object_p.clear();
        in_object_pt.clear();
        in_object_eta.clear();
        in_object_phi.clear();
      }
      for (size_t j = 0; j < y3_length; j++) {
        in_object_p.push_back(0.);
        in_object_pt.push_back(0.);
        in_object_eta.push_back(0.);
        in_object_phi.push_back(0.);
      }
      double theta_y3_1st = 0;
      for (size_t k =0; k<y3_length; k++) {
        in_object_p[k] = sqrt(pow(in_object_px[k],2) + pow(in_object_py[k],2) + pow(in_object_pz[k],2));
        in_object_pt[k] = sqrt(pow(in_object_px[k],2) + pow(in_object_py[k],2));

        //calculates the pseudorapidity to prevent a division by zero
        if (irap == 0) {
          if (fabs(in_object_pz[k]) > 1E-5) {
            theta_y3_1st = atan(in_object_pt[k]/(in_object_pz[k]));
          } else {
            theta_y3_1st = M_PI/2;
          }
          if (theta_y3_1st<0.) theta_y3_1st = theta_y3_1st + M_PI;
          in_object_eta[k] = - log(tan(0.5*theta_y3_1st));
        }
        //calculates the real rapidity
        if (irap == 1) {
          in_object_eta[k]=0.5*log((in_object_e[k]+in_object_pz[k])/(in_object_e[k]-in_object_pz[k]));
        }
        in_object_phi[k] = atan2(in_object_py[k], in_object_px[k]);
      }

      //the three-jet resolution
      //threshold y3
      double y3 = 0;

      //vector which will be filled with the
      //minimum of the distances
      double max_dmin_temp=0;

      double max_dmin = 0;

      //distance input object k, beam
      double distance_jB = 0;
      double distance_jB_min = 0;
      //distance of input object k to l
      double distance_jk = 0;
      double distance_jk_min = 0;
      //as we search the minimum of the distances
      //give them values which are for sure higher
      //than those we evaluate first in the for-loups

      size_t index_jB = 0;
      size_t index_j_jk = 0;
      size_t index_k_jk = 0;

      //to decide later if the minmum is a jB or jk
      int decide_jB = -1;

      vector<double> input_pt, input_px, input_py, input_pz;
      vector<double> input_p, input_e, input_phi, input_eta;

      if (!input_pt.empty()) {
        input_pt.clear();
        input_px.clear();
        input_px.clear();
        input_pz.clear();
        input_p.clear();
        input_e.clear();
        input_phi.clear();
        input_eta.clear();
      }

      for (size_t j = 0; j < y3_length; j++){
        input_pt.push_back(in_object_pt[j]);
        input_px.push_back(in_object_px[j]);
        input_py.push_back(in_object_py[j]);
        input_pz.push_back(in_object_pz[j]);
        input_p.push_back(in_object_p[j]);
        input_e.push_back(in_object_e[j]);
        input_phi.push_back(in_object_phi[j]);
        input_eta.push_back(in_object_eta[j]);
      }
      if (y3_length<3) {
        return -1;
      } else {
        size_t rest = y3_length;
        for (size_t i = 0; i<y3_length; i++) {
          //make the minima at the initialization step
          //of each looping bigger than the first values
          distance_jB_min = 0.36*pow(input_pt[0],2) + 10;
          //DELTA PHIs wanted not the pure difference
          distance_jk_min = min(pow(input_pt[1], 2), pow(input_pt[0], 2)) *
            (pow(input_eta[1]-input_eta[0], 2) +
             pow(_delta_phi(input_phi[1], input_phi[0]), 2)) + 10;
          //do the procedure only until we have only 2 objects left anymore
          if (rest > 2) {
            for (size_t j=0; j<rest;j++) {
              //calculate the distance between object j and the beam
              distance_jB = 0.36*pow(input_pt[j], 2);
              if(distance_jB < distance_jB_min){
                distance_jB_min = distance_jB;
                index_jB = j;
              }
              if (j > 0) {
                for(size_t k=0; k<j;k++){
                  //calculate the distance in delta eta and delta phi between object i and object j
                  distance_jk = min(pow(input_pt[j], 2),pow(input_pt[k], 2))*
                    (pow(input_eta[j]-input_eta[k], 2)+
                     pow(_delta_phi(input_phi[j],input_phi[k]), 2));
                  if (distance_jk<distance_jk_min) {
                    distance_jk_min = distance_jk;
                    index_j_jk = j;
                    index_k_jk =k;
                  }
                }
              }
            }
            //decide if the minimum is from a jB or jk combination
            if (distance_jk_min<distance_jB_min) {
              max_dmin_temp = max(distance_jk_min,max_dmin_temp);
              decide_jB = 0;
            } else {
              max_dmin_temp = max(distance_jB_min,max_dmin_temp);
              decide_jB=1;
            }
            //if we have only three jets left calculate
            //the maxima of the dmin's
            //if the minimum is a jB eliminate the input object
            if (decide_jB == 1) {
              //if index_jB is the last one nothing is to do
              if (index_jB != rest-1) {
                for (size_t i=index_jB; i<rest-1;i++) {
                  input_pt[i]=input_pt[i+1];
                  input_phi[i]=input_phi[i+1];
                  input_eta[i]=input_eta[i+1];
                  input_px[i]=input_px[i+1];
                  input_py[i]=input_py[i+1];
                  input_pz[i]=input_pz[i+1];
                  input_e[i]=input_e[i+1];
                }
              }
            }
            //if the minimum is a jk combine both input objects
            if(decide_jB==0) {
              input_px[index_k_jk] = input_px[index_k_jk]+input_px[index_j_jk];
              input_py[index_k_jk] = input_py[index_k_jk]+input_py[index_j_jk];
              input_pz[index_k_jk] = input_pz[index_k_jk]+input_pz[index_j_jk];
              input_e[index_k_jk] = input_e[index_k_jk]+input_e[index_j_jk];
              input_p[index_k_jk] = sqrt(pow(input_px[index_k_jk], 2)+
                                         pow(input_py[index_k_jk], 2)+
                                         pow(input_pz[index_k_jk], 2));
              //calculate the pt, eta and phi of the new combined momenta k_jk
              input_pt[index_k_jk] = sqrt(pow(input_px[index_k_jk], 2)+
                                          pow(input_py[index_k_jk], 2));
              //in the case of pseudorapidity
              if (irap == 0) {
                double theta_new =0;
                if (fabs(input_pz[index_k_jk]) > 1E-5){
                  theta_new = atan(input_pt[index_k_jk]/(input_pz[index_k_jk]));
                } else {
                  theta_new = M_PI/2;
                }
                if (theta_new < 0) {
                  theta_new = theta_new + M_PI;
                }
                input_eta[index_k_jk] = - log(tan(0.5*theta_new));
              }
              //in the real rapidity y is wanted
              if (irap == 1) {
                input_eta[index_k_jk] = 0.5 * log((input_e[index_k_jk]+
                                                   input_pz[index_k_jk]) /
                                                  (input_e[index_k_jk] -
                                                   input_pz[index_k_jk]));
              }
              input_phi[index_k_jk] = atan2(input_py[index_k_jk], input_px[index_k_jk]);
              if (index_j_jk != rest-1) {
                for (size_t i = index_j_jk; i<rest-1;i++) {
                  input_pt[i] = input_pt[i+1];
                  input_phi[i] = input_phi[i+1];
                  input_eta[i] = input_eta[i+1];
                  input_px[i] = input_px[i+1];
                  input_py[i] = input_py[i+1];
                  input_pz[i] = input_pz[i+1];
                  input_e[i] = input_e[i+1];
                }
              }
            }
          }
          if (rest == 3) max_dmin = max_dmin_temp;
          rest = rest-1;
        }
      }

      double et2 = 0;
      et2 = input_pt[0] + input_pt[1];
      y3 = max_dmin/pow(et2,2);

      return y3;
    }


    vector<double> EventShape::_thrust(const vector<double>& input_px, const vector<double>& input_py) {

      double thrustmax_calc = 0;
      double temp_calc = 0;
      size_t length_thrust_calc = 0;
      vector<double> thrust_values, thrust_axis_calc;
      vector<double> p_thrust_max_calc, p_dec_1_calc,  p_dec_2_calc, p_pt_beam_calc;

      if (!thrust_values.empty()){
        thrust_values.clear();
        thrust_axis_calc.clear();
        p_thrust_max_calc.clear();
        p_dec_1_calc.clear();
        p_dec_2_calc.clear();
        p_pt_beam_calc.clear();
      }

      for (size_t j = 0; j < 3; j++){
        p_pt_beam_calc.push_back(0.);
        p_dec_1_calc.push_back(0.);
        p_dec_2_calc.push_back(0.);
        p_thrust_max_calc.push_back(0.);
        thrust_axis_calc.push_back(0.);
      }

      for (size_t j = 0; j < 4; j++) {
        thrust_values.push_back(0.);
      }

      length_thrust_calc = input_px.size();
      if (input_py.size() != length_thrust_calc) {
        /// @todo Change to exception or assert
        cout<<"ERROR in thrust calculation!!! Size of input vectors differs. Change that please!"<<endl;
        return thrust_values;
      }

      double pt_sum_calc =0;
      for(size_t k=0;k<length_thrust_calc;k++){
        pt_sum_calc+=sqrt(pow(input_px[k],2)+pow(input_py[k],2));
        for(size_t j = 0; j < 3; j++){
          p_thrust_max_calc[j]=0;
        }
        //get a vector perpendicular to the beam axis and
        //perpendicular to the momentum of particle k
        //per default beam axis b = (0,0,1)
        p_pt_beam_calc[0] = input_py[k]*1;
        p_pt_beam_calc[1] = - input_px[k]*1;
        p_pt_beam_calc[2] = 0.; // GMA p_pt_beam_calc[3] = 0.;
        for(size_t i=0;i<length_thrust_calc;i++){
          if(i!=k){
            if((input_px[i]*p_pt_beam_calc[0]+input_py[i]*p_pt_beam_calc[1])>=0){
              p_thrust_max_calc[0]= p_thrust_max_calc[0] + input_px[i];
              p_thrust_max_calc[1]= p_thrust_max_calc[1] + input_py[i];
            }else{
              p_thrust_max_calc[0]= p_thrust_max_calc[0] - input_px[i];
              p_thrust_max_calc[1]= p_thrust_max_calc[1] - input_py[i];
            }
          }
        }
        p_dec_1_calc[0] = p_thrust_max_calc[0] + input_px[k];
        p_dec_1_calc[1] = p_thrust_max_calc[1] + input_py[k];
        p_dec_1_calc[2] = 0;
        p_dec_2_calc[0] = p_thrust_max_calc[0] - input_px[k];
        p_dec_2_calc[1] = p_thrust_max_calc[1] - input_py[k];
        p_dec_2_calc[2] = 0;
        temp_calc = pow(p_dec_1_calc[0], 2) + pow(p_dec_1_calc[1], 2);

        if (temp_calc>thrustmax_calc) {
          thrustmax_calc =temp_calc;
          for (size_t i=0; i<3; i++) {
            thrust_axis_calc[i] = p_dec_1_calc[i]/sqrt(thrustmax_calc);
          }
        }
        temp_calc = pow(p_dec_2_calc[0], 2)+pow(p_dec_2_calc[1], 2);
        if (temp_calc > thrustmax_calc) {
          thrustmax_calc =temp_calc;
          for (size_t i=0; i<3; i++) {
            thrust_axis_calc[i] = p_dec_2_calc[i]/sqrt(thrustmax_calc);
          }
        }
      }
      for (size_t j = 0; j < 3; j++) thrust_values[j] = thrust_axis_calc[j];
      const double thrust_calc = sqrt(thrustmax_calc)/pt_sum_calc;

      // the variable which gets returned is not the thrust but tau=1-thrust
      thrust_values[3] = 1 - thrust_calc;
      if (thrust_values[3] < 1e-20) thrust_values[3] = 1e-20;

      return thrust_values;
    }


    double EventShape::_delta_phi(double phi1, double phi2) {
      double dphi = fabs(phi2 - phi1);
      if (dphi > M_PI) dphi = 2*M_PI - dphi;
      return dphi;
    }


    // Returns the scalar product between two 4 momenta
    double EventShape::_lorentz_sp(const vector<double>& a, const vector<double>& b) {
      size_t dim = (size_t) a.size();
      if (a.size()!=b.size()) {
        cout<<"ERROR!!! Dimension of input vectors are different! Change that please!"<<endl;
        return 0;
      } else {
        double l_dot_product=a[dim-1]*b[dim-1];
        for(size_t i=0; i<dim-1;i++){
          l_dot_product-=a[i]*b[i];
        }
        return l_dot_product;
      }
    }


  }


}
