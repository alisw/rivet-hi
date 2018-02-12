// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/BinnedHistogram.hh"

namespace Rivet {


  /// @brief ATLAS azimuthal decorrelation with jet veto analysis
  /// @author James Robinson <james.robinson@cern.ch>
  class ATLAS_2014_I1307243 : public Analysis {
  public:

    /// Constructor
    ATLAS_2014_I1307243()
      : Analysis("ATLAS_2014_I1307243"),
        _fiducialRegions{2010, 2011},
        _vetoScale{20*GeV, 30*GeV},
        _yFiducial{4.4, 2.4},
        _gapCategories{"inclusive", "gap"},
        _dy_max(8),
        _nEventsInAcceptance(0),
        _sumOfAcceptedWeights(0.)
    {   }


    /// Book histograms and initialise projections before the run
    void init() {

      /// Initialise and register projections here
      FinalState fs;
      FastJets fastJets(fs, FastJets::ANTIKT, 0.6);
      fastJets.useInvisibles(true);
      declare(fastJets, "AntiKt6JetsWithInvisibles");


      /// Book histograms
      foreach (const string& gapCategory, _gapCategories ) {
        const int gapCategoryOffset = (gapCategory == "inclusive") ? 0 : 1;

        // Temporary inclusive and gap histograms
        _h_tmp_events_dy[gapCategory] = bookHisto1D(1, 1, 1);
        _h_tmp_events_dy[gapCategory]->setPath("/TMP/" + toString(gapCategory) + "_events_dy");
        _h_tmp_events_pTbar[gapCategory] = bookHisto1D(2, 1, 1);
        _h_tmp_events_pTbar[gapCategory]->setPath("/TMP/" + toString(gapCategory) + "_events_pTbar");

        // Azimuthal moment histograms
        _h_profiled_cosDeltaPhi_dy[gapCategory]       = bookProfile1D( 5+4*gapCategoryOffset, 1, 1);
        _h_profiled_cosDeltaPhi_pTbar[gapCategory]    = bookProfile1D( 6+4*gapCategoryOffset, 1, 1);
        _h_C2C1_dy[gapCategory]                       = bookScatter2D( 7+4*gapCategoryOffset, 1, 1, false);
        _h_C2C1_pTbar[gapCategory]                    = bookScatter2D( 8+4*gapCategoryOffset, 1, 1, false);
        _h_profiled_cosTwoDeltaPhi_dy[gapCategory]    = bookProfile1D(37+2*gapCategoryOffset, 1, 1);
        _h_profiled_cosTwoDeltaPhi_pTbar[gapCategory] = bookProfile1D(38+2*gapCategoryOffset, 1, 1);

        // Gap fraction vs. Q0 and cross-section in dy slices
        for (size_t dyLow = 0; dyLow < _dy_max; ++dyLow ) {
          Histo1DPtr _h_tmp_events_Q0_single_dySlice = bookHisto1D( 29+dyLow, 1, 1);
          _h_tmp_events_Q0_single_dySlice->setPath("/TMP/" + toString(gapCategory) + "_events_dySlice_" + toString(dyLow) + "_" + toString(dyLow+1) + "_Q0");
          _h_tmp_events_Q0_dySlices[gapCategory].addHistogram( dyLow, dyLow+1, _h_tmp_events_Q0_single_dySlice );
          _h_crossSection_dphi_dySlices[gapCategory].addHistogram( dyLow, dyLow+1, bookHisto1D( 13+(_dy_max*gapCategoryOffset)+dyLow, 1, 1));
        }

      }

      // Number of jets in rapidity interval
      _h_profiled_nJets_rapidity_interval_dy    = bookProfile1D( 3, 1, 1);
      _h_profiled_nJets_rapidity_interval_pTbar = bookProfile1D( 4, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get the event weight
      const double weight( event.weight() );
      bool eventAccepted( false );

      for (int iFiducialRegion = 0; iFiducialRegion < 2; ++iFiducialRegion ) {

        // Retrieve all anti-kt R=0.6 jets above _pTMin and inside |_yFiducial|
        const Jets akt6Jets = apply<JetAlg>(event, "AntiKt6JetsWithInvisibles").jetsByPt( Cuts::absrap < _yFiducial.at(iFiducialRegion) );
        // If there are fewer than 2 jets then bail
        if ( akt6Jets.size() < 2 ) { vetoEvent; }

        // Require jets to be above {60, 50} GeV
        if ( akt6Jets.at(0).momentum().pT() < 60.*GeV || akt6Jets.at(1).momentum().pT() < 50.*GeV ) { vetoEvent; }

        // Identify gap boundaries
        double yMin( std::min( akt6Jets.at(0).momentum().rapidity(), akt6Jets.at(1).momentum().rapidity() ) );
        double yMax( std::max( akt6Jets.at(0).momentum().rapidity(), akt6Jets.at(1).momentum().rapidity() ) );

        // Determine azimuthal decorrelation quantities
        const double dy( yMax - yMin );
        const double dphi( mapAngle0ToPi( akt6Jets.at(0).momentum().phi() - akt6Jets.at(1).momentum().phi() ) );
        const double pTbar( (akt6Jets.at(0).momentum().pT() + akt6Jets.at(1).momentum().pT())/2.0 );

        // Impose minimum dy for the 2011 phase space
        if ( _fiducialRegions.at(iFiducialRegion) == 2011 && dy < 1.0 ) { vetoEvent; }

        // Construct gap candidates sample
        Jets gapCandidates;
        foreach( const Jet &j, akt6Jets ) {
          if ( inRange( j.momentum().rapidity(), yMin, yMax, OPEN, OPEN ) ) {
            gapCandidates.push_back( j );
          }
        }

        // Determine gap quantities
        unsigned int nJets_rapidity_interval( 0 );
        double maximumGapQ0( 0. );
        foreach( const Jet &jet, gapCandidates ) {
          const double pT( jet.momentum().pT() );
          if ( pT > _vetoScale.at(iFiducialRegion) ) { ++nJets_rapidity_interval; }
          if ( pT > maximumGapQ0 ) { maximumGapQ0 = pT; }
        }

        // Fill histograms
        if ( weight < 0.0 ) {
          MSG_DEBUG( "Negative weight " << weight << "found!" );
        }
        fillHistograms( _fiducialRegions.at(iFiducialRegion), dy, pTbar, dphi, nJets_rapidity_interval, maximumGapQ0, weight );
        eventAccepted = true;
      }

      // Count number of accepted events
      if ( eventAccepted ) {
        _nEventsInAcceptance++;
        _sumOfAcceptedWeights += weight;
      }
      return;
    }

    void fillHistograms( const unsigned int &fiducialRegion, const double &dy, const double &pTbar, const double &dphi, const double &nJets_rapidity_interval, const double &maximumGapQ0, const double &weight) {
      // Determine gap category
      vector<string> eventGapCategories{{"inclusive"}};
      if ( nJets_rapidity_interval == 0 ) { eventGapCategories += string("gap"); }

      // Fill histograms relevant for comparison with 2010 data
      if ( fiducialRegion == _fiducialRegions.at(0) ) {
        // Fill inclusive and gap histograms
        foreach( const string &gapCategory, eventGapCategories ) {
          _h_tmp_events_dy[gapCategory]->fill( dy, weight );
          _h_crossSection_dphi_dySlices[gapCategory].fill( dy, dphi / M_PI, weight );
          _h_profiled_cosDeltaPhi_dy[gapCategory]->fill( dy, cos(M_PI - dphi), weight );
          _h_profiled_cosTwoDeltaPhi_dy[gapCategory]->fill( dy, cos(2.0 * dphi), weight );
        }
        // Fill profiled nJets_rapidity_interval
        _h_profiled_nJets_rapidity_interval_dy->fill( dy, nJets_rapidity_interval, weight );
        // Fill Q0 histograms - can fill multiple points per event
        foreach( const YODA::HistoBin1D Q0_bin, _h_tmp_events_Q0_dySlices["inclusive"].getHistograms().at(0)->bins() ) {
          const double Q0( Q0_bin.xMid() );
          _h_tmp_events_Q0_dySlices["inclusive"].fill(dy, Q0, weight);
          if ( maximumGapQ0 <= Q0 ) { _h_tmp_events_Q0_dySlices["gap"].fill(dy, Q0, weight); }
        }

      // Fill histograms relevant for comparison with 2011 data
      } else if ( fiducialRegion == _fiducialRegions.at(1) ) {
        // Fill inclusive and gap histograms
        foreach( const string &gapCategory, eventGapCategories ) {
          _h_tmp_events_pTbar[gapCategory]->fill( pTbar, weight );
          _h_profiled_cosDeltaPhi_pTbar[gapCategory]->fill( pTbar, cos(M_PI - dphi), weight );
          _h_profiled_cosTwoDeltaPhi_pTbar[gapCategory]->fill( pTbar, cos(2.0 * dphi), weight );
        }
        // Fill profiled nJets_rapidity_interval
        _h_profiled_nJets_rapidity_interval_pTbar->fill( pTbar, nJets_rapidity_interval, weight );
      }
      return;
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      // Normalise cross-section plots to correct cross-section
      const double xs_pb( crossSection()/picobarn );
      const double nEventsGenerated( sumOfWeights() );
      double xs_norm_factor( xs_pb/nEventsGenerated );
      const double ySpan(1); // all dy spans are 1
      foreach (const string& gapCategory, _gapCategories ) {
        _h_crossSection_dphi_dySlices[gapCategory].scale(xs_norm_factor/ySpan/M_PI, this);
      }

      // Create C2/C1 scatter from profiles
      foreach (const string& gapCategory, _gapCategories ) {
        divide( _h_profiled_cosTwoDeltaPhi_dy[gapCategory], _h_profiled_cosDeltaPhi_dy[gapCategory], _h_C2C1_dy[gapCategory] );
        divide( _h_profiled_cosTwoDeltaPhi_pTbar[gapCategory], _h_profiled_cosDeltaPhi_pTbar[gapCategory], _h_C2C1_pTbar[gapCategory] );
      }

      // Fill simple gap fractions
      Scatter2DPtr h_gap_fraction_dy    = bookScatter2D( 1, 1, 1);
      Scatter2DPtr h_gap_fraction_pTbar = bookScatter2D( 2, 1, 1, false);
      Rivet::Analysis::efficiency( *_h_tmp_events_dy["gap"].get(), *_h_tmp_events_dy["inclusive"].get(), h_gap_fraction_dy );
      Rivet::Analysis::efficiency( *_h_tmp_events_pTbar["gap"].get(), *_h_tmp_events_pTbar["inclusive"].get(), h_gap_fraction_pTbar );

      // Register and fill Q0 gap fractions
      for (unsigned int dyLow = 0; dyLow < _dy_max; ++dyLow ) {
        Scatter2DPtr h_gap_fraction_Q0 = bookScatter2D( 29+dyLow, 1, 1, false);
        Rivet::Analysis::efficiency( *_h_tmp_events_Q0_dySlices["gap"].getHistograms().at(dyLow).get(), *_h_tmp_events_Q0_dySlices["inclusive"].getHistograms().at(dyLow).get(), h_gap_fraction_Q0 );
      }

      /// Print summary information
      MSG_INFO( "Cross-Section/pb     : " << xs_pb );
      MSG_INFO( "Sum of weights       : " << sumOfWeights() << "  (" << _sumOfAcceptedWeights << " accepted)" );
      MSG_INFO( "nEvents              : " << numEvents() << "  (" << _nEventsInAcceptance << " accepted)" );
      MSG_INFO( "Normalisation factor : " << xs_norm_factor );
    }


  private:

    /// Member variables
    vector<unsigned int> _fiducialRegions;
    vector<double> _vetoScale, _yFiducial;
    vector<string> _gapCategories;

    unsigned int _dy_max;
    unsigned int _nEventsInAcceptance;
    double _sumOfAcceptedWeights;

    /// Histograms
    // Number of events : gap and non-gap : necessary input for gap fractions but not saved as output
    map<string, Histo1DPtr> _h_tmp_events_dy;
    map<string, Histo1DPtr> _h_tmp_events_pTbar;
    map<string, BinnedHistogram<double> > _h_tmp_events_Q0_dySlices;

    // Number of jets in rapidity interval
    Profile1DPtr _h_profiled_nJets_rapidity_interval_dy;
    Profile1DPtr _h_profiled_nJets_rapidity_interval_pTbar;

    // Azimuthal moment histograms
    map<string, Profile1DPtr> _h_profiled_cosDeltaPhi_dy;
    map<string, Profile1DPtr> _h_profiled_cosDeltaPhi_pTbar;
    map<string, Profile1DPtr> _h_profiled_cosTwoDeltaPhi_dy;
    map<string, Profile1DPtr> _h_profiled_cosTwoDeltaPhi_pTbar;
    map<string, Scatter2DPtr> _h_C2C1_dy;
    map<string, Scatter2DPtr> _h_C2C1_pTbar;

    // Cross-section vs. #Delta#phi
    map<string, BinnedHistogram<double> > _h_crossSection_dphi_dySlices;
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2014_I1307243);

}
