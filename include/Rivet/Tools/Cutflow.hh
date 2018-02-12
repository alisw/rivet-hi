#ifndef RIVET_Cutflow_HH
#define RIVET_Cutflow_HH

#include "Rivet/Tools/Utils.hh"

namespace Rivet {


  /// A tracker of numbers & fractions of events passing sequential cuts
  struct Cutflow {

    /// @brief Default constructor
    ///
    /// Does nothing! Just to allow storage in STL containers and use as a member variable without using the init list
    Cutflow() {}

    /// Proper constructor
    Cutflow(const string& cfname, const vector<string>& cutnames)
      : name(cfname), ncuts(cutnames.size()), cuts(cutnames), counts(ncuts+1, 0)
    {  }

    /// @brief Fill the pre-cut counter
    void fillinit(double weight=1.) {
      counts[0] += weight;
    }

    /// @brief Fill the @a {icut}'th post-cut counter, starting at icut=1 for first cut
    ///
    /// @note Returns the cut result to allow 'side-effect' cut-flow filling in an if-statement
    bool fill(size_t icut, bool cutresult=true, double weight=1.) {
      if (cutresult) counts[icut] += weight;
      return cutresult;
    }

    /// @brief Fill the @a {icut}'th post-cut counter, starting at icut=1 for first cut (cutvalue=true overload)
    ///
    /// This version exists to allow calling fill(i, weight) without the weight
    /// getting cast to a bool, or having to explicitly add a 'true' middle arg.
    ///
    /// @note Returns the cut result to allow 'side-effect' cut-flow filling in an if-statement
    bool fill(size_t icut, double weight) {
      return fill(icut, true, weight);
    }

    /// @brief Fill all cut-state counters from an Ncut-element results vector
    ///
    /// This function is to be used to fill all of an event's pre- and post-cut
    /// state counters at once, including the incoming event counter. It must not be
    /// mixed with calls to the @c fill(size_t, bool) and @c fillinit() methods,
    /// or double-counting will occur.
    ///
    /// @note Returns the overall cut result to allow 'side-effect' cut-flow filling in an if-statement
    bool fill(const vector<bool>& cutresults, double weight=1.) {
      if (cutresults.size() != ncuts)
        throw RangeError("Number of filled cut results needs to match the Cutflow construction");
      counts[0] += 1;
      for (size_t i = 0; i < ncuts; ++i) {
        if (cutresults[i]) counts[i+1] += weight; else break;
      }
      return all(cutresults);
    }

    /// @todo Add a fillnext(), keeping track of current ifill

    /// @todo Add a fillhead() (or vector fillnext()?)

    /// @brief Fill the N trailing post-cut counters, when supplied with an N-element results vector
    ///
    /// The @a cutresults vector represents the boolean results of the last N cuts. This function
    /// allows mixing of cut-flow filling with higher-level analyze() function escapes such as
    /// the vetoEvent directive. The initial state (state 0) is not incremented.
    ///
    /// @note Returns the overall cut result to allow 'side-effect' cut-flow filling in an if-statement
    bool filltail(const vector<bool>& cutresults, double weight=1.) {
      if (cutresults.size() > ncuts)
        throw RangeError("Number of filled cut results needs to match the Cutflow construction");
      const size_t offset = counts.size() - cutresults.size();
      for (size_t i = 0; i < cutresults.size(); ++i) {
        if (cutresults[i]) counts[offset+i] += weight; else break;
      }
      return all(cutresults);
    }

    /// Scale the cutflow weights by the given factor
    void scale(double factor) {
      for (double& x : counts) x *= factor;
    }

    /// Create a string representation
    string str() const {
      stringstream ss;
      ss << fixed << setprecision(1) << counts[0];
      const size_t count0len = ss.str().length();
      ss.str("");
      ss << name << " cut-flow:\n";
      size_t maxnamelen = 0;
      for (const string& t : cuts)
        maxnamelen = max(t.length(), maxnamelen);
      ss << setw(maxnamelen+5) << "" << "   "
         << setw(count0len) << right << "Count" << "    "
         << setw(6) << right << "A_cumu" << "    "
         << setw(6) << right << "A_incr";
      for (size_t i = 0; i <= ncuts; ++i) {
        const int pcttot = (counts[0] == 0) ? -1 : round(100*counts[i]/double(counts[0]));
        const int pctinc = (i == 0 || counts[i-1] == 0) ? -1 : round(100*counts[i]/double(counts[i-1]));
        stringstream ss2;
        ss2 << fixed << setprecision(1) << counts[i];
        const string countstr = ss2.str(); ss2.str("");
        ss2 << fixed << setprecision(3) << pcttot << "%";
        const string pcttotstr = ss2.str(); ss2.str("");
        ss2 << fixed << setprecision(3) << pctinc << "%";
        const string pctincstr = ss2.str();
        ss << "\n"
           << setw(maxnamelen+5) << left << (i == 0 ? "" : "Pass "+cuts[i-1]) << "   "
           << setw(count0len) << right << countstr << "    "
           << setw(6) << right << (pcttot < 0 ? "- " : pcttotstr) << "    "
           << setw(6) << right << (pctinc < 0 ? "- " : pctincstr);
      }
      return ss.str();
    }

    /// Print string representation to a stream
    void print(ostream& os) const {
      os << str() << flush;
    }

    string name;
    size_t ncuts;
    vector<string> cuts;
    vector<double> counts;

  };


  /// Print a Cutflow to a stream
  inline ostream& operator << (ostream& os, const Cutflow& cf) {
    return os << cf.str();
  }



  /// A container for several Cutflow objects, with some convenient batch access
  struct Cutflows {

    /// Do-nothing default constructor
    Cutflows() {  }

    /// Populating constructor
    Cutflows(const vector<Cutflow>& cutflows) : cfs(cutflows) {  }

    /// Append a provided Cutflow to the list
    void addCutflow(const Cutflow& cf) {
      cfs.push_back(cf);
    }

    /// Append a newly constructed Cutflow to the list
    void addCutflow(const string& cfname, const vector<string>& cutnames) {
      cfs.push_back(Cutflow(cfname, cutnames));
    }

    /// Access the @a i'th Cutflow
    Cutflow& operator [] (size_t i) { return cfs[i]; }
    /// Access the @a i'th Cutflow (const)
    const Cutflow& operator [] (size_t i) const { return cfs[i]; }

    /// Access the Cutflow whose name is @a name
    Cutflow& operator [] (const string& name) {
      for (Cutflow& cf : cfs)
        if (cf.name == name) return cf;
      throw UserError("Requested cut-flow name '" + name + "' does not exist");
    }
    /// Access the @a i'th Cutflow (const)
    const Cutflow& operator [] (const string& name) const {
      for (const Cutflow& cf : cfs)
        if (cf.name == name) return cf;
      throw UserError("Requested cut-flow name '" + name + "' does not exist");
    }

    /// Fill the pre-cuts state counter for all contained {Cutflow}s
    void fillinit(double weight=1.) {
      for (Cutflow& cf : cfs) cf.fillinit(weight);
    }

    /// @brief Fill the @a {icut}'th post-cut counter, starting at icut=1 for first cut, with the same result for all {Cutflow}s
    bool fill(size_t icut, bool cutresult=true, double weight=1.) {
      for (Cutflow& cf : cfs) cf.fill(icut, cutresult, weight);
      return cutresult;
    }

    /// @brief Fill the @a {icut}'th post-cut counter, starting at icut=1 for first cut, with the same result for all {Cutflow}s (cutresult=true overload)
    ///
    /// This version exists to allow calling fill(i, weight) without the weight
    /// getting cast to a bool, or having to explicitly add a 'true' middle arg.
    ///
    /// @note Returns the cut result to allow 'side-effect' cut-flow filling in an if-statement
    bool fill(size_t icut, double weight) {
      return fill(icut, true, weight);
    }

    /// @todo Add a fillnext(), keeping track of current ifill

    /// @todo Add a fillhead() (or vector fillnext()?)

    /// Scale the contained {Cutflow}s by the given factor
    void scale(double factor) {
      for (Cutflow& cf : cfs) cf.scale(factor);
    }

    /// Create a string representation
    string str() const {
      stringstream ss;
      for (const Cutflow& cf : cfs)
        ss << cf << "\n\n";
      return ss.str();
    }

    /// Print string representation to a stream
    void print(ostream& os) const {
      os << str() << flush;
    }

    vector<Cutflow> cfs;

  };

  /// Print a Cutflows to a stream
  inline ostream& operator << (ostream& os, const Cutflows& cfs) {
    return os << cfs.str();
  }


}

#endif
