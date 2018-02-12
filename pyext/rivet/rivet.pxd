from libcpp.map cimport map
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp cimport bool
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr

ctypedef int PdgId
ctypedef pair[PdgId,PdgId] PdgIdPair

cdef extern from "Rivet/AnalysisHandler.hh" namespace "Rivet":
    cdef cppclass AnalysisHandler:
        void setIgnoreBeams(bool)
        AnalysisHandler& addAnalysis(string)
        vector[string] analysisNames() const
        # Analysis* analysis(string)
        void writeData(string&)
        void readData(string&)
        double crossSection()
        void finalize()

cdef extern from "Rivet/Run.hh" namespace "Rivet":
    cdef cppclass Run:
        Run(AnalysisHandler)
        Run& setCrossSection(double) # For chaining?
        Run& setListAnalyses(bool)
        bool init(string, double) except + # $2=1.0
        bool openFile(string, double) except + # $2=1.0
        bool readEvent() except +
        bool skipEvent() except +
        bool processEvent() except +
        bool finalize() except +

cdef extern from "Rivet/Analysis.hh" namespace "Rivet":
    cdef cppclass Analysis:
        vector[PdgIdPair]& requiredBeams()
        vector[pair[double, double]] requiredEnergies()
        vector[string] authors()
        vector[string] references()
        vector[string] keywords()
        string name()
        string bibTeX()
        string bibKey()
        string collider()
        string description()
        string experiment()
        string inspireId()
        string spiresId()
        string runInfo()
        string status()
        string summary()
        string year()
        string luminosityfb()

# Might need to translate the following errors, although I believe 'what' is now
# preserved. But often, we need the exception class name.
#Error
#RangeError
#LogicError
#PidError
#InfoError
#WeightError
#UserError

cdef extern from "Rivet/AnalysisLoader.hh":
    vector[string] AnalysisLoader_analysisNames "Rivet::AnalysisLoader::analysisNames" ()
    unique_ptr[Analysis] AnalysisLoader_getAnalysis "Rivet::AnalysisLoader::getAnalysis" (string)

cdef extern from "Rivet/Tools/RivetPaths.hh" namespace "Rivet":
    vector[string] getAnalysisLibPaths()
    void setAnalysisLibPaths(vector[string])
    void addAnalysisLibPath(string)

    vector[string] getAnalysisDataPaths()
    void setAnalysisDataPaths(vector[string])
    void addAnalysisDataPath(string)
    string findAnalysisDataFile(string)

    vector[string] getAnalysisRefPaths()
    string findAnalysisRefFile(string)

    vector[string] getAnalysisInfoPaths()
    string findAnalysisInfoFile(string)

    vector[string] getAnalysisPlotPaths()
    string findAnalysisPlotFile(string)

cdef extern from "Rivet/Rivet.hh" namespace "Rivet":
    string version()

cdef extern from "Rivet/Tools/Logging.hh":
    void setLogLevel "Rivet::Log::setLevel" (string, int)
