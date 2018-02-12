# distutils: language = c++

cimport rivet as c
from cython.operator cimport dereference as deref
# Need to be careful with memory management -- perhaps use the base object that
# we used in YODA?

cdef extern from "<utility>" namespace "std" nogil:
    cdef c.unique_ptr[c.Analysis] move(c.unique_ptr[c.Analysis])

cdef class AnalysisHandler:
    cdef c.AnalysisHandler *_ptr

    def __cinit__(self):
        self._ptr = new c.AnalysisHandler()

    def __del__(self):
        del self._ptr

    def setIgnoreBeams(self, ignore=True):
        self._ptr.setIgnoreBeams(ignore)

    def addAnalysis(self, name):
        self._ptr.addAnalysis(name.encode('utf-8'))
        return self

    def analysisNames(self):
        anames = self._ptr.analysisNames()
        return [ a.decode('utf-8') for a in anames ]

    # def analysis(self, aname):
    #     cdef c.Analysis* ptr = self._ptr.analysis(aname)
    #     cdef Analysis pyobj = Analysis.__new__(Analysis)
    #     if not ptr:
    #         return None
    #     pyobj._ptr = ptr
    #     return pyobj

    def readData(self, name):
        self._ptr.readData(name.encode('utf-8'))

    def writeData(self, name):
        self._ptr.writeData(name.encode('utf-8'))

    def crossSection(self):
        return self._ptr.crossSection()

    def finalize(self):
        self._ptr.finalize()


cdef class Run:
    cdef c.Run *_ptr

    def __cinit__(self, AnalysisHandler h):
        self._ptr = new c.Run(h._ptr[0])

    def __del__(self):
        del self._ptr

    def setCrossSection(self, double x):
        self._ptr.setCrossSection(x)
        return self

    def setListAnalyses(self, choice):
        self._ptr.setListAnalyses(choice)
        return self

    def init(self, name, weight=1.0):
        return self._ptr.init(name.encode('utf-8'), weight)

    def openFile(self, name, weight=1.0):
        return self._ptr.openFile(name.encode('utf-8'), weight)

    def readEvent(self):
        return self._ptr.readEvent()

    def skipEvent(self):
        return self._ptr.skipEvent()

    def processEvent(self):
        return self._ptr.processEvent()

    def finalize(self):
        return self._ptr.finalize()


cdef class Analysis:
    cdef c.unique_ptr[c.Analysis] _ptr

    def __init__(self):
        raise RuntimeError('This class cannot be instantiated')

    def requiredBeams(self):
        return deref(self._ptr).requiredBeams()

    def requiredEnergies(self):
        return deref(self._ptr).requiredEnergies()

    def keywords(self):
        kws = deref(self._ptr).keywords()
        return [ k.decode('utf-8') for k in kws ]

    def authors(self):
        auths = deref(self._ptr).authors()
        return [ a.decode('utf-8') for a in auths ]

    def bibKey(self):
        return deref(self._ptr).bibKey().decode('utf-8')

    def name(self):
        return deref(self._ptr).name().decode('utf-8')

    def bibTeX(self):
        return deref(self._ptr).bibTeX().decode('utf-8')

    def references(self):
        refs = deref(self._ptr).references()
        return [ r.decode('utf-8') for r  in refs ]

    def collider(self):
        return deref(self._ptr).collider().decode('utf-8')

    def description(self):
        return deref(self._ptr).description().decode('utf-8')

    def experiment(self):
        return deref(self._ptr).experiment().decode('utf-8')

    def inspireId(self):
        return deref(self._ptr).inspireId().decode('utf-8')

    def spiresId(self):
        return deref(self._ptr).spiresId().decode('utf-8')

    def runInfo(self):
        return deref(self._ptr).runInfo().decode('utf-8')

    def status(self):
        return deref(self._ptr).status().decode('utf-8')

    def summary(self):
        return deref(self._ptr).summary().decode('utf-8')

    def year(self):
        return deref(self._ptr).year().decode('utf-8')

    def luminosityfb(self):
        return deref(self._ptr).luminosityfb().decode('utf-8')

#cdef object
LEVELS = dict(TRACE = 0, DEBUG = 10, INFO = 20,
              WARN = 30, WARNING = 30, ERROR = 40,
              CRITICAL = 50, ALWAYS = 50)


cdef class AnalysisLoader:
    @staticmethod
    def analysisNames():
        names = c.AnalysisLoader_analysisNames()
        return [ n.decode('utf-8') for n in names ]


    @staticmethod
    def getAnalysis(name):
        name = name.encode('utf-8')
        cdef c.unique_ptr[c.Analysis] ptr = c.AnalysisLoader_getAnalysis(name)
        cdef Analysis pyobj = Analysis.__new__(Analysis)
        if not ptr:
            return None
        pyobj._ptr = move(ptr)
        # Create python object
        return pyobj


def getAnalysisLibPaths():
    ps = c.getAnalysisLibPaths()
    return [ p.decode('utf-8') for p in ps ]

def setAnalysisLibPaths(xs):
    bs = [ x.encode('utf-8') for x in xs ]
    c.setAnalysisLibPaths(bs)

def addAnalysisLibPath(path):
    c.addAnalysisLibPath(path.encode('utf-8'))


def setAnalysisDataPaths(xs):
    bs = [ x.encode('utf-8') for x in xs ]
    c.setAnalysisDataPaths(bs)

def addAnalysisDataPath(path):
    c.addAnalysisDataPath(path.encode('utf-8'))

def getAnalysisDataPaths():
    ps = c.getAnalysisDataPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisDataFile(q):
    f = c.findAnalysisDataFile(q.encode('utf-8'))
    return f.decode('utf-8')

def getAnalysisRefPaths():
    ps = c.getAnalysisRefPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisRefFile(q):
    f = c.findAnalysisRefFile(q.encode('utf-8'))
    return f.decode('utf-8')


def getAnalysisInfoPaths():
    ps = c.getAnalysisInfoPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisInfoFile(q):
    f = c.findAnalysisInfoFile(q.encode('utf-8'))
    return f.decode('utf-8')

def getAnalysisPlotPaths():
    ps = c.getAnalysisPlotPaths()
    return [ p.decode('utf-8') for p in ps ]

def findAnalysisPlotFile(q):
    f = c.findAnalysisPlotFile(q.encode('utf-8'))
    return f.decode('utf-8')

def version():
    return c.version().decode('utf-8')

def setLogLevel(name, level):
    c.setLogLevel(name.encode('utf-8'), level)
