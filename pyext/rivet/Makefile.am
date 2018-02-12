EXTRA_DIST = \
	__init__.py \
    spiresbib.py util.py \
    plotinfo.py aopaths.py \
    core.pyx rivet.pxd \
    core.cpp

MAINTAINERCLEANFILES = core.cpp
BUILT_SOURCES = core.cpp

if WITH_CYTHON
RIVETINCLUDE = $(top_srcdir)/include/Rivet/
core.cpp: core.pyx rivet.pxd $(RIVETINCLUDE)/Analysis.hh $(RIVETINCLUDE)/AnalysisHandler.hh $(RIVETINCLUDE)/AnalysisLoader.hh $(RIVETINCLUDE)/Run.hh
	cython $(srcdir)/core.pyx --cplus \
	   -I $(srcdir) -I $(srcdir)/include \
	   -I $(builddir) -I $(builddir)/include \
	   -o $@
else
core.cpp:
	@echo "Not (re)generating core.cpp since Cython is not installed"
endif

