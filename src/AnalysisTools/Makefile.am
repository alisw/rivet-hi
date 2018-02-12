noinst_LTLIBRARIES = libRivetAnalysisTools.la

libRivetAnalysisTools_la_SOURCES = \
    MC_ParticleAnalysis.cc \
    MC_JetAnalysis.cc \
    MC_JetSplittings.cc

libRivetAnalysisTools_la_LDFLAGS = $(AM_LDFLAGS) -L$(FASTJETLIBPATH)
libRivetAnalysisTools_la_LIBS = $(FASTJETCONFIGLIBADD)