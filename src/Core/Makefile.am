SUBDIRS = yamlcpp

EXTRA_DIST = zstr

noinst_LTLIBRARIES  = libRivetCore.la

libRivetCore_la_SOURCES = \
  Run.cc Event.cc Jet.cc Particle.cc \
  ProjectionApplier.cc Projection.cc \
  Analysis.cc AnalysisLoader.cc AnalysisInfo.cc \
  AnalysisHandler.cc ProjectionHandler.cc

libRivetCore_la_CPPFLAGS = -I$(srcdir)/yamlcpp -DYAMLCPP_API=3 -DYAML_NAMESPACE=RIVET_YAML $(AM_CPPFLAGS)

libRivetCore_la_LIBADD = $(builddir)/yamlcpp/librivet-yaml-cpp.la
