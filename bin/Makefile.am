bin_SCRIPTS = rivet-config rivet-buildplugin
dist_bin_SCRIPTS = make-plots
EXTRA_DIST =

RIVETPROGS = \
    rivet \
    rivet-mkanalysis \
	rivet-findid rivet-which \
    rivet-cmphistos rivet-mkhtml
if ENABLE_PYEXT
dist_bin_SCRIPTS += $(RIVETPROGS)
else
EXTRA_DIST += $(RIVETPROGS)
endif

## bash completion
if ENABLE_PYEXT
dist_pkgdata_DATA = rivet-completion
bashcomp_dir = $(DESTDIR)$(prefix)/etc/bash_completion.d
install-data-local:
	if [[ -d "$(bashcomp_dir)" && -w "$(bashcomp_dir)" ]]; then \
      $(install_sh_DATA) rivet-completion $(bashcomp_dir)/; fi
uninstall-local:
	rm -f $(bashcomp_dir)/rivet-completion
else
EXTRA_DIST += rivet-completion
endif

## No-Python Rivet program
noinst_PROGRAMS = rivet-nopy
rivet_nopy_SOURCES = rivet-nopy.cc
rivet_nopy_CPPFLAGS = -I$(top_srcdir)/include $(AM_CPPFLAGS)
rivet_nopy_LDFLAGS = -L$(HEPMCLIBPATH)
rivet_nopy_LDADD = $(top_builddir)/src/libRivet.la -lHepMC
