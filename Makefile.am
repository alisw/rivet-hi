ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src pyext data include bin analyses doc test

EXTRA_DIST = GUIDELINES

DISTCHECK_CONFIGURE_FLAGS = --enable-debug --with-yoda=$(YODAPATH) --with-hepmc=$(HEPMCPATH) --with-fastjet=$(FASTJETPATH)

pkgconfigdir = $(libdir)/pkgconfig
dist_pkgconfig_DATA = rivet.pc

doc:
	cd doc && $(MAKE) doc

.PHONY : doc dox pyclean

clean-local:
	@rm -rf a.out


## Doxygen
if WITH_DOXYGEN
dox:
	@echo "Running Doxygen..."
	@$(DOXYGEN) $(top_builddir)/Doxyfile || true
install-dox:
	$(install_sh_DATA) doxy/html/* $(DESTDIR)$(pkgdatadir)/doc/html
mostlyclean-local:
	@rm -rf $(top_builddir)/doxy
uninstall-local:
	@rm -rf $(DESTDIR)$(pkgdatadir)/doc; \
	test ! -d $(DESTDIR)$(pkgdatadir) || rmdir $(DESTDIR)$(pkgdatadir) || true
endif


## Upload to HepForge
RSH = rsync
DEST = login.hepforge.org:rivet/downloads/
upload: dist
	$(RSH) $(DIST_ARCHIVES) $(DEST)
