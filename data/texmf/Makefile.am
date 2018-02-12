install-data-local:
	for i in `find $(srcdir) -name \*.cfg -o -name \*.pro -o -name \*.sty -o -name \*.tex -o -name \*.cnf` ; \
	do \
	$(install_sh_DATA) $$i $(DESTDIR)$(pkgdatadir)/texmf/$${i#$(srcdir)/}; \
	done

uninstall-local:
	rm -rf $(DESTDIR)$(pkgdatadir)/texmf

EXTRA_DIST = dvips tex cnf
