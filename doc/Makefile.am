dist_noinst_SCRIPTS = mk-analysis-html mk-analysis-latex mk-analysis-json

LATEXSUBSOURCES = preamble.tex intro.tex gettingstarted.tex writinganalyses.tex agilerunmc.tex acknowledgements.tex

DOCSOURCES = $(LATEXSUBSOURCES) \
  compare-histos.txt slip.eps thinker.eps cone.eps bend.eps \
  h-physrev3.bst hepnames.sty JHEP3.cls JHEP.bst make-plots.txt rivet-manual.tex \
  heppennames.sty hepnicenames.sty hepparticles.sty maybemath.sty \
  hepunits.sty underscore.sty microtype.sty

EXTRA_DIST = $(DOCSOURCES) refs.bib

DOCS =
DOCDIRS =


if ENABLE_PYEXT
DOCS += analyses.html
DOCDIRS += analyses
analyses.html: $(top_srcdir)/analyses mk-analysis-html
	LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(LD_LIBRARY_PATH) \
    DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(DYLD_LIBRARY_PATH) \
    PYTHONPATH=$(YODA_PYTHONPATH):$(RIVET_PYTHONPATH):$(PYTHONPATH) \
    $(PYTHON) mk-analysis-html
endif


if ENABLE_PYEXT
DOCS += analyses.json
analyses.json: $(top_srcdir)/analyses mk-analysis-json
	LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(LD_LIBRARY_PATH) \
    DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(DYLD_LIBRARY_PATH) \
    PYTHONPATH=$(YODA_PYTHONPATH):$(RIVET_PYTHONPATH):$(PYTHONPATH) \
    $(PYTHON) mk-analysis-json
endif


if WITH_ASCIIDOC
DOCS += compare-histos.html make-plots.html
EXTRA_DIST += compare-histos.html make-plots.html

compare-histos.html: compare-histos.txt
	asciidoc -a toc compare-histos.txt

make-plots.html: make-plots.txt
	asciidoc -a toc make-plots.txt
endif


if ENABLE_PDFMANUAL
if WITH_PDFLATEX
DOCS += rivet-manual.pdf

if ENABLE_PYEXT
analyses.bib analyses.tex: $(top_srcdir)/analyses mk-analysis-latex
	LD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(LD_LIBRARY_PATH) \
    DYLD_LIBRARY_PATH=$(top_builddir)/src/.libs:$(YAML_CPPLIBPATH):$(FASTJETLIBPATH):$(HEPMCLIBPATH):$(YODALIBPATH):$(DYLD_LIBRARY_PATH) \
    PYTHONPATH=$(YODA_PYTHONPATH):$(RIVET_PYTHONPATH):$(PYTHONPATH) \
    $(PYTHON) mk-analysis-latex
else
analyses.bib analyses.tex:
	> analyses.tex
	> analyses.bib
endif

LATEX = pdflatex --interaction=batchmode
RERUN = "(There were undefined references|Rerun to get (cross-references|the bars) right)"
RERUNBIB = "No file.*\.bbl|Citation.*undefined"
MAKEIDX = "^[^%]*\\makeindex"

rivet-manual.pdf : rivet-manual.tex $(LATEXSUBSOURCES) analyses.tex refs.bib analyses.bib
	$(LATEX) $<; true
	egrep -c $(RERUNBIB) rivet-manual.log && (bibtex rivet-manual && cp rivet-manual.toc rivet-manual.toc.bak && $(LATEX) $<); true
	for i in `seq 5`; do if egrep $(RERUN) rivet-manual.log; then echo "LaTeX re-run $i"; cp rivet-manual.toc rivet-manual.toc.bak; $(LATEX) $<; else break; fi; done; true
	if cmp -s rivet-manual.toc rivet-manual.toc.bak; then true; else $(LATEX) $<; true; fi
	rm -f rivet-manual.toc.bak; true

paper.pdf : paper.tex $(LATEXSUBSOURCES) analyses.tex refs.bib analyses.bib
	$(LATEX) $<; true
	egrep -c $(RERUNBIB) paper.log && (bibtex paper && cp paper.toc paper.toc.bak && $(LATEX) $<); true
	for i in `seq 5`; do if egrep $(RERUN) paper.log; then echo "LaTeX re-run $i"; cp paper.toc paper.toc.bak; $(LATEX) $<; else break; fi; done; true
	if cmp -s paper.toc paper.toc.bak; then true; else $(LATEX) $<; true; fi
	rm -f paper.toc.bak; true

endif
endif


################


.PHONY = all doc upload arxivtar

doc: $(DOCS)
	@true

RSH = rsync
DEST = login.hepforge.org:rivet/public_html/
upload: $(DOCS) $(DOCDIRS)
	$(RSH) -r $? $(DEST)

arxivtar: $(DOCSOURCES) rivet-manual.bbl
	for i in *.png; do convert $$i $${i/.png/.eps}; done
	tar czf rivet-manual.tar.gz \
rivet-manual.tex $(LATEXSUBSOURCES) analyses.tex refs.bib analyses.bib rivetversion.sty rivet-manual.bbl \
hepnames.sty hepnicenames.sty hepparticles.sty heppennames.sty hepunits.sty maybemath.sty microtype.sty underscore.sty \
bend.eps cone.eps thinker.eps slip.eps \
h-physrev3.bst JHEP3.cls JHEP.bst


## Install!
pkgdata_DATA = $(DOCS)

install-data-local: doc
	@$(NORMAL_INSTALL)
	test "$(DOCDIRS)" && cp -r $(DOCDIRS) $(DESTDIR)$(pkgdatadir)/

uninstall-local:
	@$(NORMAL_UNINSTALL)
	rm -rf $(DESTDIR)$(pkgdatadir)/analyses


mostlyclean-local:
	rm -rf *.aux *.log *.toc

clean-local:
	rm -rf $(DOCS) $(DOCDIRS)
