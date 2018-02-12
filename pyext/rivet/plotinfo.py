from __future__ import print_function
import os, re
from .util import texpand

class PlotParser(object):
    """
    Reads Rivet's .plot files and determines which attributes to apply to each histo path.
    """

    pat_begin_block = re.compile(r'^(#*\s*)?BEGIN (\w+) ?(\S+)?')
    pat_end_block =   re.compile(r'^(#*\s*)?END (\w+)')
    pat_comment = re.compile(r'^\s*#|^\s*$')
    pat_property = re.compile(r'^(\w+?)\s*=\s*(.*)$')
    pat_path_property  = re.compile(r'^(\S+?)::(\w+?)=(.*)$')
    pat_paths = {}

    def __init__(self, plotpaths=None, addfiles=[]):
        """
        Parameters
        ----------
        plotpaths : list of str, optional
            The directories to search for .plot files.
            The default is to call the rivet.getAnalysisPlotPaths() function to get
            the directory where the .plot files can be found. (Usually equivalent to calling :command:`rivet-config --datadir`)

        Raises
        ------
        ValueError
            If `plotpaths` is not specified and calling
            :command:`rivet-config` fails.
        """
        self.addfiles = addfiles

        self.plotpaths = plotpaths
        if not self.plotpaths:
            try:
                import rivet
                self.plotpaths = rivet.getAnalysisPlotPaths()
            except Exception as e:
                sys.stderr.write("Failed to load Rivet analysis plot paths: %s\n" % e)
                raise ValueError("No plot paths given and the rivet module could not be loaded!")


    def getSection(self, section, hpath):
        """Get a section for a histogram from a .plot file.

        Parameters
        ----------
        section : ('PLOT'|'SPECIAL'|'HISTOGRAM')
            The section that should be extracted.
        hpath : str
            The histogram path, i.e. /AnalysisID/HistogramID .

        TODO:
         * Caching! The result of the lookup is not cached so every call requires a file to be searched for and opened.
        """
        if section not in ['PLOT', 'SPECIAL', 'HISTOGRAM']:
            raise ValueError("Can't parse section \'%s\'" % section)

        ## Decompose the histo path and remove the /REF prefix if necessary
        from rivet.aopaths import AOPath
        try:
            aop = AOPath(hpath)
        except:
            print("Found analysis object with non-standard path structure:", hpath, "... skipping")
            return None

        ## Assemble the list of headers from any matching plotinfo paths and additional style files
        plotfile = aop.basepathparts()[0] + ".plot"
        ret = {'PLOT': {}, 'SPECIAL': None, 'HISTOGRAM': {}}
        for pidir in self.plotpaths:
            plotpath = os.path.join(pidir, plotfile)
            self._readHeadersFromFile(plotpath, ret, section, aop.basepath())
            ## Only read from the first file we find, otherwise erroneous attributes
            ## in dodgy plotinfo files can't be overridden by user
            if ret[section]: #< neatly excludes both empty dicts and None, used as null defaults above
                break

        ## Also look for further attributes in any user-specified files
        for extrafile in self.addfiles:
            self._readHeadersFromFile(extrafile, ret, section, hpath)
        return ret[section]


    def _readHeadersFromFile(self, plotfile, ret, section, hpath):
        """Get a section for a histogram from a .plot file."""
        if not os.access(plotfile, os.R_OK):
            return
        startreading = False
        f = open(plotfile)
        msec = None
        for line in f:
            m = self.pat_begin_block.match(line)
            if m:
                tag, pathpat = m.group(2,3)
                # pathpat could be a regex
                if pathpat not in self.pat_paths:
                    self.pat_paths[pathpat] = re.compile(pathpat)
                if tag == section:
                    m2 = self.pat_paths[pathpat].match(hpath)
                    if m2:
                        msec = m2
                        startreading = True
                        if section in ['SPECIAL']:
                            ret[section] = ''
                        continue
            if not startreading:
                continue
            if self.isEndMarker(line, section):
                startreading = False
                continue
            elif self.isComment(line):
                continue
            if section in ['PLOT', 'HISTOGRAM']:
                vm = self.pat_property.match(line)
                if vm:
                    prop, value = vm.group(1,2)
                    if msec:
                        oldval = value
                        try:
                            ## First escape backslashes *except* regex groups, then expand regex groups from path match
                            #print("\n", value)
                            value = value.encode("string-escape")
                            #print(value)
                            value = re.sub("(\\\\)(\\d)", "\\2", value) #< r-strings actually made this harder, since the \) is still treated as an escape!
                            #print(value)
                            value = msec.expand(value)
                            #print(value)
                        except Exception as e:
                            #print(e)
                            value = oldval #< roll back escapes if it goes wrong
                    ret[section][prop] = texpand(value) #< expand TeX shorthands
            elif section in ['SPECIAL']:
                ret[section] += line
        f.close()


    def getHeaders(self, hpath):
        """Get the plot headers for histogram hpath.

        This returns the PLOT section.

        Parameters
        ----------
        hpath : str
            The histogram path, i.e. /AnalysisID/HistogramID .

        Returns
        -------
        plot_section : dict
            The dictionary usually contains the 'Title', 'XLabel' and
            'YLabel' properties of the respective plot.

        See also
        --------
        :meth:`getSection`
        """
        return self.getSection('PLOT', hpath)
    ## Alias
    getPlot = getHeaders


    def getSpecial(self, hpath):
        """Get a SPECIAL section for histogram hpath.

        The SPECIAL section is only available in a few analyses.

        Parameters
        ----------
        hpath : str
            Histogram path. Must have the form /AnalysisID/HistogramID .

        See also
        --------
        :meth:`getSection`
        """
        return self.getSection('SPECIAL', hpath)


    def getHistogramOptions(self, hpath):
        """Get a HISTOGRAM section for histogram hpath.

        The HISTOGRAM section is only available in a few analyses.

        Parameters
        ----------
        hpath : str
            Histogram path. Must have the form /AnalysisID/HistogramID .

        See also
        --------
        :meth:`getSection`
        """
        return self.getSection('HISTOGRAM', hpath)


    def isEndMarker(self, line, blockname):
        m = self.pat_end_block.match(line)
        return m and m.group(2) == blockname


    def isComment(self, line):
        return self.pat_comment.match(line) is not None


    def updateHistoHeaders(self, hist):
        headers = self.getHeaders(hist.histopath)
        if "Title" in headers:
            hist.title = headers["Title"]
        if "XLabel" in headers:
            hist.xlabel = headers["XLabel"]
        if "YLabel" in headers:
            hist.ylabel = headers["YLabel"]


def mkStdPlotParser(dirs=None, addfiles=[]):
    """
    Make a PlotParser with the standard Rivet .plot locations automatically added to
    the manually set plot info dirs and additional files.
    """
    if dirs is None:
        dirs = []
    from .core import getAnalysisPlotPaths
    dirs += getAnalysisPlotPaths()
    seen = set()
    dirs = [d for d in dirs if d not in seen and not seen.add(d)]
    return PlotParser(dirs, addfiles)
