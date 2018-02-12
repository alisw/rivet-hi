#! /usr/bin/env python

import logging, re
try:
    from urllib.request import urlopen
except ImportError:
    from urllib2 import urlopen

usage = """%prog <spiresid> [<spiresid2> ...]

Given Inspire and SPIRES paper IDs, fetch the corresponding BibTeX db entry from
the SPIRES Web interface and write it to stdout. Prefix the code with I or S
appropriately.
"""

def fetch_bibtex(iscode, refid):
    if iscode.upper() == "I":
        url = "http://inspire-hep.net/record/%s/export/hx" % str(refid)
        logging.debug("Downloading Inspire BibTeX from %s" % url)
    elif iscode.upper() == "S":
        url = "http://inspire-hep.net/search?p=find+key+%s&of=hx" % str(refid)
        logging.debug("Downloading SPIRES BibTeX from %s" % url)
    hreq = urlopen(url)
    bibtexhtml = hreq.read()
    hreq.close()
    #logging.debug(bibtexhtml)
    return bibtexhtml


def extract_bibtex(html):
    ## Extract BibTeX block from HTML
    re_spiresbibtex = re.compile(r'<pre>(.*?)</pre>', re.MULTILINE | re.DOTALL)
    m = re_spiresbibtex.search(html)
    if m is None:
        return None, None
    bib = m.group(1).strip()

    ## Get BibTeX key
    re_bibtexkey = re.compile(r'^@.+?{(.+?),$', re.MULTILINE)
    m = re_bibtexkey.search(bib)
    if m is None:
        return None, bib
    key = m.group(1)

    ## Return key and BibTeX
    return key, bib


def get_bibtex_from_repo(iscode, refid):
    html = fetch_bibtex(iscode, refid)
    key, bibtex = extract_bibtex(html)
    return key, bibtex


def get_bibtexs_from_repos(iscodes_refids):
    bibdb = {}
    for iscode, refid in iscodes_refids:
        key, bibtex = get_bibtex_from_repo(iscode, refid)
        if key and bibtex:
            bibdb[refid] = (key, bibtex)
    return bibdb


if __name__ == '__main__':
    ## Parse command line options
    from optparse import OptionParser
    parser = OptionParser(usage=usage)
    opts, args = parser.parse_args()

    ## Make individual bibinfo files
    for arg in args:
        iscode = arg[0]
        refid = arg[1:]
        key, bibtex = get_bibtex_from_repo(iscode, refid)
        import sys
        f = sys.stdout
        f.write("BibKey: %s\n" % key)
        f.write("BibTeX: '%s'\n" % bibtex)

    # ## Build ref db
    # bibdb = get_bibtexs_from_spires(args)
    # for sid, (key, bibtex) in bibdb.iteritems():
    #     print key, "=>\n", bibtex

    # ## Pickle ref db9151176
    # import cPickle as pickle
    # fpkl = open("spiresbib.pkl", "w")repo
    # pickle.dump(bibdb)
    # fpkl.close()
