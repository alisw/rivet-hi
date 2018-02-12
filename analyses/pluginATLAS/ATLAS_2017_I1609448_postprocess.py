
import yoda, rivet, os
from math import sqrt

inFile  = 'OUTPUT.yoda'

hists = yoda.read( inFile )
tags = sorted(hists.keys())

# probably more fancy than it needs to be ...
def getRivetRefData(anas=None):
    "Find all Rivet reference data files"
    refhistos = {}
    for var in ("RIVET_ANALYSIS_PATH", "RIVET_DATA_PATH", "RIVET_REF_PATH", "RIVET_INFO_PATH", "RIVET_PLOT_PATH"):
        if var in os.environ:
            abspaths = map(os.path.abspath, os.environ[var].split(":"))
            os.environ[var] = ":".join(abspaths)
    rivet_data_dirs = rivet.getAnalysisRefPaths()
    dirlist = [ ]
    for d in rivet_data_dirs:
        import glob
        if anas is None:
          dirlist.append(glob.glob(os.path.join(d, '*.yoda')))
        else:
          for a in anas:
            dirlist.append(glob.glob(os.path.join(d, a+'*.yoda')))
    for filelist in dirlist:
        for infile in filelist:
            analysisobjects = yoda.read(infile)
            for path, ao in analysisobjects.iteritems():
                if not any(x in path for x in ['d01', 'd02', 'd03', 'd04']):  continue
                aop = rivet.AOPath(ao.path)
                if aop.isref():
                    ao.path = aop.basepath(keepref=False)
                    refhistos[ao.path] = ao
    return refhistos

# get hold of relevant objects in reference data files
refhistos = getRivetRefData(['ATLAS_2017_I1609448'])

def constructRmiss(hist):
    '''This recreates the constructRmiss function from the routine.'''
    rtn = yoda.mkScatter(hist)
    path = hist.annotation('Path').replace('_d', 'd')
    numer = refhistos[path.replace('y02', 'y03')]
    denom = refhistos[path.replace('y02', 'y04')]
    rmiss = refhistos[path]
    for i in range(rtn.numPoints):
        newy = (numer.points[i].y + rtn.points[i].y) / denom.points[i].y if denom.points[i].y else 0.0
        # ratio error (Rmiss = SM_num/SM_denom + BSM/SM_denom ~ Rmiss_SM + BSM/SM_denom
        rel_hist_err = rtn.points[i].yErrs[0] / denom.points[i].y if denom.points[i].y else 0.0
        newey = sqrt(rmiss.points[i].yErrs[0] ** 2 + rel_hist_err ** 2)
        rtn.points[i].y = newy
        rtn.points[i].yErrs = (newey, newey)
    return rtn

# this is where the magic happens
f = open('%s_processed.yoda' % inFile[:-5], 'w')
for h in tags:
    if 'y02' in h:
        outName = h.replace('_d', 'd')
        outName = outName.replace('y02', 'y01')
        rmiss = constructRmiss(hists[h])
        rmiss.setAnnotation('Path', outName)
        yoda.writeYODA(rmiss, f)
f.close()

