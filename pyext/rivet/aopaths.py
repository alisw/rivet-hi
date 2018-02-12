def isRefPath(path):
    return path.startswith("/REF")

def isRefAO(ao):
    return int(ao.annotation("IsRef")) == 1 or isRefPath(ao.path)

def isTmpPath(path):
    return "/_" in path #< match *any* underscore-prefixed path component

def isTmpAO(ao):
    return isTmpPath(ao.path)


class AOPath(object):
    """
    Object representation of analysis object path structures.

    TODO: move to YODA?
    """
    import re
    re_aopath = re.compile(r"^(/[^\[\]\@\#]+)(\[[A-Za-z\d\._]+\])?(#\d+|@[\d\.]+)?$")

    def __init__(self, path):
        import os
        self.origpath = path
        m = self.re_aopath.match(path)
        if not m:
            raise Exception("Supplied path '%s' does not meet required structure" % path)
        self._basepath = m.group(1)
        self._varid = m.group(2).lstrip("[").rstrip("]") if m.group(2) else None
        self._binid = int(m.group(3).lstrip("#")) if m.group(3) else None
        self._isref = isRefPath(self._basepath)

    def basepath(self, keepref=False):
        "Main 'Unix-like' part of the AO path, optionally including a /REF prefix"
        p = self._basepath.rstrip("/")
        if not keepref and p.startswith("/REF"):
            p = p[4:]
        return p

    def varpath(self, keepref=False, defaultvarid=None):
        "The basepath, plus any bracketed variation identifier"
        p = self.basepath(keepref)
        if self.varid(defaultvarid) is not None:
            p += "[%s]" % str(self.varid(defaultvarid))
        return p

    def binpath(self, keepref=False, defaultbinid=None, defaultvarid=None):
        "The varpath, plus any #-prefixed bin number identifier"
        p = self.varpath(keepref, defaultvarid)
        if self.binid(defaultbinid) is not None:
            p += "#%d" % self.binid(defaultbinid)
        return p

    def basepathparts(self, keepref=False):
        "List of basepath components, split by forward slashes"
        return self.basepath(keepref).strip("/").split("/")

    # TODO: basepathhead, basepathtail

    def dirname(self, keepref=False):
        "The non-final (i.e. dir-like) part of the basepath"
        return os.path.dirname(self.basepath(keepref))

    def dirnameparts(self, keepref=False):
        "List of dirname components, split by forward slashes"
        return self.dirname(keepref).strip("/").split("/")

    def basename(self):
        "The final (i.e. file-like) part of the basepath"
        return os.path.basename(self._basepath)

    def varid(self, default=None):
        "The variation identifier (without brackets) if there is one, otherwise None"
        return self._varid if self._varid is not None else default

    def binid(self, default=None):
        "The bin identifier (without #) if there is one, otherwise None"
        return self._binid if self._binid is not None else default

    def isref(self):
        "Is there a /REF prefix in the original path?"
        return self._isref

    def istmp(self):
        "Do any basepath components start with an underscore, used to hide them from plotting?"
        return isTmpPath(self.basepath())
