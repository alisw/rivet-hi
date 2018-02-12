"Python utility functions for use by Rivet scripts (and anyone else who wants to)"


def check_python_version(req_version=(2,6,0)):
    "Enforce the Rivet scripts' minimal Python version requirement"
    import sys
    if sys.version_info[:3] < req_version:
        sys.stderr.write( "Python version >= %s is required... exiting\n" % ".".join(req_version) )
        sys.exit(1)


def set_process_name(name):
    "Try to rename the process on Linux so it doesn't appear as 'python <scriptpath>'"
    try:
        ## Try to use this: https://code.google.com/p/py-setproctitle/
        import setproctitle
        setproctitle.setproctitle(name)
    except:
        try:
            ## Fall back to a by-hand thing that doesn't work for me...
            import ctypes
            libc = ctypes.cdll.LoadLibrary("libc.so.6")
            libc.prctl(15, name, 0, 0, 0)
        except:
            ## And then give up ;-)
            pass


def import_ET():
    "Try to import the ElementTree XML parser, which has many historical import signatures"
    ET = None
    try:
        import xml.etree.cElementTree as ET
    except ImportError:
        try:
            import cElementTree as ET
        except ImportError:
            try:
                import xml.etree.ElementTree as ET
            except:
                raise ImportError("Can't load the ElementTree XML parser (any of three historical ways)")
    return ET


def htmlify(s, para=False):
    """Modify LaTeX text strings from analysis metadata for inclusion
    in MathJax-enabled web page source code."""
    if not s:
        return s
    t = s.replace("&", "&amp;")\
        .replace("<","&lt;")\
        .replace(">","&gt;")\
        .replace(r"~", " ")
    t = t.replace(r"\pT", r"p_\perp")\
        .replace(r"\degree", r"^\circ")\
        .replace(r"\MeV", r"\text{MeV}")\
        .replace(r"\GeV", r"\text{GeV}")\
        .replace(r"\TeV", r"\text{TeV}")
    # t = t.replace(r"\;", " ")\
    #     .replace(r"\,", " ")\
    #     .replace(r"\!", "")
    if para:
        t = t.replace("\n\n", "</p><p>")
    return t


def texify(s):
    "Insert required TeX escapes"
    if not s:
        return s
    t = s \
        .replace(r"&", r"\&") \
        .replace(r"\\&", r"\&") \
        .replace(r"#", r"\#") \
        # .replace(r"_", r"\_") \
        # .replace(r"^", r"") \
    return t


def texpand(s):
    "Expand some physics-specific TeX macros."
    if not s:
        return s
    t = s \
        .replace(r"\kT", r"\ensuremath{k_\perp}\xspace") \
        .replace(r"\kt", r"\ensuremath{k_\mathrm{T}}\xspace") \
        .replace(r"\pT", r"\ensuremath{p_\perp}\xspace") \
        .replace(r"\pt", r"\ensuremath{p_\mathrm{T}}\xspace") \
        .replace(r"\sqrts", r"\ensuremath{\sqrt{s}}\xspace") \
        .replace(r"\sqrtS", r"\ensuremath{\sqrt{s}}\xspace") \
        .replace(r"\MeV", r"\text{M\eV}\xspace") \
        .replace(r"\GeV", r"\text{G\eV}\xspace") \
        .replace(r"\TeV", r"\text{T\eV}\xspace") \
        .replace(r"\eV", r"\text{e\kern-0.15ex{}V}\xspace")
    return t


def detex(tex):
    """Use pandoc (if available) to modify LaTeX text strings from
    analysis metadata for use as plain text, e.g. as printed to the terminal.

    The argument can either be a string or an iterable of strings.

    TODO: Replace \gamma, \mu, \tau, \\Upsilon, \rho, \psi, \pi, \eta, \Delta, \Omega, \omega -> no-\ form?
    TODO: Replace e^+- -> e+-?
    """
    if not tex:
        return tex
    from distutils.spawn import find_executable
    if not find_executable("pandoc"):
        return tex

    try:
        tex_is_str = type(tex) in (unicode,str)
    except NameError: # for py3
        tex_is_str = type(tex) is str

    texheader = r"""
    \newcommand{\text}[1]{#1}
    \newcommand{\ensuremath}[1]{#1}
    \newcommand{\emph}[1]{_#1_}
    \newcommand{\textrm}[1]{#1}
    \newcommand{\textit}[1]{_#1_}
    \newcommand{\textbf}[1]{*#1*}
    \newcommand{\mathrm}[1]{#1}
    \newcommand{\mathit}[1]{_#1_}
    \newcommand{\mathbf}[1]{*#1*}
    \newcommand{\bm}[1]{*#1*}
    \newcommand{\frac}[2]{#1/#2}
    \newcommand{\sqrt}[1]{sqrt(#1)}
    \newcommand{\hat}[1]{#1hat}
    \newcommand{\bar}[1]{#1bar}
    \newcommand{\d}[1]{d#1}
    \newcommand{\degree}{^\circ }
    \newcommand{\infty}{oo }
    \newcommand{\exp}{exp }
    \newcommand{\log}{log }
    \newcommand{\ln}{ln }
    \newcommand{\sin}{sin }
    \newcommand{\cos}{cos }
    \newcommand{\tan}{tan }
    \newcommand{\sinh}{sinh }
    \newcommand{\cosh}{cosh }
    \newcommand{\tanh}{tanh }
    \newcommand{\ell}{l}
    \newcommand{\varphi}{\phi}
    \newcommand{\varepsilon}{\epsilon}
    \newcommand{\sim}{~}
    \newcommand{\lesssim}{<~ }
    \newcommand{\gtrsim}{>~ }
    \newcommand{\neq}{!= }
    \newcommand{\ge}{>= }
    \newcommand{\gg}{>> }
    \newcommand{\le}{<= }
    \newcommand{\ll}{<< }
    \newcommand{\pm}{+- }
    \newcommand{\mp}{-+ }
    \newcommand{\times}{x }
    \newcommand{\cdot}{. }
    \newcommand{\dots}{... }
    \newcommand{\ldots}{... }
    \newcommand{\langle}{<}
    \newcommand{\rangle}{>}
    \newcommand{\gets}{<- }
    \newcommand{\to}{-> }
    \newcommand{\leftarrow}{<- }
    \newcommand{\rightarrow}{-> }
    \newcommand{\leftrightarrow}{<-> }
    \newcommand{\Leftarrow}{<= }
    \newcommand{\Rightarrow}{=> }
    \newcommand{\Leftrightarrow}{ }
    \newcommand{\left}{}
    \newcommand{\right}{}
    \newcommand{\!}{}
    \newcommand{\/}{}
    \newcommand{\rm}{}
    \newcommand{\it}{}
    \newcommand{\,}{ }
    \newcommand{\;}{ }
    \newcommand{\ }{ }
    \newcommand{\unit}[2]{#1 #2}
    \newcommand{\bar}[1]{#1bar}
    \newcommand{\pT}{pT }
    \newcommand{\perp}{T}
    \newcommand{\ast}{*}
    \newcommand{\MeV}{MeV }
    \newcommand{\GeV}{GeV }
    \newcommand{\TeV}{TeV }
    """
    import subprocess
    pandoc_cmd = ['pandoc','-f','latex','-t','plain','--wrap=none']
    ### check we have the right wrap option ###
    x = subprocess.Popen(pandoc_cmd,stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = x.communicate(b' ')
    x = x.wait()
    if x != 0:
        pandoc_cmd[-1] = '--no-wrap'
    p = subprocess.Popen(pandoc_cmd,stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    texbody = tex if tex_is_str else "@@".join(tex)
    # texbody = texbody.replace("$", "")
    # pandoc reads in UTF-8, need to encode / decode correctly
    plain, err = p.communicate((texheader + texbody).replace("\n", "").encode('utf-8'))
    ret = p.wait()
    if ret != 0:
        return tex
    # pandoc sends UTF-8, need to decode
    plain = plain.decode('utf-8')
    plain = plain.replace("\n", "")
    plains = plain.replace(r"\&", "&").split("@@")
    if tex_is_str:
        assert len(plains) == 1
        return plains[0] if plains[0] else tex
    else:
        return plains if plains else tex

# print detex(r"Foo \! $\int \text{bar} \d{x} \sim \; \frac{1}{3} \neq \emph{foo}$ \to \gg bar")
# print detex([r"Foo \! $\int \text{bar} \d{x} \sim", r"\frac{1}{3} \neq \emph{foo}$ \to \gg bar"])
