
import os
import numpy
import re
import codecs
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# Local code
from support import sst2pydoc as sst
from support import palvers

erfa_c = (
    "a2af.c", "a2tf.c", "ab.c", "af2a.c", "anp.c", "anpm.c",
    "apcg.c", "apcg13.c", "apci.c", "apci13.c", "apco.c",
    "apco13.c", "apcs.c", "apcs13.c", "aper.c", "aper13.c",
    "apio.c", "apio13.c", "atci13.c", "atciq.c", "atciqn.c",
    "atciqz.c", "atco13.c", "atic13.c", "aticq.c", "aticqn.c",
    "atio13.c", "atioq.c", "atoc13.c", "atoi13.c", "atoiq.c",
    "bi00.c", "bp00.c", "bp06.c", "bpn2xy.c", "c2i00a.c",
    "c2i00b.c", "c2i06a.c", "c2ibpn.c", "c2ixy.c", "c2ixys.c",
    "c2s.c", "c2t00a.c", "c2t00b.c", "c2t06a.c", "c2tcio.c",
    "c2teqx.c", "c2tpe.c", "c2txy.c", "cal2jd.c", "cp.c", "cpv.c",
    "cr.c", "d2dtf.c", "d2tf.c", "dat.c", "dtdb.c", "dtf2d.c",
    "ee00.c", "ee00a.c", "ee00b.c", "ee06a.c", "eect00.c",
    "eform.c", "eo06a.c", "eors.c", "epb.c", "epb2jd.c", "epj.c",
    "epj2jd.c", "epv00.c", "eqeq94.c", "era00.c", "fad03.c",
    "fae03.c", "faf03.c", "faju03.c", "fal03.c", "falp03.c",
    "fama03.c", "fame03.c", "fane03.c", "faom03.c", "fapa03.c",
    "fasa03.c", "faur03.c", "fave03.c", "fk52h.c", "fk5hip.c",
    "fk5hz.c", "fw2m.c", "fw2xy.c", "gc2gd.c", "gc2gde.c",
    "gd2gc.c", "gd2gce.c", "gmst00.c", "gmst06.c", "gmst82.c",
    "gst00a.c", "gst00b.c", "gst06.c", "gst06a.c", "gst94.c",
    "h2fk5.c", "hfk5z.c", "ir.c", "jd2cal.c", "jdcalf.c", "ld.c",
    "ldn.c", "ldsun.c", "num00a.c", "num00b.c", "num06a.c",
    "numat.c", "nut00a.c", "nut00b.c", "nut06a.c", "nut80.c",
    "nutm80.c", "obl06.c", "obl80.c", "p06e.c", "p2pv.c", "p2s.c",
    "pap.c", "pas.c", "pb06.c", "pdp.c", "pfw06.c", "plan94.c",
    "pm.c", "pmat00.c", "pmat06.c", "pmat76.c", "pmp.c", "pmpx.c",
    "pmsafe.c", "pn.c", "pn00.c", "pn00a.c", "pn00b.c", "pn06.c",
    "pn06a.c", "pnm00a.c", "pnm00b.c", "pnm06a.c", "pnm80.c",
    "pom00.c", "ppp.c", "ppsp.c", "pr00.c", "prec76.c", "pv2p.c",
    "pv2s.c", "pvdpv.c", "pvm.c", "pvmpv.c", "pvppv.c",
    "pvstar.c", "pvtob.c", "pvu.c", "pvup.c", "pvxpv.c", "pxp.c",
    "refco.c", "rm2v.c", "rv2m.c", "rx.c", "rxp.c", "rxpv.c",
    "rxr.c", "ry.c", "rz.c", "s00.c", "s00a.c", "s00b.c", "s06.c",
    "s06a.c", "s2c.c", "s2p.c", "s2pv.c", "s2xpv.c", "sepp.c",
    "seps.c", "sp00.c", "starpm.c", "starpv.c", "sxp.c", "sxpv.c",
    "taitt.c", "taiut1.c", "taiutc.c", "tcbtdb.c", "tcgtt.c",
    "tdbtcb.c", "tdbtt.c", "tf2a.c", "tf2d.c", "tr.c", "trxp.c",
    "trxpv.c", "tttai.c", "tttcg.c", "tttdb.c", "ttut1.c",
    "ut1tai.c", "ut1tt.c", "ut1utc.c", "utctai.c", "utcut1.c",
    "xy06.c", "xys00a.c", "xys00b.c", "xys06a.c", "zp.c", "zpv.c",
    "zr.c",
)

pal_c = (
    "pal1Atms.c", "pal1Atmt.c", "palAddet.c", "palAirmas.c",
    "palAltaz.c", "palAmp.c", "palAmpqk.c", "palAop.c", "palAoppa.c",
    "palAoppat.c", "palAopqk.c", "palAtmdsp.c", "palCaldj.c",
    "palDafin.c", "palDat.c", "palDe2h.c", "palDeuler.c",
    "palDfltin.c", "palDh2e.c", "palDjcal.c", "palDmat.c",
    "palDmoon.c", "palDrange.c", "palDs2tp.c", "palDt.c",
    "palDtp2s.c", "palDtps2c.c", "palDtt.c", "palEcleq.c", "palEcmat.c",
    "palEl2ue.c", "palEpco.c", "palEpv.c", "palEqecl.c",
    "palEqgal.c", "palEtrms.c", "palEvp.c", "palFk45z.c",
    "palFk524.c", "palFk54z.c", "palGaleq.c", "palGalsup.c",
    "palGe50.c", "palGeoc.c", "palIntin.c", "palMap.c",
    "palMappa.c", "palMapqk.c", "palMapqkz.c", "palNut.c",
    "palNutc.c", "palOap.c", "palOapqk.c", "palObs.c",
    "palOne2One.c", "palPa.c", "palPcd.c", "palPertel.c", "palPertue.c",
    "palPlanel.c", "palPlanet.c", "palPlante.c", "palPlantu.c",
    "palPm.c", "palPolmo.c", "palPrebn.c", "palPrec.c", "palPreces.c",
    "palPrenut.c", "palPv2el.c", "palPv2ue.c", "palPvobs.c",
    "palRdplan.c", "palRefco.c", "palRefro.c", "palRefv.c",
    "palRefz.c", "palRverot.c", "palRvgalc.c", "palRvlg.c",
    "palRvlsrd.c", "palRvlsrk.c", "palSubet.c", "palSupgal.c",
    "palUe2el.c", "palUe2pv.c", "palUnpcd.c"
)

# Build up source file list
sources = ["cpal.pxd", "pal.pyx"]

# Sort out path to the C files
for cfile in erfa_c:
    sources.append(os.path.join('cextern', 'erfa', 'src', cfile))

# Also read in prologues
palprologs = {}
for cfile in pal_c:
    palfile = os.path.join('cextern', 'pal', cfile)
    sources.append(palfile)
    prologs = sst.read_prologs(palfile)
    palprologs.update(prologs)

# And the version of the C PAL library
(verstring, major, minor, patch) = palvers.read_pal_version()

# Generate pal.pyx if required
infile = "pal.pyx.in"
outfile = "pal.pyx"
create_outfile = True
if os.path.exists(outfile) and \
        os.path.getmtime(infile) < os.path.getmtime(outfile):
    create_outfile = False

if create_outfile:
    paldoc_re = re.compile(r"@(pal.*)@")
    outfh = codecs.open(outfile, "w", "utf8")
    with codecs.open(infile, 'r', 'utf8') as file:
        for line in file.readlines():
            match_paldoc = paldoc_re.search(line)
            if match_paldoc is not None:
                lines = []
                funcname = match_paldoc.group(1)
                palpyname = funcname[3:].lower()
                if funcname in palprologs:
                    info = palprologs[funcname]
                    lines.append(info['purpose'])
                    lines += ("Arguments", "---------", info['arguments'])
                    if "returned value" in info:
                        lines += ("", "Returned Value", "--------------",
                                  info['returned value'])
                    if "notes" in info:
                        lines += ("", "Notes", "-----",
                                  info['notes'])
                    line = "\n".join(lines)
                else:
                    continue
            elif line.count("@"):
                if line.count("@VERSTRING@"):
                    line = line.replace("@VERSTRING@", '"'+verstring+'"')
                if line.count("@MAJOR_VERSION"):
                    line = line.replace("@MAJOR_VERSION@", str(major))
                if line.count("@MINOR_VERSION@"):
                    line = line.replace("@MINOR_VERSION@", str(minor))
                if line.count("@PATCH_VERSION@"):
                    line = line.replace("@PATCH_VERSION@", str(patch))
            outfh.write(line)

    outfh.close()

# Description
with open('README.rst') as file:
    long_description = file.read()

setup(
    name="palpy",
    version="1.8.0",
    author="Tim Jenness",
    author_email="tim.jenness@gmail.com",
    license="GPL",
    url='https://github.com/Starlink/palpy',
    description="PAL -- A Positional Astronomy Library",
    long_description=long_description,
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension(
        name="palpy",
        sources=sources,
        include_dirs=['cextern/erfa/src', 'cextern/pal', numpy.get_include()],
        language="c"
    )],
    requires=[
        'numpy',
        'Cython'
    ],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Programming Language :: Cython',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: C',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering :: Astronomy'
    ]
)
