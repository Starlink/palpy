
import os
import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

sofa_c = (
    "a2af.c", "a2tf.c", "af2a.c", "anp.c", "anpm.c", "bi00.c",
    "bp00.c", "bp06.c", "bpn2xy.c", "c2i00a.c", "c2i00b.c",
    "c2i06a.c", "c2ibpn.c", "c2ixy.c", "c2ixys.c", "c2s.c",
    "c2t00a.c", "c2t00b.c", "c2t06a.c", "c2tcio.c", "c2teqx.c",
    "c2tpe.c", "c2txy.c", "cal2jd.c", "cp.c", "cpv.c", "cr.c",
    "d2dtf.c", "d2tf.c", "dat.c", "dtdb.c", "dtf2d.c", "ee00.c",
    "ee00a.c", "ee00b.c", "ee06a.c", "eect00.c", "eform.c",
    "eo06a.c", "eors.c", "epb.c", "epb2jd.c", "epj.c",
    "epj2jd.c", "epv00.c", "eqeq94.c", "era00.c", "fad03.c",
    "fae03.c", "faf03.c", "faju03.c", "fal03.c", "falp03.c",
    "fama03.c", "fame03.c", "fane03.c", "faom03.c", "fapa03.c",
    "fasa03.c", "faur03.c", "fave03.c", "fk52h.c", "fk5hip.c",
    "fk5hz.c", "fw2m.c", "fw2xy.c", "gc2gd.c", "gc2gde.c",
    "gd2gc.c", "gd2gce.c", "gmst00.c", "gmst06.c", "gmst82.c",
    "gst00a.c", "gst00b.c", "gst06.c", "gst06a.c", "gst94.c",
    "h2fk5.c", "hfk5z.c", "ir.c", "jd2cal.c", "jdcalf.c",
    "num00a.c", "num00b.c", "num06a.c", "numat.c", "nut00a.c",
    "nut00b.c", "nut06a.c", "nut80.c", "nutm80.c", "obl06.c",
    "obl80.c", "p06e.c", "p2pv.c", "p2s.c", "pap.c", "pas.c",
    "pb06.c", "pdp.c", "pfw06.c", "plan94.c", "pm.c",
    "pmat00.c", "pmat06.c", "pmat76.c", "pmp.c", "pn.c",
    "pn00.c", "pn00a.c", "pn00b.c", "pn06.c", "pn06a.c",
    "pnm00a.c", "pnm00b.c", "pnm06a.c", "pnm80.c", "pom00.c",
    "ppp.c", "ppsp.c", "pr00.c", "prec76.c", "pv2p.c", "pv2s.c",
    "pvdpv.c", "pvm.c", "pvmpv.c", "pvppv.c", "pvstar.c",
    "pvu.c", "pvup.c", "pvxpv.c", "pxp.c", "rm2v.c", "rv2m.c",
    "rx.c", "rxp.c", "rxpv.c", "rxr.c", "ry.c", "rz.c", "s00.c",
    "s00a.c", "s00b.c", "s06.c", "s06a.c", "s2c.c", "s2p.c",
    "s2pv.c", "s2xpv.c", "sepp.c", "seps.c",
    "sp00.c", "starpm.c", "starpv.c", "sxp.c",
    "sxpv.c", "taitt.c", "taiut1.c", "taiutc.c",
    "tcbtdb.c", "tcgtt.c", "tdbtcb.c", "tdbtt.c", "tf2a.c",
    "tf2d.c", "tr.c", "trxp.c", "trxpv.c", "tttai.c", "tttcg.c",
    "tttdb.c", "ttut1.c", "ut1tai.c", "ut1tt.c", "ut1utc.c",
    "utctai.c", "utcut1.c", "xy06.c", "xys00a.c", "xys00b.c",
    "xys06a.c", "zp.c", "zpv.c", "zr.c"
)

pal_c = (
    "palAddet.c", "palAirmas.c",
    "palAmp.c", "palAmpqk.c", "palAop.c", "palAoppa.c",
    "palAoppat.c", "palAopqk.c", "palCaldj.c", "palDafin.c",
    "palDat.c", "palDe2h.c", "palDeuler.c", "palDfltin.c",
    "palDh2e.c", "palDjcal.c", "palDmat.c", "palDmoon.c",
    "palDrange.c", "palDs2tp.c", "palDt.c", "palDtp2s.c",
    "palDtps2c.c", "palDtt.c", "palEcmat.c", "palEl2ue.c",
    "palEpco.c", "palEpv.c", "palEqecl.c", "palEqgal.c",
    "palEtrms.c", "palEvp.c", "palFk45z.c", "palFk524.c",
    "palFk54z.c", "palGaleq.c", "palGalsup.c", "palGe50.c",
    "palGeoc.c", "palIntin.c", "palMap.c", "palMappa.c",
    "palMapqk.c", "palMapqkz.c", "palNut.c", "palNutc.c",
    "palOap.c", "palOapqk.c", "palObs.c", "palOne2One.c",
    "palPa.c", "palPertel.c", "palPertue.c", "palPlanel.c",
    "palPlanet.c", "palPlante.c", "palPlantu.c", "palPm.c",
    "palPrebn.c", "palPrec.c", "palPreces.c", "palPrenut.c",
    "palPv2el.c", "palPv2ue.c", "palPvobs.c", "palRdplan.c",
    "palRefco.c", "palRefro.c", "palRefz.c", "palRverot.c",
    "palRvgalc.c", "palRvlg.c", "palRvlsrd.c", "palRvlsrk.c",
    "palSubet.c", "palSupgal.c", "palUe2el.c",
    "palUe2pv.c",
    "pal1Atms.c", "pal1Atmt.c"
)

# Build up source file list
sources = [ "cpal.pxd", "pal.pyx" ]

# Sort out path to the C files
for cfile in sofa_c:
    sources.append( os.path.join( 'cextern', 'sofa', 'src', cfile ) )

for cfile in pal_c:
    sources.append( os.path.join( 'cextern', 'pal', cfile ) )


setup(
    name = "palpy",
    version = "1.0",
    author = "Tim Jenness",
    author_email = "tim.jenness@gmail.com",
    url='https://github.com/Starlink/palpy',
    description = "PAL -- A Positional Astronomy Library",
    cmdclass = { 'build_ext': build_ext },
    ext_modules = [ Extension(
        name="palpy",
        sources=sources,
        include_dirs=['cextern/sofa/src', 'cextern/pal', numpy.get_include() ],
        language="c")
    ],
    requires=[
        'numpy',
        'Cython'
        ],
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Programming Language :: Python',
        'Programming Language :: C',
        'Topic :: Scientific/Engineering :: Astronomy'
        ]
)

