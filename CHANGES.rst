1.8.0
-----

- Upgrade PAL to v0.9.4
- Upgrade ERFA to v1.3.0

1.7.0
-----

- Add test for fk524 (thanks to Scott Daniel)
- Addition of many vectorized routines (thanks to Scott Daniel)
- Upgrade PAL to v0.9.3
- Upgrade ERFA to v1.2.0

1.6.0
-----

- Add palDr2tf
- Added vectorized forms for altaz, pa, gmsta, gmst, eqeqx, refz, galeq,
  dsep, eqgal (Thanks to Scott Daniel for these).

1.5.0
-----

- Add palPcd and palUnpcd (PAL v0.9.0)

1.4.1
-----

Fix a build problem on some systems. A semi-colon had crept in to
cpal.pxd

1.4
---

- Include PAL v0.8.0
- Add ecleq (Ecliptic to J2000)

1.3
---

- Include PAL v0.7.0.
- Add polmo (polar motion).

1.2
---

- Include PAL v0.6+

- Include ERFA v1.1.1

- Submodules now use explicit https URLs to make them more robust
  (especially when forking the repo) (thanks to Srini Chandrasekharan).

- New interfaces added: altaz, dmxv (thanks to Srini Chandrasekharan
  for those 2), dav2m, deuler, dmxm, dimxv, dm2av, dvn and dvxv.

- Some routines now have "Vector" forms to allow numpy arrays to be
  given as arguments and to return numpy arrays with results. This
  provides a signficant performance enhancement over looping in
  Python code: aopqkVector, de2hVector, ds2tpVector, mapqkVector,
  mapqkzVector, pmVector (Thanks to Scott Daniel for these).

- Explicitly include GPLv3 text in distribution.

- New routine palvers() can be used to obtain the version
  information of the underyling C PAL library.

- Python docstring documentation is now included for all
  palpy routines.


1.1
---

- Includes PAL v0.5.0

- Now ships with ERFA v1.1 rather than SOFA. This clarifies
  the license of the included code.

- New interfaces added for refraction and moving sources. The new
  functions are: epv, evp, el2ue, pertel, pertue, planel, planet,
  plante, plantu, pv2el, pv2ue, rdplan, ue2el, ue2pv, atmdsp,
  refcoq, refro, refv, aopqk, prenut, refco and refz.
  Thanks to Scott Daniel for adding the final four in that list.

1.0
---

First release of palpy.
