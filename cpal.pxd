
cdef extern from "palmac.h":
    cdef double PAL__D2PI
    cdef double PAL__DAS2R
    cdef double PAL__DD2R
    cdef double PAL__DH2R
    cdef double PAL__DPI
    cdef double PAL__DPIBY2
    cdef double PAL__DR2AS
    cdef double PAL__DR2D
    cdef double PAL__DR2H
    cdef double PAL__DR2S
    cdef double PAL__DS2R

cdef extern from "pal.h":
    void palAddet ( double rm, double dm, double eq, double *rc, double *dc )
    double palAirmas ( double zd )
    void palAltaz ( double ha, double dec, double phi,
		double *az, double *azd, double *azdd,
		double *el, double *eld, double *eldd,
		double *pa, double *pad, double *padd )
    void palAmp ( double ra, double da, double date, double eq,
              double *rm, double *dm )
    void palAmpqk ( double ra, double da, double amprms[21],
                double *rm, double *dm )
    void palAop ( double rap, double dap, double date, double dut,
              double elongm, double phim, double hm, double xp,
              double yp, double tdk, double pmb, double rh,
              double wl, double tlr,
              double *aob, double *zob, double *hob,
              double *dob, double *rob )
    void palAoppa ( double date, double dut, double elongm, double phim,
                double hm, double xp, double yp, double tdk, double pmb,
                double rh, double wl, double tlr, double aoprms[14] )
    void palAoppat ( double date, double aoprms[14] )
    void palAopqk ( double rap, double dap, double aoprms[14],
                double *aob, double *zob, double *hob,
                double *dob, double *rob )
    void palAtmdsp ( double tdk, double pmb, double rh, double wl1,
                double a1, double b1, double wl2, double *a2, double *b2 )
    void palCaldj ( int iy, int im, int id, double *djm, int *j )
    void palCldj ( int iy, int im, int id, double *djm, int *j )
    void palDaf2r ( int ideg, int iamin, double asec, double *rad, int *j )
    void palDafin ( char *string, int *iptr, double *a, int *j )
    double palDat ( double dju )
    void palDav2m ( double axvec[3], double rmat[3][3] )
    double palDbear ( double a1, double b1, double a2, double b2 )
    void palDcc2s ( double v[3], double *a, double *b )
    void palDcs2c ( double a, double b, double v[3] )
    void palDd2tf ( int ndp, double days, char *sign, int ihmsf[4] )
    void palDe2h ( double ha, double dec, double phi,
                   double *az, double *el )
    void palDeuler ( const char *order, double phi, double theta, double psi,
                     double rmat[3][3] )
        #    void palDfltin ( const char string[], int *nstrt, double *dreslt, int *jflag )
    void palDh2e ( double az, double el, double phi, double *ha, double *dec)
    void palDimxv ( double dm[3][3], double va[3], double vb[3] )
    void palDjcal ( int ndp, double djm, int iymdf[4], int *j )
    void palDjcl ( double djm, int *iy, int *im, int *id, double *fd, int *j )
    void palDm2av ( double rmat[3][3], double axvec[3] )
    void palDmat ( int n, double *a, double *y, double *d, int *jf, int *iw )
    void palDmoon ( double date, double pv[6] )
    void palDmxm ( double a[3][3], double b[3][3], double c[3][3] )
    void palDmxv ( double dm[3][3], double va[3], double vb[3] )
    double palDpav ( double v1[3], double v2[3] )
    void palDr2af ( int ndp, double angle, char *sign, int idmsf[4] )
    void palDr2tf ( int ndp, double angle, char *sign, int ihmsf[4] )
    double palDrange ( double angle )
    double palDranrm ( double angle )
    void palDs2tp ( double ra, double dec, double raz, double decz,
                double *xi, double *eta, int *j )
    double palDsep ( double a1, double b1, double a2, double b2 )
    double palDsepv ( double v1[3], double v2[3] )
    double palDt ( double epoch )
    void palDtf2d ( int ihour, int imin, double sec, double *days, int *j )
    void palDtf2r ( int ihour, int imin, double sec, double *rad, int *j )
    void palDtp2s ( double xi, double eta, double raz, double decz,
                    double *ra, double *dec )
    void palDtps2c ( double xi, double eta, double ra, double dec,
                     double *raz1, double *decz1,
                     double *raz2, double *decz2, int *n )
    double palDtt ( double dju )
    double palDvdv ( double va[3], double vb[3] )
    void palDvn ( double v[3], double uv[3], double *vm )
    void palDvxv ( double va[3], double vb[3], double vc[3] )
    void palEcleq ( double dl, double db, double date, double *dr, double *dd )
    void palEcmat ( double date, double rmat[3][3] )
    void palEl2ue ( double date, int jform, double epoch, double orbinc,
                    double anode, double perih, double aorq, double e,
                    double aorl, double dm, double u[13], int *jstat )
    double palEpb ( double date )
    double palEpb2d ( double epb )
    double palEpco ( char k0, char k, double e )
    double palEpj ( double date )
    double palEpj2d ( double epj )
    void palEpv( double date, double ph[3], double vh[3],
                 double pb[3], double vb[3] )
    void palEqecl ( double dr, double dd, double date, double *dl, double *db )
    double palEqeqx ( double date )
    void palEqgal ( double dr, double dd, double *dl, double *db )
    void palEtrms ( double ep, double ev[3] )
    void palEvp ( double date, double deqx,
                  double dvb[3], double dpb[3],
                  double dvh[3], double dph[3] )
    void palFk45z ( double r1950, double d1950, double bepoch,
                    double *r2000, double *d2000 )
    void palFk524 ( double r2000, double d2000, double dr2000,
                double dd2000, double p2000, double v2000,
                double *r1950, double *d1950, double *dr1950,
                double *dd1950, double *p1950, double *v1950 )
    void palFk54z ( double r2000, double d2000, double bepoch,
                double *r1950, double *d1950,
                double *dr1950, double *dd1950 )
    void palFk5hz ( double r5, double d5, double epoch,
                double *rh, double *dh )
    void palGaleq ( double dl, double db, double *dr, double *dd )
    void palGalsup ( double dl, double db, double *dsl, double *dsb )
    void palGe50 ( double dl, double db, double *dr, double *dd )
    void palGeoc ( double p, double h, double *r, double *z )
    double palGmst ( double ut1 )
    double palGmsta ( double date, double ut1 )
    void palHfk5z ( double rh, double dh, double epoch,
                double *r5, double *d5, double *dr5, double *dd5 )
    void palIntin ( char *string, int *nstrt, long *ireslt, int *jflag )
    void palMap ( double rm, double dm, double pr, double pd,
              double px, double rv, double eq, double date,
              double *ra, double *da )
    void palMappa ( double eq, double date, double amprms[21] )
    void palMapqk ( double rm, double dm, double pr, double pd,
                double px, double rv, double amprms[21],
                double *ra, double *da )
    void palMapqkz ( double rm, double dm, double amprms[21],
                 double *ra, double *da )
    void palNut ( double date, double rmatn[3][3] )
    void palNutc ( double date, double *dpsi, double *deps, double *eps0 )
    void palOap ( char *type, double ob1, double ob2, double date,
             double dut, double elongm, double phim, double hm,
             double xp, double yp, double tdk, double pmb,
             double rh, double wl, double tlr,
             double *rap, double *dap )
    void palOapqk (char *type, double ob1, double ob2, double aoprms[14],
               double *rap, double *dap )
    int palObs( size_t n, char * c,
                char * ident, size_t identlen,
                char * name, size_t namelen,
                double * w, double * p, double * h )
    double palPa ( double ha, double dec, double phi )
    void palPcd( double disco, double *x, double *y )
    void palPertel (int jform, double date0, double date1,
                double epoch0, double orbi0, double anode0,
                double perih0, double aorq0, double e0, double am0,
                double *epoch1, double *orbi1, double *anode1,
                double *perih1, double *aorq1, double *e1, double *am1,
                int *jstat )
    void palPertue ( double date, double u[13], int *jstat )
    void palPlanel ( double date, int jform, double epoch, double orbinc,
                 double anode, double perih, double aorq,  double e,
                 double aorl, double dm, double pv[6], int *jstat )
    void palPlanet ( double date, int np, double pv[6], int *j )
    void palPlante ( double date, double elong, double phi, int jform,
                 double epoch, double orbinc, double anode, double perih,
                 double aorq, double e, double aorl, double dm,
                 double *ra, double *dec, double *r, int *jstat )
    void palPlantu ( double date, double elong, double phi, double u[13],
               double *ra, double *dec, double *r, int *jstat )
    void palPolmo( double elongm, double phim, double xp, double yp,
                   double *elong, double *phi, double *daz )
    void palPm ( double r0, double d0, double pr, double pd,
             double px, double rv, double ep0, double ep1,
             double *r1, double *d1 )
    void palPrebn ( double bep0, double bep1, double rmatp[3][3] )
    void palPrec ( double ep0, double ep1, double rmatp[3][3] )
    void palPreces ( char sys[3], double ep0, double ep1,
                     double *ra, double *dc )
    void palPrenut ( double epoch, double date, double rmatpn[3][3] )
    void palPv2el ( const double pv[6], double date, double pmass, int jformr,
                int *jform, double *epoch, double *orbinc,
                double *anode, double *perih, double *aorq, double *e,
                double *aorl, double *dm, int *jstat )
    void palPv2ue ( const double pv[6], double date, double pmass,
                double u[13], int *jstat )
    void palPvobs ( double p, double h, double stl, double pv[6] )
    void palRdplan ( double date, int np, double elong, double phi,
                 double *ra, double *dec, double *diam )
    void palRefco ( double hm, double tdk, double pmb, double rh,
                double wl, double phi, double tlr, double eps,
                double *refa, double *refb )
    void palRefcoq ( double tdk, double pmb, double rh, double wl,
                double *refa, double *refb )
    void palRefro ( double zobs, double hm, double tdk, double pmb,
                double rh, double wl, double phi, double tlr, double eps,
                double *ref )
    void palRefv ( double vu[3], double refa, double refb, double vr[3] )
    void palRefz ( double zu, double refa, double refb, double *zr )
    double palRverot ( double phi, double ra, double da, double st )
    double palRvgalc ( double r2000, double d2000 )
    double palRvlg ( double r2000, double d2000 )
    double palRvlsrd ( double r2000, double d2000 )
    double palRvlsrk ( double r2000, double d2000 )
    void palSubet ( double rc, double dc, double eq,
                double *rm, double *dm )
    void palSupgal ( double dsl, double dsb, double *dl, double *db )
    void palUe2el ( const double u[13], int jformr,
                int *jform, double *epoch, double *orbinc,
                double *anode, double *perih, double *aorq, double *e,
                double *aorl, double *dm, int *jstat )
    void palUe2pv ( double date, double u[13], double pv[], int *jstat )
    void palUnpcd( double disco, double *x, double *y )
