# Test script for PAL interface

# Copyright (C) 2012 Tim Jenness and Science and Technology
# Facilities Council.

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation; either version 3 of
# the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import unittest
import palpy as pal
import numpy as np

class TestPAL(unittest.TestCase) :

    def test_addet(self):
        rm = 2.0
        dm = -1.0
        eq = 1975.0
        ( r1, d1 ) = pal.addet( rm, dm, eq )
        self.assertAlmostEqual( r1 - rm, 2.983864874295250e-6, 11 )
        self.assertAlmostEqual( d1 - dm, 2.379650804185118e-7, 11 )

        ( r2, d2 ) = pal.subet( r1, d1, eq )
        self.assertAlmostEqual( r2 - rm, 0.0, 11 )
        self.assertAlmostEqual( d2 - dm, 0.0, 11 )

    def test_afin(self):
        (d, i) = pal.dafin( "12 34 56.7 |", 1 )
        self.assertEqual( i, 12 )
        self.assertAlmostEqual( d, 0.2196045986911432, 12 )

        (d, i) = pal.dafin( "45 00 00.000 ", 1 )
        self.assertEqual( i, 14 )
        self.assertAlmostEqual( d, pal.DPI / 4.0, 12 )

        (d, i) = pal.dafin( "30.2 < just decimal degrees", 1 )
        self.assertEqual( i, 6 )
        self.assertAlmostEqual( d, 30.2 * pal.DD2R, 12 )

        (d, i) = pal.dafin( "  30 23.6 < decimal armin", 1 )
        self.assertEqual( i, 11 )
        self.assertAlmostEqual( d, (30 + 23.6/60)*pal.DD2R, 12 )

        (d, i) = pal.dafin( " offset into string: 45.0 <<<", 22 )
        self.assertEqual( i, 27 )
        self.assertAlmostEqual( d, pal.DPI / 4.0, 12 )

        self.assertRaises( ValueError, pal.dafin, " not a sexagesimal string ", 1 )
        self.assertRaises( ValueError, pal.dafin, " 750.4 21 0.0 bad degrees ", 1 )
        self.assertRaises( ValueError, pal.dafin, " 45 30.5 0.0 bad arcminutes ", 1 )
        self.assertRaises( ValueError, pal.dafin, " 45 -30 0.0 bad arcminutes ", 1 )
        self.assertRaises( ValueError, pal.dafin, " 45 72 0.0 too many arcminutes ", 1 )
        self.assertRaises( ValueError, pal.dafin, " 45 43 85.0 too many arcseconds ", 1 )


    def test_airmass(self):
        self.assertAlmostEqual( pal.airmas( 1.2354 ),
                                3.015698990074724, 11 );

    def test_altaz(self):
        (az, azd, azdd, el, eld, eldd, pa, pad, padd) = pal.altaz( 0.7, -0.7, -0.65 )
        self.assertAlmostEqual( az, 4.400560746660174, 12 )
        self.assertAlmostEqual( azd, -0.2015438937145421, 12 )
        self.assertAlmostEqual( azdd, -0.4381266949668748, 13 )
        self.assertAlmostEqual( el, 1.026646506651396, 12 )
        self.assertAlmostEqual( eld, -0.7576920683826450, 13 )
        self.assertAlmostEqual( eldd, 0.04922465406857453, 14 )
        self.assertAlmostEqual( pa, 1.707639969653937, 12 )
        self.assertAlmostEqual( pad, 0.4717832355365627, 13 )
        self.assertAlmostEqual( padd, -0.2957914128185515, 13 )

    def test_altazVector(self):
        np.random.seed(32)
        phi = 0.5
        haIn = np.random.sample(20)*2.0*np.pi
        decIn = (np.random.sample(20)-0.5)*np.pi
        azC = np.zeros(20, dtype=np.float64)
        azdC = np.zeros(20, dtype=np.float64)
        azddC = np.zeros(20, dtype=np.float64)
        elC = np.zeros(20, dtype=np.float64)
        eldC = np.zeros(20, dtype=np.float64)
        elddC = np.zeros(20, dtype=np.float64)
        paC = np.zeros(20, dtype=np.float64)
        padC = np.zeros(20, dtype=np.float64)
        paddC = np.zeros(20, dtype=np.float64)
        for i in range(20):
            az, azd, azdd, el, eld, eldd, pa, pad, padd = pal.altaz(haIn[i], decIn[i], phi)
            azC[i] = az
            azdC[i] = azd
            azddC[i] = azdd
            elC[i] = el
            eldC[i] = eld
            elddC[i] = eldd
            paC[i] = pa
            padC[i] = pad
            paddC[i] = padd
        azT, azdT, azddT, elT, eldT, elddT, paT, padT, paddT = pal.altazVector(haIn, decIn, phi)

        for i in range(20):
            self.assertAlmostEqual(azC[i], azT[i], 12)
            self.assertAlmostEqual(azdC[i], azdT[i], 12)
            self.assertAlmostEqual(azddC[i], azddT[i], 12)
            self.assertAlmostEqual(elC[i], elT[i], 12)
            self.assertAlmostEqual(eldC[i], eldT[i], 12)
            self.assertAlmostEqual(elddC[i], elddT[i], 12)
            self.assertAlmostEqual(paC[i], paT[i], 12)
            self.assertAlmostEqual(padC[i], padT[i], 12)
            self.assertAlmostEqual(paddC[i], paddT[i], 12)

    def test_amp(self):
        (rm, dm) = pal.amp( 2.345, -1.234, 50100., 1990. )
        self.assertAlmostEqual( rm, 2.344472180027961, 6 )
        self.assertAlmostEqual( dm, -1.233573099847705, 7 )
        (rm, dm ) = pal.amp( 1.234, -0.567, 55927., 2010. )
        self.assertAlmostEqual( rm, 1.2335120411026936349, 12 )
        self.assertAlmostEqual( dm, -0.56702908706930343907, 12 )

    def test_ampqk(self):
        amprms = pal.mappa( 2010.0, 55927.0 )
        (rm, dm) = pal.ampqk( 1.234, -0.567, amprms )
        self.assertAlmostEqual( rm, 1.2335120411026936349, 11 )
        self.assertAlmostEqual( dm, -0.56702908706930343907, 11 )

    def test_aopqkVector(self):
        date = 51000.1
        dut = 25.0
        elongm = 2.1
        phim = 0.5
        hm = 3000.0
        xp = -0.5e-6
        yp = 1.0e-6
        tdk = 280.0
        pmb = 550.0
        rh = 0.6
        tlr = 0.006
        wl=0.45
        obsrms = pal.aoppa(date,dut,elongm,phim,hm,xp,yp,tdk,pmb,rh,wl,tlr)
        np.random.seed(32)
        nTests = 100
        raIn = np.random.sample(nTests)*2.0*np.pi
        decIn = (np.random.sample(nTests)-0.5)*np.pi
        azControl = None
        zeControl = None
        haControl = None
        dControl = None
        rControl = None
        for (rr,dd) in zip(raIn,decIn):
            az,ze,ha,d,r = pal.aopqk(rr,dd,obsrms)
            if azControl is None:
                azControl = np.array([az])
                zeControl = np.array([ze])
                haControl = np.array([ha])
                dControl = np.array([d])
                rControl = np.array([r])
            else:
                azControl = np.append(azControl,az)
                zeControl = np.append(zeControl,ze)
                haControl = np.append(haControl,ha)
                dControl = np.append(dControl,d)
                rControl = np.append(rControl,r)

        azTest,zeTest,haTest,dTest,rTest = pal.aopqkVector(raIn,decIn,obsrms)
        for (a1,z1,h1,d1,r1,a2,z2,h2,d2,r2) in \
            zip(azControl,zeControl,haControl,dControl,rControl,
                azTest,zeTest,haTest,dTest,rTest):

            self.assertAlmostEqual(a1,a2,12)
            self.assertAlmostEqual(z1,z2,12)
            self.assertAlmostEqual(h1,h2,12)
            self.assertAlmostEqual(d1,d2,12)
            self.assertAlmostEqual(r1,r2,12)

    def test_aop(self):
        dap = -0.1234
        date = 51000.1
        dut = 25.0
        elongm = 2.1
        phim = 0.5
        hm = 3000.0
        xp = -0.5e-6
        yp = 1.0e-6
        tdk = 280.0
        pmb = 550.0
        rh = 0.6
        tlr = 0.006

        aopres = [
            [
                1.812817787123283034,
                1.393860816635714034,
                -1.297808009092456683,
                -0.122967060534561,
                2.699270287872084
                ],
                [
                    2.019928026670621442,
                    1.101316172427482466,
                    -0.9432923558497740862,
                    -0.1232144708194224,
                    2.344754634629428
                    ],
        [
            2.019928026670621442,
            1.101267532198003760,
            -0.9432533138143315937,
            -0.1231850665614878,
            2.344715592593984
            ]
        ]
        aoptol = [ [ 10, 7, 7, 8, 7 ],
                   [ 10, 10, 10, 10, 10],
                   [ 10, 10, 10, 10, 10]
                   ]

        for i in range(len(aopres)):
            # Not very pythonic
            if i == 0:
                rap = 2.7
                wl = 0.45
            elif i==1:
                rap = 2.345
            else:
                wl = 1.0e6

            result = pal.aop( rap, dap, date, dut, elongm, phim, hm, xp, yp,
                              tdk, pmb, rh, wl, tlr )

            for j in range(len(result)):
                self.assertAlmostEqual( result[j], aopres[i][j], aoptol[i][j] )

        date = 48000.3
        wl = 0.45

        aoprms = pal.aoppa( date, dut, elongm, phim, hm, xp, yp, tdk, pmb,
                            rh, wl, tlr )

        aoppares = [ 0.4999993892136306,  0.4794250025886467,  0.8775828547167932,
                     1.363180872136126e-6, 3000., 280., 550., 0.6, 0.45, 0.006,
                     0.0001562803328459898, -1.792293660141e-7, 2.101874231495843,
                     7.601916802079765 ]
        aoppatol = [ 13, 13, 13, 13, 10, 11, 11, 13, 13, 15, 13, 13, 12, 8]
        self.assertEqual( len(aoprms), len(aoppares) )
        for i in range(len(aoprms)):
            self.assertAlmostEqual( aoprms[i], aoppares[i], aoppatol[i] )

        (rap, dap) = pal.oap( "r", 1.6, -1.01, date, dut, elongm, phim,
                              hm, xp, yp, tdk, pmb, rh, wl, tlr )
        self.assertAlmostEqual( rap, 1.601197569844787, 10 )
        self.assertAlmostEqual( dap, -1.012528566544262, 10 )

        (rap, dap) = pal.oap( "h", -1.234, 2.34, date, dut, elongm, phim,
                              hm, xp, yp, tdk, pmb, rh, wl, tlr )
        self.assertAlmostEqual( rap, 5.693087688154886463, 10 )
        self.assertAlmostEqual( dap, 0.8010281167405444, 10 )

        (rap, dap) = pal.oap( "a", 6.1, 1.1, date, dut, elongm, phim,
                              hm, xp, yp, tdk, pmb, rh, wl, tlr )
        self.assertAlmostEqual( rap, 5.894305175192448940, 10 )
        self.assertAlmostEqual( dap, 1.406150707974922, 10 )

        (rap, dap) = pal.oapqk( "r", 2.1, -0.345, aoprms )
        self.assertAlmostEqual( rap, 2.10023962776202, 10 )
        self.assertAlmostEqual( dap, -0.3452428692888919, 10 )

        (rap, dap) = pal.oapqk( "h", -0.01, 1.03, aoprms )
        self.assertAlmostEqual( rap, 1.328731933634564995, 10 )
        self.assertAlmostEqual( dap, 1.030091538647746, 10 )

        (rap, dap) = pal.oapqk( "a", 4.321, 0.987, aoprms )
        self.assertAlmostEqual( rap, 0.4375507112075065923, 10 )
        self.assertAlmostEqual( dap, -0.01520898480744436, 10 )

        aoprms = pal.aoppat( date + pal.DS2R, aoprms )
        self.assertAlmostEqual( aoprms[13], 7.602374979243502, 8 )

    def test_bear(self):
        a1 = 1.234
        b1 = -0.123
        a2 = 2.345
        b2 = 0.789
        self.assertAlmostEqual( pal.dbear( a1, b1, a2, b2),
                                0.7045970341781791, 12 )
        d1 = pal.dcs2c( a1, b1 )
        d2 = pal.dcs2c( a2, b2 )
        self.assertAlmostEqual( pal.dpav( d1, d2 ), 0.7045970341781791, 12 )

    def test_caldj(self):
        djm = pal.caldj( 1999, 12, 31 )
        self.assertEqual( djm, 51543 )

        self.assertRaises( ValueError, pal.caldj, -5000, 12, 31 )
        self.assertRaises( ValueError, pal.caldj, 1970, 13, 1 )
        self.assertRaises( ValueError, pal.caldj, 1970, 1, 32 )

    def test_caf2r(self):
        dr = pal.daf2r( 76, 54, 32.1 )
        self.assertAlmostEqual( dr, 1.342313819975276, 12 )

    def test_cc2s(self):
        (da, db) = pal.dcc2s( np.array( [100., -50., 25. ] ) )
        self.assertAlmostEqual( da, -0.4636476090008061, 12 )
        self.assertAlmostEqual( db, 0.2199879773954594, 12 )

    def test_cd2tf(self):
        ( sign, hours, minutes, seconds, fraction ) = pal.dd2tf( 4, -0.987654321 )
        self.assertEqual( sign, "-" )
        self.assertEqual( hours, 23 )
        self.assertEqual( minutes, 42 )
        self.assertEqual( seconds, 13 )
        self.assertEqual( fraction, 3333 )

    def test_cldj(self):
        d = pal.cldj( 1899, 12, 31 )
        self.assertEqual( d, 15019 )
        self.assertRaises( ValueError, pal.cldj, -5000, 12, 31 )
        self.assertRaises( ValueError, pal.cldj, 1970, 13, 1 )
        self.assertRaises( ValueError, pal.cldj, 1970, 1, 32 )

    def test_cr2af(self):
        (sign, deg, min, sec, f) = pal.dr2af( 4, 2.345 )
        self.assertEqual( sign, "+" )
        self.assertEqual( deg, 134 )
        self.assertEqual( min, 21 )
        self.assertEqual( sec, 30 )
        self.assertEqual( f, 9706 )

    def test_cr2tf(self):
        (sign, hr, min, sec, f) = pal.dr2tf( 4, -3.01234 )
        self.assertEqual( sign, "-" )
        self.assertEqual( hr, 11 )
        self.assertEqual( min, 30 )
        self.assertEqual( sec, 22 )
        self.assertEqual( f, 6484 )

    def test_ctf2d(self):
        dd = pal.dtf2d( 23, 56, 59.1 )
        self.assertAlmostEqual( dd, 0.99790625, 12 )
        self.assertRaises( ValueError, pal.dtf2d, 24, 1, 32 )
        self.assertRaises( ValueError, pal.dtf2d, 23, -1, 32 )
        self.assertRaises( ValueError, pal.dtf2d, 23, 1, 60 )

    def test_ctf2r(self):
        dr = pal.dtf2r( 23, 56, 59.1 )
        self.assertAlmostEqual( dr, 6.270029887942679, 12 )
        self.assertRaises( ValueError, pal.dtf2r, 24, 1, 32 )
        self.assertRaises( ValueError, pal.dtf2r, 23, -1, 32 )
        self.assertRaises( ValueError, pal.dtf2r, 23, 1, 60 )

    def test_dat(self):
        self.assertEqual( pal.dat( 43900 ), 18 )
        self.assertAlmostEqual( pal.dtt( 40404 ), 39.709746, 12 )
        self.assertAlmostEqual( pal.dt( 500 ), 4686.7, 10 )
        self.assertAlmostEqual( pal.dt( 1400 ), 408, 11 )
        self.assertAlmostEqual( pal.dt( 1950 ), 27.99145626 )

    def test_djcal(self):
        djm = 50123.9999
        ( y, m, d, f ) = pal.djcal( 4, djm )
        self.assertEqual( y, 1996 )
        self.assertEqual( m, 2 )
        self.assertEqual( d, 10 )
        self.assertEqual( f, 9999 )

        ( y, m, d, f ) = pal.djcl( djm )
        self.assertEqual( y, 1996 )
        self.assertEqual( m, 2 )
        self.assertEqual( d, 10 )
        self.assertAlmostEqual( f, 0.9999, 7 )

    def test_dmat(self):
        da = np.array(
           [ [ 2.22,     1.6578,     1.380522 ],
            [  1.6578,   1.380522,   1.22548578 ],
            [  1.380522, 1.22548578, 1.1356276122 ] ]
        )
        dv = np.array( [ 2.28625, 1.7128825, 1.429432225 ] )
        (da, dv, dd) = pal.dmat( da, dv )
        self.assertAlmostEqual( dd, 0.003658344147359863, 12 )
        self.assertAlmostEqual( dv[0], 1.002346480763383, 12 )
        self.assertAlmostEqual( dv[1], 0.03285594016974583489, 12 )
        self.assertAlmostEqual( dv[2], 0.004760688414885247309, 12 )

        da = np.array( [ [ 0., 1. ], [ 0., 1. ] ] )
        dv = np.array( [ 1., 1. ] )
        self.assertRaises( ArithmeticError, pal.dmat, da, dv )

    def test_e2h(self):
        dp = -0.7
        hin = -0.3
        din = -1.1
        (da, de) = pal.de2h( hin, din, dp )
        self.assertAlmostEqual( da, 2.820087515852369, 12 )
        self.assertAlmostEqual( de, 1.132711866443304, 12 )

        (dh, dd) = pal.dh2e( da, de, dp )
        self.assertAlmostEqual( dh, hin, 12 )
        self.assertAlmostEqual( dd, din, 12 )

    def test_de2hVector(self):
        nTests = 100
        phi = 0.35
        np.random.seed(32)
        raIn = np.random.random_sample(nTests)*np.pi*2.0
        decIn = (np.random.random_sample(nTests)-0.5)*np.pi
        azControl = None
        elControl = None
        for (rr,dd) in zip (raIn,decIn):
            az, el = pal.de2h(rr, dd, phi)
            if azControl is None:
                azControl = np.array([az])
                elControl = np.array([el])
            else:
                azControl = np.append(azControl,az)
                elControl = np.append(elControl,el)

        azTest, elTest = pal.de2hVector(raIn, decIn, phi)
        for (a1,e1,a2,e2) in zip (azControl,elControl,azTest,elTest):
            self.assertAlmostEqual(a1,a2,12)
            self.assertAlmostEqual(e1,e2,12)

    def test_ecleq(self):
        (dr,dd) = pal.ecleq( 1.234, -0.123, 43210.0)
        self.assertAlmostEqual( dr, 1.229910118208851, 5 )
        self.assertAlmostEqual( dd, 0.2638461400411088, 5 )

    def test_ecmat(self):
        expected = np.array( [
            [ 1.0,                    0.0,                   0.0 ],
            [ 0.0, 0.91749307789883549624, 0.3977517467060596168 ],
            [ 0.0, -0.3977517467060596168, 0.91749307789883549624 ] ] );

        rmat = pal.ecmat( 55966.46 )
        np.testing.assert_array_almost_equal( rmat, expected, decimal=12 )

    def test_epb(self):
        self.assertAlmostEqual( pal.epb( 45123 ), 1982.419793168669, 8)

    def test_epb2d(self):
        self.assertAlmostEqual( pal.epb2d(1975.5), 42595.5995279655, 7)

    def test_epco(self):
        self.assertAlmostEqual( pal.epco("B","J", 2000), 2000.001277513665, 7 )
        self.assertAlmostEqual( pal.epco("J","B", 1950), 1949.999790442300, 7 )
        self.assertAlmostEqual( pal.epco("J","j", 2000), 2000, 7 )

    def test_epj(self):
        self.assertAlmostEqual( pal.epj(42999), 1976.603696098563, 7)

    def test_epj2d(self):
        self.assertAlmostEqual( pal.epj2d(2010.077), 55225.124250, 6 )

    def test_eqecl(self):
        (dl, db) = pal.eqecl( 0.789, -0.123, 46555 )
        self.assertAlmostEqual( dl, 0.7036566430349022, 6 )
        self.assertAlmostEqual( db, -0.4036047164116848, 6 )

    def test_eqeqx(self):
        self.assertAlmostEqual( pal.eqeqx( 53736 ), -0.8834195072043790156e-5, 15 )

    def test_eqeqxVector(self):
        np.random.seed(32)
        dateIn = 53000.0 + np.random.sample(20)*5000.0
        eqControl = np.zeros(20, dtype=np.float64)
        for i, d in enumerate(dateIn):
            eqControl[i] = pal.eqeqx(d)
        eqTest = pal.eqeqxVector(dateIn)
        for et, ec in zip(eqTest, eqControl):
            self.assertAlmostEqual(et, ec, 15)

    def test_eqgal(self):
        (dl, db) = pal.eqgal( 5.67, -1.23 )
        self.assertAlmostEqual( dl, 5.612270780904526, 12 )
        self.assertAlmostEqual( db, -0.6800521449061520, 12 )

    def test_eqgalVector(self):
        np.random.seed(32)
        raIn = np.random.sample(10)*2.0*np.pi
        decIn = (np.random.sample(10)-0.5)*np.pi
        dlControl = None
        dbControl = None
        for (ra, dec) in zip(raIn, decIn):
            dl, db = pal.eqgal(ra, dec)
            if dlControl is None:
                dlControl = np.array([dl])
                dbControl = np.array([db])
            else:
                np.append(dlControl, dl)
                np.append(dbControl, db)
        dlTest, dbTest = pal.eqgalVector(raIn, decIn)
        for (dlt, dbt, dlc, dbc) in zip(dlTest, dbTest, dlControl, dbControl):
            self.assertAlmostEqual(dlt, dlc, 12)
            self.assertAlmostEqual(dbt, dbc, 12)

    def test_etrms(self):
        ev = pal.etrms( 1976.9 )
        self.assertAlmostEqual( ev[0], -1.621617102537041e-6, 18 )
        self.assertAlmostEqual( ev[1], -3.310070088507914e-7, 18 )
        self.assertAlmostEqual( ev[2], -1.435296627515719e-7, 18 )

    def test_evp(self):
        vbex = np.array([ 1.6957348127008098514e-07,
                         -9.1093446116039685966e-08,
                         -3.9528532243991863036e-08 ])
        pbex = np.array([-0.49771075259730546136,
                         -0.80273812396332311359,
                         -0.34851593942866060383])
        vhex = np.array([ 1.6964379181455713805e-07,
                         -9.1147224045727438391e-08,
                         -3.9553158272334222497e-08])
        phex = np.array([-0.50169124421419830639,
                         -0.80650980174901798492,
                          -0.34997162028527262212])

        vbex2 = np.array([-0.0109187426811683,
                          -0.0124652546173285,
                          -0.0054047731809662])
        pbex2 = np.array([-0.7714104440491060,
                          +0.5598412061824225,
                          +0.2425996277722475])
        vhex2 = np.array([-0.0109189182414732,
                          -0.0124718726844084,
                          -0.0054075694180650])
        phex2 = np.array([-0.7757238809297653,
                          +0.5598052241363390,
                          +0.2426998466481708])

        (dvb, dpb, dvh, dph) = pal.evp( 2010.0, 2012.0 )
        np.testing.assert_allclose( dvb, vbex, atol=1e-12 )
        np.testing.assert_allclose( dpb, pbex, atol=1e-12 )
        np.testing.assert_allclose( dvh, vhex, atol=1e-12 )
        np.testing.assert_allclose( dph, phex, atol=1e-12 )

        (dph, dvh, dpb, dvb) = pal.epv( 53411.52501161 )
        np.testing.assert_allclose( dvb, vbex2, atol=1e-12 )
        np.testing.assert_allclose( dpb, pbex2, atol=1e-12 )
        np.testing.assert_allclose( dvh, vhex2, atol=1e-12 )
        np.testing.assert_allclose( dph, phex2, atol=1e-12 )

    def test_fk45z(self):
        (r2000, d2000) = pal.fk45z( 1.2, -0.3, 1960 )
        self.assertAlmostEqual( r2000, 1.2097812228966762227, 11 )
        self.assertAlmostEqual( d2000, -0.29826111711331398935, 12 )

    def test_fk52h(self):
        inr5 = 1.234
        ind5 = -0.987
        epoch = 1980
        (rh, dh) = pal.fk5hz( inr5, ind5, epoch )
        self.assertAlmostEqual( rh, 1.234000136713611301, 13 )
        self.assertAlmostEqual( dh, -0.9869999702020807601, 13 )
        (r5, d5, dr5, dd5) = pal.hfk5z( rh, dh, epoch )
        self.assertAlmostEqual( r5, inr5, 13 )
        self.assertAlmostEqual( d5, ind5, 13 )
        self.assertAlmostEqual( dr5, 0.000000006822074, 13 )
        self.assertAlmostEqual( dd5, -0.000000002334012, 13 )

    def test_fk54z(self):
        (r1950, d1950, dr1950, dd1950) = pal.fk54z( 1.2, -0.3, 1960 )
        self.assertAlmostEqual( r1950, 1.1902221805755279771, 12 )
        self.assertAlmostEqual( d1950, -0.30178317645793828472, 12 )
        self.assertAlmostEqual( dr1950, -1.7830874775952945507e-08, 12 )
        self.assertAlmostEqual( dd1950, 7.196059425334821089e-09, 12 )

    def test_flotin(self):
        pass

    def test_galeq(self):
        (dr, dd) = pal.galeq( 5.67, -1.23 )
        self.assertAlmostEqual( dr, 0.04729270418071426, 12 )
        self.assertAlmostEqual( dd, -0.7834003666745548, 12 )

    def test_galeqVector(self):
        np.random.seed(32)
        dlIn = np.random.sample(10)*2.0*np.pi
        dbIn = (np.random.sample(10)-0.5)*np.pi
        drControl = np.zeros(10, dtype=np.float64)
        ddControl = np.zeros(10, dtype=np.float64)
        for (i, cc) in enumerate(zip(dlIn, dbIn)):
            dr, dd = pal.galeq(cc[0], cc[1])
            drControl[i] = dr
            ddControl[i] = dd
        drTest, ddTest = pal.galeqVector(dlIn, dbIn)
        for (drt, ddt, drc, ddc) in zip(drTest, ddTest, drControl, ddControl):
            self.assertAlmostEqual(drt, drc, 12)
            self.assertAlmostEqual(ddt, ddc, 12)

    def test_galsup(self):
        (dsl, dsb) = pal.galsup( 6.1, -1.4 )
        self.assertAlmostEqual( dsl, 4.567933268859171, 12 )
        self.assertAlmostEqual( dsb, -0.01862369899731829, 12 )

    def test_geoc(self):
        lat = 19.822838905884 * pal.DD2R
        alt = 4120.0
        (r, z) = pal.geoc( lat, alt )
        self.assertAlmostEqual( r, 4.01502667039618e-05, 10 )
        self.assertAlmostEqual( z, 1.43762411970295e-05, 10 )

    def test_ge50(self):
        (dr, dd) = pal.ge50( 6.1, -1.55 )
        self.assertAlmostEqual( dr, 0.1966825219934508, 12 )
        self.assertAlmostEqual( dd, -0.4924752701678960, 12 )

    def test_gmst(self):
        self.assertAlmostEqual( pal.gmst( 53736 ), 1.754174971870091203, 12 )
        self.assertAlmostEqual( pal.gmsta( 53736, 0 ), 1.754174971870091203, 12 )

    def test_gmstVector(self):
        np.random.seed(32)
        dateIn = 53000.0 + np.random.sample(20)*5000.0
        gmControl = np.zeros(20, dtype=np.float64)
        for i, d in enumerate(dateIn):
             gm = pal.gmst(d)
             gmControl[i] = gm
        gmTest = pal.gmstVector(dateIn)
        for gt, gc in zip(gmTest, gmControl):
            self.assertAlmostEqual(gt, gc, 12)

    def test_gmstaVector(self):
        np.random.seed(32)
        dateIn = 53000 + np.random.random_integers(0, 1000, 20)
        utIn = np.random.sample(20)
        gmControl = np.zeros(20, dtype=np.float64)
        for i, d in enumerate(zip(dateIn, utIn)):
            dd = pal.gmsta(d[0], d[1])
            gmControl[i] = dd
        gmTest = pal.gmstaVector(dateIn, utIn)
        for gt, gc in zip(gmTest, gmControl):
            self.assertAlmostEqual(gt, gc, 12)

    def test_intin(self):
        s = "  -12345, , -0  2000  +     "
        i = 1
        (n, i, sign) = pal.intin( s, i )
        self.assertEqual( i, 10 )
        self.assertEqual( n, -12345 )
        self.assertLess( sign, 0 )

        (n, i, sign) = pal.intin( s, i )
        self.assertEqual( i, 12 )
        self.assertIsNone( n )

        (n, i, sign) = pal.intin( s, i )
        self.assertEqual( i, 17 )
        self.assertEqual( n, 0 )
        self.assertLess( sign, 0 )

        (n, i, sign) = pal.intin( s, i )
        self.assertEqual( i, 23 )
        self.assertEqual( n, 2000 )
        self.assertGreater( sign, 0 )

        (n, i, sign) = pal.intin( s, i )
        self.assertEqual( i, 29 )
        self.assertIsNone( n )
        self.assertIsNone( sign )

    def test_map(self):
        (ra, da) = pal.map( 6.123, -0.999, 1.23e-5, -0.987e-5,
                            0.123, 32.1, 1999, 43210.9 )
        self.assertAlmostEqual( ra, 6.117130429775647, 6 )
        self.assertAlmostEqual( da, -1.000880769038632, 7 )

    def test_mappa(self):
        expected = np.array( [1.9986310746064646082,
                          -0.1728200754134739392,
                          0.88745394651412767839,
                          0.38472374350184274094,
                          -0.17245634725219796679,
                          0.90374808622520386159,
                          0.3917884696321610738,
                          2.0075929387510784968e-08,
                          -9.9464149073251757597e-05,
                          -1.6125306981057062306e-05,
                          -6.9897255793245634435e-06,
                          0.99999999489900059935,
                          0.99999983777998024959,
                          -0.00052248206600935195865,
                          -0.00022683144398381763045,
                          0.00052248547063364874764,
                          0.99999986339269864022,
                          1.4950491424992534218e-05,
                          0.00022682360163333854623,
                          -1.5069005133483779417e-05,
                          0.99999997416198904698 ] )
        amprms = pal.mappa( 2010.0, 55927 )
        np.testing.assert_array_almost_equal( amprms, expected, decimal=12 )

    def test_mapqkz(self):
        amprms = pal.mappa( 2010, 55927 )
        (ra, da) = pal.mapqkz( 1.234, -0.567, amprms )
        self.assertAlmostEqual( ra, 1.2344879748414849807, 12 )
        self.assertAlmostEqual( da, -0.56697099554368701746, 12 )

        (ra, da) = pal.mapqk( 1.234, -0.567, 0., 0., 0., 0., amprms )
        self.assertAlmostEqual( ra, 1.2344879748414849807, 7 )
        self.assertAlmostEqual( da, -0.56697099554368701746, 7 )

    def test_mapqkzVector(self):
        amprms = pal.mappa(2010, 55927)
        np.random.seed(32)
        nTests = 100
        raIn = np.random.sample(nTests)*2.0*np.pi
        decIn = (np.random.sample(nTests)-0.5)*np.pi
        rControl=None
        dControl=None
        for (rr,dd) in zip(raIn,decIn):
            r,d = pal.mapqkz(rr, dd, amprms)
            if rControl is None:
                rControl = np.array([r])
                dControl = np.array([d])
            else:
                rControl = np.append(rControl,r)
                dControl = np.append(dControl,d)

        rTest,dTest = pal.mapqkzVector(raIn,decIn,amprms)
        for (r1,d1,r2,d2) in zip(rControl,dControl,rTest,dTest):
            self.assertAlmostEqual(r1,r2,12)
            self.assertAlmostEqual(d1,d2,12)

    def test_mapqkVector(self):
        amprms = pal.mappa(2010, 55927)
        np.random.seed(32)
        nTests = 100
        raIn = np.random.sample(nTests)*2.0*np.pi
        decIn = (np.random.sample(nTests)-0.5)*np.pi
        pmr = (np.random.sample(nTests)-0.5)*0.01
        pmd = (np.random.sample(nTests)-0.5)*0.01
        px = 0.00045+np.random.sample(nTests)*0.001
        rv = 200.0*np.random.sample(nTests)
        rControl=None
        dControl=None
        for (rr,dd,pr,pd,x,v) in zip(raIn,decIn,pmr,pmd,px,rv):
            r,d = pal.mapqk(rr, dd, pr, pd, x, v, amprms)
            if rControl is None:
                rControl = np.array([r])
                dControl = np.array([d])
            else:
                rControl = np.append(rControl,r)
                dControl = np.append(dControl,d)

        rTest,dTest = pal.mapqkVector(raIn,decIn,pmr,pmd,px,rv,amprms)
        for (r1,d1,r2,d2) in zip(rControl,dControl,rTest,dTest):
            self.assertAlmostEqual(r1,r2,12)
            self.assertAlmostEqual(d1,d2,12)

    def test_moon(self):
        expected = np.array( [
             0.00229161514616454,
             0.000973912029208393,
             0.000669931538978146,
             -3.44709700068209e-09,
             5.44477533462392e-09,
             2.11785724844417e-09
             ] )
        pv = pal.dmoon( 48634.4687174074 )
        np.testing.assert_array_almost_equal( pv, expected, decimal=12 )

    def test_nut(self):
        expected = np.array( [
            [  9.999999969492166e-1, 7.166577986249302e-5,  3.107382973077677e-5 ],
            [ -7.166503970900504e-5, 9.999999971483732e-1, -2.381965032461830e-5 ],
            [ -3.107553669598237e-5, 2.381742334472628e-5,  9.999999992335206818e-1 ]
        ] )

        rmatn = pal.nut( 46012.32 )
        np.testing.assert_array_almost_equal( rmatn, expected, decimal=3 )

        (dpsi, deps, eps0) = pal.nutc( 54388.0 )
        self.assertAlmostEqual( eps0, 0.4090749229387258204, 14 )

        (dpsi, deps, eps0) = pal.nutc( 53736.0 )
        self.assertAlmostEqual( dpsi, -0.9630912025820308797e-5, 13 )
        self.assertAlmostEqual( deps, 0.4063238496887249798e-4, 13 )

    def test_obs(self):
        obsdata = pal.obs()
        mmt = obsdata["MMT"]
        self.assertEqual( mmt["name"], "MMT 6.5m, Mt Hopkins" )
        self.assertAlmostEqual( mmt["long"], 1.935300584055477, 8 )
        self.assertAlmostEqual( mmt["lat"], 0.5530735081550342238, 10 )
        self.assertAlmostEqual( mmt["height"], 2608, 10 )

        self.assertEqual( len(obsdata), 85 )

    def test_pa(self):
        self.assertAlmostEqual( pal.pa( -1.567, 1.5123, 0.987 ),
                                -1.486288540423851, 12 )
        self.assertAlmostEqual( pal.pa( 0, 0.789, 0.789 ), 0, 12 )

    def test_paVector(self):
        np.random.seed(32)
        haIn = np.random.sample(20)*2.0*np.pi
        decIn = (np.random.sample(20)-0.5)*np.pi
        phi = 0.3
        paControl = np.zeros(20, dtype=np.float64)
        for i in range(20):
            paControl[i] = pal.pa(haIn[i], decIn[i], phi)
        paTest = pal.paVector(haIn, decIn, phi)
        for pt, pc in zip(paTest, paControl):
            self.assertAlmostEqual(pt, pc, 12)

    def test_pcd(self):
        disco = 178.585
        REFX = 0.0123
        REFY = -0.00987
        (x, y) = pal.pcd( disco, REFX, REFY )
        self.assertAlmostEqual( x, 0.01284630845735895, 14 )
        self.assertAlmostEqual( y, -0.01030837922553926, 14 )

        (x, y) = pal.unpcd( disco, x, y )
        self.assertAlmostEqual( x, REFX, 14 )
        self.assertAlmostEqual( y, REFY, 14 )

        # Round trip
        (x, y) = pal.pcd( -disco, REFX, REFY )
        (x, y) = pal.unpcd( -disco, x, y )
        self.assertAlmostEqual( x, REFX, 14 )
        self.assertAlmostEqual( y, REFY, 14 )

    def test_planet(self):
        # palEl2ue
        u = pal.el2ue( 50000, 1, 49000, 0.1, 2, 0.2,
                       3, 0.05, 3, 0.003312 )
        expectedue1 = np.array([1.000878908362435284,  -0.3336263027874777288,  50000.,
                2.840425801310305210,   0.1264380368035014224, -0.2287711835229143197,
                -0.01301062595106185195, 0.5657102158104651697,  0.2189745287281794885,
                2.852427310959998500,  -0.01552349065435120900,
                50000., 0.0])
        np.testing.assert_allclose( u, expectedue1, atol=1e-12 )

        # palPertel
        (epoch, orbinc, anode, perih, aorq, e, aorl) = pal.pertel(2, 43000., 43200., 43000.,
                                                                  0.2, 3, 4, 5, 0.02, 6)
        self.assertAlmostEqual( epoch, 43200, 10 )
        self.assertAlmostEqual( orbinc, 0.1995661466545422381, 7 )
        self.assertAlmostEqual( anode, 2.998052737821591215, 7 )
        self.assertAlmostEqual( perih, 4.009516448441143636, 6 )
        self.assertAlmostEqual( aorq, 5.014216294790922323, 7 )
        self.assertAlmostEqual( e, 0.02281386258309823607, 7 )
        self.assertAlmostEqual( aorl, 0.01735248648779583748, 6 )

        # palPertue
        unew = pal.pertue( 50100, u )
        expectedue3 = np.array([
            1.000000000000000,
            -0.3329769417028020949,
            50100.,
            2.638884303608524597,
            1.070994304747824305,
            0.1544112080167568589,
            -0.2188240619161439344,
            0.5207557453451906385,
            0.2217782439275216936,
            2.852118859689216658,
            0.01452010174371893229,
            50100., 0. ])
        np.testing.assert_allclose( unew, expectedue3, atol=1e-12 )

        # palPlanel
        pv = pal.planel( 50600, 2, 50500, 0.1, 3, 5, 2, 0.3, 4, 0 )
        expectedpv2 = np.array ([1.947628959288897677,
            -1.013736058752235271,
            -0.3536409947732733647,
            2.742247411571786194e-8,
            1.170467244079075911e-7,
            3.709878268217564005e-8])
        np.testing.assert_allclose( pv, expectedpv2, atol=1e-12 )

        # palPlanet
        self.assertRaises( ValueError, pal.planet, 1e6, 0 )
        self.assertRaises( ValueError, pal.planet, 1e6, 9 )

        pv = pal.planet( -320000, 3 )
        self.assertAlmostEqual( pv[0], 0.9308038666827242603, 11)
        self.assertAlmostEqual( pv[1], 0.3258319040252137618, 11)
        self.assertAlmostEqual( pv[2], 0.1422794544477122021, 11)
        self.assertAlmostEqual( pv[3], -7.441503423889371696e-8, 17)
        self.assertAlmostEqual( pv[4], 1.699734557528650689e-7, 17)
        self.assertAlmostEqual( pv[5], 7.415505123001430864e-8, 17)

        pv = pal.planet( 43999.9, 1 )
        self.assertAlmostEqual( pv[0], 0.2945293959257422246, 11)
        self.assertAlmostEqual( pv[1], -0.2452204176601052181, 11)
        self.assertAlmostEqual( pv[2], -0.1615427700571978643, 11)
        self.assertAlmostEqual( pv[3], 1.636421147459047057e-7, 18)
        self.assertAlmostEqual( pv[4], 2.252949422574889753e-7, 18)
        self.assertAlmostEqual( pv[5], 1.033542799062371839e-7, 18)

        # palPlante
        (ra, dec, r) = pal.plante( 50600., -1.23, 0.456, 2, 50500.,
                                   0.1, 3., 5., 2., 0.3, 4., 0.0 )
        self.assertAlmostEqual( ra, 6.222958101333794007, 6 )
        self.assertAlmostEqual( dec, 0.01142220305739771601, 6 )
        self.assertAlmostEqual( r, 2.288902494080167624, 8 )

        # palPlantu
        u = np.array( [1.0005, -0.3, 55000., 2.8, 0.1, -0.2,
                       -0.01, 0.5, 0.22, 2.8, -0.015, 55001., 0.0] )
        (ra, dec, r) = pal.plantu( 55001., -1.23, 0.456, u )
        self.assertAlmostEqual( ra, 0.3531814831241686647, 6 )
        self.assertAlmostEqual( dec, 0.06940344580567131328, 6 )
        self.assertAlmostEqual( r, 3.031687170873274464, 8 )

        # palPv2el
        pv = np.array( [ 0.3, -0.2, 0.1, -0.9e-7, 0.8e-7, -0.7e-7 ] )
        (jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm ) = pal.pv2el( pv, 50000, 0.00006, 1 )
        self.assertEqual( jform, 1 )
        self.assertAlmostEqual( epoch, 50000, 10 )
        self.assertAlmostEqual( orbinc, 1.52099895268912, 12 )
        self.assertAlmostEqual( anode, 2.720503180538650, 12 )
        self.assertAlmostEqual( perih, 2.194081512031836, 12 )
        self.assertAlmostEqual( aorq, 0.2059371035373771, 12 )
        self.assertAlmostEqual( e, 0.9866822985810528, 12 )
        self.assertAlmostEqual( aorl, 0.2012758344836794, 12 )
        self.assertAlmostEqual( dm, 0.184074050795182, 12 )

        # palPv2ue
        expectedue2 = np.array( [1.00006, -4.856142884511782, 50000., 0.3, -0.2,
                                0.1,  -0.4520378601821727,  0.4018114312730424,
                                -.3515850023639121, 0.3741657386773941,
                                -0.2511321445456515, 50000., 0.] )
        u = pal.pv2ue( pv, 50000., 0.00006 )
        np.testing.assert_allclose( u, expectedue2, atol=1e-12 )

        # Planets
        (ra, dec, diam) = pal.rdplan( 40999.9, 0, 0.1, -0.9 )
        self.assertAlmostEqual( ra, 5.772270359389275837, 6 )
        self.assertAlmostEqual( dec, -0.2089207338795416192, 7 )
        self.assertAlmostEqual( diam, 9.415338935229717875e-3, 10 )

        (ra, dec, diam) = pal.rdplan( 41999.9, 1, 1.1, -0.9 )
        self.assertAlmostEqual( ra, 3.866363420052936653, 6 )
        self.assertAlmostEqual( dec, -0.2594430577550113130, 7 )
        self.assertAlmostEqual( diam, 4.638468996795023071e-5, 14 )

        (ra, dec, diam) = pal.rdplan( 42999.9, 2, 2.1, 0.9 )
        self.assertAlmostEqual( ra, 2.695383203184077378, 6 )
        self.assertAlmostEqual( dec, 0.2124044506294805126, 7 )
        self.assertAlmostEqual( diam, 4.892222838681000389e-5, 14 )

        (ra, dec, diam) = pal.rdplan( 43999.9, 3, 3.1, 0.9 )
        self.assertAlmostEqual( ra, 2.908326678461540165, 6 )
        self.assertAlmostEqual( dec, 0.08729783126905579385, 7 )
        self.assertAlmostEqual( diam, 8.581305866034962476e-3, 7 )

        (ra, dec, diam) = pal.rdplan( 44999.9, 4, -0.1, 1.1 )
        self.assertAlmostEqual( ra, 3.429840787472851721, 6 )
        self.assertAlmostEqual( dec, -0.06979851055261161013, 7 )
        self.assertAlmostEqual( diam, 4.540536678439300199e-5, 14 )

        (ra, dec, diam) = pal.rdplan( 45999.9, 5, -1.1, 0.1 )
        self.assertAlmostEqual( ra, 4.864669466449422548, 6 )
        self.assertAlmostEqual( dec, -0.4077714497908953354, 7 )
        self.assertAlmostEqual( diam, 1.727945579027815576e-4, 14 )

        (ra, dec, diam) = pal.rdplan( 46999.9, 6, -2.1, -0.1 )
        self.assertAlmostEqual( ra, 4.432929829176388766, 6 )
        self.assertAlmostEqual( dec, -0.3682820877854730530, 7 )
        self.assertAlmostEqual( diam, 8.670829016099083311e-5, 14 )

        (ra, dec, diam) = pal.rdplan( 47999.9, 7, -3.1, -1.1 )
        self.assertAlmostEqual( ra, 4.894972492286818487, 6 )
        self.assertAlmostEqual( dec, -0.4084068901053653125, 7 )
        self.assertAlmostEqual( diam, 1.793916783975974163e-5, 14 )

        (ra, dec, diam) = pal.rdplan( 48999.9, 8, 0, 0 )
        self.assertAlmostEqual( ra, 5.066050284760144000, 6 )
        self.assertAlmostEqual( dec, -0.3744690779683850609, 7 )
        self.assertAlmostEqual( diam, 1.062210086082700563e-5, 14 )

        # palUe2el
        (jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm ) = pal.ue2el( u, 1 )
        self.assertEqual( jform, 1 )
        self.assertAlmostEqual( epoch, 50000.00, 10 )
        self.assertAlmostEqual( orbinc, 1.520998952689120, 12 )
        self.assertAlmostEqual( anode, 2.720503180538650, 12 )
        self.assertAlmostEqual( perih, 2.194081512031836, 12 )
        self.assertAlmostEqual( aorq, 0.2059371035373771, 12 )
        self.assertAlmostEqual( e, 0.9866822985810528, 12 )
        self.assertAlmostEqual( aorl, 0.2012758344836794, 12 )

        # palUe2pv
        (u2, pv) = pal.ue2pv( 50010., u )

        # Update final two elements of the test UE array
        expectedue2[11] = 50010.;
        expectedue2[12] = 0.7194308220038886856;
        np.testing.assert_allclose( u2, expectedue2, atol=1e-12 )

        expectedpv = np.array( [
            0.07944764084631667011, -0.04118141077419014775,
            0.002915180702063625400, -0.6890132370721108608e-6,
            0.4326690733487621457e-6, -0.1763249096254134306e-6 ] )
        np.testing.assert_allclose( pv, expectedpv, atol=1e-12 )

    def test_pm(self):
        (ra, dec) = pal.pm( 5.43, -0.87, -0.33e-5, 0.77e-5, 0.7,
                              50.3*365.2422/365.25, 1899, 1943 )
        self.assertAlmostEqual( ra, 5.429855087793875, 10 )
        self.assertAlmostEqual( dec, -0.8696617307805072, 10 )

        (ra, dec) = pal.pm( 0.01686756, -1.093989828, -1.78323516e-5,
                            2.336024047e-6, 0.74723, -21.6,
                            pal.epj( 50083.0 ), pal.epj( 53736.0 ) )
        self.assertAlmostEqual( ra, 0.01668919069414242368, 13 )
        self.assertAlmostEqual( dec, -1.093966454217127879, 13 )

    def test_pmVector(self):
        ep0 = 52000.0
        ep1 = 53510.0
        np.random.seed(32)
        nTests = 100
        raIn = np.random.sample(nTests)*2.0*np.pi
        decIn = (np.random.sample(nTests)-0.5)*np.pi
        pmr = 0.01*np.random.sample(nTests)
        pmd = 0.01*np.random.sample(nTests)
        px = 0.00045 + 0.001*np.random.sample(nTests)
        rv = 1000.0*np.random.sample(nTests)
        rControl = None
        dControl = None
        for (rr,dd,pr,pd,x,v) in zip(raIn,decIn,pmr,pmd,px,rv):
            r, d = pal.pm(rr,dd,pr,pd,x,v,ep0,ep1)
            if rControl is None:
                rControl = np.array([r])
                dControl = np.array([d])
            else:
                rControl = np.append(rControl,r)
                dControl = np.append(dControl,d)

        rTest, dTest = pal.pmVector(raIn,decIn,pmr,pmd,px,rv,ep0,ep1)
        for (r1,d1,r2,d2) in zip(rControl,dControl,rTest,dTest):
            self.assertAlmostEqual(r1,r2,12)
            self.assertAlmostEqual(d1,d2,12)

    def test_polmo(self):
        (elong, phi, daz) = pal.polmo( 0.7, -0.5, 1.0e-6, -2.0e-6)

        self.assertAlmostEqual( elong, 0.7000004837322044, 12 )
        self.assertAlmostEqual( phi, -0.4999979467222241, 12 )
        self.assertAlmostEqual( daz, 1.008982781275728e-6, 12 )

    def test_prebn(self):
        expected = np.array( [
            [ 9.999257613786738e-1, -1.117444640880939e-2, -4.858341150654265e-3 ],
            [ 1.117444639746558e-2,  9.999375635561940e-1, -2.714797892626396e-5 ],
            [ 4.858341176745641e-3, -2.714330927085065e-5,  9.999881978224798e-1 ]
        ] )
        rmatp = pal.prebn( 1925, 1975 )
        np.testing.assert_array_almost_equal( rmatp, expected, 12 )

    def test_prec(self):
        expected = np.array( [
            [ 0.9999856154510, -0.0049192906204,    -0.0021376320580 ],
            [  0.0049192906805,  0.9999879002027,    -5.2297405698747e-06 ],
            [ 0.0021376319197, -5.2859681191735e-06, 0.9999977152483 ] ] )
        rmat = pal.prec( 1990, 2012 )
        np.testing.assert_array_almost_equal( rmat, expected, 12 )

    def test_preces(self):
        (ra, dc) = pal.preces( "FK4", 1925, 1950, 6.28, -1.123 )
        self.assertAlmostEqual( ra, 0.002403604864728447, 12 )
        self.assertAlmostEqual( dc, -1.120570643322045, 12 )

        (ra, dec) = pal.preces( "FK5", 2050, 1990, 0.0123, 1.0987 )
        self.assertAlmostEqual( ra, 6.282003602708382, 5 )
        self.assertAlmostEqual( dc, -1.120570643322045, 6 )

    def test_pvobs(self):
        expected = np.array( [ -4.7683600138836167813e-06,
                           1.0419056712717953176e-05,
                           4.099831053320363277e-05,
                          -7.5976959740661272483e-10,
                          -3.4771429582640930371e-10,
                           0.0 ] )
        pv = pal.pvobs( 1.3, 10000, 2 )
        np.testing.assert_array_almost_equal( pv, expected, decimal=12 )

    def test_range(self):
        self.assertAlmostEqual( pal.drange( -4 ), 2.283185307179586, 12 )

    def test_ranorm(self):
        self.assertAlmostEqual( pal.dranrm( -0.1 ), 6.183185307179587, 12 )

    def test_ref(self):
        self.assertAlmostEqual( pal.refro(1.4, 3456.7, 280, 678.9, 0.9, 0.55,
                                -0.3, 0.006, 1e-9 ),
                                0.00106715763018568, 12 )
        self.assertAlmostEqual( pal.refro(1.4, 3456.7, 280, 678.9, 0.9, 1000,
                                -0.3, 0.006, 1e-9),
                                0.001296416185295403, 12 )

        (refa, refb) = pal.refcoq( 275.9, 709.3, 0.9, 101 )
        self.assertAlmostEqual( refa, 2.324736903790639e-4, 12 )
        self.assertAlmostEqual( refb, -2.442884551059e-7, 15 )

        (refa, refb) = pal.refco( 2111.1, 275.9, 709.3, 0.9, 101,
                                -1.03, 0.0067, 1e-12 )
        self.assertAlmostEqual( refa, 2.324673985217244e-4, 12 )
        self.assertAlmostEqual( refb, -2.265040682496e-7, 15 )

        (refa, refb) = pal.refcoq( 275.9, 709.3, 0.9, 0.77 )
        self.assertAlmostEqual( refa, 2.007406521596588e-4, 12 )
        self.assertAlmostEqual( refb, -2.264210092590e-7, 15 )

        (refa, refb) = pal.refco( 2111.1, 275.9, 709.3, 0.9, 0.77,
                                -1.03, 0.0067, 1e-12 )
        self.assertAlmostEqual( refa, 2.007202720084551e-4, 12 )
        self.assertAlmostEqual( refb, -2.223037748876e-7, 15 )

        (refa2, refb2) = pal.atmdsp( 275.9, 709.3, 0.9, 0.77,
                                     refa, refb, 0.5 )
        self.assertAlmostEqual( refa2, 2.034523658888048e-4, 12 )
        self.assertAlmostEqual( refb2, -2.250855362179e-7, 15 )

        vu = pal.dcs2c( 0.345, 0.456 )
        vr = pal.refv( vu, refa, refb )
        self.assertAlmostEqual( vr[0], 0.8447487047790478, 12 )
        self.assertAlmostEqual( vr[1], 0.3035794890562339, 12 )
        self.assertAlmostEqual( vr[2], 0.4407256738589851, 12 )

        vu = pal.dcs2c( 3.7, 0.03 )
        vr = pal.refv( vu, refa, refb )
        self.assertAlmostEqual( vr[0], -0.8476187691681673, 12 )
        self.assertAlmostEqual( vr[1], -0.5295354802804889, 12 )
        self.assertAlmostEqual( vr[2], 0.0322914582168426, 12 )

        zr = pal.refz( 0.567, refa, refb )
        self.assertAlmostEqual( zr, 0.566872285910534, 12 )

        zr = pal.refz( 1.55, refa, refb )
        self.assertAlmostEqual( zr, 1.545697350690958, 12 )

        np.random.seed(32)
        zuIn = np.random.sample(20)*1.4
        zrControl = np.zeros(20, dtype=np.float64)
        for i, zu in enumerate(zuIn):
            zr = pal.refz(zu, refa, refb)
            zrControl[i] = zr
        zrTest = pal.refzVector(zuIn, refa, refb)
        for zt, zc in zip(zrTest, zrControl):
            self.assertAlmostEqual(zt, zc, 12)

    def test_refc(self): # This is the SOFA test
        (refa, refb) = pal.refcoq( 10.0+273.15, 800.0, 0.9, 0.4)
        self.assertAlmostEqual( refa, 0.2264949956241415009e-3, 15 )
        self.assertAlmostEqual( refb, -0.2598658261729343970e-6, 18 )

    def test_rv(self):
        self.assertAlmostEqual( pal.rverot( -0.777, 5.67, -0.3, 3.19 ),
                                 -0.1948098355075913, 6 )
        self.assertAlmostEqual( pal.rvgalc( 1.11, -0.99 ),
                                158.9630759840254, 3 )
        self.assertAlmostEqual( pal.rvlg( 3.97, 1.09 ),
                                -197.818762175363, 3 )
        self.assertAlmostEqual( pal.rvlsrd( 6.01, 0.1 ),
                                -4.082811335150567, 4 )
        self.assertAlmostEqual( pal.rvlsrk( 6.01, 0.1 ),
                                -5.925180579830265, 4 )

    def test_rvgalc(self):
        self.assertAlmostEqual( pal.rvgalc(2.7, -1.0), 213.98084425751144977, 12 )

    def test_rvlg(self):
        self.assertAlmostEqual( pal.rvlg(2.7, -1.0), 291.79205281252404802, 12 )

    def test_rvlsrd(self):
        self.assertAlmostEqual( pal.rvlsrd(2.7, -1.0), 9.620674692097630043, 12 )

    def test_rvlsrk(self):
        self.assertAlmostEqual( pal.rvlsrk(2.7, -1.0), 12.556356851411955233, 12 )

    def test_sep(self):
        d1 = np.array( [ 1.0, 0.1, 0.2 ] )
        d2 = np.array( [ -3.0, 1e-3, 0.2 ] )
        (ad1, bd1) = pal.dcc2s( d1 )
        (ad2, bd2) = pal.dcc2s( d2 )

        self.assertAlmostEqual( pal.dsep( ad1, bd1, ad2, bd2 ),
                                2.8603919190246608, 7 )
        self.assertAlmostEqual( pal.dsepv( d1, d2 ),
                                2.8603919190246608, 7 )

    def test_dsepVector(self):
        np.random.seed(32)
        ra1 = np.random.sample(20)*2.0*np.pi
        dec1 = np.random.sample(20)*2.0*np.pi
        ra2 = np.random.sample(20)*2.0*np.pi
        dec2 = np.random.sample(20)*2.0*np.pi
        ddControl = np.zeros(20, dtype=np.float64)
        for (i, rr) in enumerate(zip(ra1, dec1, ra2, dec2)):
            dd = pal.dsep(rr[0], rr[1], rr[2], rr[3])
            ddControl[i] = dd
        ddTest = pal.dsepVector(ra1, dec1, ra2, dec2)
        for (ddc, ddt) in zip (ddTest, ddControl):
            self.assertAlmostEqual(ddc, ddt, 12)

    def test_supgal(self):
        (dl, db) = pal.supgal( 6.1, -1.4 )
        self.assertAlmostEqual( dl, 3.798775860769474, 12 )
        self.assertAlmostEqual( db,  -0.1397070490669407, 12 )

    def test_tp(self):
        dr0 = 3.1
        dd0 = -0.9
        dr1 = dr0 + 0.2
        dd1 = dd0 - 0.1
        (dx, dy) = pal.ds2tp( dr1, dd1, dr0, dd0 )
        self.assertAlmostEqual( dx, 0.1086112301590404, 12 )
        self.assertAlmostEqual( dy, -0.1095506200711452, 12 )

        (dr2, dd2) = pal.dtp2s( dx, dy, dr0, dd0 )
        self.assertAlmostEqual( dr2 - dr1, 0.0, 12 )
        self.assertAlmostEqual( dd2 - dd1, 0.0, 12 )

        (dr01, dd01, dr02, dd02) = pal.dtps2c( dx, dy, dr2, dd2 )
        self.assertAlmostEqual( dr01, dr0, 12 )
        self.assertAlmostEqual( dd01, dd0 )
        self.assertIsNone( dr02 )
        self.assertIsNone( dd02 )

    def test_ds2tpVector(self):
        raz = 0.3
        decz = 0.1
        nTests = 100
        np.random.seed(32)
        raIn = raz+(np.random.sample(nTests)-0.5)*0.2
        decIn = decz+(np.random.sample(nTests)-0.5)*0.2
        xiControl = None
        etaControl = None
        for (rr,dd) in zip(raIn,decIn):
            xi, eta = pal.ds2tp(rr, dd, raz, decz)
            if xiControl is None:
                xiControl = np.array([xi])
                etaControl = np.array([eta])
            else:
                xiControl = np.append(xiControl,xi)
                etaControl = np.append(etaControl,eta)

        xiTest, etaTest = pal.ds2tpVector(raIn, decIn, raz, decz)
        for (x1,e1,x2,e2) in zip(xiControl,etaControl,xiTest,etaTest):
            self.assertAlmostEqual(x1,x2,12)
            self.assertAlmostEqual(e1,e2,12)

    def test_vecmat(self):
        # Not everything is implemented here
        dav = np.array( [ -0.123, 0.0987, 0.0654 ] )
        dav2m_expected = np.array( [
            [  0.9930075842721269,  0.05902743090199868, -0.1022335560329612 ],
            [ -0.07113807138648245, 0.9903204657727545,  -0.1191836812279541 ],
            [  0.09420887631983825, 0.1256229973879967,   0.9875948309655174 ],
        ] )

        deuler_expected = np.array( [
            [ -0.1681574770810878,  0.1981362273264315,  0.9656423242187410 ],
            [ -0.2285369373983370,  0.9450659587140423, -0.2337117924378156 ],
            [ -0.9589024617479674, -0.2599853247796050, -0.1136384607117296 ] ] )

        dmxm_expected = np.array( [
            [ -0.09010460088585805,  0.3075993402463796,  0.9472400998581048 ],
            [ -0.3161868071070688,   0.8930686362478707, -0.3200848543149236 ],
            [ -0.9444083141897035,  -0.3283459407855694,  0.01678926022795169 ] ] )

        drm1 = pal.dav2m( dav )
        np.testing.assert_array_almost_equal( drm1, dav2m_expected, decimal=12 )

        drm2 = pal.deuler( "YZY", 2.345, -0.333, 2.222 )
        np.testing.assert_array_almost_equal( drm2, deuler_expected, decimal=12 )

        drm = pal.dmxm( drm2, drm1 )
        np.testing.assert_array_almost_equal( drm, dmxm_expected, decimal=12 )

        dv1 = pal.dcs2c( 3.0123, -0.999 )
        dcs2c_expected = np.array( [-0.5366267667260525, 0.06977111097651444, -0.8409302618566215] )
        np.testing.assert_array_almost_equal( dv1, dcs2c_expected, decimal=12 )

        dv2 = pal.dmxv( drm1, dv1 )
        dv3 = pal.dmxv( drm2, dv2 )
        dmxv_expected = np.array( [-0.7267487768696160, 0.5011537352639822, 0.4697671220397141] )
        np.testing.assert_array_almost_equal( dv3, dmxv_expected, decimal=12 )

        dv4 = pal.dimxv( drm, dv3 )
        dimxv_expected = np.array( [ -0.5366267667260526, 0.06977111097651445, -0.8409302618566215 ] )
        np.testing.assert_array_almost_equal( dv4, dimxv_expected, decimal=12 )

        dv5 = pal.dm2av( drm )
        dm2av_expected = np.array( [ 0.006889040510209034, -1.577473205461961, 0.5201843672856759 ] )
        np.testing.assert_array_almost_equal( dv5, dm2av_expected, decimal=12 )

        dv5 *= 1000.0

        (dv6, dvm) = pal.dvn( dv5 )
        dvn_expected = np.array( [0.004147420704640065, -0.9496888606842218, 0.3131674740355448] )
        np.testing.assert_array_almost_equal( dv6, dvn_expected, decimal=12 )
        self.assertAlmostEqual( dvm, 1661.042127339937, 9 )

        dv7 = pal.dvxv( dv6, dv1 )
        dvxv_expected = np.array( [ 0.7767720597123304, -0.1645663574562769, -0.5093390925544726 ] )
        np.testing.assert_array_almost_equal( dv7, dvxv_expected, decimal=12 )

if __name__ == '__main__':
    unittest.main()
