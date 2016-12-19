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

from __future__ import with_statement
import unittest
import palpy as pal
import numpy as np


class TestPal(unittest.TestCase):

    def test_addet(self):
        rm = 2.0
        dm = -1.0
        eq = 1975.0
        (r1, d1) = pal.addet(rm, dm, eq)
        self.assertAlmostEqual(r1 - rm, 2.983864874295250e-6, 11)
        self.assertAlmostEqual(d1 - dm, 2.379650804185118e-7, 11)

        (r2, d2) = pal.subet(r1, d1, eq)
        self.assertAlmostEqual(r2 - rm, 0.0, 11)
        self.assertAlmostEqual(d2 - dm, 0.0, 11)

    def test_afin(self):
        (d, i) = pal.dafin("12 34 56.7 |", 1)
        self.assertEqual(i, 12)
        self.assertAlmostEqual(d, 0.2196045986911432, 12)

        (d, i) = pal.dafin("45 00 00.000 ", 1)
        self.assertEqual(i, 14)
        self.assertAlmostEqual(d, pal.DPI / 4.0, 12)

        (d, i) = pal.dafin("30.2 < just decimal degrees", 1)
        self.assertEqual(i, 6)
        self.assertAlmostEqual(d, 30.2 * pal.DD2R, 12)

        (d, i) = pal.dafin("  30 23.6 < decimal armin", 1)
        self.assertEqual(i, 11)
        self.assertAlmostEqual(d, (30 + 23.6/60)*pal.DD2R, 12)

        (d, i) = pal.dafin(" offset into string: 45.0 <<<", 22)
        self.assertEqual(i, 27)
        self.assertAlmostEqual(d, pal.DPI / 4.0, 12)

        self.assertRaises(ValueError, pal.dafin, " not a sexagesimal string ", 1)
        self.assertRaises(ValueError, pal.dafin, " 750.4 21 0.0 bad degrees ", 1)
        self.assertRaises(ValueError, pal.dafin, " 45 30.5 0.0 bad arcminutes ", 1)
        self.assertRaises(ValueError, pal.dafin, " 45 -30 0.0 bad arcminutes ", 1)
        self.assertRaises(ValueError, pal.dafin, " 45 72 0.0 too many arcminutes ", 1)
        self.assertRaises(ValueError, pal.dafin, " 45 43 85.0 too many arcseconds ", 1)

    def test_airmass(self):
        self.assertAlmostEqual(pal.airmas(1.2354),
                               3.015698990074724, 11)

    def test_airmass_vector(self):
        """
        Test that airmassVector gives results consistent with airmass
        """
        np.random.seed(145)
        n_samples = 1000
        zd = np.random.random_sample(n_samples)*0.5*np.pi
        test_am = pal.airmasVector(zd)
        for ii in range(n_samples):
            control_am = pal.airmas(zd[ii])
            self.assertEqual(control_am, test_am[ii])
            self.assertFalse(np.isnan(test_am[ii]))

    def test_altaz(self):
        (az, azd, azdd, el, eld, eldd, pa, pad, padd) = pal.altaz(0.7, -0.7, -0.65)
        self.assertAlmostEqual(az, 4.400560746660174, 12)
        self.assertAlmostEqual(azd, -0.2015438937145421, 12)
        self.assertAlmostEqual(azdd, -0.4381266949668748, 13)
        self.assertAlmostEqual(el, 1.026646506651396, 12)
        self.assertAlmostEqual(eld, -0.7576920683826450, 13)
        self.assertAlmostEqual(eldd, 0.04922465406857453, 14)
        self.assertAlmostEqual(pa, 1.707639969653937, 12)
        self.assertAlmostEqual(pad, 0.4717832355365627, 13)
        self.assertAlmostEqual(padd, -0.2957914128185515, 13)

    def test_altaz_vector(self):
        np.random.seed(32)
        phi = 0.5
        ha_in = np.random.sample(20)*2.0*np.pi
        dec_in = (np.random.sample(20)-0.5)*np.pi
        az_c = np.zeros(20, dtype=np.float64)
        azd_c = np.zeros(20, dtype=np.float64)
        azdd_c = np.zeros(20, dtype=np.float64)
        el_c = np.zeros(20, dtype=np.float64)
        eld_c = np.zeros(20, dtype=np.float64)
        eldd_c = np.zeros(20, dtype=np.float64)
        pa_c = np.zeros(20, dtype=np.float64)
        pad_c = np.zeros(20, dtype=np.float64)
        padd_c = np.zeros(20, dtype=np.float64)
        for i in range(20):
            az, azd, azdd, el, eld, eldd, pa, pad, padd = pal.altaz(ha_in[i], dec_in[i], phi)
            az_c[i] = az
            azd_c[i] = azd
            azdd_c[i] = azdd
            el_c[i] = el
            eld_c[i] = eld
            eldd_c[i] = eldd
            pa_c[i] = pa
            pad_c[i] = pad
            padd_c[i] = padd
        azT, azdT, azddT, elT, eldT, elddT, paT, padT, paddT = pal.altazVector(ha_in, dec_in, phi)

        for i in range(20):
            self.assertAlmostEqual(az_c[i], azT[i], 12)
            self.assertAlmostEqual(azd_c[i], azdT[i], 12)
            self.assertAlmostEqual(azdd_c[i], azddT[i], 12)
            self.assertAlmostEqual(el_c[i], elT[i], 12)
            self.assertAlmostEqual(eld_c[i], eldT[i], 12)
            self.assertAlmostEqual(eldd_c[i], elddT[i], 12)
            self.assertAlmostEqual(pa_c[i], paT[i], 12)
            self.assertAlmostEqual(pad_c[i], padT[i], 12)
            self.assertAlmostEqual(padd_c[i], paddT[i], 12)

        # test that an exception is raised if input arrays have
        # different lengths
        self.assertRaises(ValueError, pal.altazVector, ha_in[:10], dec_in, phi)

    def test_amp(self):
        (rm, dm) = pal.amp(2.345, -1.234, 50100., 1990.)
        self.assertAlmostEqual(rm, 2.344472180027961, 6)
        self.assertAlmostEqual(dm, -1.233573099847705, 7)
        (rm, dm) = pal.amp(1.234, -0.567, 55927., 2010.)
        self.assertAlmostEqual(rm, 1.233512033578303857, 12)
        self.assertAlmostEqual(dm, -0.56702909748530827549, 12)

    def test_ampqk(self):
        amprms = pal.mappa(2010.0, 55927.0)
        (rm, dm) = pal.ampqk(1.234, -0.567, amprms)
        self.assertAlmostEqual(rm, 1.233512033578303857, 11)
        self.assertAlmostEqual(dm, -0.56702909748530827549, 11)

    def test_ampqk_vector(self):
        """
        Test that ampqkVector produces results consistent with ampqk
        """
        np.random.seed(144)
        n_samples = 200
        amprms = pal.mappa(2010.0, 55927.0)
        ra_in = np.random.random_sample(n_samples)*2.0*np.pi
        dec_in = (np.random.random_sample(n_samples)-0.5)*np.pi

        testRa, testDec = pal.ampqkVector(ra_in, dec_in, amprms)

        for ii in range(n_samples):
            controlRa, controlDec = pal.ampqk(ra_in[ii], dec_in[ii], amprms)
            self.assertEqual(controlRa, testRa[ii])
            self.assertEqual(controlDec, testDec[ii])
            self.assertFalse(np.isnan(testRa[ii]))
            self.assertFalse(np.isnan(testDec[ii]))

        # test that ampqkVector and mapqkzVector invert each other
        ra_roundtrip, dec_roundtrip = pal.mapqkzVector(testRa, testDec, amprms)
        pal.DR2AS = 3600.0*np.degrees(1.0)
        distance = pal.DR2AS*pal.dsepVector(ra_roundtrip, dec_roundtrip,
                                            ra_in, dec_in)
        np.testing.assert_array_almost_equal(distance, np.zeros(n_samples), 9)

        # test that exceptions are raised when input arrays are not of the same
        # length
        self.assertRaises(ValueError, pal.ampqkVector, ra_in[:17], dec_in, amprms)

    def test_aopqk_vector(self):
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
        wl = 0.45
        obsrms = pal.aoppa(date, dut, elongm, phim, hm, xp, yp, tdk, pmb, rh, wl, tlr)
        np.random.seed(32)
        n_tests = 100
        ra_in = np.random.sample(n_tests)*2.0*np.pi
        dec_in = (np.random.sample(n_tests)-0.5)*np.pi
        az_control = None
        ze_control = None
        ha_control = None
        d_control = None
        r_control = None
        for (rr, dd) in zip(ra_in, dec_in):
            az, ze, ha, d, r = pal.aopqk(rr, dd, obsrms)
            if az_control is None:
                az_control = np.array([az])
                ze_control = np.array([ze])
                ha_control = np.array([ha])
                d_control = np.array([d])
                r_control = np.array([r])
            else:
                az_control = np.append(az_control, az)
                ze_control = np.append(ze_control, ze)
                ha_control = np.append(ha_control, ha)
                d_control = np.append(d_control, d)
                r_control = np.append(r_control, r)

        azTest, zeTest, haTest, dTest, rTest = pal.aopqkVector(ra_in, dec_in, obsrms)
        for (a1, z1, h1, d1, r1, a2, z2, h2, d2, r2) in \
            zip(az_control, ze_control, ha_control, d_control, r_control,
                azTest, zeTest, haTest, dTest, rTest):

            self.assertAlmostEqual(a1, a2, 12)
            self.assertAlmostEqual(z1, z2, 12)
            self.assertAlmostEqual(h1, h2, 12)
            self.assertAlmostEqual(d1, d2, 12)
            self.assertAlmostEqual(r1, r2, 12)

        # test that an exception is raised if input arrays have
        # different lengths
        self.assertRaises(ValueError, pal.aopqkVector, ra_in[:6], dec_in, obsrms)

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
        aoptol = [[10, 7, 7, 8, 7],
                  [10, 10, 10, 10, 10],
                  [10, 10, 10, 10, 10]]

        for i in range(len(aopres)):
            # Not very pythonic
            if i == 0:
                rap = 2.7
                wl = 0.45
            elif i == 1:
                rap = 2.345
            else:
                wl = 1.0e6

            result = pal.aop(rap, dap, date, dut, elongm, phim, hm, xp, yp,
                             tdk, pmb, rh, wl, tlr)

            for j in range(len(result)):
                self.assertAlmostEqual(result[j], aopres[i][j], aoptol[i][j])

        date = 48000.3
        wl = 0.45

        aoprms = pal.aoppa(date, dut, elongm, phim, hm, xp, yp, tdk, pmb,
                           rh, wl, tlr)

        aoppares = [0.4999993892136306, 0.4794250025886467, 0.8775828547167932,
                    1.363180872136126e-6, 3000., 280., 550., 0.6, 0.45, 0.006,
                    0.0001562803328459898, -1.792293660141e-7, 2.101874231495843,
                    7.601916802079765]
        aoppatol = [13, 13, 13, 13, 10, 11, 11, 13, 13, 15, 13, 13, 12, 8]
        self.assertEqual(len(aoprms), len(aoppares))
        for i in range(len(aoprms)):
            self.assertAlmostEqual(aoprms[i], aoppares[i], aoppatol[i])

        (rap, dap) = pal.oap("r", 1.6, -1.01, date, dut, elongm, phim,
                             hm, xp, yp, tdk, pmb, rh, wl, tlr)
        self.assertAlmostEqual(rap, 1.601197569844787, 10)
        self.assertAlmostEqual(dap, -1.012528566544262, 10)

        (rap, dap) = pal.oap("h", -1.234, 2.34, date, dut, elongm, phim,
                             hm, xp, yp, tdk, pmb, rh, wl, tlr)
        self.assertAlmostEqual(rap, 5.693087688154886463, 10)
        self.assertAlmostEqual(dap, 0.8010281167405444, 10)

        (rap, dap) = pal.oap("a", 6.1, 1.1, date, dut, elongm, phim,
                             hm, xp, yp, tdk, pmb, rh, wl, tlr)
        self.assertAlmostEqual(rap, 5.894305175192448940, 10)
        self.assertAlmostEqual(dap, 1.406150707974922, 10)

        (rap, dap) = pal.oapqk("r", 2.1, -0.345, aoprms)
        self.assertAlmostEqual(rap, 2.10023962776202, 10)
        self.assertAlmostEqual(dap, -0.3452428692888919, 10)

        (rap, dap) = pal.oapqk("h", -0.01, 1.03, aoprms)
        self.assertAlmostEqual(rap, 1.328731933634564995, 10)
        self.assertAlmostEqual(dap, 1.030091538647746, 10)

        (rap, dap) = pal.oapqk("a", 4.321, 0.987, aoprms)
        self.assertAlmostEqual(rap, 0.4375507112075065923, 10)
        self.assertAlmostEqual(dap, -0.01520898480744436, 10)

        aoprms = pal.aoppat(date + pal.DS2R, aoprms)
        self.assertAlmostEqual(aoprms[13], 7.602374979243502, 8)

    def test_oapqk_vector(self):
        """
        Test that oapqkVector gives results consistent with oapqk
        """

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
        wl = 0.45
        aoprms = pal.aoppa(date, dut, elongm, phim, hm, xp, yp, tdk, pmb,
                           rh, wl, tlr)

        np.random.seed(133)
        n_samples = 200
        ob1 = np.random.random_sample(n_samples) * 2.0 * np.pi
        ob2 = np.random.random_sample(n_samples) * 0.5 * np.pi
        # restrict ob2 to 0 < ob2 < pi/2 because we will also be testing az-zenith distance
        # coordinate pairs

        for typeFlag in ['r', 'a', 'h']:
            testRa, testDec = pal.oapqkVector(typeFlag, ob1, ob2, aoprms)
            for ii in range(n_samples):
                controlRa, controlDec = pal.oapqk(typeFlag, ob1[ii], ob2[ii], aoprms)
                self.assertEqual(testRa[ii], controlRa)
                self.assertEqual(testDec[ii], controlDec)

        # verify that pal.aopqkVector and pal.oapqkVector invert each other;
        # we limit our test raApparent, decApparent values to be within 75 degrees
        # of zenith, because of the limited accuracy of these rout_ines at
        # large zenith distance (even though we approximate the sky as a flat
        # surface and demand that all of the sample points be within 50 degrees
        # of zenith in this gross approximation, heuristically, this causes
        # the actual maximum zenith distance to be 74.1 degrees).
        raApCenter, decApCenter = pal.oapqk('h', 0.0, 0.0, aoprms)

        rr = np.random.random_sample(n_samples)*np.radians(50.0)
        theta = np.random.random_sample(n_samples)*2.0*np.pi

        ra_ap_list = raApCenter + rr*np.cos(theta)
        dec_ap_list = decApCenter + rr*np.sin(theta)

        azList, zdList, haList, \
            decObList, raObList = pal.aopqkVector(ra_ap_list, dec_ap_list, aoprms)

        testRa, testDec = pal.oapqkVector('r', raObList, decObList, aoprms)
        np.testing.assert_array_almost_equal(testRa, ra_ap_list, 12)
        np.testing.assert_array_almost_equal(testDec, dec_ap_list, 12)

        testRa, testDec = pal.oapqkVector('h', haList, decObList, aoprms)
        np.testing.assert_array_almost_equal(testRa, ra_ap_list, 12)
        np.testing.assert_array_almost_equal(testDec, dec_ap_list, 12)

        testRa, testDec = pal.oapqkVector('a', azList, zdList, aoprms)
        np.testing.assert_array_almost_equal(testRa, ra_ap_list, 12)
        np.testing.assert_array_almost_equal(testDec, dec_ap_list, 12)

        # test that an exception is thrown if the input arrays do not
        # have the same length
        self.assertRaises(ValueError, pal.oapqkVector, 'a',
                          azList, zdList[:17], aoprms)

    def test_bear(self):
        a1 = 1.234
        b1 = -0.123
        a2 = 2.345
        b2 = 0.789
        self.assertAlmostEqual(pal.dbear(a1, b1, a2, b2),
                               0.7045970341781791, 12)
        d1 = pal.dcs2c(a1, b1)
        d2 = pal.dcs2c(a2, b2)
        self.assertAlmostEqual(pal.dpav(d1, d2), 0.7045970341781791, 12)

    def test_dbear_vector(self):
        """
        Test that dbearVector gives the same
        results as dbear
        """
        np.random.seed(122)
        n_samples = 100
        a1_in = np.random.random_sample(n_samples)*2.0*np.pi
        b1_in = (np.random.random_sample(n_samples)-0.5)*np.pi
        a2_in = np.random.random_sample(n_samples)*2.0*np.pi
        b2_in = (np.random.random_sample(n_samples)-0.5)*np.pi

        # test case where a2, b2 have the same number of elements
        # as a1, b1
        b_test = pal.dbearVector(a1_in, b1_in, a2_in, b2_in)
        for i in range(len(a1_in)):
            b_control = pal.dbear(a1_in[i], b1_in[i], a2_in[i], b2_in[i])
            self.assertEqual(b_test[i], b_control)

        # test case where a2, b2 have only one element
        b_test = pal.dbearVector(a1_in, b1_in, a2_in[7:8], b2_in[7:8])
        for i in range(len(a1_in)):
            b_control = pal.dbear(a1_in[i], b1_in[i], a2_in[7], b2_in[7])
            self.assertEqual(b_test[i], b_control)

        # test that exceptions are raied when the numpy arrays are of
        # incorrect size
        self.assertRaises(ValueError, pal.dbearVector, a1_in, b1_in[:9], a2_in, b2_in)

        self.assertRaises(ValueError, pal.dbearVector, a1_in, b1_in, a2_in, b2_in[:8])

        self.assertRaises(ValueError, pal.dbearVector, a1_in, b1_in, a2_in[:9], b2_in[:9])

    def test_dpav_vector(self):
        """
        Test that dpavVector is consistent with dpav
        """
        np.random.seed(127)
        n_samples = 200
        phi = np.random.random_sample(n_samples)*2.0*np.pi
        theta = (np.random.random_sample(n_samples)-0.5)*np.pi

        v1 = np.array([[np.cos(th), np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)]
                       for th, ph in zip(theta, phi)])

        phi = np.random.random_sample(n_samples)*2.0*np.pi
        theta = (np.random.random_sample(n_samples)-0.5)*np.pi

        v2 = np.array([[np.cos(th), np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph)]
                       for th, ph in zip(theta, phi)])

        test_pa = pal.dpavVector(v1, v2)
        for ii in range(n_samples):
            pa_control = pal.dpav(v1[ii], v2[ii])
            self.assertEqual(pa_control, test_pa[ii])
            self.assertFalse(np.isnan(pa_control))

        # test the case where we only feed one point in as v2
        test_pa = pal.dpavVector(v1, v2[4])
        for ii in range(n_samples):
            pa_control = pal.dpav(v1[ii], v2[4])
            self.assertEqual(pa_control, test_pa[ii])
            self.assertFalse(np.isnan(pa_control))

        # test that exceptions are raised when they should be
        self.assertRaises(ValueError, pal.dpavVector, v1, v2[:19])

        v3 = np.random.random_sample((n_samples, 6))
        self.assertRaises(ValueError, pal.dpavVector, v3, v2)

        self.assertRaises(ValueError, pal.dpavVector, v1, v3)

        self.assertRaises(ValueError, pal.dpavVector, v1,
                          np.random.random_sample(9))

    def test_caldj(self):
        djm = pal.caldj(1999, 12, 31)
        self.assertEqual(djm, 51543)

        self.assertRaises(ValueError, pal.caldj, -5000, 12, 31)
        self.assertRaises(ValueError, pal.caldj, 1970, 13, 1)
        self.assertRaises(ValueError, pal.caldj, 1970, 1, 32)

    def test_caldj_vector(self):
        """
        Test that caldjVector gives results consistent with cadj
        """
        np.random.seed(143)
        n_samples = 200
        iy = np.random.randint(1000, 10001, n_samples)
        im = np.random.randint(1, 12, n_samples)
        iday = np.random.randint(1, 21, n_samples)

        test_mjd = pal.caldjVector(iy, im, iday)

        for ii in range(n_samples):
            control_mjd = pal.caldj(iy[ii], im[ii], iday[ii])
            self.assertEqual(control_mjd, test_mjd[ii])

        # test that caldjVector and djcalVector invert each other
        iyTest, imTest, idTest, fracTest = pal.djcalVector(5, test_mjd)
        np.testing.assert_array_equal(iyTest, iy)
        np.testing.assert_array_equal(imTest, im)
        np.testing.assert_array_equal(idTest, iday)

        # test that bad values are handled properly
        iy[5] = -4800
        im[9] = 13
        im[11] = 0
        iday[17] = 33

        test_mjd = pal.caldjVector(iy, im, iday)
        for ii in range(n_samples):
            if ii in (5, 9, 11, 17):
                self.assertTrue(np.isnan(test_mjd[ii]))
                self.assertRaises(ValueError, pal.caldj, iy[ii], im[ii],
                                  iday[ii])
            else:
                control_mjd = pal.caldj(iy[ii], im[ii], iday[ii])
                self.assertEqual(control_mjd, test_mjd[ii])

        # test that exceptions are raised when the input arrays are
        # of different lengths
        self.assertRaises(ValueError, pal.caldjVector, iy, im[:11], iday)

        self.assertRaises(ValueError, pal.caldjVector, iy, im, iday[:8])

    def test_daf2r(self):
        dr = pal.daf2r(76, 54, 32.1)
        self.assertAlmostEqual(dr, 1.342313819975276, 12)

    def test_daf2r_vector(self):
        """
        Test that daf2rVector returns the same results as daf2r
        """
        np.random.seed(123)
        n_samples = 100
        deg = np.random.randint(0, 360, n_samples)
        i_min = np.random.randint(0, 60, n_samples)
        asec = np.random.random_sample(n_samples)*60.0

        radian_test = pal.daf2rVector(deg, i_min, asec)
        for ii in range(len(deg)):
            radian_control = pal.daf2r(deg[ii], i_min[ii], asec[ii])
            self.assertEqual(radian_control, radian_test[ii])

        # test that bad values get set to np.NaN
        deg[9] = 360
        i_min[17] = 60
        asec[21] = 60.0

        radian_test = pal.daf2rVector(deg, i_min, asec)
        for ii, rad in enumerate(radian_test):
            if ii in (9, 17, 21):
                self.assertTrue(np.isnan(rad))
            else:
                self.assertFalse(np.isnan(rad))

    def test_cc2s(self):
        (da, db) = pal.dcc2s(np.array([100., -50., 25.]))
        self.assertAlmostEqual(da, -0.4636476090008061, 12)
        self.assertAlmostEqual(db, 0.2199879773954594, 12)

    def test_dcc2s_vector(self):
        """
        Test that dcc2sVector returns the same results as dcc2s
        """
        np.random.seed(124)
        input_data = np.random.random_sample((20, 3)) * 10.0
        aTest, b_test = pal.dcc2sVector(input_data)
        for ix in range(input_data.shape[0]):
            aControl, b_control = pal.dcc2s(input_data[ix])
            self.assertEqual(aControl, aTest[ix])
            self.assertEqual(b_control, b_test[ix])

        # test that an exception is raised if you don't pass in
        # 3-D cartesian points
        dummy_data = np.random.random_sample((20, 5))*10.0
        self.assertRaises(ValueError, pal.dcc2sVector, dummy_data)

    def test_dcs2c_vector(self):
        """
        Test that dcs2cVector returns the same results as dcs2c
        """
        np.random.seed(125)
        n_samples = 100
        ra = np.random.random_sample(n_samples)*2.0*np.pi
        dec = (np.random.random_sample(n_samples)-0.5)*np.pi
        v_test = pal.dcs2cVector(ra, dec)
        for ii in range(n_samples):
            v_control = pal.dcs2c(ra[ii], dec[ii])
            np.testing.assert_array_equal(v_control, v_test[ii])

        # test that dcs2cVector and dcc2sVector do, in fact, invert one another
        aTest, b_test = pal.dcc2sVector(v_test)
        np.testing.assert_array_almost_equal(np.cos(ra), np.cos(aTest), 12)
        np.testing.assert_array_almost_equal(np.sin(ra), np.sin(aTest), 12)
        np.testing.assert_array_almost_equal(np.cos(dec), np.cos(b_test), 12)
        np.testing.assert_array_almost_equal(np.sin(dec), np.sin(b_test), 12)

        # test that an exception is raised if you pass in arrays of different lengths
        self.assertRaises(ValueError, pal.dcs2cVector, ra, dec[:16])

    def test_cd2tf(self):
        (sign, hours, minutes, seconds, fraction) = pal.dd2tf(4, -0.987654321)
        self.assertEqual(sign, "-")
        self.assertEqual(hours, 23)
        self.assertEqual(minutes, 42)
        self.assertEqual(seconds, 13)
        self.assertEqual(fraction, 3333)

    def test_dd2tf_vector(self):
        """
        Test that dd2tfVector gives the same results as dd2tf
        """
        np.random.seed(126)
        n_samples = 100
        days_list = (np.random.sample(n_samples)-0.5)*1200.0

        for ndp in [2, 3, 4, 5]:
            testSign, testIh, testIm, testIs, testFrac = pal.dd2tfVector(ndp, days_list)
            for ix, days in enumerate(days_list):
                controlSign, controlIh,\
                    controlIm, controlIs,\
                    controlFrac = pal.dd2tf(ndp, days)

                self.assertEqual(controlSign, testSign[ix])
                self.assertEqual(controlIm, testIm[ix])
                self.assertEqual(controlIs, testIs[ix])
                self.assertEqual(controlFrac, testFrac[ix])

    def test_cldj(self):
        d = pal.cldj(1899, 12, 31)
        self.assertEqual(d, 15019)
        self.assertRaises(ValueError, pal.cldj, -5000, 12, 31)
        self.assertRaises(ValueError, pal.cldj, 1970, 13, 1)
        self.assertRaises(ValueError, pal.cldj, 1970, 1, 32)

    def test_dr2af(self):
        (sign, deg, min, sec, f) = pal.dr2af(4, 2.345)
        self.assertEqual(sign, "+")
        self.assertEqual(deg, 134)
        self.assertEqual(min, 21)
        self.assertEqual(sec, 30)
        self.assertEqual(f, 9706)

    def test_dr2af_vector(self):
        """
        Test that dr2afVector produces the same results as
        dr2af
        """
        np.random.seed(128)
        n_samples = 200
        angle_list = (np.random.random_sample(n_samples)-0.5)*4.0*np.pi
        for npd in [2, 3, 4, 5]:
            testSign, testDeg, \
                testMin, testSec, testFrac = pal.dr2afVector(npd, angle_list)

            for ii in range(n_samples):
                controlSign, controlDeg, \
                    controlMin, controlSec, controlFrac = pal.dr2af(npd, angle_list[ii])

                self.assertEqual(controlSign, testSign[ii])
                self.assertEqual(controlDeg, testDeg[ii])
                self.assertEqual(controlMin, testMin[ii])
                self.assertEqual(controlSec, testSec[ii])
                self.assertEqual(controlFrac, testFrac[ii])

    def test_dr2tf(self):
        (sign, hr, min, sec, f) = pal.dr2tf(4, -3.01234)
        self.assertEqual(sign, "-")
        self.assertEqual(hr, 11)
        self.assertEqual(min, 30)
        self.assertEqual(sec, 22)
        self.assertEqual(f, 6484)

    def test_dr2tf_vector(self):
        """
        Test that dr2tfVector produces the same results as
        dr2tf
        """
        np.random.seed(128)
        n_samples = 200
        angle_list = (np.random.random_sample(n_samples)-0.5)*4.0*np.pi
        for npd in [2, 3, 4, 5]:
            testSign, testHr, \
                testMin, testSec, testFrac = pal.dr2tfVector(npd, angle_list)

            for ii in range(n_samples):
                controlSign, controlHr, \
                    controlMin, controlSec, controlFrac = pal.dr2tf(npd, angle_list[ii])

                self.assertEqual(controlSign, testSign[ii])
                self.assertEqual(controlHr, testHr[ii])
                self.assertEqual(controlMin, testMin[ii])
                self.assertEqual(controlSec, testSec[ii])
                self.assertEqual(controlFrac, testFrac[ii])

    def test_dtf2d(self):
        dd = pal.dtf2d(23, 56, 59.1)
        self.assertAlmostEqual(dd, 0.99790625, 12)
        self.assertRaises(ValueError, pal.dtf2d, 24, 1, 32)
        self.assertRaises(ValueError, pal.dtf2d, 23, -1, 32)
        self.assertRaises(ValueError, pal.dtf2d, 23, 1, 60)

    def test_dtf2d_vector(self):
        """
        Test that dtf2dVector gives results consistent with
        dtf2d
        """
        np.random.seed(131)
        n_samples = 100
        i_hour = np.random.randint(0, 24, n_samples)
        i_min = np.random.randint(0, 60, n_samples)
        sec = np.random.random_sample(n_samples)*60.0

        test_days = pal.dtf2dVector(i_hour, i_min, sec)
        for ii in range(n_samples):
            control_days = pal.dtf2d(i_hour[ii], i_min[ii], sec[ii])
            self.assertEqual(test_days[ii], control_days)
            self.assertFalse(np.isnan(test_days[ii]))

        # test that NaN's are produced when bad inputs are given
        i_hour[5] = 24
        i_min[7] = 60
        sec[19] = 60.001
        test_days = pal.dtf2dVector(i_hour, i_min, sec)
        for ii in range(n_samples):
            if ii not in (5, 7, 19):
                control_days = pal.dtf2d(i_hour[ii], i_min[ii], sec[ii])
                self.assertEqual(test_days[ii], control_days)
                self.assertFalse(np.isnan(test_days[ii]))
            else:
                self.assertTrue(np.isnan(test_days[ii]))

        # test that exceptions are raised when you pass in
        # mis-matched input arrays
        self.assertRaises(ValueError, pal.dtf2dVector, i_hour, i_min[:19], sec)

        self.assertRaises(ValueError, pal.dtf2dVector, i_hour, i_min, sec[:19])

    def test_dtf2r(self):
        dr = pal.dtf2r(23, 56, 59.1)
        self.assertAlmostEqual(dr, 6.270029887942679, 12)
        self.assertRaises(ValueError, pal.dtf2r, 24, 1, 32)
        self.assertRaises(ValueError, pal.dtf2r, 23, -1, 32)
        self.assertRaises(ValueError, pal.dtf2r, 23, 1, 60)

    def test_dtf2r_vector(self):
        """
        Test that dtf2rVector gives results consistent with
        dtf2r
        """
        np.random.seed(131)
        n_samples = 100
        i_hour = np.random.randint(0, 24, n_samples)
        i_min = np.random.randint(0, 60, n_samples)
        sec = np.random.random_sample(n_samples)*60.0

        test_rad = pal.dtf2rVector(i_hour, i_min, sec)
        for ii in range(n_samples):
            control_rad = pal.dtf2r(i_hour[ii], i_min[ii], sec[ii])
            self.assertEqual(test_rad[ii], control_rad)
            self.assertFalse(np.isnan(test_rad[ii]))

        # test that NaN's are produced when bad inputs are given
        i_hour[5] = 24
        i_min[7] = 60
        sec[19] = 60.001
        test_rad = pal.dtf2rVector(i_hour, i_min, sec)
        for ii in range(n_samples):
            if ii not in (5, 7, 19):
                control_rad = pal.dtf2r(i_hour[ii], i_min[ii], sec[ii])
                self.assertEqual(test_rad[ii], control_rad)
                self.assertFalse(np.isnan(test_rad[ii]))
            else:
                self.assertTrue(np.isnan(test_rad[ii]))

        # test that exceptions are raised when you pass in
        # mis-matched input arrays
        self.assertRaises(ValueError, pal.dtf2rVector, i_hour, i_min[:19], sec)

        self.assertRaises(ValueError, pal.dtf2rVector, i_hour, i_min, sec[:19])

    def test_dat(self):
        self.assertEqual(pal.dat(43900), 18)
        self.assertAlmostEqual(pal.dtt(40404), 39.709746, 12)
        self.assertAlmostEqual(pal.dt(500), 4686.7, 10)
        self.assertAlmostEqual(pal.dt(1400), 408, 11)
        self.assertAlmostEqual(pal.dt(1950), 27.99145626)

    def test_djcal(self):
        djm = 50123.9999
        (y, m, d, f) = pal.djcal(4, djm)
        self.assertEqual(y, 1996)
        self.assertEqual(m, 2)
        self.assertEqual(d, 10)
        self.assertEqual(f, 9999)

        (y, m, d, f) = pal.djcl(djm)
        self.assertEqual(y, 1996)
        self.assertEqual(m, 2)
        self.assertEqual(d, 10)
        self.assertAlmostEqual(f, 0.9999, 7)

    def test_djcal_vector(self):
        """
        Test that djcalVector returns results consistent with djcal
        """
        np.random.seed(142)
        n_samples = 200
        mjd = (np.random.random_sample(n_samples)-0.5)*100000.0

        for ndp in [2, 3, 4, 5]:

            testY, testM, \
                testD, testFrac = pal.djcalVector(ndp, mjd)

            for ii in range(n_samples):
                controlY, controlM, \
                    controlD, controlFrac = pal.djcal(ndp, mjd[ii])

                self.assertEqual(controlY, testY[ii])
                self.assertEqual(controlM, testM[ii])
                self.assertEqual(controlD, testD[ii])
                self.assertEqual(controlFrac, testFrac[ii])

        ndp = 4
        # test that unacceptable dates result in -1 being placed in all
        # of the output slots
        mjd[4] = -2468570
        mjd[5] = -2468571
        mjd[9] = 1.0e8
        mjd[10] = 1.0e9
        testY, testM, testD, testFrac = pal.djcalVector(ndp, mjd)
        for ii in range(n_samples):
            if ii in (5, 10):
                self.assertEqual(testY[ii], -1)
                self.assertEqual(testM[ii], -1)
                self.assertEqual(testD[ii], -1)
                self.assertEqual(testFrac[ii], -1)
                self.assertRaises(ValueError, pal.djcal, ndp, mjd[ii])
            else:
                controlY, controlM, \
                    controlD, controlFrac = pal.djcal(ndp, mjd[ii])
                self.assertEqual(controlY, testY[ii])
                self.assertEqual(controlM, testM[ii])
                self.assertEqual(controlD, testD[ii])
                self.assertEqual(controlFrac, testFrac[ii])

    def test_dmat(self):
        da = np.array(
            [[2.22, 1.6578, 1.380522],
             [1.6578, 1.380522, 1.22548578],
             [1.380522, 1.22548578, 1.1356276122]])
        dv = np.array([2.28625, 1.7128825, 1.429432225])
        (da, dv, dd) = pal.dmat(da, dv)
        self.assertAlmostEqual(dd, 0.003658344147359863, 12)
        self.assertAlmostEqual(dv[0], 1.002346480763383, 12)
        self.assertAlmostEqual(dv[1], 0.03285594016974583489, 12)
        self.assertAlmostEqual(dv[2], 0.004760688414885247309, 12)

        da = np.array([[0., 1.], [0., 1.]])
        dv = np.array([1., 1.])
        self.assertRaises(ArithmeticError, pal.dmat, da, dv)

    def test_e2h(self):
        dp = -0.7
        hin = -0.3
        din = -1.1
        (da, de) = pal.de2h(hin, din, dp)
        self.assertAlmostEqual(da, 2.820087515852369, 12)
        self.assertAlmostEqual(de, 1.132711866443304, 12)

        (dh, dd) = pal.dh2e(da, de, dp)
        self.assertAlmostEqual(dh, hin, 12)
        self.assertAlmostEqual(dd, din, 12)

    def test_de2h_vector(self):
        n_tests = 100
        phi = 0.35
        np.random.seed(32)
        ha_in = np.random.random_sample(n_tests)*np.pi*2.0
        dec_in = (np.random.random_sample(n_tests)-0.5)*np.pi
        az_control = None
        el_control = None
        for (ha, dd) in zip(ha_in, dec_in):
            az, el = pal.de2h(ha, dd, phi)
            if az_control is None:
                az_control = np.array([az])
                el_control = np.array([el])
            else:
                az_control = np.append(az_control, az)
                el_control = np.append(el_control, el)

        azTest, elTest = pal.de2hVector(ha_in, dec_in, phi)
        for (a1, e1, a2, e2) in zip(az_control, el_control, azTest, elTest):
            self.assertAlmostEqual(a1, a2, 12)
            self.assertAlmostEqual(e1, e2, 12)

        # test that an exception is raised if inputs are not
        # of the same length
        self.assertRaises(ValueError, pal.de2hVector, ha_in[:10], dec_in, phi)

    def test_dh2e_vector(self):
        """
        Test that dh2eVector gives results consistent with dh2e
        """
        np.random.seed(142)
        n_samples = 200
        phi = 1.432
        az = np.random.random_sample(n_samples)*2.0*np.pi
        el = (np.random.random_sample(n_samples)-0.5)*np.pi

        testHa, testDec = pal.dh2eVector(az, el, phi)

        for ii in range(n_samples):
            controlHa, controlDec = pal.dh2e(az[ii], el[ii], phi)
            self.assertEqual(controlHa, testHa[ii])
            self.assertEqual(controlDec, testDec[ii])
            self.assertFalse(np.isnan(testHa[ii]))
            self.assertFalse(np.isnan(testDec[ii]))

        # test that dh2eVector and de2hVector invert each other
        testAz, testEl = pal.de2hVector(testHa, testDec, phi)
        pal.DR2ASs = 3600.0*np.degrees(1.0)
        distance = pal.DR2ASs*pal.dsepVector(testAz, testEl, az, el)
        np.testing.assert_array_almost_equal(distance,
                                             np.zeros(n_samples), 9)

        # test that an exception is raised when the input arrays
        # are of different lengths
        self.assertRaises(ValueError, pal.dh2eVector, az[:40], el, phi)

    def test_ecleq(self):
        (dr, dd) = pal.ecleq(1.234, -0.123, 43210.0)
        self.assertAlmostEqual(dr, 1.229910118208851, 5)
        self.assertAlmostEqual(dd, 0.2638461400411088, 5)

    def test_ecleq_vector(self):
        """
        Test that ecleqVector produces results consistent with
        ecleq
        """
        mjd = 58734.2
        np.random.seed(138)
        n_samples = 200
        dl = np.random.random_sample(n_samples)*2.0*np.pi
        db = (np.random.random_sample(n_samples)-0.5)*np.pi

        testRa, testDec = pal.ecleqVector(dl, db, mjd)

        for ii in range(n_samples):
            controlRa, controlDec = pal.ecleq(dl[ii], db[ii], mjd)
            self.assertEqual(controlRa, testRa[ii])
            self.assertEqual(controlDec, testDec[ii])

        # test that ecleqVector and eqeclVector invert
        # one another
        testDl, testDb = pal.eqeclVector(testRa, testDec, mjd)
        pal.DR2AS = 3600.0*np.degrees(1.0)
        distance = pal.DR2AS*pal.dsepVector(testDl, testDb, dl, db)
        np.testing.assert_array_almost_equal(distance,
                                             np.zeros(len(distance)), 4)

        # test that an exception is raised if input arrays are
        # of different lenghts
        self.assertRaises(ValueError, pal.ecleqVector, dl[:4], db, mjd)

    def test_ecmat(self):
        expected = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 0.91749307789883549624, 0.3977517467060596168],
            [0.0, -0.3977517467060596168, 0.91749307789883549624]])

        rmat = pal.ecmat(55966.46)
        np.testing.assert_array_almost_equal(rmat, expected, decimal=12)

    def test_epb(self):
        self.assertAlmostEqual(pal.epb(45123), 1982.419793168669, 8)

    def test_epb2d(self):
        self.assertAlmostEqual(pal.epb2d(1975.5), 42595.5995279655, 7)

    def test_epco(self):
        self.assertAlmostEqual(pal.epco("B", "J", 2000), 2000.001277513665, 7)
        self.assertAlmostEqual(pal.epco("J", "B", 1950), 1949.999790442300, 7)
        self.assertAlmostEqual(pal.epco("J", "j", 2000), 2000, 7)

    def test_epj(self):
        self.assertAlmostEqual(pal.epj(42999), 1976.603696098563, 7)

    def test_epj_vector(self):
        """
        Test that epjVector returns results consistent with epj
        """
        np.random.seed(45738)
        n_samples = 300
        date = 43000.0 + np.random.random_sample(n_samples)*10000.0
        test_epj = pal.epjVector(date)
        for ii in range(n_samples):
            control_epj = pal.epj(date[ii])
            self.assertEqual(control_epj, test_epj[ii])
            self.assertFalse(np.isnan(test_epj[ii]))

    def test_epj2d(self):
        self.assertAlmostEqual(pal.epj2d(2010.077), 55225.124250, 6)

    def test_epj2d_vector(self):
        """
        Test that epj2dVector returns results consistent with epj
        """
        np.random.seed(45367)
        n_samples = 300
        epj = 2000.0 + np.random.random_sample(n_samples)*50.0
        test_mjd = pal.epj2dVector(epj)
        for ii in range(n_samples):
            control_mjd = pal.epj2d(epj[ii])
            self.assertEqual(control_mjd, test_mjd[ii])
            self.assertFalse(np.isnan(test_mjd[ii]))

    def test_eqecl(self):
        (dl, db) = pal.eqecl(0.789, -0.123, 46555)
        self.assertAlmostEqual(dl, 0.7036566430349022, 6)
        self.assertAlmostEqual(db, -0.4036047164116848, 6)

    def test_eqecl_vector(self):
        """
        Test that eqeclVector produces results consistent with
        eqecl
        """
        mjd = 53000.0
        np.random.seed(137)
        n_samples = 200
        ra = np.random.random_sample(n_samples)*2.0*np.pi
        dec = (np.random.random_sample(n_samples)-0.5)*np.pi

        testDb, testDl = pal.eqeclVector(ra, dec, mjd)
        for ii in range(n_samples):
            controlDb, controlDl = pal.eqecl(ra[ii], dec[ii], mjd)
            self.assertEqual(controlDb, testDb[ii])
            self.assertEqual(controlDl, testDl[ii])

        # test that an exception is raised if inpu arrays are of
        # different lengths
        self.assertRaises(ValueError, pal.eqeclVector, ra, dec[:8], mjd)

    def test_eqeqx(self):
        self.assertAlmostEqual(pal.eqeqx(53736), -0.8834195072043790156e-5, 15)

    def test_eqeqx_vector(self):
        np.random.seed(32)
        date_in = 53000.0 + np.random.sample(20)*5000.0
        eq_control = np.zeros(20, dtype=np.float64)
        for i, d in enumerate(date_in):
            eq_control[i] = pal.eqeqx(d)
        eq_test = pal.eqeqxVector(date_in)
        for et, ec in zip(eq_test, eq_control):
            self.assertAlmostEqual(et, ec, 15)

    def test_eqgal(self):
        (dl, db) = pal.eqgal(5.67, -1.23)
        self.assertAlmostEqual(dl, 5.612270780904526, 12)
        self.assertAlmostEqual(db, -0.6800521449061520, 12)

    def test_eqgal_vector(self):
        np.random.seed(32)
        ra_in = np.random.sample(10)*2.0*np.pi
        dec_in = (np.random.sample(10)-0.5)*np.pi
        dl_control = None
        db_control = None
        for (ra, dec) in zip(ra_in, dec_in):
            dl, db = pal.eqgal(ra, dec)
            if dl_control is None:
                dl_control = np.array([dl])
                db_control = np.array([db])
            else:
                np.append(dl_control, dl)
                np.append(db_control, db)
        dlTest, db_test = pal.eqgalVector(ra_in, dec_in)
        for (dlt, dbt, dlc, dbc) in zip(dlTest, db_test, dl_control, db_control):
            self.assertAlmostEqual(dlt, dlc, 12)
            self.assertAlmostEqual(dbt, dbc, 12)

        # test that an exception is raised when the input
        # arrays are of different lengths
        self.assertRaises(ValueError, pal.eqgalVector, ra_in[:2], dec_in)

    def test_etrms(self):
        ev = pal.etrms(1976.9)
        self.assertAlmostEqual(ev[0], -1.621617102537041e-6, 18)
        self.assertAlmostEqual(ev[1], -3.310070088507914e-7, 18)
        self.assertAlmostEqual(ev[2], -1.435296627515719e-7, 18)

    def test_evp(self):
        vbex = np.array([1.6957348127008098514e-07,
                         -9.1093446116039685966e-08,
                         -3.9528532243991863036e-08])
        pbex = np.array([-0.49771075259730546136,
                         -0.80273812396332311359,
                         -0.34851593942866060383])
        vhex = np.array([1.6964379181455713805e-07,
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

        (dvb, dpb, dvh, dph) = pal.evp(2010.0, 2012.0)
        np.testing.assert_allclose(dvb, vbex, atol=1e-12)
        np.testing.assert_allclose(dpb, pbex, atol=1e-12)
        np.testing.assert_allclose(dvh, vhex, atol=1e-12)
        np.testing.assert_allclose(dph, phex, atol=1e-12)

        (dph, dvh, dpb, dvb) = pal.epv(53411.52501161)
        np.testing.assert_allclose(dvb, vbex2, atol=1e-12)
        np.testing.assert_allclose(dpb, pbex2, atol=1e-12)
        np.testing.assert_allclose(dvh, vhex2, atol=1e-12)
        np.testing.assert_allclose(dph, phex2, atol=1e-12)

    def test_fk524(self):
        """
        Test that fk524 gives results consistent with data published
        in:

        'Explanatory Supplement to The Astronomical Almanac'
        Seidelmann, P. Kenneth (1992)
        University Science Books

        Table 3.58.1
        """

        # fk4 ra
        hr_in = np.array([0, 3, 6, 14, 21, 1, 20, 11, 14])
        min_in = np.array([17, 17, 11, 36, 4, 48, 15, 50, 54])
        sec_in = np.array([28.774, 55.847, 43.975, 11.250, 39.935,
                           48.784, 3.004, 6.172, 59.224])

        fk4_ra = pal.dtf2rVector(hr_in, min_in, sec_in)

        # fk4 dec
        sgn = np.array([-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0])
        deg_in = np.array([65, 43, 74, 60, 38, 89, 89, 38, 0])
        amin_in = np.array([10, 15, 44, 37, 29, 1, 8, 4, 1])
        asec_in = np.array([6.70, 35.74, 12.46, 48.85, 59.10, 43.74, 18.48,
                            39.15, 58.08])

        fk4_dec = sgn*pal.daf2rVector(deg_in, amin_in, asec_in)

        # fk4 mura
        sgn = np.array([1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        hr_in = np.zeros(9, dtype=np.int)
        min_in = np.zeros(9, dtype=np.int)
        sec_in = np.array([27.141, 27.827, 3.105, 49.042,
                           35.227, 18.107, 11.702, 33.873, 0.411])

        fk4_mura = 0.01*sgn*pal.dtf2rVector(hr_in, min_in, sec_in)

        # fk4 mudec
        fk4_mudec = 0.01*pal.DAS2R*np.array([116.74, 74.76, -21.12, 71.20,
                                             318.47, -0.43, -0.09, -580.57,
                                             -2.73])

        fk4_px = np.array([0.134, 0.156, 0.115, 0.751, 0.292, 0.000,
                           0.000, 0.116, 0.000])

        fk4_vr = np.array([8.70, 86.80, 35.00, -22.20, -64.00, 0.00,
                           0.00, -98.30, 0.00])

        # fk5 ra
        hr_in = np.array([0, 3, 6, 14, 21, 2, 21, 11, 14])
        min_in = np.array([20, 19, 10, 39, 6, 31, 8, 52, 57])
        sec_in = np.array([4.3100, 55.6785, 14.5196, 36.1869,
                           54.5901, 49.8131, 46.0652, 58.7461,
                           33.2650])

        fk5_ra = pal.dtf2rVector(hr_in, min_in, sec_in)

        # fk5 dec
        sgn = np.array([-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0])
        deg_in = np.array([64, 43, 74, 60, 38, 89, 88, 37, 0])
        amin_in = np.array([52, 4, 45, 50, 44, 15, 57, 43, 10])
        asec_in = np.array([29.332, 10.830, 11.036, 7.393, 44.969,
                            50.661, 23.667, 7.456, 3.240])

        fk5_dec = sgn*pal.daf2rVector(deg_in, amin_in, asec_in)

        # fk5 mura
        sgn = np.array([1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        hr_in = np.zeros(9, dtype=np.int)
        min_in = np.zeros(9, dtype=np.int)
        sec_in = np.array([26.8649, 27.7694, 3.1310, 49.5060,
                           35.3528, 21.7272, 8.4469, 33.7156,
                           0.4273])

        fk5_mura = 0.01*sgn*pal.dtf2rVector(hr_in, min_in, sec_in)

        # fk5 mudec
        fk5_mudec = 0.01*pal.DAS2R*np.array([116.285, 73.050, -21.304,
                                             69.934, 320.206, -1.571,
                                             0.171, -581.216,
                                             -2.402])

        fk5_px = np.array([0.1340, 0.1559, 0.1150, 0.7516,
                           0.2923, 0.0000, 0.0000,
                           0.1161, 0.0000])

        fk5_vr = np.array([8.74, 86.87, 35.00, -22.18, -63.89, 0.00,
                           0.00, -97.81, 0.00])

        for ra5, dec5, mura5, mudec5, px5, vr5, \
            ra4, dec4, mura4, mudec4, px4, vr4 in \
            zip(fk5_ra, fk5_dec, fk5_mura, fk5_mudec, fk5_px, fk5_vr,
                fk4_ra, fk4_dec, fk4_mura, fk4_mudec, fk4_px, fk4_vr):

            ra, dec, mura, \
                mudec, px, vr = pal.fk524(ra5, dec5, mura5, mudec5, px5, vr5)

            # dpos is the angular separation (in arcsec) between
            # the result of pal.fk524 and the true fk4 coordinates
            # of the sample point
            dpos = pal.DR2AS*pal.dsep(ra4, dec4, ra, dec)

            dmura = pal.DR2AS*np.abs(mura-mura4)
            dmudec = pal.DR2AS*np.abs(mudec-mudec4)
            dpx = np.abs(px-px4)
            dvr = np.abs(vr-vr4)

            self.assertLess(dpos, 0.001)
            self.assertLess(dmura, 0.01)
            self.assertLess(dmudec, 0.01)
            self.assertLess(dpx, 0.001)
            self.assertLess(dvr, 0.01)

    def test_fk524_vectors(self):
        """
        Test that fk524Vector returns results consistent with
        fk524
        """
        np.random.seed(135)
        n_samples = 200
        ra5 = np.random.random_sample(n_samples)*2.0*np.pi
        dec5 = (np.random.random_sample(n_samples)-0.5)*np.pi
        mura5 = 0.01*(np.random.random_sample(n_samples)-0.5)*200.0*np.radians(1.0/3600.0)
        mudec5 = 0.01*(np.random.random_sample(n_samples)-0.5)*200.0*np.radians(1.0/3600.0)
        px5 = np.random.random_sample(n_samples)
        vr5 = (np.random.random_sample(n_samples)-0.5)*200.0

        testRa, testDec, \
            testMura, testMudec, \
            testPx, testVr = pal.fk524Vector(ra5, dec5, mura5, mudec5, px5, vr5)

        for ii in range(n_samples):
            controlRa, controlDec, \
                controlMura, controlMudec, \
                controlPx, controlVr = pal.fk524(ra5[ii], dec5[ii], mura5[ii],
                                                 mudec5[ii], px5[ii], vr5[ii])

            self.assertEqual(controlRa, testRa[ii])
            self.assertEqual(controlDec, testDec[ii])
            self.assertEqual(controlMura, testMura[ii])
            self.assertEqual(controlMudec, testMudec[ii])
            self.assertEqual(controlPx, testPx[ii])
            self.assertEqual(controlVr, testVr[ii])

        # test that exceptions are raised when the input arrays are
        # of different lengths
        self.assertRaises(ValueError, pal.fk524Vector, ra5, dec5[:6], mura5,
                          mudec5, px5, vr5)

        self.assertRaises(ValueError, pal.fk524Vector, ra5, dec5, mura5[:8],
                          mudec5, px5, vr5)

        self.assertRaises(ValueError, pal.fk524Vector, ra5, dec5, mura5,
                          mudec5[:6], px5, vr5)

        self.assertRaises(ValueError, pal.fk524Vector, ra5, dec5, mura5,
                          mudec5, px5[:9], vr5)

        self.assertRaises(ValueError, pal.fk524Vector, ra5, dec5, mura5,
                          mudec5, px5, vr5[:9])

    def test_fk524_original(self):
        r1950, d1950, dr1950, dd1950, p1950, v1950 = pal.fk524(4.567, -1.23, -3e-5, 8e-6, 0.29, -35.0)
        self.assertAlmostEqual(r1950, 4.543778603272084, 12)
        self.assertAlmostEqual(d1950, -1.229642790187574, 12)
        self.assertAlmostEqual(dr1950, -2.957873121769244e-5, 17)
        self.assertAlmostEqual(dd1950, 8.117725309659079e-6, 17)
        self.assertAlmostEqual(p1950, 0.2898494999992917, 12)
        self.assertAlmostEqual(v1950, -35.026862824252680, 11)

    def test_fk45z(self):
        (r2000, d2000) = pal.fk45z(1.2, -0.3, 1960)
        self.assertAlmostEqual(r2000, 1.2097812228966762227, 11)
        self.assertAlmostEqual(d2000, -0.29826111711331398935, 12)

    def test_fk45z_vector(self):
        """
        Test that fk45zVector produces results consistent with fk45z
        """
        epoch = 1960
        np.random.seed(136)
        n_samples = 200
        r1950 = np.random.random_sample(n_samples)*2.0*np.pi
        d1950 = (np.random.random_sample(n_samples)-0.5)*np.pi

        testR2000, testD2000 = pal.fk45zVector(r1950, d1950, epoch)

        for ii in range(n_samples):
            controlR2000, controlD2000 = pal.fk45z(r1950[ii], d1950[ii], epoch)
            self.assertEqual(controlR2000, testR2000[ii])
            self.assertEqual(controlD2000, testD2000[ii])

        # test that fk45zVector and fk54zVector invert each other
        testR1950, testD1950, \
            testDR1950, testDD1950 = pal.fk54zVector(testR2000, testD2000, epoch)

        distance = pal.dsepVector(testR1950, testD1950, r1950, d1950)
        np.testing.assert_array_almost_equal(distance*pal.DR2ASs,
                                             np.zeros(len(distance)), 4)

        # test that an exception is raised if the input arrays are of different
        # lengths
        self.assertRaises(ValueError, pal.fk45zVector, r1950, d1950[:8], epoch)

    def test_fk52h(self):
        inr5 = 1.234
        ind5 = -0.987
        epoch = 1980
        (rh, dh) = pal.fk5hz(inr5, ind5, epoch)
        self.assertAlmostEqual(rh, 1.234000136713611301, 13)
        self.assertAlmostEqual(dh, -0.9869999702020807601, 13)
        (r5, d5, dr5, dd5) = pal.hfk5z(rh, dh, epoch)
        self.assertAlmostEqual(r5, inr5, 13)
        self.assertAlmostEqual(d5, ind5, 13)
        self.assertAlmostEqual(dr5, 0.000000006822074, 13)
        self.assertAlmostEqual(dd5, -0.000000002334012, 13)

    def test_fk5hz_vector(self):
        """
        Test that fk5hzVector returns results consistent with
        fk5hz
        """
        np.random.seed(132)
        n_samples = 200
        ra_list = np.random.random_sample(n_samples)*2.0*np.pi
        dec_list = (np.random.random_sample(n_samples)-0.5)*np.pi
        testRa, testDec = pal.fk5hzVector(ra_list, dec_list, 2000.0)

        for ii in range(n_samples):
            controlRa, controlDec = pal.fk5hz(ra_list[ii], dec_list[ii], 2000.0)
            self.assertEqual(controlRa, testRa[ii])
            self.assertEqual(controlDec, testDec[ii])

        # test that an exception is raised if the input ra_list and dec_list
        # are of different sizes
        self.assertRaises(ValueError, pal.fk5hzVector, ra_list, dec_list[:24],
                          2000.0)

    def test_hkf5z_vector(self):
        """
        Test that hkf5zVector produces results consistent with
        hkf5z
        """
        np.random.seed(133)
        n_samples = 200
        ra_list = np.random.random_sample(n_samples)*2.0*np.pi
        dec_list = (np.random.random_sample(n_samples)-0.5)*np.pi
        testRa, testDec, testDr, testDd = pal.hfk5zVector(ra_list, dec_list, 2000.0)
        for ii in range(n_samples):
            controlRa, controlDec, \
                controlDr, controlDd = pal.hfk5z(ra_list[ii], dec_list[ii], 2000.0)
            self.assertEqual(testRa[ii], controlRa)
            self.assertEqual(testDec[ii], controlDec)
            self.assertEqual(testDr[ii], controlDr)
            self.assertEqual(testDd[ii], controlDd)

        # test hfk5zVector anf fk5hzVector invert each other
        ra5, dec5 = pal.fk5hzVector(testRa, testDec, 2000.0)
        np.testing.assert_array_almost_equal(ra5, ra_list, 12)
        np.testing.assert_array_almost_equal(dec5, dec_list, 12)

        # test that an exception is raised if ra_list and dec_list are
        # of different lengths
        self.assertRaises(ValueError, pal.hfk5zVector, ra_list, dec_list[:77],
                          2000.0)

    def test_fk54z(self):
        (r1950, d1950, dr1950, dd1950) = pal.fk54z(1.2, -0.3, 1960)
        self.assertAlmostEqual(r1950, 1.1902221805755279771, 12)
        self.assertAlmostEqual(d1950, -0.30178317645793828472, 12)
        self.assertAlmostEqual(dr1950, -1.7830874775952945507e-08, 12)
        self.assertAlmostEqual(dd1950, 7.196059425334821089e-09, 12)

    def test_fk54z_vector(self):
        """
        Test that fk54zVector returns results consistent with fk54z
        """
        epoch = 1960
        np.random.seed(136)
        n_samples = 200
        r2000 = np.random.random_sample(n_samples)*2.0*np.pi
        d2000 = (np.random.random_sample(n_samples)-0.5)*np.pi
        testR1950, testD1950, \
            testDr1950, testDd1950 = pal.fk54zVector(r2000, d2000, epoch)

        for ii in range(n_samples):
            controlR1950, controlD1950, \
                controlDr1950, controlDd1950 = pal.fk54z(r2000[ii], d2000[ii], epoch)

            self.assertEqual(controlR1950, testR1950[ii])
            self.assertEqual(controlD1950, testD1950[ii])
            self.assertEqual(controlDr1950, testDr1950[ii])
            self.assertEqual(controlDd1950, testDd1950[ii])

        # test that an exception is raised if input arrays are of
        # different lengths
        self.assertRaises(ValueError, pal.fk54zVector, r2000[:11], d2000, epoch)

    def test_flotin(self):
        pass

    def test_galeq(self):
        (dr, dd) = pal.galeq(5.67, -1.23)
        self.assertAlmostEqual(dr, 0.04729270418071426, 12)
        self.assertAlmostEqual(dd, -0.7834003666745548, 12)

    def test_galeq_vector(self):
        np.random.seed(32)
        dl_in = np.random.sample(10)*2.0*np.pi
        db_in = (np.random.sample(10)-0.5)*np.pi
        dr_control = np.zeros(10, dtype=np.float64)
        dd_control = np.zeros(10, dtype=np.float64)
        for (i, cc) in enumerate(zip(dl_in, db_in)):
            dr, dd = pal.galeq(cc[0], cc[1])
            dr_control[i] = dr
            dd_control[i] = dd
        drTest, dd_test = pal.galeqVector(dl_in, db_in)
        for (drt, ddt, drc, ddc) in zip(drTest, dd_test, dr_control, dd_control):
            self.assertAlmostEqual(drt, drc, 12)
            self.assertAlmostEqual(ddt, ddc, 12)

        # test that an exception is raise if input arrays
        # are of different lengths
        self.assertRaises(ValueError, pal.galeqVector, dl_in[:3], db_in)

    def test_galsup(self):
        (dsl, dsb) = pal.galsup(6.1, -1.4)
        self.assertAlmostEqual(dsl, 4.567933268859171, 12)
        self.assertAlmostEqual(dsb, -0.01862369899731829, 12)

    def test_galsup_vector(self):
        """
        Test that galsupVector gives results consistent with galsup
        """
        np.random.seed(134)
        n_samples = 200

        ll_list = np.random.random_sample(n_samples)*2.0*np.pi
        bb_list = (np.random.random_sample(n_samples)-0.5)*np.pi

        testSl, testSb = pal.galsupVector(ll_list, bb_list)

        for ii in range(n_samples):
            controlSl, controlSb = pal.galsup(ll_list[ii], bb_list[ii])
            self.assertEqual(controlSl, testSl[ii])
            self.assertEqual(controlSb, testSb[ii])

        # test that an exception is raised if you pass in arrays with
        # different lengths
        self.assertRaises(ValueError, pal.galsupVector, ll_list[:55], bb_list)

    def test_geoc(self):
        lat = 19.822838905884 * pal.DD2R
        alt = 4120.0
        (r, z) = pal.geoc(lat, alt)
        self.assertAlmostEqual(r, 4.01502667039618e-05, 10)
        self.assertAlmostEqual(z, 1.43762411970295e-05, 10)

    def test_ge50(self):
        (dr, dd) = pal.ge50(6.1, -1.55)
        self.assertAlmostEqual(dr, 0.1966825219934508, 12)
        self.assertAlmostEqual(dd, -0.4924752701678960, 12)

    def test_ge50_vector(self):
        """
        Test that ge50Vector returns results consistent with ge50
        """
        np.random.seed(133)
        n_samples = 200
        ll_list = np.random.random_sample(n_samples)*2.0*np.pi
        bb_list = (np.random.random_sample(n_samples)-0.5)*np.pi
        testRa, testDec = pal.ge50Vector(ll_list, bb_list)
        for ii in range(n_samples):
            controlRa, controlDec = pal.ge50(ll_list[ii], bb_list[ii])
            self.assertEqual(controlRa, testRa[ii])
            self.assertEqual(controlDec, testDec[ii])

        # test that an exception is raised when you pass in different
        # numbers of longitudes as you do latitudes
        self.assertRaises(ValueError, pal.ge50Vector, ll_list, bb_list[:45])

    def test_gmst(self):
        self.assertAlmostEqual(pal.gmst(53736), 1.754174971870091203, 12)
        self.assertAlmostEqual(pal.gmsta(53736, 0), 1.754174971870091203, 12)

    def test_gmst_vector(self):
        np.random.seed(32)
        date_in = 53000.0 + np.random.sample(20)*5000.0
        gm_control = np.zeros(20, dtype=np.float64)
        for i, d in enumerate(date_in):
            gm = pal.gmst(d)
            gm_control[i] = gm
        gm_test = pal.gmstVector(date_in)
        for gt, gc in zip(gm_test, gm_control):
            self.assertAlmostEqual(gt, gc, 12)

    def test_gmsta_vector(self):
        np.random.seed(32)
        date_in = 53000 + np.random.randint(0, 1001, 20)
        ut_in = np.random.sample(20)
        gm_control = np.zeros(20, dtype=np.float64)
        for i, d in enumerate(zip(date_in, ut_in)):
            dd = pal.gmsta(d[0], d[1])
            gm_control[i] = dd
        gm_test = pal.gmstaVector(date_in, ut_in)
        for gt, gc in zip(gm_test, gm_control):
            self.assertAlmostEqual(gt, gc, 12)

        # test that an exception is raised if input arrays are
        # of different lengths
        self.assertRaises(ValueError, pal.gmstaVector, date_in[:8], ut_in)

    def test_intin(self):
        s = "  -12345, , -0  2000  +     "
        i = 1
        (n, i, sign) = pal.intin(s, i)
        self.assertEqual(i, 10)
        self.assertEqual(n, -12345)
        self.assertLess(sign, 0)

        (n, i, sign) = pal.intin(s, i)
        self.assertEqual(i, 12)
        self.assertIsNone(n)

        (n, i, sign) = pal.intin(s, i)
        self.assertEqual(i, 17)
        self.assertEqual(n, 0)
        self.assertLess(sign, 0)

        (n, i, sign) = pal.intin(s, i)
        self.assertEqual(i, 23)
        self.assertEqual(n, 2000)
        self.assertGreater(sign, 0)

        (n, i, sign) = pal.intin(s, i)
        self.assertEqual(i, 29)
        self.assertIsNone(n)
        self.assertIsNone(sign)

    def test_map(self):
        (ra, da) = pal.map(6.123, -0.999, 1.23e-5, -0.987e-5,
                           0.123, 32.1, 1999, 43210.9)
        self.assertAlmostEqual(ra, 6.117130429775647, 6)
        self.assertAlmostEqual(da, -1.000880769038632, 7)

    def test_mappa(self):
        expected = np.array([1.9986310746064646082,
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
                            0.99999997416198904698])
        amprms = pal.mappa(2010.0, 55927)
        np.testing.assert_array_almost_equal(amprms, expected, decimal=12)

    def test_mapqk(self):
        """
        Test mapqk by taking the geocentric apparent positions of Arcturus
        as downloaded from aa.usno.mil/data/docs/geocentric.php and trying
        to calculate it from Arcturus' mean position, proper motion, parallax,
        and radial velocity.

        ra_0 and dec_0 are the mean position; ra_app and dec_app are the
        geocentric apparent position.

        The data below represents the position of Arcturus on
        JD (UT) 2457000.375 as reported by
        http://aa.usno.navy.mil/data/docs/geocentric.php
        """

        ra_0, i = pal.dafin("14 15 39.67207 ", 1)
        ra_0 *= 15.0
        dec_0, i = pal.dafin("19 10 56.673", 1)

        pm_ra = -1.0939*pal.DAS2R
        pm_ra /= np.cos(dec_0)
        pm_dec = -2.00006*pal.DAS2R
        v_rad = -5.19
        px = 0.08883*pal.DAS2R

        # time is the TDB MJD calculated from a JD of 2457000.375 with astropy.time
        amprms = pal.mappa(2000.0, 56999.87537249177)

        ra_test, dec_test = pal.mapqk(ra_0, dec_0, pm_ra, pm_dec, px, v_rad, amprms)

        ra_app, i = pal.dafin("14 16 19.59", 1)
        ra_app *= 15.0
        dec_app, i = pal.dafin("19 6 19.56", 1)

        # find the angular distance from the known mean position
        # to the calculated mean postion
        dd = pal.dsep(ra_test, dec_test, ra_app, dec_app)
        dd *= pal.DR2AS
        self.assertAlmostEqual(dd, 0.0, 0)

    def test_mapqkz(self):
        """
        Run inputs through mapqk with zero proper motion, parallax
        and radial velocity.  Then run the same inputs through mapqkz.
        Verify that the results are the same.
        """
        amprms = pal.mappa(2010, 55927)
        (ra_c, da_c) = pal.mapqk(1.234, -0.567, 0.0, 0.0, 0.0, 0.0, amprms)
        (ra, da) = pal.mapqkz(1.234, -0.567, amprms)

        self.assertAlmostEqual(ra, ra_c, 12)
        self.assertAlmostEqual(da, da_c, 12)

    def test_mapqkz_vector(self):
        amprms = pal.mappa(2010, 55927)
        np.random.seed(32)
        n_tests = 100
        ra_in = np.random.sample(n_tests)*2.0*np.pi
        dec_in = (np.random.sample(n_tests)-0.5)*np.pi
        r_control = None
        d_control = None
        for (rr, dd) in zip(ra_in, dec_in):
            r, d = pal.mapqkz(rr, dd, amprms)
            if r_control is None:
                r_control = np.array([r])
                d_control = np.array([d])
            else:
                r_control = np.append(r_control, r)
                d_control = np.append(d_control, d)

        rTest, dTest = pal.mapqkzVector(ra_in, dec_in, amprms)
        for (r1, d1, r2, d2) in zip(r_control, d_control, rTest, dTest):
            self.assertAlmostEqual(r1, r2, 12)
            self.assertAlmostEqual(d1, d2, 12)

        # test that an exception is raised if the input arrays
        # are of inconsistent lengths
        self.assertRaises(ValueError, pal.mapqkzVector, ra_in[:4], dec_in, amprms)

    def test_mapqk_vector(self):
        amprms = pal.mappa(2010, 55927)
        np.random.seed(32)
        n_tests = 100
        ra_in = np.random.sample(n_tests)*2.0*np.pi
        dec_in = (np.random.sample(n_tests)-0.5)*np.pi
        pmr = (np.random.sample(n_tests)-0.5)*0.01
        pmd = (np.random.sample(n_tests)-0.5)*0.01
        px = 0.00045+np.random.sample(n_tests)*0.001
        rv = 200.0*np.random.sample(n_tests)
        r_control = None
        d_control = None
        for (rr, dd, pr, pd, x, v) in zip(ra_in, dec_in, pmr, pmd, px, rv):
            r, d = pal.mapqk(rr, dd, pr, pd, x, v, amprms)
            if r_control is None:
                r_control = np.array([r])
                d_control = np.array([d])
            else:
                r_control = np.append(r_control, r)
                d_control = np.append(d_control, d)

        rTest, dTest = pal.mapqkVector(ra_in, dec_in, pmr, pmd, px, rv, amprms)
        for (r1, d1, r2, d2) in zip(r_control, d_control, rTest, dTest):
            self.assertAlmostEqual(r1, r2, 12)
            self.assertAlmostEqual(d1, d2, 12)

        # test that an exception is raised if the input arrays
        # are of inconsistent shapes
        self.assertRaises(ValueError, pal.mapqkVector, ra_in, dec_in[:10], pmr,
                          pmd, px, rv, amprms)

        self.assertRaises(ValueError, pal.mapqkVector, ra_in, dec_in, pmr[:10],
                          pmd, px, rv, amprms)

        self.assertRaises(ValueError, pal.mapqkVector, ra_in, dec_in, pmr,
                          pmd[:10], px, rv, amprms)

        self.assertRaises(ValueError, pal.mapqkVector, ra_in, dec_in, pmr,
                          pmd, px[:10], rv, amprms)

        self.assertRaises(ValueError, pal.mapqkVector, ra_in, dec_in, pmr,
                          pmd, px, rv[:10], amprms)

    def test_moon(self):
        expected = np.array([
            0.00229161514616454,
            0.000973912029208393,
            0.000669931538978146,
            -3.44709700068209e-09,
            5.44477533462392e-09,
            2.11785724844417e-09
        ])
        pv = pal.dmoon(48634.4687174074)
        np.testing.assert_array_almost_equal(pv, expected, decimal=12)

    def test_dmoon_vector(self):
        """
        Test that dmoonVector returns results consistent with dmoon
        """
        np.random.seed(141)
        n_samples = 1000
        date = np.random.random_sample(n_samples)*10000.0 + 43000.0
        testX, testY, testZ, testXd, testYd, testZd = pal.dmoonVector(date)

        for ii in range(n_samples):
            x, y, z, xd, yd, zd = pal.dmoon(date[ii])
            self.assertEqual(x, testX[ii])
            self.assertEqual(y, testY[ii])
            self.assertEqual(z, testZ[ii])
            self.assertEqual(xd, testXd[ii])
            self.assertEqual(yd, testYd[ii])
            self.assertEqual(zd, testZd[ii])

    def test_nut(self):
        expected = np.array([
            [9.999999969492166e-1, 7.166577986249302e-5, 3.107382973077677e-5],
            [-7.166503970900504e-5, 9.999999971483732e-1, -2.381965032461830e-5],
            [-3.107553669598237e-5, 2.381742334472628e-5, 9.999999992335206818e-1]
        ])

        rmatn = pal.nut(46012.32)
        np.testing.assert_array_almost_equal(rmatn, expected, decimal=3)

        (dpsi, deps, eps0) = pal.nutc(54388.0)
        self.assertAlmostEqual(eps0, 0.4090749229387258204, 14)

        (dpsi, deps, eps0) = pal.nutc(53736.0)
        self.assertAlmostEqual(dpsi, -0.9630912025820308797e-5, 13)
        self.assertAlmostEqual(deps, 0.4063238496887249798e-4, 13)

    def test_obs(self):
        obsdata = pal.obs()
        mmt = obsdata["MMT"]
        self.assertEqual(mmt["name"], "MMT 6.5m, Mt Hopkins")
        self.assertAlmostEqual(mmt["long"], 1.935300584055477, 8)
        self.assertAlmostEqual(mmt["lat"], 0.5530735081550342238, 10)
        self.assertAlmostEqual(mmt["height"], 2608, 10)

        self.assertEqual(len(obsdata), 85)

    def test_pa(self):
        self.assertAlmostEqual(pal.pa(-1.567, 1.5123, 0.987),
                               -1.486288540423851, 12)
        self.assertAlmostEqual(pal.pa(0, 0.789, 0.789), 0, 12)

    def test_pa_vector(self):
        np.random.seed(32)
        ha_in = np.random.sample(20)*2.0*np.pi
        dec_in = (np.random.sample(20)-0.5)*np.pi
        phi = 0.3
        pa_control = np.zeros(20, dtype=np.float64)
        for i in range(20):
            pa_control[i] = pal.pa(ha_in[i], dec_in[i], phi)
        pa_test = pal.paVector(ha_in, dec_in, phi)
        for pt, pc in zip(pa_test, pa_control):
            self.assertAlmostEqual(pt, pc, 12)

        # test that an exception is raised if input arrays
        # are of inconsistent length
        self.assertRaises(ValueError, pal.paVector, ha_in[:4], dec_in, phi)

    def test_pcd(self):
        disco = 178.585
        refx = 0.0123
        refy = -0.00987
        (x, y) = pal.pcd(disco, refx, refy)
        self.assertAlmostEqual(x, 0.01284630845735895, 14)
        self.assertAlmostEqual(y, -0.01030837922553926, 14)

        (x, y) = pal.unpcd(disco, x, y)
        self.assertAlmostEqual(x, refx, 14)
        self.assertAlmostEqual(y, refy, 14)

        # Round trip
        (x, y) = pal.pcd(-disco, refx, refy)
        (x, y) = pal.unpcd(-disco, x, y)
        self.assertAlmostEqual(x, refx, 14)
        self.assertAlmostEqual(y, refy, 14)

    def test_pcd_vector(self):
        """
        Test that pcdVector returns the same results as running
        pcd on each element of the vectors
        """
        np.random.seed(120)
        x_in = 2.0*(np.random.random_sample(100)-0.5)
        y_in = 2.0*(np.random.random_sample(100)-0.5)
        disco = 191.0

        xTestList, yTestList = pal.pcdVector(disco, x_in, y_in)
        for xTest, yTest, xx, yy in zip(xTestList, yTestList, x_in, y_in):
            xControl, yControl = pal.pcd(disco, xx, yy)
            self.assertEqual(xControl, xTest)
            self.assertEqual(yControl, yTest)

        # make sure that x and y were changed by applying pcd
        # (also make sure that x_in and y_in were not changed
        # by pcdVector)
        with self.assertRaises(AssertionError):
            np.testing.assert_array_almost_equal(xTestList, x_in, 12)

        with self.assertRaises(AssertionError):
            np.testing.assert_array_almost_equal(yTestList, y_in, 12)

        # make sure that an exceptoion is raised if you pass in
        # mismatched x and y vectors
        self.assertRaises(ValueError, pal.pcdVector, disco, x_in, y_in[:10])

    def test_unpcd_vector(self):
        """
        Test that unpcdVector produces the same results as running unpcd
        on each element of the input vectors
        """
        np.random.seed(121)
        disco = 132.0
        x_in = 2.0*(np.random.random_sample(120)-0.5)
        y_in = 2.0*(np.random.random_sample(120)-0.5)

        xTestList, yTestList = pal.unpcdVector(disco, x_in, y_in)

        # make sure that unpcdVector and unpcd give the same results
        for xTest, yTest, xx, yy in zip(xTestList, yTestList, x_in, y_in):
            xControl, yControl = pal.unpcd(disco, xx, yy)
            self.assertEqual(xTest, xControl)
            self.assertEqual(yTest, yControl)

        # make sure that xTestList and yTestList are different from x_in and y_in
        with self.assertRaises(AssertionError):
            np.testing.assert_array_almost_equal(xTestList, x_in, 12)

        with self.assertRaises(AssertionError):
            np.testing.assert_array_almost_equal(yTestList, y_in, 12)

        # use pcdVector to transform back into the distorted coordinates and
        # verify that those are equivalent to x_in and y_in
        xDistorted, yDistorted = pal.pcdVector(disco, xTestList, yTestList)
        np.testing.assert_array_almost_equal(xDistorted, x_in, 12)
        np.testing.assert_array_almost_equal(yDistorted, y_in, 12)

        # test that an exception is raised if you pass in different numbers
        # of x and y coordinates
        self.assertRaises(ValueError, pal.unpcdVector, disco, x_in, y_in[:30])

    def test_planet(self):
        # palEl2ue
        u = pal.el2ue(50000, 1, 49000, 0.1, 2, 0.2,
                      3, 0.05, 3, 0.003312)
        expectedue1 = np.array([1.000878908362435284, -0.3336263027874777288, 50000.,
                               2.840425801310305210, 0.1264380368035014224, -0.2287711835229143197,
                               -0.01301062595106185195, 0.5657102158104651697, 0.2189745287281794885,
                               2.852427310959998500, -0.01552349065435120900,
                               50000., 0.0])
        np.testing.assert_allclose(u, expectedue1, atol=1e-12)

        # palPertel
        (epoch, orbinc, anode, perih, aorq, e, aorl) = pal.pertel(2, 43000., 43200., 43000.,
                                                                  0.2, 3, 4, 5, 0.02, 6)
        self.assertAlmostEqual(epoch, 43200, 10)
        self.assertAlmostEqual(orbinc, 0.1995661466545422381, 7)
        self.assertAlmostEqual(anode, 2.998052737821591215, 7)
        self.assertAlmostEqual(perih, 4.009516448441143636, 6)
        self.assertAlmostEqual(aorq, 5.014216294790922323, 7)
        self.assertAlmostEqual(e, 0.02281386258309823607, 7)
        self.assertAlmostEqual(aorl, 0.01735248648779583748, 6)

        # palPertue
        unew = pal.pertue(50100, u)
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
            50100., 0.])
        np.testing.assert_allclose(unew, expectedue3, atol=1e-12)

        # palPlanel
        pv = pal.planel(50600, 2, 50500, 0.1, 3, 5, 2, 0.3, 4, 0)
        expectedpv2 = np.array([1.947628959288897677,
                               -1.013736058752235271,
                               -0.3536409947732733647,
                               2.742247411571786194e-8,
                               1.170467244079075911e-7,
                               3.709878268217564005e-8])
        np.testing.assert_allclose(pv, expectedpv2, atol=1e-12)

        # palPlanet
        self.assertRaises(ValueError, pal.planet, 1e6, 0)
        self.assertRaises(ValueError, pal.planet, 1e6, 9)

        pv = pal.planet(-320000, 3)
        self.assertAlmostEqual(pv[0], 0.9308038666827242603, 11)
        self.assertAlmostEqual(pv[1], 0.3258319040252137618, 11)
        self.assertAlmostEqual(pv[2], 0.1422794544477122021, 11)
        self.assertAlmostEqual(pv[3], -7.441503423889371696e-8, 17)
        self.assertAlmostEqual(pv[4], 1.699734557528650689e-7, 17)
        self.assertAlmostEqual(pv[5], 7.415505123001430864e-8, 17)

        pv = pal.planet(43999.9, 1)
        self.assertAlmostEqual(pv[0], 0.2945293959257422246, 11)
        self.assertAlmostEqual(pv[1], -0.2452204176601052181, 11)
        self.assertAlmostEqual(pv[2], -0.1615427700571978643, 11)
        self.assertAlmostEqual(pv[3], 1.636421147459047057e-7, 18)
        self.assertAlmostEqual(pv[4], 2.252949422574889753e-7, 18)
        self.assertAlmostEqual(pv[5], 1.033542799062371839e-7, 18)

        # palPlante
        (ra, dec, r) = pal.plante(50600., -1.23, 0.456, 2, 50500.,
                                  0.1, 3., 5., 2., 0.3, 4., 0.0)
        self.assertAlmostEqual(ra, 6.222958101333794007, 6)
        self.assertAlmostEqual(dec, 0.01142220305739771601, 6)
        self.assertAlmostEqual(r, 2.288902494080167624, 8)

        # palPlantu
        u = np.array([1.0005, -0.3, 55000., 2.8, 0.1, -0.2,
                     -0.01, 0.5, 0.22, 2.8, -0.015, 55001., 0.0])
        (ra, dec, r) = pal.plantu(55001., -1.23, 0.456, u)
        self.assertAlmostEqual(ra, 0.3531814831241686647, 6)
        self.assertAlmostEqual(dec, 0.06940344580567131328, 6)
        self.assertAlmostEqual(r, 3.031687170873274464, 8)

        # palPv2el
        pv = np.array([0.3, -0.2, 0.1, -0.9e-7, 0.8e-7, -0.7e-7])
        (jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm) = pal.pv2el(pv, 50000, 0.00006, 1)
        self.assertEqual(jform, 1)
        self.assertAlmostEqual(epoch, 50000, 10)
        self.assertAlmostEqual(orbinc, 1.52099895268912, 12)
        self.assertAlmostEqual(anode, 2.720503180538650, 12)
        self.assertAlmostEqual(perih, 2.194081512031836, 12)
        self.assertAlmostEqual(aorq, 0.2059371035373771, 12)
        self.assertAlmostEqual(e, 0.9866822985810528, 12)
        self.assertAlmostEqual(aorl, 0.2012758344836794, 12)
        self.assertAlmostEqual(dm, 0.184074050795182, 12)

        # palPv2ue
        expectedue2 = np.array([1.00006, -4.856142884511782, 50000., 0.3, -0.2,
                                0.1, -0.4520378601821727, 0.4018114312730424,
                                -.3515850023639121, 0.3741657386773941,
                                -0.2511321445456515, 50000., 0.])
        u = pal.pv2ue(pv, 50000., 0.00006)
        np.testing.assert_allclose(u, expectedue2, atol=1e-12)

        # Planets
        (ra, dec, diam) = pal.rdplan(40999.9, 0, 0.1, -0.9)
        self.assertAlmostEqual(ra, 5.772270359389275837, 6)
        self.assertAlmostEqual(dec, -0.2089207338795416192, 7)
        self.assertAlmostEqual(diam, 9.415338935229717875e-3, 10)

        (ra, dec, diam) = pal.rdplan(41999.9, 1, 1.1, -0.9)
        self.assertAlmostEqual(ra, 3.866363420052936653, 6)
        self.assertAlmostEqual(dec, -0.2594430577550113130, 7)
        self.assertAlmostEqual(diam, 4.638468996795023071e-5, 14)

        (ra, dec, diam) = pal.rdplan(42999.9, 2, 2.1, 0.9)
        self.assertAlmostEqual(ra, 2.695383203184077378, 6)
        self.assertAlmostEqual(dec, 0.2124044506294805126, 7)
        self.assertAlmostEqual(diam, 4.892222838681000389e-5, 14)

        (ra, dec, diam) = pal.rdplan(43999.9, 3, 3.1, 0.9)
        self.assertAlmostEqual(ra, 2.908326678461540165, 6)
        self.assertAlmostEqual(dec, 0.08729783126905579385, 7)
        self.assertAlmostEqual(diam, 8.581305866034962476e-3, 7)

        (ra, dec, diam) = pal.rdplan(44999.9, 4, -0.1, 1.1)
        self.assertAlmostEqual(ra, 3.429840787472851721, 6)
        self.assertAlmostEqual(dec, -0.06979851055261161013, 7)
        self.assertAlmostEqual(diam, 4.540536678439300199e-5, 14)

        (ra, dec, diam) = pal.rdplan(45999.9, 5, -1.1, 0.1)
        self.assertAlmostEqual(ra, 4.864669466449422548, 6)
        self.assertAlmostEqual(dec, -0.4077714497908953354, 7)
        self.assertAlmostEqual(diam, 1.727945579027815576e-4, 14)

        (ra, dec, diam) = pal.rdplan(46999.9, 6, -2.1, -0.1)
        self.assertAlmostEqual(ra, 4.432929829176388766, 6)
        self.assertAlmostEqual(dec, -0.3682820877854730530, 7)
        self.assertAlmostEqual(diam, 8.670829016099083311e-5, 14)

        (ra, dec, diam) = pal.rdplan(47999.9, 7, -3.1, -1.1)
        self.assertAlmostEqual(ra, 4.894972492286818487, 6)
        self.assertAlmostEqual(dec, -0.4084068901053653125, 7)
        self.assertAlmostEqual(diam, 1.793916783975974163e-5, 14)

        (ra, dec, diam) = pal.rdplan(48999.9, 8, 0, 0)
        self.assertAlmostEqual(ra, 5.066050284760144000, 6)
        self.assertAlmostEqual(dec, -0.3744690779683850609, 7)
        self.assertAlmostEqual(diam, 1.062210086082700563e-5, 14)

        # palUe2el
        (jform, epoch, orbinc, anode, perih, aorq, e, aorl, dm) = pal.ue2el(u, 1)
        self.assertEqual(jform, 1)
        self.assertAlmostEqual(epoch, 50000.00, 10)
        self.assertAlmostEqual(orbinc, 1.520998952689120, 12)
        self.assertAlmostEqual(anode, 2.720503180538650, 12)
        self.assertAlmostEqual(perih, 2.194081512031836, 12)
        self.assertAlmostEqual(aorq, 0.2059371035373771, 12)
        self.assertAlmostEqual(e, 0.9866822985810528, 12)
        self.assertAlmostEqual(aorl, 0.2012758344836794, 12)

        # palUe2pv
        (u2, pv) = pal.ue2pv(50010., u)

        # Update final two elements of the test UE array
        expectedue2[11] = 50010.
        expectedue2[12] = 0.7194308220038886856
        np.testing.assert_allclose(u2, expectedue2, atol=1e-12)

        expectedpv = np.array([
            0.07944764084631667011, -0.04118141077419014775,
            0.002915180702063625400, -0.6890132370721108608e-6,
            0.4326690733487621457e-6, -0.1763249096254134306e-6])
        np.testing.assert_allclose(pv, expectedpv, atol=1e-12)

    def test_pm(self):
        (ra, dec) = pal.pm(5.43, -0.87, -0.33e-5, 0.77e-5, 0.7,
                           50.3*365.2422/365.25, 1899, 1943)
        self.assertAlmostEqual(ra, 5.429855087793875, 10)
        self.assertAlmostEqual(dec, -0.8696617307805072, 10)

        (ra, dec) = pal.pm(0.01686756, -1.093989828, -1.78323516e-5,
                           2.336024047e-6, 0.74723, -21.6,
                           pal.epj(50083.0), pal.epj(53736.0))
        self.assertAlmostEqual(ra, 0.01668919069414242368, 13)
        self.assertAlmostEqual(dec, -1.093966454217127879, 13)

    def test_pm_vector(self):
        ep0 = 52000.0
        ep1 = 53510.0
        np.random.seed(32)
        n_tests = 100
        ra_in = np.random.sample(n_tests)*2.0*np.pi
        dec_in = (np.random.sample(n_tests)-0.5)*np.pi
        pmr = 0.01*np.random.sample(n_tests)
        pmd = 0.01*np.random.sample(n_tests)
        px = 0.00045 + 0.001*np.random.sample(n_tests)
        rv = 1000.0*np.random.sample(n_tests)
        r_control = None
        d_control = None
        for (rr, dd, pr, pd, x, v) in zip(ra_in, dec_in, pmr, pmd, px, rv):
            r, d = pal.pm(rr, dd, pr, pd, x, v, ep0, ep1)
            if r_control is None:
                r_control = np.array([r])
                d_control = np.array([d])
            else:
                r_control = np.append(r_control, r)
                d_control = np.append(d_control, d)

        rTest, dTest = pal.pmVector(ra_in, dec_in, pmr, pmd, px, rv, ep0, ep1)
        for (r1, d1, r2, d2) in zip(r_control, d_control, rTest, dTest):
            self.assertAlmostEqual(r1, r2, 12)
            self.assertAlmostEqual(d1, d2, 12)

        # test that an exception is raised if input arrays are of
        # inconsistent length
        self.assertRaises(ValueError, pal.pmVector, ra_in, dec_in[:3], pmr,
                          pmd, px, rv, ep0, ep1)

        self.assertRaises(ValueError, pal.pmVector, ra_in, dec_in, pmr[:3],
                          pmd, px, rv, ep0, ep1)

        self.assertRaises(ValueError, pal.pmVector, ra_in, dec_in, pmr,
                          pmd[:3], px, rv, ep0, ep1)

        self.assertRaises(ValueError, pal.pmVector, ra_in, dec_in, pmr,
                          pmd, px[:3], rv, ep0, ep1)

        self.assertRaises(ValueError, pal.pmVector, ra_in, dec_in, pmr,
                          pmd, px, rv[:3], ep0, ep1)

    def test_polmo(self):
        (elong, phi, daz) = pal.polmo(0.7, -0.5, 1.0e-6, -2.0e-6)

        self.assertAlmostEqual(elong, 0.7000004837322044, 12)
        self.assertAlmostEqual(phi, -0.4999979467222241, 12)
        self.assertAlmostEqual(daz, 1.008982781275728e-6, 12)

    def test_prebn(self):
        expected = np.array([
            [9.999257613786738e-1, -1.117444640880939e-2, -4.858341150654265e-3],
            [1.117444639746558e-2, 9.999375635561940e-1, -2.714797892626396e-5],
            [4.858341176745641e-3, -2.714330927085065e-5, 9.999881978224798e-1]
        ])
        rmatp = pal.prebn(1925, 1975)
        np.testing.assert_array_almost_equal(rmatp, expected, 12)

    def test_prec(self):
        expected = np.array([
            [0.9999856154510, -0.0049192906204, -0.0021376320580],
            [0.0049192906805, 0.9999879002027, -5.2297405698747e-06],
            [0.0021376319197, -5.2859681191735e-06, 0.9999977152483]])
        rmat = pal.prec(1990, 2012)
        np.testing.assert_array_almost_equal(rmat, expected, 12)

    def test_preces(self):
        (ra, dc) = pal.preces("FK4", 1925, 1950, 6.28, -1.123)
        self.assertAlmostEqual(ra, 0.002403604864728447, 12)
        self.assertAlmostEqual(dc, -1.120570643322045, 12)

        (ra, dec) = pal.preces("FK5", 2050, 1990, 0.0123, 1.0987)
        self.assertAlmostEqual(ra, 6.282003602708382, 5)
        self.assertAlmostEqual(dc, -1.120570643322045, 6)

    def test_pvobs(self):
        expected = np.array([-4.7683600138836167813e-06,
                            1.0419056712717953176e-05,
                            4.099831053320363277e-05,
                            -7.5976959740661272483e-10,
                            -3.4771429582640930371e-10,
                            0.0])
        pv = pal.pvobs(1.3, 10000, 2)
        np.testing.assert_array_almost_equal(pv, expected, decimal=12)

    def test_range(self):
        self.assertAlmostEqual(pal.drange(-4), 2.283185307179586, 12)

    def test_drange_vector(self):
        """
        Test that drangeVector returns results consistent with drange
        """
        np.random.seed(140)
        n_samples = 1000
        angle_in = np.random.random_sample(n_samples)*10.0*np.pi+2.0*np.pi
        angle_in = np.append(angle_in, np.random.random_sample(n_samples)*(-10.0)*np.pi)

        test_angle = pal.drangeVector(angle_in)

        for ii in range(len(angle_in)):
            control_angle = pal.drange(angle_in[ii])
            self.assertEqual(control_angle, test_angle[ii])
            self.assertAlmostEqual(np.sin(angle_in[ii]), np.sin(test_angle[ii]), 10)
            self.assertAlmostEqual(np.cos(angle_in[ii]), np.cos(test_angle[ii]), 10)

    def test_ranorm(self):
        self.assertAlmostEqual(pal.dranrm(-0.1), 6.183185307179587, 12)

    def test_dranrm_vector(self):
        """
        Test that dranrmVector returns the same results as dranrm
        """
        np.random.seed(74310)
        n_samples = 100
        angle_in = (np.random.random_sample(100) - 0.5) * 4.0 * np.pi
        angle_in = np.array([aa + 2.0*np.pi if aa > 0.0 else aa for aa in angle_in])

        # make sure we created input angles that are outside the
        # 0 to 2pi range
        for aa in angle_in:
            if aa > 0.0 and aa < 2.0 * np.pi:
                raise RuntimeError("did not create correct inputs for test_dranrmVector")

        test_angle = pal.dranrmVector(angle_in)

        # verify that the resulting angles are equivalent to the input
        # angles normalized to the 0-2pi ranges
        self.assertGreater(test_angle.min(), 0.0)
        self.assertLess(test_angle.max(), 2.0 * np.pi)
        np.testing.assert_array_almost_equal(np.cos(test_angle), np.cos(angle_in), 12)
        np.testing.assert_array_almost_equal(np.sin(test_angle), np.sin(angle_in), 12)

        for ii in range(n_samples):
            control_angle = pal.dranrm(angle_in[ii])
            self.assertEqual(control_angle, test_angle[ii])

    def test_ref(self):
        self.assertAlmostEqual(pal.refro(1.4, 3456.7, 280, 678.9, 0.9, 0.55,
                               -0.3, 0.006, 1e-9),
                               0.00106715763018568, 12)
        self.assertAlmostEqual(pal.refro(1.4, 3456.7, 280, 678.9, 0.9, 1000,
                               -0.3, 0.006, 1e-9),
                               0.001296416185295403, 12)

        (refa, refb) = pal.refcoq(275.9, 709.3, 0.9, 101)
        self.assertAlmostEqual(refa, 2.324736903790639e-4, 12)
        self.assertAlmostEqual(refb, -2.442884551059e-7, 15)

        (refa, refb) = pal.refco(2111.1, 275.9, 709.3, 0.9, 101,
                                 -1.03, 0.0067, 1e-12)
        self.assertAlmostEqual(refa, 2.324673985217244e-4, 12)
        self.assertAlmostEqual(refb, -2.265040682496e-7, 15)

        (refa, refb) = pal.refcoq(275.9, 709.3, 0.9, 0.77)
        self.assertAlmostEqual(refa, 2.007406521596588e-4, 12)
        self.assertAlmostEqual(refb, -2.264210092590e-7, 15)

        (refa, refb) = pal.refco(2111.1, 275.9, 709.3, 0.9, 0.77,
                                 -1.03, 0.0067, 1e-12)
        self.assertAlmostEqual(refa, 2.007202720084551e-4, 12)
        self.assertAlmostEqual(refb, -2.223037748876e-7, 15)

        (refa2, refb2) = pal.atmdsp(275.9, 709.3, 0.9, 0.77,
                                    refa, refb, 0.5)
        self.assertAlmostEqual(refa2, 2.034523658888048e-4, 12)
        self.assertAlmostEqual(refb2, -2.250855362179e-7, 15)

        vu = pal.dcs2c(0.345, 0.456)
        vr = pal.refv(vu, refa, refb)
        self.assertAlmostEqual(vr[0], 0.8447487047790478, 12)
        self.assertAlmostEqual(vr[1], 0.3035794890562339, 12)
        self.assertAlmostEqual(vr[2], 0.4407256738589851, 12)

        vu = pal.dcs2c(3.7, 0.03)
        vr = pal.refv(vu, refa, refb)
        self.assertAlmostEqual(vr[0], -0.8476187691681673, 12)
        self.assertAlmostEqual(vr[1], -0.5295354802804889, 12)
        self.assertAlmostEqual(vr[2], 0.0322914582168426, 12)

        zr = pal.refz(0.567, refa, refb)
        self.assertAlmostEqual(zr, 0.566872285910534, 12)

        zr = pal.refz(1.55, refa, refb)
        self.assertAlmostEqual(zr, 1.545697350690958, 12)

        np.random.seed(32)
        zu_int = np.random.sample(20)*1.4
        zr_control = np.zeros(20, dtype=np.float64)
        for i, zu in enumerate(zu_int):
            zr = pal.refz(zu, refa, refb)
            zr_control[i] = zr
        zr_test = pal.refzVector(zu_int, refa, refb)
        for zt, zc in zip(zr_test, zr_control):
            self.assertAlmostEqual(zt, zc, 12)

    def test_refv_vector(self):
        """
        Test that refvVector gives the same results as iterating over
        vectors with refv
        """
        np.random.seed(118)
        input_data = np.random.random_sample((4, 3)) * 5.0
        ref_a = 0.02
        ref_b = 0.05
        output_data = pal.refvVector(input_data, ref_a, ref_b)
        for inVector, outVector in zip(input_data, output_data):
            control_vector = pal.refv(inVector, ref_a, ref_b)
            np.testing.assert_array_equal(outVector, control_vector)

    def test_refro_vector(self):
        """
        Test that refroVector returns the same results as just calling
        refro on each element of the input vector individually.
        """
        np.random.seed(119)
        hm = 150.0
        tdk = 289.0
        pmb = 670.0
        rh = 0.7
        wl = 0.6
        phi = 0.3
        tlr = 0.008
        eps = 1.0e-10
        zobs = np.random.random_sample(20)*np.pi*0.5
        test_output = pal.refroVector(zobs, hm, tdk, pmb, rh, wl, phi, tlr, eps)
        for zz, test in zip(zobs, test_output):
            control = pal.refro(zz, hm, tdk, pmb, rh, wl, phi, tlr, eps)
            self.assertEqual(test, control)

    def test_refc(self):  # This is the SOFA test
        (refa, refb) = pal.refcoq(10.0+273.15, 800.0, 0.9, 0.4)
        self.assertAlmostEqual(refa, 0.2264949956241415009e-3, 15)
        self.assertAlmostEqual(refb, -0.2598658261729343970e-6, 18)

    def test_rv(self):
        self.assertAlmostEqual(pal.rverot(-0.777, 5.67, -0.3, 3.19),
                               -0.1948098355075913, 6)
        self.assertAlmostEqual(pal.rvgalc(1.11, -0.99),
                               158.9630759840254, 3)
        self.assertAlmostEqual(pal.rvlg(3.97, 1.09),
                               -197.818762175363, 3)
        self.assertAlmostEqual(pal.rvlsrd(6.01, 0.1),
                               -4.082811335150567, 4)
        self.assertAlmostEqual(pal.rvlsrk(6.01, 0.1),
                               -5.925180579830265, 4)

    def test_rvgalc(self):
        self.assertAlmostEqual(pal.rvgalc(2.7, -1.0), 213.98084425751144977, 12)

    def test_rvlg(self):
        self.assertAlmostEqual(pal.rvlg(2.7, -1.0), 291.79205281252404802, 12)

    def test_rvlsrd(self):
        self.assertAlmostEqual(pal.rvlsrd(2.7, -1.0), 9.620674692097630043, 12)

    def test_rvlsrk(self):
        self.assertAlmostEqual(pal.rvlsrk(2.7, -1.0), 12.556356851411955233, 12)

    def test_sep(self):
        d1 = np.array([1.0, 0.1, 0.2])
        d2 = np.array([-3.0, 1e-3, 0.2])
        (ad1, bd1) = pal.dcc2s(d1)
        (ad2, bd2) = pal.dcc2s(d2)

        self.assertAlmostEqual(pal.dsep(ad1, bd1, ad2, bd2),
                               2.8603919190246608, 7)
        self.assertAlmostEqual(pal.dsepv(d1, d2),
                               2.8603919190246608, 7)

    def test_dsep_vector(self):
        np.random.seed(32)
        ra1 = np.random.sample(20)*2.0*np.pi
        dec1 = np.random.sample(20)*2.0*np.pi
        ra2 = np.random.sample(20)*2.0*np.pi
        dec2 = np.random.sample(20)*2.0*np.pi
        dd_control = np.zeros(20, dtype=np.float64)
        for (i, rr) in enumerate(zip(ra1, dec1, ra2, dec2)):
            dd = pal.dsep(rr[0], rr[1], rr[2], rr[3])
            dd_control[i] = dd
        dd_test = pal.dsepVector(ra1, dec1, ra2, dec2)
        for (ddc, ddt) in zip(dd_test, dd_control):
            self.assertAlmostEqual(ddc, ddt, 12)

        # test that an exception is raised if input arrays
        # have different lengths
        self.assertRaises(ValueError, pal.dsepVector, ra1, dec1[:5], ra2, dec2)

        self.assertRaises(ValueError, pal.dsepVector, ra1, dec1, ra2[:5], dec2)

        self.assertRaises(ValueError, pal.dsepVector, ra1, dec1, ra2, dec2[:5])

    def test_dsepv_vector(self):
        """
        Test that dsepvVector produces results that agree
        with dsepv
        """
        np.random.seed(130)
        n_samples = 100
        v1 = (np.random.random_sample((n_samples, 3))-0.5)*100.0
        v2 = (np.random.random_sample((n_samples, 3))-0.5)*100.0

        test_sep = pal.dsepvVector(v1, v2)
        for ii in range(n_samples):
            control_sep = pal.dsepv(v1[ii], v2[ii])
            self.assertEqual(control_sep, test_sep[ii])
            self.assertFalse(np.isnan(test_sep[ii]))

        test_sep = pal.dsepvVector(v1, v2[6])
        for ii in range(n_samples):
            control_sep = pal.dsepv(v1[ii], v2[6])
            self.assertEqual(control_sep, test_sep[ii])
            self.assertFalse(np.isnan(test_sep[ii]))

        # test that exceptions are raised when they ought to be
        self.assertRaises(ValueError, pal.dsepvVector, v1, v2[0:19])

        v3 = np.random.random_sample((n_samples, 2))
        self.assertRaises(ValueError, pal.dsepvVector, v3, v2)

        self.assertRaises(ValueError, pal.dsepvVector, v1, v3)

        self.assertRaises(ValueError, pal.dsepvVector, v1,
                          np.random.random_sample(4))

    def test_supgal(self):
        (dl, db) = pal.supgal(6.1, -1.4)
        self.assertAlmostEqual(dl, 3.798775860769474, 12)
        self.assertAlmostEqual(db, -0.1397070490669407, 12)

    def test_supgal_vector(self):
        """
        Test that supgalVector returns results consistent with supgal
        """

        np.random.seed(134)
        n_samples = 200
        dsl_list = np.random.random_sample(n_samples)*2.0*np.pi
        dsb_list = (np.random.random_sample(n_samples)-0.5)*np.pi

        testDl, testDb = pal.supgalVector(dsl_list, dsb_list)

        for ii in range(n_samples):
            controlDl, controlDb = pal.supgal(dsl_list[ii], dsb_list[ii])
            self.assertEqual(controlDl, testDl[ii])
            self.assertEqual(controlDb, testDb[ii])

        # test that supgalVector and galsupVector invert each other
        testDsl, testDsb = pal.galsupVector(testDl, testDb)
        np.testing.assert_array_almost_equal(testDsl, dsl_list, 9)
        np.testing.assert_array_almost_equal(testDsb, dsb_list, 9)

        # test that an exception is raised if the input arrays are of
        # different length
        self.assertRaises(ValueError, pal.supgalVector, dsl_list[:16], dsb_list)

    def test_tp(self):
        dr0 = 3.1
        dd0 = -0.9
        dr1 = dr0 + 0.2
        dd1 = dd0 - 0.1
        (dx, dy) = pal.ds2tp(dr1, dd1, dr0, dd0)
        self.assertAlmostEqual(dx, 0.1086112301590404, 12)
        self.assertAlmostEqual(dy, -0.1095506200711452, 12)

        (dr2, dd2) = pal.dtp2s(dx, dy, dr0, dd0)
        self.assertAlmostEqual(dr2 - dr1, 0.0, 12)
        self.assertAlmostEqual(dd2 - dd1, 0.0, 12)

        (dr01, dd01, dr02, dd02) = pal.dtps2c(dx, dy, dr2, dd2)
        self.assertAlmostEqual(dr01, dr0, 12)
        self.assertAlmostEqual(dd01, dd0)
        self.assertIsNone(dr02)
        self.assertIsNone(dd02)

    def test_ds2tp_vector(self):
        raz = 0.3
        decz = 0.1
        n_tests = 100
        np.random.seed(32)
        ra_in = raz + (np.random.sample(n_tests) - 0.5) * 0.2
        dec_in = decz + (np.random.sample(n_tests) - 0.5) * 0.2
        xi_control = None
        eta_control = None
        for (rr, dd) in zip(ra_in, dec_in):
            xi, eta = pal.ds2tp(rr, dd, raz, decz)
            if xi_control is None:
                xi_control = np.array([xi])
                eta_control = np.array([eta])
            else:
                xi_control = np.append(xi_control, xi)
                eta_control = np.append(eta_control, eta)

        xiTest, etaTest = pal.ds2tpVector(ra_in, dec_in, raz, decz)
        for (x1, e1, x2, e2) in zip(xi_control, eta_control, xiTest, etaTest):
            self.assertAlmostEqual(x1, x2, 12)
            self.assertAlmostEqual(e1, e2, 12)

        # test that bad input values get mapped to NaNs
        ra_in[5] = raz + np.pi
        dec_in[7] = decz + np.pi
        xiTest, etaTest = pal.ds2tpVector(ra_in, dec_in, raz, decz)
        for ii in range(len(ra_in)):
            if ii in (5, 7):
                self.assertTrue(np.isnan(xiTest[ii]))
                self.assertTrue(np.isnan(etaTest[ii]))
                self.assertRaises(ValueError, pal.ds2tp, ra_in[ii], dec_in[ii], raz, decz)
            else:
                self.assertEqual(xiTest[ii], xi_control[ii])
                self.assertEqual(etaTest[ii], eta_control[ii])

        # test that an exception is raised if input arrays
        # are not of the same shape
        self.assertRaises(ValueError, pal.ds2tpVector, ra_in[:3], dec_in,
                          raz, decz)

    def test_dtp2s_vector(self):
        """
        Test that dtp2sVector produces results consistent with
        dtp2s
        """

        np.random.seed(139)
        n_samples = 200
        raz = 0.9
        decz = -0.3

        xi = (np.random.random_sample(n_samples)-0.5)*10.0
        eta = (np.random.random_sample(n_samples)-0.5)*10.0
        testRa, testDec = pal.dtp2sVector(xi, eta, raz, decz)

        for ii in range(n_samples):
            controlRa, controlDec = pal.dtp2s(xi[ii], eta[ii], raz, decz)
            self.assertEqual(testRa[ii], controlRa)
            self.assertEqual(testDec[ii], controlDec)
            self.assertFalse(np.isnan(testRa[ii]))
            self.assertFalse(np.isnan(testDec[ii]))

        # test that dtp2sVector and ds2tpVector invert each other
        testXi, testEta = pal.ds2tpVector(testRa, testDec, raz, decz)
        np.testing.assert_array_almost_equal(testXi, xi, 12)
        np.testing.assert_array_almost_equal(testEta, eta, 12)

        # test that an exception is raised if input arrays are of
        # different sizes
        self.assertRaises(ValueError, pal.dtp2sVector, xi, eta[:20], raz, decz)

    def test_vecmat(self):
        # Not everything is implemented here
        dav = np.array([-0.123, 0.0987, 0.0654])
        dav2m_expected = np.array([
            [0.9930075842721269, 0.05902743090199868, -0.1022335560329612],
            [-0.07113807138648245, 0.9903204657727545, -0.1191836812279541],
            [0.09420887631983825, 0.1256229973879967, 0.9875948309655174],
        ])

        deuler_expected = np.array([
            [-0.1681574770810878, 0.1981362273264315, 0.9656423242187410],
            [-0.2285369373983370, 0.9450659587140423, -0.2337117924378156],
            [-0.9589024617479674, -0.2599853247796050, -0.1136384607117296]])

        dmxm_expected = np.array([
            [-0.09010460088585805, 0.3075993402463796, 0.9472400998581048],
            [-0.3161868071070688, 0.8930686362478707, -0.3200848543149236],
            [-0.9444083141897035, -0.3283459407855694, 0.01678926022795169]])

        drm1 = pal.dav2m(dav)
        np.testing.assert_array_almost_equal(drm1, dav2m_expected, decimal=12)

        drm2 = pal.deuler("YZY", 2.345, -0.333, 2.222)
        np.testing.assert_array_almost_equal(drm2, deuler_expected, decimal=12)

        drm = pal.dmxm(drm2, drm1)
        np.testing.assert_array_almost_equal(drm, dmxm_expected, decimal=12)

        dv1 = pal.dcs2c(3.0123, -0.999)
        dcs2c_expected = np.array([-0.5366267667260525, 0.06977111097651444, -0.8409302618566215])
        np.testing.assert_array_almost_equal(dv1, dcs2c_expected, decimal=12)

        dv2 = pal.dmxv(drm1, dv1)
        dv3 = pal.dmxv(drm2, dv2)
        dmxv_expected = np.array([-0.7267487768696160, 0.5011537352639822, 0.4697671220397141])
        np.testing.assert_array_almost_equal(dv3, dmxv_expected, decimal=12)

        dv4 = pal.dimxv(drm, dv3)
        dimxv_expected = np.array([-0.5366267667260526, 0.06977111097651445, -0.8409302618566215])
        np.testing.assert_array_almost_equal(dv4, dimxv_expected, decimal=12)

        dv5 = pal.dm2av(drm)
        dm2av_expected = np.array([0.006889040510209034, -1.577473205461961, 0.5201843672856759])
        np.testing.assert_array_almost_equal(dv5, dm2av_expected, decimal=12)

        dv5 *= 1000.0

        (dv6, dvm) = pal.dvn(dv5)
        dvn_expected = np.array([0.004147420704640065, -0.9496888606842218, 0.3131674740355448])
        np.testing.assert_array_almost_equal(dv6, dvn_expected, decimal=12)
        self.assertAlmostEqual(dvm, 1661.042127339937, 9)

        dv7 = pal.dvxv(dv6, dv1)
        dvxv_expected = np.array([0.7767720597123304, -0.1645663574562769, -0.5093390925544726])
        np.testing.assert_array_almost_equal(dv7, dvxv_expected, decimal=12)

if __name__ == '__main__':
    unittest.main()
