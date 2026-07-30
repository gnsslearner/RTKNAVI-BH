#include "rtklib.h"
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#ifdef IERS_MODEL

static void st1idiu(const double *xsta, const double *xsun, const double *xmon, double fac2sun, double fac2mon, double *xcorsta)
{
    /* Initialized data */

    const double  dhi = -.0025;
    const double  dli = -7e-4;

    double de, dn, dr, rsta, rmon, rsun;
    double cosla, demon, sinla, dnmon, desun, drmon, dnsun, drsun,
        cosphi, sinphi, cos2phi;

    /* Compute the normalized position vector of the IGS station. */
    rsta = norm(xsta, 3);
    sinphi = xsta[2] / rsta;
    cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
    cos2phi = cosphi * cosphi - sinphi * sinphi;
    sinla = xsta[1] / cosphi / rsta;
    cosla = xsta[0] / cosphi / rsta;
    /* Compute the normalized position vector of the Moon. */
    rmon = norm(xmon, 3);
    /* Compute the normalized position vector of the Sun. */
    rsun = norm(xsun, 3);
    /* Computing 2nd power */
    drsun = -3.0 *dhi *  sinphi * cosphi * fac2sun * xsun[2] * (xsun[0] * sinla - xsun[1] * cosla) / (rsun *  rsun);
    /* Computing 2nd power */
    drmon = -3.0*dhi  * sinphi * cosphi * fac2mon * xmon[2] * (xmon[0] * sinla - xmon[1] * cosla) / (rmon * rmon);
    /* Computing 2nd power */
    dnsun = -3.0 * dli * cos2phi * fac2sun * xsun[2] * (xsun[0] * sinla - xsun[1] * cosla) / (rsun* rsun);
    /* Computing 2nd power */
    dnmon = -3.0 *dli *  cos2phi * fac2mon * xmon[2] * (xmon[0] * sinla - xmon[1] * cosla) / (rmon * rmon);
    /* Computing 2nd power */
    desun = -3.0 *dli *  sinphi * fac2sun * xsun[2] * (xsun[0] * cosla + xsun[1] * sinla) / (rsun * rsun);
    /* Computing 2nd power */
    demon = -3.0*dli * sinphi * fac2mon * xmon[2] * (xmon[0] * cosla + xmon[1] * sinla) / (rmon * rmon);
    dr = drsun + drmon;
    dn = dnsun + dnmon;
    de = desun + demon;
    /*  Compute the corrections for the station. */
    xcorsta[0] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
    xcorsta[1] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
    xcorsta[2] = dr * sinphi + dn * cosphi;
}

static void st1isem(const double *xsta, const double *xsun, const double *xmon, double fac2sun, double fac2mon, double *xcorsta)
{
    /* Initialized data */

    const double  dhi = -0.0022;
    const double  dli = -7e-4;

    double costwola, sintwola, de, dn, dr, rsta, rmon, rsun;
    double cosla, demon, sinla, dnmon, desun, drmon, dnsun, drsun,
        cosphi, sinphi;

    /* Compute the normalized position vector of the IGS station. */
    rsta = norm(xsta, 3);
    sinphi = xsta[2] / rsta;
    cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
    sinla = xsta[1] / cosphi / rsta;
    cosla = xsta[0] / cosphi / rsta;
    costwola = cosla * cosla - sinla * sinla;
    sintwola = 2.0 * cosla * sinla;
    /* Compute the normalized position vector of the Moon. */
    rmon = norm(xmon, 3);
    /* Compute the normalized position vector of the Sun. */
    rsun = norm(xsun, 3);

    drsun = -0.75 *dhi *  (cosphi * cosphi) * fac2sun * ((xsun[0] * xsun[0] - xsun[1] * xsun[1])
        * sintwola - 2.0 *xsun[0] * xsun[1] * costwola) / (rsun *	rsun);

    drmon = -0.75 *dhi *  (cosphi * cosphi) * fac2mon * ((xmon[0] * xmon[0] - xmon[1] * xmon[1])
        * sintwola - 2.0 *xmon[0] * xmon[1] * costwola) / (rmon *rmon);

    dnsun = 1.5 *dli *  sinphi * cosphi * fac2sun * ((xsun[0] * xsun[0] - xsun[1] * xsun[1])
        * sintwola - 2.0 *xsun[0] * xsun[1] * costwola) / (rsun *rsun);

    dnmon = 1.5 *dli * sinphi * cosphi * fac2mon * ((xmon[0] * xmon[0] - xmon[1] * xmon[1])
        * sintwola - 2.0 * xmon[0] * xmon[1] * costwola) / (rmon *rmon);

    desun = -1.5 * dli *  cosphi * fac2sun * ((xsun[0] * xsun[0] - xsun[1] * xsun[1]) *
        costwola + 2.0 *xsun[0] * xsun[1] * sintwola) / (rsun * rsun);

    demon = -1.5 *dli *  cosphi * fac2mon * ((xmon[0] * xmon[0] - xmon[1] * xmon[1]) *
        costwola + 2.0 * xmon[0] * xmon[1] * sintwola) / (rmon *rmon);
    dr = drsun + drmon;
    dn = dnsun + dnmon;
    de = desun + demon;
    xcorsta[0] = dr * cosla * cosphi - de * sinla - dn * sinphi * cosla;
    xcorsta[1] = dr * sinla * cosphi + de * cosla - dn * sinphi * sinla;
    xcorsta[2] = dr * sinphi + dn * cosphi;

}

static void st1l1(const double *xsta, const double *xsun, const double *xmon, double fac2sun, double fac2mon, double *xcorsta)
{
    /* Initialized data */
    const double l1d = 0.0012;
    const double l1sd = 0.0024;


    double costwola, sintwola, l1, de, dn, rsta, rmon, rsun;
    double cosla, demon, sinla, dnmon, desun, dnsun, cosphi, sinphi;

    /* Compute the normalized position vector of the IGS station. */
    rsta = norm(xsta, 3);
    sinphi = xsta[2] / rsta;
    cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
    sinla = xsta[1] / cosphi / rsta;
    cosla = xsta[0] / cosphi / rsta;

    /* Compute the normalized position vector of the Moon. */
    rmon = norm(xmon, 3);
    /* Compute the normalized position vector of the Sun. */
    rsun = norm(xsun, 3);
    /* Compute the station corrections for the diurnal band. */
    l1 = l1d;
    dnsun = -l1 * (sinphi * sinphi) * fac2sun * xsun[2] * (xsun[0] * cosla + xsun[1] * sinla) / (rsun * rsun);
    dnmon = -l1 * (sinphi * sinphi) * fac2mon * xmon[2] * (xmon[0] * cosla + xmon[1] * sinla) / (rmon * rmon);

    desun = l1 * sinphi * (cosphi * cosphi - sinphi * sinphi) * fac2sun * xsun[2] * (
        xsun[0] * sinla - xsun[1] * cosla) / (rsun * rsun);
    demon = l1 * sinphi * (cosphi * cosphi - sinphi * sinphi) * fac2mon * xmon[2] * (
        xmon[0] * sinla - xmon[1] * cosla) / (rmon * rmon);
    de = (desun + demon) * 3.0;
    dn = (dnsun + dnmon) * 3.0;
    xcorsta[0] = -de * sinla - dn * sinphi * cosla;
    xcorsta[1] = de * cosla - dn * sinphi * sinla;
    xcorsta[2] = dn * cosphi;

    /* Compute the station corrections for the semi-diurnal band. */
    l1 = l1sd;
    costwola = cosla * cosla - sinla * sinla;
    sintwola = 2.0*cosla * sinla;

    dnsun = -l1 / 2.0 * sinphi * cosphi * fac2sun * ((xsun[0] * xsun[0] - xsun[1] * xsun[1])
        * costwola + 2.0 *xsun[0] * xsun[1] * sintwola) / (rsun *rsun);

    dnmon = -l1 / 2.0 * sinphi * cosphi * fac2mon * ((xmon[0] * xmon[0] - xmon[1] * xmon[1])
        * costwola + 2.0 *xmon[0] * xmon[1] * sintwola) / (rmon *rmon);

    desun = -l1 / 2.0 * (sinphi *sinphi) * cosphi * fac2sun * ((xsun[0] * xsun[0] - xsun[1] * xsun[1])
        * sintwola - 2.0 *xsun[0] * xsun[1] * costwola) / (rsun *rsun);

    demon = -l1 / 2.0 * (sinphi * sinphi) * cosphi * fac2mon * ((xmon[0] * xmon[0] - xmon[1] * xmon[1])
        * sintwola - 2.0 *xmon[0] * xmon[1] * costwola) / (rmon *rmon);
    de = (desun + demon) * 3.0;
    dn = (dnsun + dnmon) * 3.0;
    xcorsta[0] = xcorsta[0] - de * sinla - dn * sinphi * cosla;
    xcorsta[1] = xcorsta[1] + de * cosla - dn * sinphi * sinla;
    xcorsta[2] += dn * cosphi;
}

static void step2diu(const double *xsta, double fhr, double t, double *xcorsta)
{
    const double datdi[279]	/* was [9][31] */ = { -3., 0., 2., 0., 0.,
        -.01, 0., 0., 0., -3., 2., 0., 0., 0., -.01, 0., 0., 0., -2., 0., 1., -1., 0., -.02,
        0., 0., 0., -2., 0., 1., 0., 0., -.08, 0., -.01, .01, -2., 2., -1., 0., 0., -.02,
        0., 0., 0., -1., 0., 0., -1., 0., -.1, 0., 0., 0., -1., 0., 0., 0., 0., -.51, 0.,
        -.02, .03, -1., 2., 0., 0., 0., .01, 0., 0., 0., 0., -2., 1., 0., 0., .01, 0., 0.,
        0., 0., 0., -1., 0., 0., .02, 0., 0., 0., 0., 0., 1., 0., 0., .06, 0., 0., 0., 0., 0.,
        1., 1., 0., .01, 0., 0., 0., 0., 2., -1., 0., 0., .01, 0., 0., 0., 1., -3., 0., 0.,
        1., -.06, 0., 0., 0., 1., -2., 0., -1., 0., .01, 0., 0., 0., 1., -2., 0., 0., 0.,
        -1.23, -.07, .06, .01, 1., -1., 0., 0., -1., .02, 0., 0., 0., 1., -1., 0., 0., 1.,
        .04, 0., 0., 0., 1., 0., 0., -1., 0., -.22, .01, .01, 0., 1., 0., 0., 0., 0., 12.,
        -.8, -.67, -.03, 1., 0., 0., 1., 0., 1.73, -.12, -.1, 0., 1., 0., 0., 2., 0., -.04,
        0., 0., 0., 1., 1., 0., 0., -1., -.5, -.01, .03, 0., 1., 1., 0., 0., 1., .01, 0., 0.,
        0., 0., 1., 0., 1., -1., -.01, 0., 0., 0., 1., 2., -2., 0., 0., -.01, 0., 0., 0., 1.,
        2., 0., 0., 0., -.11, .01, .01, 0., 2., -2., 1., 0., 0., -.01, 0., 0., 0., 2., 0.,
        -1., 0., 0., -.02, 0., 0., 0., 3., 0., 0., 0., 0., 0., 0., 0., 0., 3., 0., 0., 1., 0.,
        0., 0., 0., 0. };

    const double c_b2 = 360.0;
    const double deg2rad = .017453292519943295;

    double h;
    int i, j;
    double p, s, de, dn, dr, pr, ps, zla, tau, zns, rsta, cosla,
        sinla, thetaf, cosphi, sinphi;


    /*  Compute the phase angles in degrees. */
    s = ((t * 1.85139e-6 - .0014663889) * t + 481267.88194) * t +
        218.31664563;
    tau = fhr * 15. + 280.4606184 + ((t * -2.58e-8 + 3.8793e-4) * t +
        36000.7700536) * t + (-s);
    pr = (((t * 7e-9 + 2.1e-8) * t + 3.08889e-4) * t + 1.396971278) * t;
    s += pr;
    h = (((t * -6.54e-9 + 2e-8) * t + 3.0322222e-4) * t + 36000.7697489)
        * t + 280.46645;
    p = (((t * 5.263e-8 - 1.24991e-5) * t - .01032172222) * t +
        4069.01363525) * t + 83.35324312;
    zns = (((t * 1.65e-8 - 2.13944e-6) * t - .00207561111) * t +
        1934.13626197) * t + 234.95544499;
    ps = (((t * -3.34e-9 - 1.778e-8) * t + 4.5688889e-4) * t +
        1.71945766667) * t + 282.93734098;
    /* Reduce angles to between the range 0 and 360. */
    s = fmod(s, c_b2);
    tau = fmod(tau, c_b2);
    h = fmod(h, c_b2);
    p = fmod(p, c_b2);
    zns = fmod(zns, c_b2);
    ps = fmod(ps, c_b2);

    rsta = norm(xsta, 3);
    sinphi = xsta[2] / rsta;
    cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
    cosla = xsta[0] / cosphi / rsta;
    sinla = xsta[1] / cosphi / rsta;
    zla = atan2(xsta[1], xsta[0]);
    for (i = 0; i < 3; i++)
    {
        /* Initialize. */
        xcorsta[i] = 0.0;
    }
    for (j = 1; j <= 31; ++j)
    {
        /* Convert from degrees to radians. */
        thetaf = (tau + datdi[j * 9 - 9] * s + datdi[j * 9 - 8] * h + datdi[
            j * 9 - 7] * p + datdi[j * 9 - 6] * zns + datdi[j * 9 - 5] *
                ps) * deg2rad;
        dr = datdi[j * 9 - 4] * 2. * sinphi * cosphi * sin(thetaf + zla) +
            datdi[j * 9 - 3] * 2. * sinphi * cosphi * cos(thetaf + zla);

        dn = datdi[j * 9 - 2] * (cosphi * cosphi - sinphi * sinphi) * sin(thetaf +
            zla) + datdi[j * 9 - 1] * (cosphi * cosphi - sinphi * sinphi) * cos(
                thetaf + zla);
        /*      DE=DATDI(8,J)*SINPHI*COS(THETAF+ZLA)+ */
        /*     Modified 20 June 2007 */
        de = datdi[j * 9 - 2] * sinphi * cos(thetaf + zla) - datdi[j * 9 - 1]
            * sinphi * sin(thetaf + zla);
        xcorsta[0] = xcorsta[0] + dr * cosla * cosphi - de * sinla - dn *
            sinphi * cosla;
        xcorsta[1] = xcorsta[1] + dr * sinla * cosphi + de * cosla - dn *
            sinphi * sinla;
        xcorsta[2] = xcorsta[2] + dr * sinphi + dn * cosphi;
    }
    for (i = 0; i < 3; i++)
    {
        xcorsta[i] /= 1e3;
    }
}

static void step2lon(const double *xsta, double t, double *xcorsta)
{
    const  double datdi[45]	/* was [9][5] */ = { 0., 0., 0., 1., 0., .47, .23,
        .16, .07, 0., 2., 0., 0., 0., -.2, -.12, -.11, -.05, 1., 0., -1., 0., 0., -.11,
        -.08, -.09, -.04, 2., 0., 0., 0., 0., -.13, -.11, -.15, -.07, 2., 0., 0., 1., 0.,
        -.05, -.05, -.06, -.03 };
    const double c_b2 = 360.0;
    const double deg2rad = .017453292519943295;

    double h;
    int i, j;
    double p, s, de, dn, dr, pr, ps, zns, rsta, cosla, sinla,
        thetaf, cosphi, dn_tot, sinphi, dr_tot;

    s = ((t * 1.85139e-6 - .0014663889) * t + 481267.88194) * t +
        218.31664563;
    pr = (((t * 7e-9 + 2.1e-8) * t + 3.08889e-4) * t + 1.396971278) * t;
    s += pr;
    h = (((t * -6.54e-9 + 2e-8) * t + 3.0322222e-4) * t + 36000.7697489)
        * t + 280.46645;
    p = (((t * 5.263e-8 - 1.24991e-5) * t - .01032172222) * t +
        4069.01363525) * t + 83.35324312;
    zns = (((t * 1.65e-8 - 2.13944e-6) * t - .00207561111) * t +
        1934.13626197) * t + 234.95544499;
    ps = (((t * -3.34e-9 - 1.778e-8) * t + 4.5688889e-4) * t +
        1.71945766667) * t + 282.93734098;

    rsta = norm(xsta, 3);
    sinphi = xsta[2] / rsta;

    cosphi = sqrt(xsta[0] * xsta[0] + xsta[1] * xsta[1]) / rsta;
    cosla = xsta[0] / cosphi / rsta;
    sinla = xsta[1] / cosphi / rsta;
    /* Reduce angles to between the range 0 and 360. */
    s = fmod(s, c_b2);
    /*      TAU = DMOD(TAU,360D0) */
    h = fmod(h, c_b2);
    p = fmod(p, c_b2);
    zns = fmod(zns, c_b2);
    ps = fmod(ps, c_b2);
    dr_tot = 0.;
    dn_tot = 0.;
    for (i = 0; i < 3; i++)
    {
        xcorsta[i] = 0.;
    }
    for (j = 1; j <= 5; j++)
    {
        thetaf = (datdi[j * 9 - 9] * s + datdi[j * 9 - 8] * h + datdi[j * 9
            - 7] * p + datdi[j * 9 - 6] * zns + datdi[j * 9 - 5] * ps) *
            deg2rad;

        dr = datdi[j * 9 - 4] * (sinphi * sinphi * 3.0 - 1.0) / 2.0 * cos(thetaf) +
            datdi[j * 9 - 2] * (sinphi * sinphi * 3.0 - 1.0) / 2.0 * sin(thetaf);
        dn = datdi[j * 9 - 3] * (cosphi * sinphi * 2.0) * cos(thetaf) + datdi[
            j * 9 - 1] * (cosphi * sinphi * 2.0) * sin(thetaf);
            de = 0.0;
            dr_tot += dr;
            dn_tot += dn;
            xcorsta[0] = xcorsta[0] + dr * cosla * cosphi - de * sinla - dn *
                sinphi * cosla;
            xcorsta[1] = xcorsta[1] + dr * sinla * cosphi + de * cosla - dn *
                sinphi * sinla;
            xcorsta[2] = xcorsta[2] + dr * sinphi + dn * cosphi;
    }
    for (i = 0; i < 3; i++)
    {
        xcorsta[i] /= 1e3;
    }

}

extern void dehanttideinel(const double *rr, const int *year, const int *mon,const int *day,const double *fhr,
    const double *rs, const double *rm,double *drt)
{
    const double h20 = 0.6078;
    const double l20 = 0.0847;
    const double h3 = 0.292;
    const double l3 = 0.015;

    int i;
    double dxtide[3], xcorsta[3];

    /* SCALAR PRODUCT OF STATION VECTOR WITH SUN/MOON VECTOR */
    double rsta, rsun, rmon, scs, scm, scsun, scmon;

    rsta = norm(rr, 3);
    rsun = norm(rs, 3);
    rmon = norm(rm, 3);
    scs = dot(rr, rs, 3);
    scm = dot(rr, rm, 3);
    scsun = scs / rsta / rsun;
    scmon = scm / rsta / rmon;

    /* COMPUTATION OF NEW H2 AND L2 */
    double cosphi, h2, l2;
    cosphi = sqrt(rr[0] * rr[0] + rr[1] * rr[1]) / rsta;
    h2 = h20 - (1.0 - cosphi * cosphi * 1.5) * 6e-4;
    l2 = l20 + (1.0 - cosphi * cosphi * 1.5) * 2e-4;

    /* P2 term */
    double p2sun, p2mon;
    p2sun = (h2 / 2.0 - l2) * 3.0 * (scsun * scsun) - h2 / 2.0;
    p2mon = (h2 / 2.0 - l2) * 3.0 * (scmon * scmon) - h2 / 2.0;

    /* P3 term */
    double p3sun, p3mon;
    p3sun = (h3 - l3 * 3.0) * 2.5 * (scsun * scsun * scsun) + (l3 - h3) * 1.5 *scsun;
    p3mon = (h3 - l3 * 3.0) * 2.5 * (scmon * scmon * scmon) + (l3 - h3) * 1.5 *scmon;

    /* TERM IN DIRECTION OF SUN/MOON VECTOR */
    double x2sun, x2mon, x3sun, x3mon;
    x2sun = l2 * 3.0 * scsun;
    x2mon = l2 * 3.0 * scmon;
    x3sun = l3 * 3.0 / 2.0 * (scsun * scsun * 5.0 - 1.0);
    x3mon = l3 * 3.0 / 2.0 * (scmon * scmon * 5.0 - 1.0);

    /* FACTORS FOR SUN/MOON USING IAU CURRENT BEST ESTIMATES*/
    const double mass_ratio_sun = 332946.0482;
    const double mass_ratio_moon = .0123000371;
    const double re = 6378136.6;

    double fac2sun, fac2mon, fac3sun, fac3mon;
    fac2sun = mass_ratio_sun * re * pow(re / rsun, 3);
    fac2mon = mass_ratio_moon * re * pow(re / rmon, 3);
    fac3sun = fac2sun * (re / rsun);
    fac3mon = fac2mon * (re / rmon);

    /* TOTAL DISPLACEMENT */

    for (i = 0; i < 3; i++)
    {
        dxtide[i] = fac2sun * (x2sun * rs[i] / rsun + p2sun * rr[i] / rsta)
            + fac2mon * (x2mon * rm[i] / rmon + p2mon * rr[i] / rsta)
            + fac3sun * (x3sun * rs[i] / rsun + p3sun * rr[i] / rsta)
            + fac3mon * (x3mon * rm[i] / rmon + p3mon * rr[i] / rsta);
    }

    /* FIRST, FOR THE DIURNAL BAND */
    st1idiu(rr, rs, rm, fac2sun, fac2mon, xcorsta);
    for (i = 0; i < 3; i++)
    {
        dxtide[i] += xcorsta[i];
    }
    /* SECOND, FOR THE SEMI-DIURNAL BAND */
    st1isem(rr, rs, rm, fac2sun, fac2mon, xcorsta);
    for (i = 0; i < 3; i++)
    {
        dxtide[i] += xcorsta[i];
    }

    /* CORRECTIONS FOR THE LATITUDE DEPENDENCE OF LOVE NUMBERS (PART L^(1) ) */
    st1l1(rr, rs, rm, fac2sun, fac2mon, xcorsta);
    for (i = 0; i < 3; i++)
    {
        dxtide[i] += xcorsta[i];
    }

    gtime_t tt0,tgps0,tutc,tgps;
    double ep[6], t;
    const double ep0[] = { 2000, 1, 1, 12, 0, 0 };
    tt0=epoch2time(ep0);
    tgps0= timeadd(tt0, -19.0 - 32.184);
    ep[0] = (double)*year; ep[1] = (double)*mon; ep[2] = (double)*day;
    ep[3] = ep[4] = ep[5] = 0;
    tutc = epoch2time(ep);
    tutc = timeadd(tutc,*fhr * 3600);
    tgps = utc2gpst(tutc);
    t = timediff(tgps, tgps0) / 86400.0 / 36525.0;

    /*  SECOND, WE CAN CALL THE SUBROUTINE STEP2DIU, FOR THE DIURNAL BAND CORRECTIONS */
    step2diu(rr, *fhr, t, xcorsta);
    for (i = 0; i < 3; i++)
    {
        dxtide[i] += xcorsta[i];
    }

    /*  CORRECTIONS FOR THE LONG-PERIOD BAND */
    step2lon(rr, t, xcorsta);
    for (i = 0; i < 3; i++)
    {
        dxtide[i] += xcorsta[i];
    }

    drt[0] = dxtide[0];
    drt[1] = dxtide[1];
    drt[2] = dxtide[2];
}


/* -----------------------------------------------------------------------------
* Name       : gmf
* Description: compute Global Mapping function (tropospheric mapping function)
*-----------------------------------------------------------------------------*/
extern void gmf(const double *mjd, const double *lat, const double *lon, const double *hgt, 
    const double *zd,double *gmfh, double *gmfw)
{
    const double ah_mean[55] =
    {
        +1.2517e+02, +8.503e-01, +6.936e-02, -6.760e+00, +1.771e-01,
        +1.130e-02, +5.963e-01, +1.808e-02, +2.801e-03, -1.414e-03,
        -1.212e+00, +9.300e-02, +3.683e-03, +1.095e-03, +4.671e-05,
        +3.959e-01, -3.867e-02, +5.413e-03, -5.289e-04, +3.229e-04,
        +2.067e-05, +3.000e-01, +2.031e-02, +5.900e-03, +4.573e-04,
        -7.619e-05, +2.327e-06, +3.845e-06, +1.182e-01, +1.158e-02,
        +5.445e-03, +6.219e-05, +4.204e-06, -2.093e-06, +1.540e-07,
        -4.280e-08, -4.751e-01, -3.490e-02, +1.758e-03, +4.019e-04,
        -2.799e-06, -1.287e-06, +5.468e-07, +7.580e-08, -6.300e-09,
        -1.160e-01, +8.301e-03, +8.771e-04, +9.955e-05, -1.718e-06,
        -2.012e-06, +1.170e-08, +1.790e-08, -1.300e-09, +1.000e-10
    };

    const double bh_mean[55] =
    {
        +0.000e+00, +0.000e+00, +3.249e-02, +0.000e+00, +3.324e-02,
        +1.850e-02, +0.000e+00, -1.115e-01, +2.519e-02, +4.923e-03,
        +0.000e+00, +2.737e-02, +1.595e-02, -7.332e-04, +1.933e-04,
        +0.000e+00, -4.796e-02, +6.381e-03, -1.599e-04, -3.685e-04,
        +1.815e-05, +0.000e+00, +7.033e-02, +2.426e-03, -1.111e-03,
        -1.357e-04, -7.828e-06, +2.547e-06, +0.000e+00, +5.779e-03,
        +3.133e-03, -5.312e-04, -2.028e-05, +2.323e-07, -9.100e-08,
        -1.650e-08, +0.000e+00, +3.688e-02, -8.638e-04, -8.514e-05,
        -2.828e-05, +5.403e-07, +4.390e-07, +1.350e-08, +1.800e-09,
        +0.000e+00, -2.736e-02, -2.977e-04, +8.113e-05, +2.329e-07,
        +8.451e-07, +4.490e-08, -8.100e-09, -1.500e-09, +2.000e-10
    };

    const double ah_amp[55] =
    {
        -2.738e-01, -2.837e+00, +1.298e-02, -3.588e-01, +2.413e-02,
        +3.427e-02, -7.624e-01, +7.272e-02, +2.160e-02, -3.385e-03,
        +4.424e-01, +3.722e-02, +2.195e-02, -1.503e-03, +2.426e-04,
        +3.013e-01, +5.762e-02, +1.019e-02, -4.476e-04, +6.790e-05,
        +3.227e-05, +3.123e-01, -3.535e-02, +4.840e-03, +3.025e-06,
        -4.363e-05, +2.854e-07, -1.286e-06, -6.725e-01, -3.730e-02,
        +8.964e-04, +1.399e-04, -3.990e-06, +7.431e-06, -2.796e-07,
        -1.601e-07, +4.068e-02, -1.352e-02, +7.282e-04, +9.594e-05,
        +2.070e-06, -9.620e-08, -2.742e-07, -6.370e-08, -6.300e-09,
        +8.625e-02, -5.971e-03, +4.705e-04, +2.335e-05, +4.226e-06,
        +2.475e-07, -8.850e-08, -3.600e-08, -2.900e-09, +0.000e+00
    };

    const double bh_amp[55] =
    {
        +0.000e+00, +0.000e+00, -1.136e-01, +0.000e+00, -1.868e-01,
        -1.399e-02, +0.000e+00, -1.043e-01, +1.175e-02, -2.240e-03,
        +0.000e+00, -3.222e-02, +1.333e-02, -2.647e-03, -2.316e-05,
        +0.000e+00, +5.339e-02, +1.107e-02, -3.116e-03, -1.079e-04,
        -1.299e-05, +0.000e+00, +4.861e-03, +8.891e-03, -6.448e-04,
        -1.279e-05, +6.358e-06, -1.417e-07, +0.000e+00, +3.041e-02,
        +1.150e-03, -8.743e-04, -2.781e-05, +6.367e-07, -1.140e-08,
        -4.200e-08, +0.000e+00, -2.982e-02, -3.000e-03, +1.394e-05,
        -3.290e-05, -1.705e-07, +7.440e-08, +2.720e-08, -6.600e-09,
        +0.000e+00, +1.236e-02, -9.981e-04, -3.792e-05, -1.355e-05,
        +1.162e-06, -1.789e-07, +1.470e-08, -2.400e-09, -4.000e-10
    };

    const double aw_mean[55] = {
        +5.640e+01, +1.555e+00, -1.011e+00, -3.975e+00, +3.171e-02,
        +1.065e-01, +6.175e-01, +1.376e-01, +4.229e-02, +3.028e-03,
        +1.688e+00, -1.692e-01, +5.478e-02, +2.473e-02, +6.059e-04,
        +2.278e+00, +6.614e-03, -3.505e-04, -6.697e-03, +8.402e-04,
        +7.033e-04, -3.236e+00, +2.184e-01, -4.611e-02, -1.613e-02,
        -1.604e-03, +5.420e-05, +7.922e-05, -2.711e-01, -4.406e-01,
        -3.376e-02, -2.801e-03, -4.090e-04, -2.056e-05, +6.894e-06,
        +2.317e-06, +1.941e+00, -2.562e-01, +1.598e-02, +5.449e-03,
        +3.544e-04, +1.148e-05, +7.503e-06, -5.667e-07, -3.660e-08,
        +8.683e-01, -5.931e-02, -1.864e-03, -1.277e-04, +2.029e-04,
        +1.269e-05, +1.629e-06, +9.660e-08, -1.015e-07, -5.000e-10
    };

    const double bw_mean[55] =
    {
        +0.000e+00, +0.000e+00, +2.592e-01, +0.000e+00, +2.974e-02,
        -5.471e-01, +0.000e+00, -5.926e-01, -1.030e-01, -1.567e-02,
        +0.000e+00, +1.710e-01, +9.025e-02, +2.689e-02, +2.243e-03,
        +0.000e+00, +3.439e-01, +2.402e-02, +5.410e-03, +1.601e-03,
        +9.669e-05, +0.000e+00, +9.502e-02, -3.063e-02, -1.055e-03,
        -1.067e-04, -1.130e-04, +2.124e-05, +0.000e+00, -3.129e-01,
        +8.463e-03, +2.253e-04, +7.413e-05, -9.376e-05, -1.606e-06,
        +2.060e-06, +0.000e+00, +2.739e-01, +1.167e-03, -2.246e-05,
        -1.287e-04, -2.438e-05, -7.561e-07, +1.158e-06, +4.950e-08,
        +0.000e+00, -1.344e-01, +5.342e-03, +3.775e-04, -6.756e-05,
        -1.686e-06, -1.184e-06, +2.768e-07, +2.730e-08, +5.700e-09
    };

    const double aw_amp[55] =
    {
        +1.023e-01, -2.695e+00, +3.417e-01, -1.405e-01, +3.175e-01,
        +2.116e-01, +3.536e+00, -1.505e-01, -1.660e-02, +2.967e-02,
        +3.819e-01, -1.695e-01, -7.444e-02, +7.409e-03, -6.262e-03,
        -1.836e+00, -1.759e-02, -6.256e-02, -2.371e-03, +7.947e-04,
        +1.501e-04, -8.603e-01, -1.360e-01, -3.629e-02, -3.706e-03,
        -2.976e-04, +1.857e-05, +3.021e-05, +2.248e+00, -1.178e-01,
        +1.255e-02, +1.134e-03, -2.161e-04, -5.817e-06, +8.836e-07,
        -1.769e-07, +7.313e-01, -1.188e-01, +1.145e-02, +1.011e-03,
        +1.083e-04, +2.570e-06, -2.140e-06, -5.710e-08, +2.000e-08,
        -1.632e+00, -6.948e-03, -3.893e-03, +8.592e-04, +7.577e-05,
        +4.539e-06, -3.852e-07, -2.213e-07, -1.370e-08, +5.800e-09
    };

    const double bw_amp[55] =
    {
        +0.000e+00, +0.000e+00, -8.865e-02, +0.000e+00, -4.309e-01,
        +6.340e-02, +0.000e+00, +1.162e-01, +6.176e-02, -4.234e-03,
        +0.000e+00, +2.530e-01, +4.017e-02, -6.204e-03, +4.977e-03,
        +0.000e+00, -1.737e-01, -5.638e-03, +1.488e-04, +4.857e-04,
        -1.809e-04, +0.000e+00, -1.514e-01, -1.685e-02, +5.333e-03,
        -7.611e-05, +2.394e-05, +8.195e-06, +0.000e+00, +9.326e-02,
        -1.275e-02, -3.071e-04, +5.374e-05, -3.391e-05, -7.436e-06,
        +6.747e-07, +0.000e+00, -8.637e-02, -3.807e-03, -6.833e-04,
        -3.861e-05, -2.268e-05, +1.454e-06, +3.860e-07, -1.068e-07,
        +0.000e+00, -2.658e-02, -1.947e-03, +7.131e-04, -3.506e-05,
        +1.885e-07, +5.792e-07, +3.990e-08, +2.000e-08, -5.700e-09
    };

    int i, j, k, m, n,ir;
    double p[100], t, ah, bh, ch, ap[55], bp[55], aw, bw, cw;
    double c0h, aha, c10h, c11h, ahm, awa, phh, awm, doy, sum1,
        dfac[20], beta, undu, a_ht, b_ht, c_ht, sine, ht_corr_coef,
        gamma, hs_km, topcon, ht_corr, tgmfh;

    /* +--------------------------------------------------------------------- */
    /*     Reference day is 28 January 1980 */
    /*     This is taken from Niell (1996) to be consistent */
    /* ---------------------------------------------------------------------- */
    doy = * mjd - 44239. - 27;
    /*     Define a parameter t */
    t = sin(*lat);
    /*     Define degree n and order m EGM */
    n = 9;
    m = 9;
    /*     Determine n!  (factorial)  moved by 1 */
    dfac[0] = 1.;
    for (i = 1; i <= 2 * n + 1; i++)
    {
        dfac[i] = dfac[i - 1] * i;
    }
    /*     Determine Legendre functions (Heiskanen and Moritz, */
    /*     Physical Geodesy, 1967, eq. 1-62) */
    for (i = 0; i <= n; i++)
    {
        for (j = 0; j <= MIN(i, m); j++)
        {
            ir = (i - j) / 2;
            sum1 = 0.;
            for (k = 0; k <= ir; k++)
            {
                sum1 += pow(-1, k) * dfac[2 * i - 2 * k] /
                    dfac[k] / dfac[i - k] / dfac[i - j - 2 * k] *
                    pow(t, i - j - 2 * k);
            }
            /*         Legendre functions moved by 1 */
            /* Computing 2nd power */
            p[i + 1 + (j + 1) * 10 - 11] = 1. / pow(2, i) * sqrt(
                pow(1 - t * t, j)) * sum1;
        }
    }
    /*     Calculate spherical harmonics */
    k = 0;
    for (i = 0; i <= n; i++)
    {
        for (j = 0; j <= i; j++)
        {
            ap[k] = p[i + 1 + (j + 1) * 10 - 11] * cos(j * (*lon));
            bp[k] = p[i + 1 + (j + 1) * 10 - 11] * sin(j * (*lon));
            k++;
        }
    }
    /*     Compute hydrostatic mapping function */
    bh = 0.0029;
    c0h = 0.062;
    if (*lat < 0.)
    {
        /* SOUTHERN HEMISPHERE */
        phh = PI;
        c11h = 0.007;
        c10h = 0.002;
    }
    else
    {
        /* NORTHERN HEMISPHERE */
        phh = 0.0;
        c11h = 0.005;
        c10h = 0.001;
    }
    ch = c0h + ((cos(doy / 365.25 * PI + phh) + 1.0) *
        c11h / 2.0 + c10h) * (1.0 - cos(*lat));
    ahm = 0.;
    aha = 0.;
    for (i = 0; i < 55; i++)
    {
        ahm += (ah_mean[i] * ap[i] + bh_mean[i] * bp[i]) * 1e-5;
        aha += (ah_amp[i] * ap[i] + bh_amp[i] * bp[i]) * 1e-5;
    }
    ah = ahm + aha * cos(doy / 365.25 * PI);
    sine = sin(PI / 2.0 - *zd);
    beta = bh / (sine + ch);
    gamma = ah / (sine + beta);
    topcon = ah / (bh / (ch + 1.0) + 1.0) + 1.0;
    tgmfh = topcon / (sine + gamma);

    /*     Height correction for hydrostatic mapping function from Niell (1996) */
    a_ht = 2.53e-5;
    b_ht = 0.00549;
    c_ht = 0.00114;
    hs_km = *hgt / 1e3;
    beta = b_ht / (sine + c_ht);
    gamma = a_ht / (sine + beta);
    topcon = a_ht / (b_ht / (c_ht + 1.0) + 1.0) + 1.0;
    ht_corr_coef = 1.0 / sine - topcon / (sine + gamma);
    ht_corr = ht_corr_coef * hs_km;
    tgmfh += ht_corr;

    /*     Compute wet mapping function */
    bw = 0.00146;
    cw = 0.04391;
    awm = 0.0;
    awa = 0.0;
    for (i = 0; i < 55; i++)
    {
        awm += (aw_mean[i] * ap[i] + bw_mean[i] * bp[i]) * 1e-5;
        awa += (aw_amp[i] * ap[i] + bw_amp[i] * bp[i]) * 1e-5;
    }
    aw = awm + awa * cos(doy / 365.25 * PI);
    beta = bw / (sine + cw);
    gamma = aw / (sine + beta);
    topcon = aw / (bw / (cw + 1.0) + 1.0) + 1.0;

    if (gmfw) *gmfw = topcon / (sine + gamma);
    if (gmfh) *gmfh = tgmfh;
}

#endif /*IERS_MODEL*/
