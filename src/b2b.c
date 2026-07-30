
#include "rtklib.h"

#define MaskBDS 63       /* Number of BDS masks*/
#define MaskGPS 37       /* Number of GPS masks*/
#define MaskGAL 37       /* Number of Galileo masks*/
#define MaskGLO 37       /* Number of GLONASS masks*/
#define MAXAGEB2b_ORB 140     /* Max age of B2b orbit corrections (seconds)*/
#define MAXAGEB2b_CBIAS 86400   /* Max age of B2b code bias corrections (seconds)*/
#define MAXAGEB2b_CLOCK 100    /* Max age of B2b clock corrections (seconds)*/
#define MU_GPS   3.9860050E14     /* gravitational constant         ref [1] */
#define MU_GLO   3.9860044E14     /* gravitational constant         ref [2] */
#define MU_GAL   3.986004418E14   /* earth gravitational constant   ref [7] */
#define MU_CMP   3.986004418E14   /* earth gravitational constant   ref [9] */
#define J2_GLO   1.0826257E-3     /* 2nd zonal harmonic of geopot   ref [2] */
#define OMGE_GLO 7.292115E-5      /* earth angular velocity (rad/s) ref [2] */
#define OMGE_GAL 7.2921151467E-5  /* earth angular velocity (rad/s) ref [7] */
#define OMGE_CMP 7.292115E-5      /* earth angular velocity (rad/s) ref [9] */
#define SIN_5 -0.0871557427476582 /* sin(-5.0 deg) */
#define COS_5  0.9961946980917456 /* cos(-5.0 deg) */
#define RTOL_KEPLER 1E-13         /* relative tolerance for Kepler equation */
#define MAX_ITER_KEPLER 30        /* max number of iteration of Kelpler */


#define SQR(x)      ((x)*(x))

#ifdef ENAB2B
/* Adjust time to be within ¡À3.5 days of reference epoch*/
extern gtime_t adjweek(gtime_t t, gtime_t t0) /*mdf*/
{
    double tt = timediff(t, t0);

    if (tt < -302400.0) return timeadd(t, 604800.0); /* Add one GPS week*/
    if (tt > 302400.0) return timeadd(t, -604800.0); /* Subtract one GPS week*/
    return t;
}


/* Convert satellite slot number to satellite number*/
extern int satSlot2Sat(int SatSlot)
{
    int sat, sys;

    if (SatSlot <= 63 && SatSlot >= 1) { /* BDS satellites (1-63)*/
        sys = SYS_CMP;
        sat = satno(sys, SatSlot);
    }
    else if (SatSlot >= 64 && SatSlot <= 100) { /* GPS satellites (64-100)*/
        sys = SYS_GPS;
        sat = satno(sys, SatSlot - MaskBDS);
    }
    else if (SatSlot >= 101 && SatSlot <= 137) { /*Galileo satellites (101-137)*/
        sys = SYS_GAL;
        sat = satno(sys, SatSlot - MaskBDS - MaskGPS);
    }
    else if (SatSlot >= 138 && SatSlot <= 174) { /* GLONASS satellites (138-174)*/
        sys = SYS_GLO;
        sat = satno(sys, SatSlot - MaskBDS - MaskGPS - MaskGLO);
    }
    else if (SatSlot >= 175 || SatSlot <= 0) { /* Invalid slot range*/
        sat = 0;
    }

    return sat;
}

/* Convert satellite number to satellite slot index*/
extern int sat2Slot(int sat)
{
    int prn = 0, sys = 0, SatSlot = 0;

    sys = satsys(sat, &prn); /* Extract system type and PRN*/

    if (sys == SYS_CMP) {
        SatSlot = prn;
    }
    else if (sys == SYS_GPS) {
        SatSlot = prn + MaskBDS;
    }
    else if (sys == SYS_GAL) {
        SatSlot = prn + MaskBDS + MaskGPS;
    }
    else if (sys == SYS_GLO) {
        SatSlot = prn + MaskBDS + MaskGPS + MaskGLO;
    }
    else {
        SatSlot = 0;
    }

    return SatSlot;
}


/* decode msg type 1 --------------------------------------------------*/
static int decodeB2bType1(const b2bmsg_t *b2bmsg, b2bsat_t *b2bsat)
{
    int i, type, pos,sod, iodssr,iodp,slot,count ;
    double ep[6];
    gtime_t msgtime, reftime;

    pos = 0;
    type = (int)getbitu(b2bmsg->msg, pos, 6);
    if (type != 1) return 0;
    pos += 6;

    sod = (int)getbitu(b2bmsg->msg, pos, 17);
    msgtime = gpst2time(b2bmsg->week, b2bmsg->tow);
    msgtime = gpst2bdt(msgtime);
    time2epoch(msgtime, ep);
    ep[3] = ep[4] = ep[5] = 0.0;
    reftime = epoch2time(ep);
    reftime = timeadd(reftime, (double)sod);
    reftime = adjweek(reftime, msgtime);
    reftime = bdt2gpst(reftime);
    pos += 17;

    pos += 4;

    iodssr = (int)getbitu(b2bmsg->msg, pos, 2);
    pos += 2;

    iodp = (int)getbitu(b2bmsg->msg, pos, 4);
    pos += 4;

    /*sat mask*/
    count = 0;
    for (i = 0; i < MaskBDS+ MaskGPS+ MaskGAL+ MaskGLO; i++)
    {
        slot = (int)getbitu(b2bmsg->msg, pos + i, 1);
        b2bsat->SatSlot[i] = slot;
        if (slot==1)
        {         
            b2bsat->b2bsats[count].sat = satSlot2Sat(i + 1);
            b2bsat->b2bsats[count].b2btype1.Iodp = iodp;
            b2bsat->b2bsats[count].b2btype1.IodSsr = iodssr;
            b2bsat->b2bsats[count].b2btype1.t0 = reftime;
            count++;
        }
    }
    b2bsat->nsat = count;

    return 1;
}

/* decode msg type 2 --------------------------------------------------*/
static int decodeB2bType2(const b2bmsg_t *b2bmsg, b2bsat_t *b2bsat)
{
    int i, type, pos,sod, iodssr, slotnum,iodn,iodcorr, cla, val,ind;
    double ep[6],deph[3];
    int sat,j;
    gtime_t msgtime, reftime;

    pos = 0;
    type = (int)getbitu(b2bmsg->msg, pos, 6);
    if (type != 2) return 0;
    pos += 6;

    sod = (int)getbitu(b2bmsg->msg, pos, 17);
    msgtime = gpst2time(b2bmsg->week, b2bmsg->tow);
    msgtime = gpst2bdt(msgtime);
    time2epoch(msgtime, ep);
    ep[3] = ep[4] = ep[5] = 0.0;
    reftime = epoch2time(ep);
    reftime = timeadd(reftime, (double)sod);
    reftime = adjweek(reftime, msgtime);
    reftime = bdt2gpst(reftime);
    pos += 17;

    pos += 4;

    iodssr = (int)getbitu(b2bmsg->msg, pos, 2);
    pos += 2;

    for (i = 0; i < 6; i++)
    {
        slotnum = (int)getbitu(b2bmsg->msg, pos, 9);
        sat = satSlot2Sat(slotnum);
        pos += 9;

        ind = -1;
        for (j = 0; j < b2bsat->nsat; j++) {
            if (sat != b2bsat->b2bsats[j].sat) continue;
            ind = j;
            break;
        }

        iodn = (int)getbitu(b2bmsg->msg, pos, 10);
        pos += 10;

        iodcorr = (int)getbitu(b2bmsg->msg, pos, 3);
        pos += 3;

        deph[0] = getbits(b2bmsg->msg, pos, 15)*0.0016;
        pos += 15;
        deph[1] = getbits(b2bmsg->msg, pos, 13)*0.0064;
        pos += 13;
        deph[2] = getbits(b2bmsg->msg, pos, 13)*0.0064;
        pos += 13;

        cla = (int)getbitu(b2bmsg->msg, pos, 3);
        pos += 3;
        val = (int)getbitu(b2bmsg->msg, pos, 3);
        pos += 3;
       
        
        if (sat > 0 && ind >= 0)
        {            
            b2bsat->b2bsats[ind].b2btype2.t0 = reftime;
            b2bsat->b2bsats[ind].b2btype2.IodSsr = iodssr;
            b2bsat->b2bsats[ind].b2btype2.IODN = iodn;
            b2bsat->b2bsats[ind].b2btype2.IodCorr = iodcorr;
            b2bsat->b2bsats[ind].b2btype2.OrbCorr[0] = deph[0];
            b2bsat->b2bsats[ind].b2btype2.OrbCorr[1] = deph[1];
            b2bsat->b2bsats[ind].b2btype2.OrbCorr[2] = deph[2];
            b2bsat->b2bsats[ind].b2btype2.UraClass = cla;
            b2bsat->b2bsats[ind].b2btype2.UraValue = val;            
        }
    }

    return 1;
}

/* decode msg type 3 --------------------------------------------------*/
static int decodeB2bType3(const b2bmsg_t *b2bmsg, b2bsat_t *b2bsat)
{
    int i, j, type, iodssr, pos, sat,  nsat,sod, slotnum, ncode, icode,ind;
    double ep[6],cbias;
    gtime_t msgtime,reftime;
  
    pos = 0;
    type = (int)getbitu(b2bmsg->msg, pos, 6);
    if (type != 3) return 0;
    pos += 6;

    sod = (int)getbitu(b2bmsg->msg, pos, 17);
    msgtime = gpst2time(b2bmsg->week, b2bmsg->tow);
    msgtime = gpst2bdt(msgtime);
    time2epoch(msgtime, ep);
    ep[3] = ep[4] = ep[5] = 0.0;
    reftime = epoch2time(ep);
    reftime = timeadd(reftime, (double)sod);
    reftime = adjweek(reftime, msgtime);
    reftime = bdt2gpst(reftime);
    pos += 17;

    pos += 4;

    iodssr = (int)getbitu(b2bmsg->msg, pos, 2);
    pos += 2;

    nsat = (int)getbitu(b2bmsg->msg, pos, 5);
    pos += 5;

    for (i = 0; i < nsat; i++)
    {
        slotnum = (int)getbitu(b2bmsg->msg, pos, 9);
        sat = satSlot2Sat(slotnum);        
        pos += 9;

        ind = -1;
        for (j = 0; j < b2bsat->nsat; j++) {
            if (sat != b2bsat->b2bsats[j].sat) continue;
            ind = j;
            break;
        }

        ncode = (int)getbitu(b2bmsg->msg, pos, 4);
        pos += 4;

        for (j = 0; j < ncode; j++)
        {
            icode = (int)getbitu(b2bmsg->msg, pos, 4);
            pos += 4;

            cbias = getbits(b2bmsg->msg, pos, 12)*0.017;
            pos += 12;

            if (sat>0&&ind >= 0)
            {                
                b2bsat->b2bsats[ind].b2btype3.t0 = reftime;
                b2bsat->b2bsats[ind].b2btype3.IodSsr = iodssr;
                b2bsat->b2bsats[ind].b2btype3.SatDCB[icode] = cbias;
            }
        }
    }

    return 1;
}

/* Get satellite number from subtype index*/
static int B2bSubtype2Sat(const int ind, b2bsat_t* b2bsat)
{
    int i, count = 0, sat = -1;

    for (i = 0; i < 255; i++) {
        if (b2bsat->SatSlot[i]) {
            count++;
        }
        if (count == ind + 1) {
            sat = satSlot2Sat(i + 1);
            return sat;
        }
    }
    return -1; /* Not found*/
}


/* decode msg type 4 --------------------------------------------------*/
static int decodeB2bType4(const b2bmsg_t *b2bmsg, b2bsat_t *b2bsat, int *mmi)
{
    int i,k, type, pos, sat,sod, iodssr,iodp, subtype, iodcorr,ind;
    double ep[6],c0;
    gtime_t msgtime,reftime;
    
    *mmi = 1;

    pos = 0;
    type = (int)getbitu(b2bmsg->msg, pos, 6);
    if (type != 4) return 0;
    pos += 6;

    sod = (int)getbitu(b2bmsg->msg, pos, 17);
    msgtime = gpst2time(b2bmsg->week, b2bmsg->tow);
    msgtime =gpst2bdt(msgtime);
    time2epoch(msgtime, ep);
    ep[3] = ep[4] = ep[5] = 0.0;
    reftime = epoch2time(ep);
    reftime = timeadd(reftime, (double)sod);
    reftime = adjweek(reftime, msgtime);
    reftime = bdt2gpst(reftime);
    pos += 17;

    pos += 4;

    iodssr = (int)getbitu(b2bmsg->msg, pos, 2);
    pos += 2;

    iodp = (int)getbitu(b2bmsg->msg, pos, 4);
    pos += 4;

    subtype = (int)getbitu(b2bmsg->msg, pos, 5);
    pos += 5;

    for (i = 0; i < 23&&i+ subtype * 23< b2bsat->nsat; i++)
    {
        sat = B2bSubtype2Sat(i + subtype*23, b2bsat);

        ind = -1;
        for (k = 0; k < b2bsat->nsat; k++) {
            if (sat != b2bsat->b2bsats[k].sat) continue;
            ind = k;
            break;
        }
       
        iodcorr = (int)getbitu(b2bmsg->msg, pos, 3);       
        pos += 3;

        c0= getbits(b2bmsg->msg, pos, 15)*0.0016;       
        pos += 15;

        if (sat>0&&ind >= 0)
        {
            b2bsat->b2bsats[ind].b2btype4.t0 = reftime;
            b2bsat->b2bsats[ind].b2btype4.Iodp = iodp;
            b2bsat->b2bsats[ind].b2btype4.IodSsr = iodssr;            
            b2bsat->b2bsats[ind].b2btype4.IodCorr = iodcorr;
            b2bsat->b2bsats[ind].b2btype4.C0 = c0;
        }
    }
    if (i + subtype * 23 >= b2bsat->nsat) *mmi = 0;

    return 1;
}

/* -----------------------------------------------------------------------------
* Name       : decodeB2bMessage
* Description: decode B2b message
* Parameters : *b2bmsg        b2bmsg_t       I       B2b message struct
*              *b2bsat        b2bsat_t       IO      B2b correction struct
* Return     : status(0:no clock updated; 13:clock updated)
* Notes      :
*              message format:
*              +------------+-----------------------------------------+----+
*              |    type    |          data message                   | 00 |
*              +------------+-----------------------------------------+----+
*              |<---- 6---->|<------------------ 456----------------->|  2 |
*-----------------------------------------------------------------------------*/
extern int decodeB2bMessage(const b2bmsg_t *b2bmsg, b2bsat_t *b2bsat)
{
    int i,j,type,stat, mmi,sync=0;
    gtime_t time;

    b2bsat->clkupdate = 0;
    
    type = (int)getbitu(b2bmsg->msg, 0, 6);

    if (b2bmsg->week == 0) return 0;
    time = gpst2time(b2bmsg->week, b2bmsg->tow);

    stat = 0; mmi = 1;
    switch (type)
    {
    case  1: stat = decodeB2bType1(b2bmsg, b2bsat);break;
    case  2: stat = decodeB2bType2(b2bmsg, b2bsat);break;
    case  3: stat = decodeB2bType3(b2bmsg, b2bsat);break;
    case  4: stat = decodeB2bType4(b2bmsg, b2bsat, &mmi);break;/*mdf*/
    case  5: stat = 0; break;
    case  6: stat = 0; break;
    case  7: stat = 0; break;
    case 63: stat = 1; break; /* null message */
    }

    sync = 0;

    if (type==4&&stat&&!mmi&&b2bsat->nsat>0)
    {
        for (i = 0; i < b2bsat->nsat; i++)
        {
            if (fabs(b2bsat->b2bsats[i].b2btype4.C0) < 26.212)
            {
                if (fabs(timediff(time, b2bsat->b2bsats[i].b2btype2.t0)) > 96 ||
                    fabs(timediff(time, b2bsat->b2bsats[i].b2btype4.t0)) > 9 ||
                    b2bsat->b2bsats[i].b2btype1.IodSsr != b2bsat->b2bsats[i].b2btype2.IodSsr ||
                    b2bsat->b2bsats[i].b2btype1.IodSsr != b2bsat->b2bsats[i].b2btype4.IodSsr ||
                    b2bsat->b2bsats[i].b2btype1.Iodp != b2bsat->b2bsats[i].b2btype4.Iodp ||
                    b2bsat->b2bsats[i].b2btype2.IodCorr != b2bsat->b2bsats[i].b2btype4.IodCorr)
                {
                    break;
                }
            }
        }
        if (i >= b2bsat->nsat)
        {
            sync = 1;
            b2bsat->clkupdate = 1;
        }       
    }
    return stat ? (sync?13:0) : 0;
}

/* Calculate variance from B2b URA value*/
static double varUraB2b(double* ura)
{
    double URA;

    URA = (pow(3, ura[0]) * (1 + 0.25 * ura[1]) - 1) * 1E-3; /* URA calculation formula*/
    return SQR(URA); /* Return variance (square of URA)*/
}



extern int satpos_B2b(gtime_t time, gtime_t teph, int sat, const nav_t* nav,
    double* rs, double* dts, double* var, int* svh)
{
    const b2bsatp_t* b2bsatp;
    b2bsat_t b2bsat = nav->b2bsat;
    eph_t* eph;
    double t1, t2, t4, er[3], ea[3], ec[3], rc[3], dant[3] = { 0 }, tk;
    int i, sys, slot, index, num, count = 0, IodCorr4 = -1, ind = -1;
    /*int maskSlot[MaskNSAT];*/

    /* only suport gps, galileo and bds */
    sys = satsys(sat, NULL);
    switch (sys) {
    case SYS_GPS:setseleph(SYS_GPS, 0); break;
    case SYS_GAL:setseleph(SYS_GAL, 0); break;
    case SYS_CMP:setseleph(SYS_CMP, 1); break;
    default:*svh = -1; return 0;
    }

    for (i = 0; i < nav->b2bsat.nsat; i++) {
        b2bsatp = nav->b2bsat.b2bsats + i;
        if (b2bsatp->sat == sat)break;
    }
    if (i >= nav->b2bsat.nsat) {

        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
        *svh = -1;
        return 0;
    }
    if (b2bsatp->b2btype4.Iodp != b2bsatp->b2btype1.Iodp) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
        *svh = -1;
        return 0;

    }
    double OrbCorr[3] = { 0 };
    OrbCorr[0] = b2bsatp->b2btype2.OrbCorr[0];
    OrbCorr[1] = b2bsatp->b2btype2.OrbCorr[1];
    OrbCorr[2] = b2bsatp->b2btype2.OrbCorr[2];
    double C0 = b2bsatp->b2btype4.C0;

    if (b2bsatp->b2btype1.IodSsr != b2bsatp->b2btype2.IodSsr &&
        b2bsatp->b2btype1.IodSsr != b2bsatp->b2btype4.IodSsr) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
        *svh = -1;
        return 0;
    }
    if (b2bsatp->b2btype2.IodCorr != b2bsatp->b2btype4.IodCorr) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
        *svh = -1;
        return 0;
    }

    if (fabs(OrbCorr[0]) > 26.2 || fabs(OrbCorr[1]) > 26.2 || fabs(OrbCorr[2]) > 26.2 || fabs(C0) > 26.2) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
        *svh = -1;
        return 0;
    }

    t2 = timediff(time, b2bsatp->b2btype2.t0);
    t4 = timediff(time, b2bsatp->b2btype4.t0);

    if (fabs(t2) > MAXAGEB2b_ORB || fabs(t4) > MAXAGEB2b_CLOCK) {
        ephpos(time, teph, sat, nav, -1, rs, dts, var, svh);
        *svh = -1;
        return 0;
    }

    /* satellite postion and clock by broadcast ephemeris */
    if (!ephpos(time, teph, sat, nav, b2bsatp->b2btype2.IODN, rs, dts, var, svh)) {
        *svh = -1;
        return 0;
    }


    /* Calculate radial, along-track, and cross-track unit vectors in ECEF */
    if (!normv3(rs + 3, ea)) {
        *svh = -1;
        return 0;
    }
    cross3(rs, rs + 3, rc);
    if (!normv3(rc, ec)) {
        *svh = -1;
        return 0;
    }
    cross3(ea, ec, er);

    for (i = 0; i < 3; i++) {
        rs[i] += -(er[i] * OrbCorr[0] + ea[i] * OrbCorr[1] + ec[i] * OrbCorr[2]);
    }

    if (fabs(C0) > 26.21) {
        C0 = 0;
    }
    dts[0] -= C0 / CLIGHT;

    /* Variance by B2b URA (User Range Accuracy) */
    double Ura[2] = { b2bsatp->b2btype2.UraClass,b2bsatp->b2btype2.UraValue };
    *var = varUraB2b(Ura);
    return 1;
}

extern double cbias_B2b(const obsd_t* obs, const nav_t* nav, int freqidx)
{
    int i,sat = obs->sat, sys;
    double t, dcb = 0.0;
    double freq[NFREQ] = { 0 };
    const b2bsatp_t *b2bsatp;
    
    if (!(sys = satsys(sat, NULL))) return 0.0;

    for (i = 0; i < nav->b2bsat.nsat; i++) {
        b2bsatp = nav->b2bsat.b2bsats + i;
        if (b2bsatp->sat == sat)break;
    }

    if (i >= nav->b2bsat.nsat) {
        return  0.0;
    }
    t = timediff(obs->time, b2bsatp->b2btype3.t0);
    if (fabs(t) > MAXAGEB2b_CBIAS) return  0.0;

    /* Check if the satellite¡¯s ionospheric correction (IodSsr) is valid*/
    if (b2bsatp->b2btype3.IodSsr != b2bsatp->b2btype2.IodSsr &&
        b2bsatp->b2btype3.IodSsr != b2bsatp->b2btype4.IodSsr &&
        b2bsatp->b2btype3.IodSsr != b2bsatp->b2btype1.IodSsr)
        return  0.0;

    /* Specific handling for GPS system to account for differential code bias (DCB)*/
    if (sys == SYS_GPS) {
        switch (obs->code[freqidx])
        {
        case CODE_L1C:dcb = b2bsatp->b2btype3.SatDCB[0]; break;/*L1C/A*/
        case CODE_L1W:
        case CODE_L1P: dcb = b2bsatp->b2btype3.SatDCB[1]; break;/*L1P*/
        case CODE_L1L: dcb = b2bsatp->b2btype3.SatDCB[4]; break;/*L1C(P)*/
        case CODE_L1X: dcb = b2bsatp->b2btype3.SatDCB[5]; break;/*L1C(D+P)*/
        case CODE_L2L: dcb = b2bsatp->b2btype3.SatDCB[7]; break;/*L2C(L)*/
        case CODE_L2X: dcb = b2bsatp->b2btype3.SatDCB[8]; break;/*L2C(M+L)*/
        case CODE_L5I: dcb = b2bsatp->b2btype3.SatDCB[11]; break;/*L5I*/
        case CODE_L5Q: dcb = b2bsatp->b2btype3.SatDCB[12]; break;/*L5Q*/
        case CODE_L5X: dcb = b2bsatp->b2btype3.SatDCB[12]; break;/*L5 I+Q*/
        default:dcb = 0.0;
            break;
        }
        return 0.0;
    }
    /* Specific handling for BeiDou (COMPASS) system*/
    else if (sys == SYS_CMP) {
        switch (obs->code[freqidx])
        {
        case CODE_L2I: dcb = b2bsatp->b2btype3.SatDCB[0]; break;/*B1I*/
        case CODE_L1D: dcb = b2bsatp->b2btype3.SatDCB[1]; break;/*B1C(D)*/
        case CODE_L1P: dcb = b2bsatp->b2btype3.SatDCB[2];; break;/*B1C(P)*/
        case CODE_L5D: dcb = b2bsatp->b2btype3.SatDCB[4]; break;/*B2a(D)*/
        case CODE_L5P: dcb = b2bsatp->b2btype3.SatDCB[5]; break;/*B2a(P)*/
        case CODE_L7I: dcb = b2bsatp->b2btype3.SatDCB[7]; break;/*B2b-I*/
        case CODE_L7Q: dcb = b2bsatp->b2btype3.SatDCB[8]; break;/*B2b-Q*/
        case CODE_L6I: dcb = b2bsatp->b2btype3.SatDCB[12]; break;/*B3I*/
        default:dcb = 0.0;
            break;
        }
        return 0.0;
    }
    else
        return 0.0;
}

#endif /*ENAB2B*/
