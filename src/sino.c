
#include "rtklib.h"

#define SINOSYNC1       0xAA    /* sino message start sync code 1 */
#define SINOSYNC2       0x44    /* sino message start sync code 2 */
#define SINOSYNC3       0x12    /* sino message start sync code 3 */
#define SINOHLEN        28      /* sino message header length (bytes) */

/* message IDs */
#define ID_RANGECMP     140     /* sino range compressed */
#define ID_RANGE        43      /* sino range measurement */
#define ID_RAWEPHEM     41      /* sino raw ephemeris */
#define ID_IONUTC       8       /* sino iono and utc data */
#define ID_GPSEPHEM     71      /* sino decoded gps/bd2 ephemeris */
#define ID_BD3EPHEM     72      /* sino decoded bd3 ephemeris */
#define ID_GLOEPHEMERIS 723     /* sino decoded glonass ephemeris */
#define ID_GALEPHEMERIS 1122    /* sino decoded galileo ephemeris */
#define ID_B2BCORR      1697    /* sino b2b correction data */
#define ID_HASCORR      1797    /* sino has correction data*/


#define WL1             0.1902936727984
#define WL2             0.2442102134246
#define MAXVAL          8388608.0
#define OFF_FRQNO       -7      /* F/W ver.3.620 */

#define SQR(x)          ((x)*(x))

/* get fields (little-endian) ------------------------------------------------*/
#define U1(p) (*((uint8_t *)(p)))
#define I1(p) (*((int8_t  *)(p)))
static uint16_t U2(uint8_t *p) {uint16_t u; memcpy(&u,p,2); return u;}
static int16_t  I2(uint8_t *p) {int16_t  u; memcpy(&u,p,2);return u;}
static uint32_t U4(uint8_t *p) {uint32_t u; memcpy(&u,p,4); return u;}
static int32_t  I4(uint8_t *p) {int32_t  i; memcpy(&i,p,4); return i;}
static float    R4(uint8_t *p) {float    r; memcpy(&r,p,4); return r;}
static double   R8(uint8_t *p) {double   r; memcpy(&r,p,8); return r;}

/* extend sign ---------------------------------------------------------------*/
static int32_t exsign(uint32_t v, int bits)
{
    return (int32_t)(v&(1<<(bits-1))?v|(~0u<<bits):v);
}
/* checksum ------------------------------------------------------------------*/
static uint8_t chksum(const uint8_t *buff, int len)
{
    uint8_t sum=0;
    int i;
    for (i=0;i<len;i++) sum^=buff[i];
    return sum;
}
/* adjust weekly rollover of GPS time ----------------------------------------*/
static gtime_t adjweek(gtime_t time, double tow)
{
    double tow_p;
    int week;
    tow_p=time2gpst(time,&week);
    if      (tow<tow_p-302400.0) tow+=604800.0;
    else if (tow>tow_p+302400.0) tow-=604800.0;
    return gpst2time(week,tow);
}
/* UTC 8-bit week -> full week -----------------------------------------------*/
static void adj_utcweek(gtime_t time, double *utc)
{
    int week;
    
    time2gpst(time,&week);
    utc[3]+=week/256*256;
    if      (utc[3]<week-127) utc[3]+=256.0;
    else if (utc[3]>week+127) utc[3]-=256.0;
    utc[5]+=utc[3]/256*256;
    if      (utc[5]<utc[3]-127) utc[5]+=256.0;
    else if (utc[5]>utc[3]+127) utc[5]-=256.0;
}
/* get observation data index ------------------------------------------------*/
static int obsindex(obs_t *obs, gtime_t time, int sat)
{
    int i,j;
    
    if (obs->n>=MAXOBS) return -1;
    for (i=0;i<obs->n;i++) {
        if (obs->data[i].sat==sat) return i;
    }
    obs->data[i].time=time;
    obs->data[i].sat=sat;
    for (j=0;j<NFREQ+NEXOBS;j++) {
        obs->data[i].L[j]=obs->data[i].P[j]=0.0;
        obs->data[i].D[j]=0.0;
        obs->data[i].SNR[j]=obs->data[i].LLI[j]=0;
        obs->data[i].code[j]=CODE_NONE;
    }
    obs->n++;
    return i;
}
/* URA value (m) to URA index ------------------------------------------------*/
static int uraindex(double value)
{
    static const double ura_eph[]={
        2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,
        3072.0,6144.0,0.0
    };
    int i;
    for (i=0;i<15;i++) if (ura_eph[i]>=value) break;
    return i;
}
/* signal type to obs code ---------------------------------------------------*/
static int sig2code(int sys, int sigtype)
{
    if (sys==SYS_GPS) {
        switch (sigtype) {
            case  0: return CODE_L1C; /* L1C/A */
            case  2: return CODE_L5C; /* L5C   */
            case  5: return CODE_L2P; /* L2P   */
            case  9: return CODE_L2W; /* L2P(Y),semi-codeless */
            case 14: return CODE_L5Q; /* L5Q    */
            case 16: return CODE_L1L; /* L1C(P) */
            case 17: return CODE_L2S; /* L2C(M) */
        }
    }
    else if (sys==SYS_GLO) {
        switch (sigtype) {
            case  0: return CODE_L1C; /* L1C/A */
            case  1: return CODE_L2C; /* L2C/A */
            case  5: return CODE_L2P; /* L2P */
            case  6: return CODE_L3Q; /* L3Q */
        }
    }
    else if (sys==SYS_GAL) {
        switch (sigtype) {
            case  1: return CODE_L1B; /* E1B  */
            case  2: return CODE_L1C; /* E1C  */
            case  6: return CODE_L6B; /* E6B  */
            case  7: return CODE_L6C; /* E6C  */
            case 12: return CODE_L5Q; /* E5aQ */
            case 17: return CODE_L7Q; /* E5bQ */
            case 20: return CODE_L8Q; /* AltBOCQ*/
        }
    }
    else if (sys==SYS_QZS) {
        switch (sigtype) {
            case  0: return CODE_L1C; /* L1C/A */
            case 14: return CODE_L5Q; /* L5Q   */
            case 16: return CODE_L1L; /* L1C(P)*/
            case 17: return CODE_L2S; /* L2C(M)*/
            case 27: return CODE_L6L; /* L6P   */
        }
    }
    else if (sys==SYS_CMP) {
        switch (sigtype) {
            case  0: return CODE_L2I; /* B1I with D1 */
            case  1: return CODE_L7I; /* B2I with D1 */
            case  2: return CODE_L6I; /* B3I with D1 */
            case  8: return CODE_L1P; /* B1C(P) */
            case 12: return CODE_L5P; /* B2a(P) */
            case 19: return CODE_L7D; /* B2b(I) */
        }
    }
    else if (sys==SYS_IRN) {
        switch (sigtype) {
            case  0: return CODE_L5A; /* L5  */
        }
    }
    else if (sys==SYS_SBS) {
        switch (sigtype) {
            case  0: return CODE_L1C; /* L1C/A */
            case  6: return CODE_L5I; /* L5I */
        }
    }
    return 0;
}
/* decode receiver tracking status ---------------------------------------------
* decode receiver tracking status
* args   : uint32_t stat I  tracking status field
*          int    *sys   O      system (SYS_???)
*          int    *code  O      signal code (CODE_L??)
*          int    *track O      tracking state
*                         0=idle                      7=freq-lock loop
*                         2=wide freq band pull-in    9=channel alignment
*                         3=narrow freq band pull-in 10=code search
*                         4=phase lock loop          11=aided phase lock loop
*          int    *plock O      phase-lock flag   (0=not locked, 1=locked)
*          int    *clock O      code-lock flag    (0=not locked, 1=locked)
*          int    *parity O     parity known flag (0=not known,  1=known)
*          int    *halfc O      phase measurement (0=half-cycle not added,
*                                                  1=added)
* return : freq-index (-1:error)
*-----------------------------------------------------------------------------*/
static int decode_track_stat(uint32_t stat, int *sys, int *code, int *track,
                             int *plock, int *clock, int *parity, int *halfc)
{
    int satsys,sigtype,idx=-1;
    
    *code=CODE_NONE;
    *track =stat&0x1F;
    *plock =(stat>>10)&1;
    *parity=(stat>>11)&1;
    *clock =(stat>>12)&1;
    satsys =(stat>>16)&7;
    *halfc =(stat>>28)&1;
    sigtype=(stat>>21)&0x1F;
    
    switch (satsys) {
        case 0: *sys=SYS_GPS; break;
        case 1: *sys=SYS_GLO; break;
        case 2: *sys=SYS_SBS; break;
        case 3: *sys=SYS_GAL; break; 
        case 4: *sys=SYS_CMP; break; 
        case 5: *sys=SYS_QZS; break; 
        case 6: *sys=SYS_IRN; break; 
        default:
            trace(2,"sino unknown system: sys=%d\n",satsys);
            return -1;
    }
    if (!(*code=sig2code(*sys,sigtype))||(idx=code2idx(*sys,*code))<0) {
        trace(2,"sino signal type error: sys=%d sigtype=%d\n",*sys,sigtype);
        return -1;
    }
    return idx;
}
/* check code priority and return freq-index ---------------------------------*/
static int checkpri(const char *opt, int sys, int code, int idx)
{
    int nex=NEXOBS;
    
    if (sys==SYS_GPS) {
        if (strstr(opt,"-GL1L")&&idx==0) return (code==CODE_L1L)?0:-1;
        if (strstr(opt,"-GL2S")&&idx==1) return (code==CODE_L2X)?1:-1;
        if (strstr(opt,"-GL2P")&&idx==1) return (code==CODE_L2P)?1:-1;
        if (code==CODE_L1L) return (nex<1)?-1:NFREQ;
        if (code==CODE_L2S) return (nex<2)?-1:NFREQ+1;
        if (code==CODE_L2P) return (nex<3)?-1:NFREQ+2;
    }
    else if (sys==SYS_GLO) {
        if (strstr(opt,"-RL2C")&&idx==1) return (code==CODE_L2C)?1:-1;
        if (code==CODE_L2C) return (nex<1)?-1:NFREQ;
    }
    else if (sys==SYS_GAL) {
        if (strstr(opt,"-EL6B")&&idx==3) return (code==CODE_L6B)?3:-1;
        if (code==CODE_L6B) return (nex<2)?-1:NFREQ;
    }
    else if (sys==SYS_QZS) {
        if (strstr(opt,"-JL1L")&&idx==0) return (code==CODE_L1L)?0:-1;
        if (strstr(opt,"-JL1Z")&&idx==0) return (code==CODE_L1Z)?0:-1;
        if (code==CODE_L1L) return (nex<1)?-1:NFREQ;
        if (code==CODE_L1Z) return (nex<2)?-1:NFREQ+1;
    }
    else if (sys==SYS_CMP) {
        if (strstr(opt,"-CL1P")&&idx==0) return (code==CODE_L1P)?0:-1;
        if (strstr(opt,"-CL7D")&&idx==0) return (code==CODE_L7D)?0:-1;
        if (code==CODE_L1P) return (nex<1)?-1:NFREQ;
        if (code==CODE_L7D) return (nex<2)?-1:NFREQ+1;
    }
    return idx<NFREQ?idx:-1;
}
/* decode RANGECMPB ----------------------------------------------------------*/
static int decode_rangecmpb(raw_t *raw)
{
    uint8_t *p=raw->buff+SINOHLEN;
    char *q;
    double psr,adr,adr_rolls,lockt,tt,dop,snr,freq,glo_bias=0.0;
    int i,index,nobs,prn,sat,sys,code,idx,track,plock,clock,parity,halfc,lli;
    
    if ((q=strstr(raw->opt,"-GLOBIAS="))) sscanf(q,"-GLOBIAS=%lf",&glo_bias);
    
    nobs=U4(p);
    if (raw->len<SINOHLEN+4+nobs*24) {
        trace(2,"oem4 rangecmpb length error: len=%d nobs=%d\n",raw->len,nobs);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype+strlen(raw->msgtype)," nobs=%d",nobs);
    }
    for (i=0,p+=4;i<nobs;i++,p+=24) {
        if ((idx=decode_track_stat(U4(p),&sys,&code,&track,&plock,&clock,
                                   &parity,&halfc))<0) {
            continue;
        }
        prn=U1(p+17);
        if (sys==SYS_GLO) prn-=37;
        if (sys == SYS_CMP) prn -= 140;
        if (sys == SYS_QZS) prn += 62;
        if (sys==SYS_SBS&&prn>=MINPRNQZS_S&&prn<=MAXPRNQZS_S&&code==CODE_L1C) {
            sys=SYS_QZS;
            prn+=10;
            code=CODE_L1Z; /* QZS L1S */
        }
        if (!(sat=satno(sys,prn))) {
            trace(3,"sino rangecmpb satellite number error: sys=%d,prn=%d\n",sys,prn);
            continue;
        }
        if (sys==SYS_GLO&&!parity) continue; /* invalid if GLO parity unknown */
        
        if ((idx=checkpri(raw->opt,sys,code,idx))<0) continue;
        
        dop=exsign(U4(p+4)&0xFFFFFFF,28)/256.0;
        psr=(U4(p+7)>>4)/128.0+U1(p+11)*2097152.0;
        
        if ((freq=sat2freq(sat,(uint8_t)code,&raw->nav))!=0.0) {
            adr=I4(p+12)/256.0;
            adr_rolls=(psr*freq/CLIGHT+adr)/MAXVAL;
            adr=-adr+MAXVAL*floor(adr_rolls+(adr_rolls<=0?-0.5:0.5));
            if (sys==SYS_GLO) adr+=glo_bias*freq/CLIGHT;
        }
        else {
            adr=1e-9;
        }
        lockt=(U4(p+18)&0x1FFFFF)/32.0; /* lock time */
        
        if (raw->tobs[sat-1][idx].time!=0) {
            tt=timediff(raw->time,raw->tobs[sat-1][idx]);
            lli=(lockt<65535.968&&lockt-raw->lockt[sat-1][idx]+0.05<=tt)?LLI_SLIP:0;
        }
        else {
            lli=0;
        }
        if (!parity) lli|=LLI_HALFC;
        if (halfc  ) lli|=LLI_HALFA;
        raw->tobs [sat-1][idx]=raw->time;
        raw->lockt[sat-1][idx]=lockt;
        raw->halfc[sat-1][idx]=halfc;
        
        snr=((U2(p+20)&0x3FF)>>5)+20.0;
        if (!clock) psr=0.0;     /* code unlock */
        if (!plock) adr=dop=0.0; /* phase unlock */
        
        if (fabs(timediff(raw->obs.data[0].time,raw->time))>1E-9) {
            raw->obs.n=0;
        }
        if ((index=obsindex(&raw->obs,raw->time,sat))>=0) {
            raw->obs.data[index].L  [idx]=adr;
            raw->obs.data[index].P  [idx]=psr;
            raw->obs.data[index].D  [idx]=(float)dop;
            raw->obs.data[index].SNR[idx]=(uint16_t)(snr/SNR_UNIT+0.5);
            raw->obs.data[index].LLI[idx]=(uint8_t)lli;
            raw->obs.data[index].code[idx]=(uint8_t)code;
        }
    }
    return 1;
}
/* decode RANGEB -------------------------------------------------------------*/
static int decode_rangeb(raw_t *raw)
{
    uint8_t *p=raw->buff+SINOHLEN;
    char *q;
    double psr,adr,dop,snr,lockt,tt,freq,glo_bias=0.0;
    int i,index,nobs,prn,sat,sys,code,idx,track,plock,clock,parity,halfc,lli;
    int gfrq;
    
    if ((q=strstr(raw->opt,"-GLOBIAS="))) sscanf(q,"-GLOBIAS=%lf",&glo_bias);
    
    nobs=U4(p);
    if (raw->len<SINOHLEN+4+nobs*44) {
        trace(2,"sino rangeb length error: len=%d nobs=%d\n",raw->len,nobs);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype+strlen(raw->msgtype)," nobs=%d",nobs);
    }
    for (i=0,p+=4;i<nobs;i++,p+=44) {
        if ((idx=decode_track_stat(U4(p+40),&sys,&code,&track,&plock,&clock,&parity,&halfc))<0) {
            continue;
        }
        prn=U2(p);
        if (sys==SYS_GLO) prn-=37;
        if (sys == SYS_CMP) prn -= 140;
        if (sys == SYS_QZS) prn += 62;
        if (sys==SYS_SBS&&prn>=MINPRNQZS_S&&prn<=MAXPRNQZS_S&&code==CODE_L1C) {
            sys=SYS_QZS;
            prn+=10;
            code=CODE_L1Z; /* QZS L1S */
        }
        if (!(sat=satno(sys,prn))) {
            trace(3,"sino rangeb satellite number error: sys=%d,prn=%d\n",sys,prn);
            continue;
        }
        if (sys==SYS_GLO&&!parity) continue;
        
        if ((idx=checkpri(raw->opt,sys,code,idx))<0) continue;
        
        gfrq =U2(p+ 2); /* GLONASS Frequency+7*/
        psr  =R8(p+ 4);
        adr  =R8(p+16);
        dop  =R4(p+28);
        snr  =R4(p+32);
        lockt=R4(p+36);
        
        if (sys==SYS_GLO) {
            freq=sat2freq(sat,(uint8_t)code,&raw->nav);
            adr-=glo_bias*freq/CLIGHT;
            if (!raw->nav.glo_fcn[prn-1]) {
                raw->nav.glo_fcn[prn-1]=gfrq; /* Frequency+7 */
            }
        }
        if (raw->tobs[sat-1][idx].time!=0) {
            tt=timediff(raw->time,raw->tobs[sat-1][idx]);
            lli=lockt-raw->lockt[sat-1][idx]+0.05<=tt?LLI_SLIP:0;
        }
        else {
            lli=0;
        }
        if (!parity) lli|=LLI_HALFC;
        if (halfc  ) lli|=LLI_HALFA;
        raw->tobs [sat-1][idx]=raw->time;
        raw->lockt[sat-1][idx]=lockt;
        raw->halfc[sat-1][idx]=halfc;
        
        if (!clock) psr=0.0;     /* code unlock */
        if (!plock) adr=dop=0.0; /* phase unlock */
        
        if (fabs(timediff(raw->obs.data[0].time,raw->time))>1E-9) {
            raw->obs.n=0;
        }
        if ((index=obsindex(&raw->obs,raw->time,sat))>=0) {
            raw->obs.data[index].L  [idx]=-adr;
            raw->obs.data[index].P  [idx]=psr;
            raw->obs.data[index].D  [idx]=(float)dop;
            raw->obs.data[index].SNR[idx]=(uint16_t)(snr/SNR_UNIT+0.5);
            raw->obs.data[index].LLI[idx]=(uint8_t)lli;
            raw->obs.data[index].code[idx]=(uint8_t)code;
        }
    }
    return 1;
}
/* decode RAWEPHEMB ----------------------------------------------------------*/
static int decode_rawephemb(raw_t *raw)
{
    eph_t eph={0};
    uint8_t *p=raw->buff+SINOHLEN,subframe[30*5]={0};
    int prn,sat;
    
    if (raw->len<SINOHLEN+102) {
        trace(2,"sino rawephemb length error: len=%d\n",raw->len);
        return -1;
    }
    prn=U4(p);
    if (!(sat=satno(SYS_GPS,prn))) {
        trace(2,"sino rawephemb satellite number error: prn=%d\n",prn);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype+strlen(raw->msgtype)," prn=%d",prn);
    }
    memcpy(subframe,p+12,30*3); /* subframe 1-3 */
    
    if (!decode_frame(subframe,&eph,NULL,NULL,NULL)) {
        trace(2,"sino rawephemb subframe error: prn=%d\n",prn);
        return -1;
    }
    if (!strstr(raw->opt,"-EPHALL")) {
        if (eph.iode==raw->nav.eph[sat-1].iode&&
            eph.iodc==raw->nav.eph[sat-1].iodc) return 0;
    }
    eph.sat=sat;
    raw->nav.eph[sat-1]=eph;
    raw->ephsat=sat;
    raw->ephset=0;
    return 2;
}


/* decode IONUTCB ------------------------------------------------------------*/
static int decode_ionutcb(raw_t *raw)
{
    uint8_t *p=raw->buff+SINOHLEN;
    int i;
    
    if (raw->len<SINOHLEN+108) {
        trace(2,"sino ionutcb length error: len=%d\n",raw->len);
        return -1;
    }
    for (i=0;i<8;i++) raw->nav.ion_gps[i]=R8(p+i*8);
    raw->nav.utc_gps[0]=R8(p+ 72); /* A0 */
    raw->nav.utc_gps[1]=R8(p+ 80); /* A1 */
    raw->nav.utc_gps[2]=U4(p+ 68); /* tot */
    raw->nav.utc_gps[3]=U4(p+ 64); /* WNt */
    raw->nav.utc_gps[4]=I4(p+ 96); /* dt_LS */
    raw->nav.utc_gps[5]=U4(p+ 88); /* WN_LSF */
    raw->nav.utc_gps[6]=U4(p+ 92); /* DN */
    raw->nav.utc_gps[7]=I4(p+100); /* dt_LSF */
    return 9;
}

/* decode BD3EPHEMB ------------------------------------------------------*/
static int decode_bd3ephemb(raw_t *raw)
{
    const double ArefMEO_BDS = 27906100;     /* Reference semi-major axis (A) for BDS MEO satellites*/
    const double ArefIGSOGEO_BDS = 42162200; /* Reference semi-major axis (A) for BDS GEO/IGSO satellites*/
    eph_t eph = { 0 };
    uint8_t *p = raw->buff + SINOHLEN;
    int prn, sat,week,toc,urai;
    double sow,delA;
    gtime_t time;

    if (raw->len < SINOHLEN + 208) {
        trace(2, "sino bd3ephemb length error: len=%d\n", raw->len);
        return -1;
    }
    prn = U1(p);   p += 1;
    p += 1;
    eph.satType= U1(p);   p += 1;
    eph.svh = U1(p) & 1; p += 1;
    urai= U1(p);   p += 1;
    eph.iode = U1(p);   p += 1; /* IODE */
    eph.iodc = U1(p);   p += 1; /* IODC*/
    p += 9;
    eph.toes = U4(p);   p += 4;
    toc = U4(p);        p += 4;
    delA = R8(p);  p += 8;
    eph.Adot = R8(p);   p += 8;
    eph.deln= R8(p);    p += 8;
    eph.ndot = R8(p);   p += 8;
    eph.M0 = R8(p);   p += 8;
    eph.e = R8(p);   p += 8;
    eph.omg = R8(p);   p += 8;
    eph.OMG0 = R8(p);   p += 8;
    eph.i0 = R8(p);     p += 8;
    eph.OMGd = R8(p);   p += 8;
    eph.idot = R8(p);   p += 8;
    eph.cuc = R8(p);   p += 8;
    eph.cus = R8(p);   p += 8;
    eph.crc = R8(p);   p += 8;
    eph.crs = R8(p);   p += 8;
    eph.cic = R8(p);   p += 8;
    eph.cis = R8(p);   p += 8;
    eph.f0 = R8(p);   p += 8;
    eph.f1 = R8(p);   p += 8;
    eph.f2 = R8(p);   p += 8;
    eph.tgd[0] = R8(p);   p += 8; /* TGD1 B1Cp (s) */
    eph.tgd[1] = R8(p);   p += 8; /* TGD2 B2ap (s) */
    eph.tgd[2] = R8(p);   p += 8; /* TGD2 b1cd (s) */

    if (!(sat = satno(SYS_CMP, prn))) {
        trace(2, "sino bd3ephemb satellite error: prn=%d\n", prn);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype + strlen(raw->msgtype), " prn=%d", prn);
    }
    eph.sat = sat;
    if(eph.satType<3)
    {
        eph.A = ArefIGSOGEO_BDS+ delA;
    }
    else
    {
        eph.A = ArefMEO_BDS + delA;
    }
    eph.sva = uraindex(urai);
    time = gpst2bdt(raw->time);
    sow=time2bdt(time, &week);
    if (eph.toes - sow > 320400) week = week - 1;
    else if (eph.toes - sow < -320400) week = week + 1;
    eph.week = week;
    eph.toe = bdt2gpst(bdt2time(eph.week, eph.toes));
    eph.toc = bdt2gpst(bdt2time(eph.week, toc));
    eph.toc = adjweek(eph.toc, eph.toes);
    eph.ttr = raw->time;
    eph.Enmt = 3;

    if (!strstr(raw->opt, "-EPHALL")) {
        if (timediff(eph.toe, raw->nav.eph[sat - 1 + MAXSAT].toe) == 0.0&&
            timediff(eph.toc, raw->nav.eph[sat - 1 + MAXSAT].toc) == 0.0) {
            return 0; /* unchanged */
        }
    }
    raw->nav.eph[sat - 1 + MAXSAT * 1] = eph;
    raw->ephsat = sat;
    raw->ephset = 1;/* 0:D1/D2,1:CNAV */
    return 2;
}

/* decode GPSEPHEMB ------------------------------------------------------*/
static int decode_gpsephemb(raw_t *raw)
{
    eph_t eph = { 0 };
    uint8_t *p = raw->buff + SINOHLEN;
    double sqrtA,toc,toe,tow,gtow;
    int prn,sat,sys,blFlag, bHealth, urai,gweek;

    if (raw->len < SINOHLEN + 224) {
        trace(2, "sino gpsephemb length error: len=%d\n", raw->len);
        return -1;
    }
    p += 2;
    blFlag= U1(p);   p += 1;
    bHealth = U1(p);   p += 1;
    eph.svh = bHealth;
    prn = U1(p);   p += 1;
    if (prn >= 1&&prn<=32)
    {
        sys = SYS_GPS;
        if (!(sat = satno(SYS_GPS, prn))) {
            trace(2, "sino gpsephemb satellite error: prn=%d\n", prn);
            return -1;
        }
        
    }
    else if (prn >= 141&&prn<=203)
    {
        prn -= 140;
        sys = SYS_CMP;
        if (!(sat = satno(SYS_CMP, prn))) {
            trace(2, "sino gpsephemb satellite error: prn=%d\n", prn);
            return -1;
        }
    }
    else
    {
        trace(2, "sino gpsephemb satellite error: prn=%d\n", prn);
        return -1;
    }
    p += 5;
    eph.iodc = I2(p);   p += 2; /* IODC */
    urai= U2(p);   p += 2;
    eph.week = U2(p);   p += 2;
    eph.iode = I4(p);   p += 4; /* IODE */
    tow = I4(p);   p += 4; /* tow */
    eph.toes = R8(p);   p += 8;
    toc = R8(p);   p += 8;
    eph.f2 = R8(p);   p += 8;
    eph.f1 = R8(p);   p += 8;
    eph.f0 = R8(p);   p += 8;
    eph.M0 = R8(p);   p += 8;
    eph.deln = R8(p);   p += 8;
    eph.e = R8(p);   p += 8;
    sqrtA = R8(p);   p += 8;
    eph.OMG0 = R8(p);   p += 8;
    eph.i0 = R8(p);   p += 8;
    eph.omg = R8(p);   p += 8;
    eph.OMGd = R8(p);   p += 8;
    eph.idot = R8(p);   p += 8;
    eph.cuc = R8(p);   p += 8;
    eph.cus = R8(p);   p += 8;
    eph.crc = R8(p);   p += 8;
    eph.crs = R8(p);   p += 8;
    eph.cic = R8(p);   p += 8;
    eph.cis = R8(p);   p += 8;
    eph.tgd[0] = R8(p);   p += 8; /* TGD1 */
    eph.tgd[1] = R8(p);   p += 8; /* TGD2 */

    if (raw->outtype) {
        sprintf(raw->msgtype + strlen(raw->msgtype), " prn=%d", prn);
    }
    eph.sat = sat;
    eph.A = SQR(sqrtA);
    eph.sva = urai;
    if (sys == SYS_CMP)
    {
        eph.toe = bdt2gpst(bdt2time(eph.week, eph.toes));
        eph.toc = bdt2gpst(bdt2time(eph.week, toc));
        eph.toc = adjweek(eph.toc, eph.toes);
        eph.ttr = bdt2gpst(bdt2time(eph.week, tow));
        eph.ttr = adjweek(eph.ttr, eph.toes);
    }   
    else
    {
        eph.week = adjgpsweek(eph.week);
        eph.toe = gpst2time(eph.week, eph.toes);
        eph.toc = gpst2time(eph.week, toc);
        eph.toc = adjweek(eph.toc, eph.toes);
        eph.ttr = gpst2time(eph.week, tow);
        eph.ttr = adjweek(eph.ttr, eph.toes);
    }

    if (!strstr(raw->opt, "-EPHALL")) {
        if (timediff(raw->nav.eph[sat - 1].toe, eph.toe) == 0.0&&
            timediff(raw->nav.eph[sat - 1].toc, eph.toc) == 0.0) return 0;
    }
    raw->nav.eph[sat - 1] = eph;
    raw->ephsat = sat;
    raw->ephset = 0;
    return 2;
}


/* decode GLOEPHEMERISB ------------------------------------------------------*/
static int decode_gloephemerisb(raw_t *raw)
{
    uint8_t *p=raw->buff+SINOHLEN;
    geph_t geph={0};
    double tow,tof,toff;
    int prn,sat,week;
    
    if (raw->len<SINOHLEN+144) {
        trace(2,"sino gloephemerisb length error: len=%d\n",raw->len);
        return -1;
    }
    prn=U2(p)-37;
    
    if (!(sat=satno(SYS_GLO,prn))) {
        trace(2,"sino gloephemerisb prn error: prn=%d\n",prn);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype+strlen(raw->msgtype)," prn=%d",prn);
    }
    geph.frq   =U2(p+  2)+OFF_FRQNO;
    week       =U2(p+  6);
    tow        =floor(U4(p+8)/1000.0+0.5); /* rounded to integer sec */
    toff       =U4(p+ 12);
    geph.iode  =U4(p+ 20)&0x7F;
    geph.svh   =(U4(p+24)<4)?0:1; /* 0:healthy,1:unhealthy */
    geph.pos[0]=R8(p+ 28);
    geph.pos[1]=R8(p+ 36);
    geph.pos[2]=R8(p+ 44);
    geph.vel[0]=R8(p+ 52);
    geph.vel[1]=R8(p+ 60);
    geph.vel[2]=R8(p+ 68);
    geph.acc[0]=R8(p+ 76);
    geph.acc[1]=R8(p+ 84);
    geph.acc[2]=R8(p+ 92);
    geph.taun  =R8(p+100);
    geph.dtaun =R8(p+108);
    geph.gamn  =R8(p+116);
    tof        =U4(p+124)-toff; /* glonasst->gpst */
    geph.age   =U4(p+136);
    geph.toe=gpst2time(week,tow);
    tof+=floor(tow/86400.0)*86400;
    if      (tof<tow-43200.0) tof+=86400.0;
    else if (tof>tow+43200.0) tof-=86400.0;
    geph.tof=gpst2time(week,tof);
    
    if (!strstr(raw->opt,"-EPHALL")) {
        if (fabs(timediff(geph.toe,raw->nav.geph[prn-1].toe))<1.0&&
            geph.svh==raw->nav.geph[prn-1].svh) return 0; /* unchanged */
    }
    geph.sat=sat;
    raw->nav.geph[prn-1]=geph;
    raw->ephsat=sat;
    raw->ephset=0;
    return 2;
}

/* decode GALEPHEMERISB ------------------------------------------------------*/
static int decode_galephemerisb(raw_t *raw)
{
    eph_t eph={0};
    uint8_t *p=raw->buff+SINOHLEN;
    double tow,sqrtA,af0_fnav,af1_fnav,af2_fnav,af0_inav,af1_inav,af2_inav,tt;
    int prn,sat,week,rcv_fnav,rcv_inav,svh_e1b,svh_e5a,svh_e5b,dvs_e1b,dvs_e5a;
    int dvs_e5b,toc_fnav,toc_inav,set,sel_eph=3; /* 1:I/NAV+2:F/NAV */
    
    if (strstr(raw->opt,"-GALINAV")) sel_eph=1;
    if (strstr(raw->opt,"-GALFNAV")) sel_eph=2;
    
    if (raw->len<SINOHLEN+220) {
        trace(2,"sino galephemrisb length error: len=%d\n",raw->len);
        return -1;
    }
    prn       =U4(p);   p+=4;
    rcv_fnav  =U4(p)&1; p+=4;
    rcv_inav  =U4(p)&1; p+=4;
    svh_e1b   =U1(p)&3; p+=1;
    svh_e5a   =U1(p)&3; p+=1;
    svh_e5b   =U1(p)&3; p+=1;
    dvs_e1b   =U1(p)&1; p+=1;
    dvs_e5a   =U1(p)&1; p+=1;
    dvs_e5b   =U1(p)&1; p+=1;
    eph.sva   =U1(p); /*;*/
    if (eph.sva <= 49) { eph.sva= eph.sva * 0.01;}
    else if (eph.sva <= 74)  {eph.sva= 0.5 + (eph.sva - 50) * 0.02;}
    else if (eph.sva <= 99) { eph.sva= 1.0 + (eph.sva - 75) * 0.04;}
    else if (eph.sva <= 125) {eph.sva= 2.0 + (eph.sva - 100) * 0.16;}

    p+=1+1; /* SISA index */
    eph.iode  =U4(p);   p+=4;   /* IODNav */
    eph.toes  =U4(p);   p+=4;
    sqrtA     =R8(p);   p+=8;
    eph.deln  =R8(p);   p+=8;
    eph.M0    =R8(p);   p+=8;
    eph.e     =R8(p);   p+=8;
    eph.omg   =R8(p);   p+=8;
    eph.cuc   =R8(p);   p+=8;
    eph.cus   =R8(p);   p+=8;
    eph.crc   =R8(p);   p+=8;
    eph.crs   =R8(p);   p+=8;
    eph.cic   =R8(p);   p+=8;
    eph.cis   =R8(p);   p+=8;
    eph.i0    =R8(p);   p+=8;
    eph.idot  =R8(p);   p+=8;
    eph.OMG0  =R8(p);   p+=8;
    eph.OMGd  =R8(p);   p+=8;
    toc_fnav  =U4(p);   p+=4;
    af0_fnav  =R8(p);   p+=8;
    af1_fnav  =R8(p);   p+=8;
    af2_fnav  =R8(p);   p+=8;
    toc_inav  =U4(p);   p+=4;
    af0_inav  =R8(p);   p+=8;
    af1_inav  =R8(p);   p+=8;
    af2_inav  =R8(p);   p+=8;
    eph.tgd[0]=R8(p);   p+=8; /* BGD: E5A-E1 (s) */
    eph.tgd[1]=R8(p);         /* BGD: E5B-E1 (s) */
    
    if (!(sat=satno(SYS_GAL,prn))) {
        trace(2,"sino galephemeris satellite error: prn=%d\n",prn);
        return -1;
    }
    if (raw->outtype) {
        sprintf(raw->msgtype+strlen(raw->msgtype)," prn=%d",prn);
    }
    set=rcv_fnav?1:0; /* 0:I/NAV,1:F/NAV */
    if (!(sel_eph&1)&&set==0) return 0;
    if (!(sel_eph&2)&&set==1) return 0;
    
    eph.sat =sat;
    eph.A   =SQR(sqrtA);
    eph.f0  =set?af0_fnav:af0_inav;
    eph.f1  =set?af1_fnav:af1_inav;
    eph.f2  =set?af2_fnav:af2_inav;
    eph.svh =((svh_e5b<<7)|(dvs_e5b<<6)|(svh_e5a<<4)|(dvs_e5a<<3)|
             (svh_e1b<<1)|dvs_e1b);
    eph.code=set?((1<<1)+(1<<8)):((1<<0)+(1<<2)+(1<<9));
    eph.iodc=eph.iode;
    tow=time2gpst(raw->time,&week);
    eph.week=week; /* gps-week = gal-week */
    eph.toe=gpst2time(eph.week,eph.toes);
    
    tt=timediff(eph.toe,raw->time);
    if      (tt<-302400.0) eph.week++;
    else if (tt> 302400.0) eph.week--;
    eph.toe=gpst2time(eph.week,eph.toes);
    eph.toc=adjweek(raw->time,set?toc_fnav:toc_inav);
    eph.ttr=raw->time;

    if (!strstr(raw->opt,"-EPHALL")) {
        if (eph.iode==raw->nav.eph[sat-1+MAXSAT*set].iode&&
            timediff(eph.toe,raw->nav.eph[sat-1+MAXSAT*set].toe)==0.0&&
            timediff(eph.toc,raw->nav.eph[sat-1+MAXSAT*set].toc)==0.0) {
            return 0; /* unchanged */
        }
    }
    raw->nav.eph[sat-1+MAXSAT*set]=eph;
    raw->ephsat=sat;
    raw->ephset=set;
    return 2;
}

#ifdef ENAB2B
/* decode Sino message -----------------------------------------*/
static int decode_b2bcorrb(raw_t *raw)
{
    b2bmsg_t b2bmsg = { 0 };
    int i,week, prn,type;
    double tow;       
    uint8_t msg[58]; 
    int stat;

    tow=time2gpst(raw->time,&week);
    b2bmsg.week = week;
    b2bmsg.tow = tow;
    if (raw->b2bmsg.week == week && fabs(raw->b2bmsg.tow - tow) < 0.5) return 0;/*mdf*/
    b2bmsg.prn = (int)getbitu(raw->buff, SINOHLEN * 8 + 4 * 8, 6);
    for (i = 0; i < 58; i++)
    {
        b2bmsg.msg[i] = (uint8_t)getbitu(raw->buff, SINOHLEN * 8 + 4 * 8 + 12 + i * 8, 8);
    }
    b2bmsg.msg[57] = (uint8_t)(b2bmsg.msg[57] & 0xfc);
    type = (int)getbitu(b2bmsg.msg, 0, 6);
    b2bmsg.type = type;
    raw->b2bmsg = b2bmsg;

    stat=decodeB2bMessage(&b2bmsg, &raw->b2bsat);

    return stat;
}
#endif /* ENAB2B*/

#ifdef ENAHAS
static int decode_hascorrb(raw_t* raw)
{
    haspage_t haspage = { 0 };
    int  week;
    double tow;
    int stat;

    tow = time2gpst(raw->time, &week);
    haspage.week = week;
    haspage.tow = tow;
    if (raw->haspage.week == week && fabs(raw->haspage.tow - tow) < 0.5) return 0;

    if(raw->len- SINOHLEN < 64)  return 0;
    memcpy(haspage.msg, raw->buff + SINOHLEN,64*sizeof(uint8_t));

    stat = decodeHasPage(&haspage,raw->haspageset,&raw->hassat);
    return stat;
}
#endif /*ENAHAS*/
/* decode Sino message -----------------------------------------*/
static int decode_sino(raw_t *raw)
{
    double tow;
    char tstr[32];
    int week,type=U2(raw->buff+4);
    
    trace(3,"decode_sino: type=%3d len=%d\n",type,raw->len);
    
    /* check crc32 */
    if (rtk_crc32(raw->buff,raw->len)!=U4(raw->buff+raw->len)) {
        trace(2,"sino crc error: type=%3d len=%d\n",type,raw->len);
        return -1;
    }
    week=U2(raw->buff+14);
    
    if (type==20||week==0) {
        trace(3,"sino time error: type=%3d week=%d\n",type,week);
        return 0;
    }
    week=adjgpsweek(week);
    tow =U4(raw->buff+16)*0.001;
    raw->time=gpst2time(week,tow);

    if (raw->outtype) {
        time2str(gpst2time(week,tow),tstr,2);
        sprintf(raw->msgtype,"SINO %4d (%4d): %s",type,raw->len,tstr);
    }
    switch (type) {
        case ID_RANGECMP       : return decode_rangecmpb       (raw);
        case ID_RANGE          : return decode_rangeb          (raw);
        case ID_RAWEPHEM       : return decode_rawephemb       (raw);
        case ID_IONUTC         : return decode_ionutcb         (raw);
        case ID_GPSEPHEM       : return decode_gpsephemb       (raw);
        case ID_BD3EPHEM       : return decode_bd3ephemb       (raw);
        case ID_GLOEPHEMERIS   : return decode_gloephemerisb   (raw);
        case ID_GALEPHEMERIS   : return decode_galephemerisb   (raw);  
#ifdef ENAB2B
        case ID_B2BCORR        : return decode_b2bcorrb        (raw);
#endif /*ENAB2B*/

#ifdef ENAHAS
        case ID_HASCORR        : return decode_hascorrb        (raw);
#endif /*ENAHAS*/
    }
    return 0;
}

/* sync header ---------------------------------------------------------------*/
static int sync_sino(uint8_t *buff, uint8_t data)
{
    buff[0]=buff[1]; buff[1]=buff[2]; buff[2]=data;
    return buff[0]==SINOSYNC1&&buff[1]== SINOSYNC2&&buff[2]== SINOSYNC3;
}

/* input Sino raw data from stream -------------------------------
* fetch next Sino raw data and input a mesasge from stream
* args   : raw_t *raw       IO  receiver raw data control struct
*          uint8_t data     I   stream data (1 byte)
* return : status (-1: error message, 0: no message, 1: input observation data,
*                  2: input ephemeris, 9: input ion/utc parameter
*                  )
*
* notes  : to specify input options for oem4, set raw->opt to the following
*          option strings separated by spaces.
*
*          -EPHALL : input all ephemerides
*          -GL1L   : select 1L for GPS L1 (default 1C)
*          -GL2S   : select 2S for GPS L2 (default 2W)
*          -GL2P   : select 2P for GPS L2 (default 2W)
*          -RL2C   : select 2C for GLO G2 (default 2P)
*          -EL6B   : select 6B for GAL E6 (default 6C)
*          -JL1L   : select 1L for QZS L1 (default 1C)
*          -JL1Z   : select 1Z for QZS L1 (default 1C)
*          -CL1P   : select 1P for BDS B1 (default 2I)
*          -CL7D   : select 7D for BDS B2 (default 7I)
*          -GALINAV: select I/NAV for Galileo ephemeris (default: all)
*          -GALFNAV: select F/NAV for Galileo ephemeris (default: all)
*          -GLOBIAS=bias: GLONASS code-phase bias (m)
*-----------------------------------------------------------------------------*/
extern int input_sino(raw_t *raw, uint8_t data)
{
    trace(5,"input_sino: data=%02x\n",data);
    
    /* synchronize frame */
    if (raw->nbyte==0) {
        if (sync_sino(raw->buff,data)) raw->nbyte=3;
        return 0;
    }
    raw->buff[raw->nbyte++]=data;
    
    if (raw->nbyte==10&&(raw->len=U2(raw->buff+8)+SINOHLEN)>MAXRAWLEN-4) {
        trace(2,"sino length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (raw->nbyte<10||raw->nbyte<raw->len+4) {
        /*printf("iptsino0\n");*/
        return 0;
    }
    raw->nbyte=0;
    
    return decode_sino(raw);
}

/* input Sino raw data from file ---------------------------------
* fetch next Sino raw data and input a message from file
* args   : raw_t  *raw      IO  receiver raw data control struct
*          FILE   *fp       I   file pointer
* return : status(-2: end of file, -1...9: same as above)
*-----------------------------------------------------------------------------*/
extern int input_sinof(raw_t *raw, FILE *fp)
{
    int i,data;
    
    trace(4,"input_sinof:\n");
    
    /* synchronize frame */
    if (raw->nbyte==0) {
        for (i=0;;i++) {
            if ((data=fgetc(fp))==EOF) return -2;
            if (sync_sino(raw->buff,(uint8_t)data)) break;
            if (i>=4096) return 0;
        }
    }
    if (fread(raw->buff+3,7,1,fp)<1) return -2;
    raw->nbyte=10;
    
    if ((raw->len=U2(raw->buff+8)+SINOHLEN)>MAXRAWLEN-4) {
        trace(2,"sino length error: len=%d\n",raw->len);
        raw->nbyte=0;
        return -1;
    }
    if (fread(raw->buff+10,raw->len-6,1,fp)<1) return -2;
    raw->nbyte=0;
    
    /* decode Sino message */
    return decode_sino(raw);
}
