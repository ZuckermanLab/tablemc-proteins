#include <cstdio>
#include <cstdlib>
#include "solvation.h"
#include "ffield.h"
#include "util.h"
#include "mc.h"
#include "rotations.h"

//this code adapted from Sundar's MRMC code
//int             numOfGBCubicSplinePoints;
#define         GB_NUM_SPLINE       500
double          gbCubicSplineParams[4*GB_NUM_SPLINE];
double          gbCubicSplineInterval;
double          gbCubicSplineIntervalInv;
void setup_gb_spline(double kappa)
{
    int i;
    double f;
    double x[GB_NUM_SPLINE],y[GB_NUM_SPLINE],h[GB_NUM_SPLINE],a[GB_NUM_SPLINE],b[GB_NUM_SPLINE],c[GB_NUM_SPLINE],d[GB_NUM_SPLINE];
    gbCubicSplineInterval=10.0/GB_NUM_SPLINE;
    gbCubicSplineIntervalInv=1/gbCubicSplineInterval;
    for(i=0;i<=(GB_NUM_SPLINE-1);i++)
    {
        x[i] = i*gbCubicSplineInterval;
        //y[i] = 1.0/(sqrt((x[i]*x[i])+exp(-(x[i]*x[i])/GB_KAPPA)));
        y[i] = 1.0/(sqrt(x[i]+exp(-x[i]/kappa))); //this way, we apply to r^2/alpha_i alpha_j)
    }
    for(i=0;i<=(GB_NUM_SPLINE-2);i++){h[i] = x[i+1] - x[i];}//setup matrix h
    //setup matrix elements of tri-diagonal matrix
    for(i=1;i<=(GB_NUM_SPLINE-2);i++)
    {
        a[i] = h[i-1];
        b[i] = 2*((h[i]+h[i-1]));
        c[i] = h[i];
        d[i] = 6*((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1]);
    }
    //run tri-diagonal matrix algorithm
    for(i=2;i<=(GB_NUM_SPLINE-2);i++)
    {
        f = a[i]/b[i-1];
        b[i] -= f*c[i-1];
        d[i] -= f*d[i-1];
    }
    //back substitution
    d[GB_NUM_SPLINE-2] = d[GB_NUM_SPLINE-2]/b[GB_NUM_SPLINE-2];//reuse d as the second derivative
    for(i=(GB_NUM_SPLINE-3);i>=1;i--){d[i] = (d[i]-c[i]*d[i+1])/b[i];}
    //setup end conditions
    d[0] = 0;
    d[GB_NUM_SPLINE-1] = 0;
    //calculate cubic spline parameters
    for(i=0;i<=(GB_NUM_SPLINE-2);i++)
    {
        gbCubicSplineParams[4*i]   = (d[i+1]-d[i])/(6*h[i]);                          // this is parameter a
        gbCubicSplineParams[4*i+1] = d[i]/2;                                          // parameter b
        gbCubicSplineParams[4*i+2] = (y[i+1]-y[i])/h[i] - (2*h[i]*d[i]+h[i]*d[i+1])/6;// c
        gbCubicSplineParams[4*i+3] = y[i];                                            // d
    }
}

double splintGB(double x)
{
    int k;
    double xs,xs2,y;
    k = (int)(x*gbCubicSplineIntervalInv);
    //justin adds this to cover our rear end if k is outside the range
    if (k>=GB_NUM_SPLINE-2) return 1/sqrt(x); //spline function goes to 1/x as x -> infinity
    xs = x - k*gbCubicSplineInterval;
    xs2 = xs*xs;
    y = gbCubicSplineParams[4*k]*xs2*xs + gbCubicSplineParams[4*k+1]*xs2 + gbCubicSplineParams[4*k+2]*xs + gbCubicSplineParams[4*k+3];
    return y;
}

//same thing but includes derivative
void splintGB_diff(const double x, double * y, double * dy)
{
    int k;
    double xs,xs2;
    k = (int)(x*gbCubicSplineIntervalInv);
    //justin adds this to cover our rear end if k is outside the range
    if (k>=GB_NUM_SPLINE-2) {
        *y=1/sqrt(x);
        *dy=-((*y)/(2*x)); //-(1/2)x^(-3/2) = x^(-1/2) * (-1/(2x))
        return;
    }
    xs = x - k*gbCubicSplineInterval;
    xs2 = xs*xs;
    *y = gbCubicSplineParams[4*k]*xs2*xs + gbCubicSplineParams[4*k+1]*xs2 + gbCubicSplineParams[4*k+2]*xs + gbCubicSplineParams[4*k+3];
    *dy = 3*gbCubicSplineParams[4*k]*xs2 + 2*gbCubicSplineParams[4*k+1]*xs + gbCubicSplineParams[4*k+2];
    //return y;
}

//formulas are from esolv1.f in tinker
double topology::dipolar_gb_energy(gb_param_info * gb_params, double * per_frag_born_radii, double * center, double * orient)
{
    int ifrag,jfrag,k,l;
    double * dipoles;
    double qi,qj,aiaj, cutoff2, egbtot, rij[3],r2,r,sum1, sum2, sum3, sum4;
    double gf2, gf, gf3, gf5, expterm, expc1, a00,a10,a20,a01,a11,gc2thru4[3],gu1[3],gux[3][3],term,term1,term2, esym, ewi, ewj,egb;
#ifdef TIMERS
    switch_timer(TIMER_GB_POLAR);
#endif
    //to keep as close as possible to the fortran, the zero elements in the "g" arrays are not used
    egbtot=0.0;
    cutoff2=gb_params->rgbmax*gb_params->rgbmax;
    //compute the dipole moment for every fragment in the frame of reference of the simulation
    dipoles=(double *) checkalloc(3*nfrag,sizeof(double));
    for (ifrag=0; ifrag<nfrag; ifrag++) rotate_vector_by_quat(&orient[4*ifrag],&fragtypes[frags[ifrag].type]->dipole[0],&dipoles[3*ifrag]);
    for (ifrag=0; ifrag<nfrag; ifrag++)
        for (jfrag=ifrag; jfrag<nfrag; jfrag++) if (fragtypes[frags[ifrag].type]->is_charged && fragtypes[frags[jfrag].type]->is_charged) {
                //need to skip if either fragment is nonpolar (no charge and no dipole)
                for (k=0; k<3; k++) rij[k]=center[3*jfrag+k]-center[3*ifrag+k];
                r2=rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
                if (r2>cutoff2) continue;
                qi=fragtypes[frags[ifrag].type]->qtot;
                qj=fragtypes[frags[jfrag].type]->qtot;
                aiaj=per_frag_born_radii[ifrag]*per_frag_born_radii[jfrag];
                expterm=exp(-r2/(gb_params->kappa*aiaj));
                gf2=1/(r2+aiaj*expterm); //eventually replace with spline
                gf=sqrt(gf2);
                gf3=gf*gf2;
                gf5=gf3*gf2;
                a00 = gf;
                a10 = -gf3;
                a20 = 3.0 * gf5;
                a01 = expc1 * a10;
                a11 = expc1 * a20;
                a00 = gb_params->c0 * a00;
                a01 = gb_params->c0 * a01;
                a10 = gb_params->c1 * a10;
                a11 = gb_params->c1 * a11;
                //gux[1]=xr*a10, etc.
                for (k=0; k<3; k++) gu1[k]=rij[k]*a10;
                //gc[2]=xr*a01, etc.
                for (k=0; k<3; k++) gc2thru4[k]=rij[k]*a01;
                //gux[2]=a10*xr2*a11, similarly for guy[3] and guz[4]
                for (k=0; k<3; k++) gux[k][k]=a10+rij[k]*rij[k]*a11;
                //gux[3] = xr*yr*a11.  First index replaces x,y,z second replaces 2 thru 4
                for (k=0; k<3; k++) for (l=k+1; l<3; l++) {
                    gux[k][l]=rij[k]*rij[l]*a11;
                    //gux[l][k]=gux[k][l];
                }
                 /*esym = ci * ck * gc(1)
     &                     - (uxi*(uxk*gux(2)+uyk*guy(2)+uzk*guz(2))
     &                       +uyi*(uxk*gux(3)+uyk*guy(3)+uzk*guz(3))
     &                       +uzi*(uxk*gux(4)+uyk*guy(4)+uzk*guz(4)))*/
                esym = qi*qj*a00;
                //it's a symmetric sum, so can double up when l!=k
                for (k=0; k<3; k++) for (l=k; l<3; l++) {
                    term=dipoles[3*ifrag+k]*dipoles[3*jfrag+l]*gux[k][l];
                    if (l!=k) term*=2.0;
                    esym-=term;
                }
                /*ewi = ci*(uxk*gc(2)+uyk*gc(3)+uzk*gc(4))
     &                -ck*(uxi*gux(1)+uyi*guy(1)+uzi*guz(1))*/
                ewi=0.0;
                for (k=0; k<3; k++) {
                    term1=qi*dipoles[3*jfrag+k]*gc2thru4[k];
                    term2=qj*dipoles[3*ifrag+k]*gu1[k];
                    ewi+=(term1-term2);
                }
                /*ewk = ci*(uxk*gux(1)+uyk*guy(1)+uzk*guz(1))
     &                -ck*(uxi*gc(2)+uyi*gc(3)+uzi*gc(4))*/
                ewj=0.0;
                for (k=0; k<3; k++) {
                    term1=qi*dipoles[3*jfrag+k]*gu1[k];
                    term2=qj*dipoles[3*ifrag+k]*gc2thru4[k];
                    ewj+=(term1-term2);
                }
                egb=esym+0.5*(ewi+ewj);
                //scale by 1/2 if it is a self energy
                if (ifrag==jfrag) egb=0.5*egb;
                egbtot+=egb;
        }
    free(dipoles);
    return egbtot;
}
//the below code is taken from the "eff.c" and "gbsa.c" codes that Ramu sent me (ultimately from NAB in AmberTools

   /*FGB taylor coefficients follow */
   /* from A to H :                 */
   /* 1/3 , 2/5 , 3/7 , 4/9 , 5/11  */
   /* 4/3 , 12/5 , 24/7 , 40/9 , 60/11 */

#define TA 0.33333333333333333333
#define TB 0.4
#define TC 0.42857142857142857143
#define TD 0.44444444444444444444
#define TDD 0.45454545454545454545

#define TE 1.33333333333333333333
#define TF 2.4
#define TG 3.42857142857142857143
#define TH 4.44444444444444444444
#define THH 5.45454545454545454545
//The solvation radii for the fragments will defined such that the Generalized Kirkwood expression equals the PB solvation energy.
#define GBOFFSET 0.0
//#define GB_OVERLAP_FACTOR 0.80
void topology::calculate_born_radii(gb_param_info * params, double * center, double * per_frag_born_radii)
{
    double rgbmax1i, rgbmax2i, rgbmaxpsmax2, rgbmax;
    double ri, ri1i, sumi, r2, dij, dij1i, dij2i, sj, sj2, uij, tmpsd, dumbo,theta, psi;
    int i,j,k;
#ifdef TIMERS
    switch_timer(TIMER_BORN_RADII);
#endif
    //int gb = params->gbtype;
    rgbmax = params->rgbmax;
    rgbmax1i = 1.0 / params->rgbmax;
    rgbmax2i = rgbmax1i * rgbmax1i;
    rgbmaxpsmax2 = params->rgbmaxpsmax2;

    //need to use the born formula for

    for (i=0; i<nfrag; i++) {
        //uncharged fragments have no defined born radius
        if (!fragtypes[frags[i].type]->is_charged) {
            per_frag_born_radii[i]=0.0;
            continue;
        }
        //ri = rborn[i] - gboffset;
        ri = fragtypes[frags[i].type]->solvation_radius - params->radii_offset;
        ri1i = 1. / ri;
        sumi = 0.0;

        for (j=i+1; j<nfrag; j++) {
            r2 = 0.0;
            for (k=0; k<3; k++) r2+=(center[3*j+k]-center[3*i+k])*(center[3*j+k]-center[3*i+k]);



         /* Select atom j from the pair list.  Non-graceful error handling. */

         /*for (k = 0; k < lpears[i] + upears[i]; k++) {

            if (pearlist[i] == NULL) {
               fprintf(nabout,
                       "NULL pair list entry in egb loop 1, taskid = %d\n",
                       mytaskid);
               fflush(nabout);
            }
            j = pearlist[i][k];*/




            if (r2 > rgbmaxpsmax2)
               continue;
            dij1i = 1.0 / sqrt(r2);
            dij = r2 * dij1i;
            //gb_OVERLAP_FACTOR was atoms[i].fs  -- use the volume radius here, since this is part of volume excluion
            sj = params->radii_scale_factor * (fragtypes[frags[j].type]->volume_radius - params->radii_offset);
            sj2 = sj * sj;

            /*
             * ---following are from the Appendix of Schaefer and Froemmel,
             * JMB 216:1045-1066, 1990;  Taylor series expansion for d>>s
             * is by Andreas Svrcek-Seiler; smooth rgbmax idea is from
             * Andreas Svrcek-Seiler and Alexey Onufriev.
             */

            if (dij > rgbmax + sj)
               continue;

            if ((dij > rgbmax - sj)) {
               uij = 1. / (dij - sj);
               sumi -= 0.125 * dij1i * (1.0 + 2.0 * dij * uij +
                                        rgbmax2i * (r2 -
                                                    4.0 * rgbmax *
                                                    dij - sj2) +
                                        2.0 * log((dij - sj) * rgbmax1i));

            } else if (dij > 4.0 * sj) {
               dij2i = dij1i * dij1i;
               tmpsd = sj2 * dij2i;
               dumbo =
                   TA + tmpsd * (TB +
                                 tmpsd * (TC +
                                          tmpsd * (TD + tmpsd * TDD)));
               sumi -= sj * tmpsd * dij2i * dumbo;

            } else if (dij > ri + sj) {
               sumi -= 0.5 * (sj / (r2 - sj2) +
                              0.5 * dij1i * log((dij - sj) / (dij + sj)));

            } else if (dij > fabs(ri - sj)) {
               theta = 0.5 * ri1i * dij1i * (r2 + ri * ri - sj2);
               uij = 1. / (dij + sj);
               sumi -= 0.25 * (ri1i * (2. - theta) - uij +
                               dij1i * log(ri * uij));

            } else if (ri < sj) {
               sumi -= 0.5 * (sj / (r2 - sj2) + 2. * ri1i +
                              0.5 * dij1i * log((sj - dij) / (sj + dij)));

            }

            /*if (gb == 7 || gb == 8) {
               if (dij < rborn[i] + rborn[j] + GBNECKCUT) {
                  mdist = dij - neckMaxPos[NeckIdx[i]][NeckIdx[j]];
                  mdist2 = mdist * mdist;
                  mdist6 = mdist2 * mdist2 * mdist2;
                  neck = neckMaxVal[NeckIdx[i]][NeckIdx[j]] /
                                    (1.0 + mdist2 + 0.3 * mdist6);
                  sumi -= gbneckscale * neck;
               }
            } // gb == 7 || gb == 8*/

         } // end loop over pairlist

         /*if (gb == 1) {

            // "standard" (HCT) effective radii:
            reff[i] = 1.0 / (ri1i + sumi);
            if (reff[i] < 0.0)
               reff[i] = 30.0;

         } else {*/

            // "gbao" formulas:

            psi = -ri * sumi;
            per_frag_born_radii[i] = 1.0 / (ri1i - tanh((params->alpha - params->beta * psi +
                                          params->gamma * psi * psi) * psi) / fragtypes[frags[i].type]->solvation_radius);
         //}

//         if (gb_debug)
#ifdef DEBUG_NON_TABULATED
            fprintf(stdout, "%d\t%15.7f\t%15.7f\n", i + 1, rborn[i],
                    reff[i]);
#endif
      }
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
}

void topology::read_gb_params(char * line, gb_param_info * params, double cutoff)
{
    int itype;
    double max_vol_radius,solv_radius,volume,vol_radius;
    char fname[255],fileline[255],fragname[MAX_FRAGMENT_NAME];
    FILE * input;
    //order: radii filename, epsp, epsw, kappa, sasa_coef, alpha, beta, gamma,  radii_scale_factor, radii_offset
    //try 1 78.5 4 0.015 0.8 0 2.91 1.0 0.0 for starters
    //the alpha, beta, gamma are GB-OBC I
    params->rgbmax = cutoff;
    sscanf(line,"%s %lg %lg %lg %lg %lg %lg %lg %lg %lg",fname,&params->epsp,&params->epsw,&params->kappa,&params->sasa_coef,
        &params->alpha,&params->beta,&params->gamma,&params->radii_scale_factor,&params->radii_offset);
    //now need to initialize c0, c1, and rgbmaxpsmax2
    //see eq. 19, Schneiders and Ponder paper
    params->c0=COUL_CONST*1.0*(params->epsp-params->epsw)/(params->epsp*(1.0*params->epsw+0.0*params->epsp));
    params->c1=COUL_CONST*2.0*(params->epsp-params->epsw)/(params->epsp*(2.0*params->epsw+1.0*params->epsp));
    max_vol_radius=0.0;
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("could not open radii file %s\n",fname);
        die();
    }
    while (!feof(input)) {
        fgets(fileline,sizeof(fileline),input);
        sscanf(fileline,"%s %lg %lg\n",fragname,&solv_radius,&volume);
        itype=fragtypebyname(fragname);
        if (itype<0) continue; /*{
            printf("undefined type %s\n",itype);
            die();
        }*/
        if (fragtypes[itype]->is_charged) {
            fragtypes[itype]->solvation_radius = solv_radius;
        } else {
            fragtypes[itype]->solvation_radius = 0.0;
        }
        fragtypes[itype]->volume = volume;
        //we are in initialization, so ok to use slow pow function.
        vol_radius=pow((3/(4*M_PI))*volume,1.0/3.0);
        fragtypes[itype]->volume_radius = vol_radius;
        if (vol_radius > max_vol_radius) max_vol_radius=vol_radius;
    }
    params->rgbmaxpsmax2 = params->rgbmax + params->radii_scale_factor*(max_vol_radius-params->radii_offset);
    params->rgbmaxpsmax2 = params->rgbmaxpsmax2 * params->rgbmaxpsmax2;
}

void topology::print_gb_info(gb_param_info * params)
{
    int itype;
    printf("Solute and solvent dielectric constant:  %.1f %.1f\n",params->epsp,params->epsw);
    printf("Kappa in Still formula:                  %.1f\n",params->kappa);
    //printf("SASA coefficient:                        %.3f kcal/(mol A^2)\n",params->sasa_coef);
    printf("SASA not yet available.\n");
    printf("Alpha, beta, gamma:                      %.2f %.2f %.2f\n",params->alpha,params->beta,params->gamma);
    printf("Radii scale factor and offset:           %.2f %.2f\n",params->radii_scale_factor,params->radii_offset);
    printf("Table of radii (type name, solvation radius (A), volume(A^3), volume radius (A):\n");
    for (itype=0; itype<nfragtypes; itype++) {
        printf("%32s %.3f %.3f %.3f\n",fragtypes[itype]->fragname,fragtypes[itype]->solvation_radius,fragtypes[itype]->volume,fragtypes[itype]->volume_radius);
    }
}

