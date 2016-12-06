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
    double fk,xs,xs2,y;
    fk=(x*gbCubicSplineIntervalInv);
    //justin adds this to cover our rear end if k is outside the range
    if (fk>=GB_NUM_SPLINE-2) return 1/sqrt(x); //spline function goes to 1/x as x -> infinity
    k=(int) fk;
    xs = x - k*gbCubicSplineInterval;
    xs2 = xs*xs;
    y = gbCubicSplineParams[4*k]*xs2*xs + gbCubicSplineParams[4*k+1]*xs2 + gbCubicSplineParams[4*k+2]*xs + gbCubicSplineParams[4*k+3];
    return y;
}

//same thing but includes derivative
void splintGB_diff(const double x, double * y, double * dy)
{
    int k;
    double fk,xs,xs2;
    fk=(x*gbCubicSplineIntervalInv);
    //justin adds this to cover our rear end if k is outside the range
    if (fk>=GB_NUM_SPLINE-2) {
        *y=1/sqrt(x);
        *dy=-((*y)/(2*x)); //-(1/2)x^(-3/2) = x^(-1/2) * (-1/(2x))
        return;
    }
    k = (int) fk;
    xs = x - k*gbCubicSplineInterval;
    xs2 = xs*xs;
    *y = gbCubicSplineParams[4*k]*xs2*xs + gbCubicSplineParams[4*k+1]*xs2 + gbCubicSplineParams[4*k+2]*xs + gbCubicSplineParams[4*k+3];
    *dy = 3*gbCubicSplineParams[4*k]*xs2 + 2*gbCubicSplineParams[4*k+1]*xs + gbCubicSplineParams[4*k+2];
    //return y;
}

//read a CHARMM long format cor file with born radii in WMAIN
//fm2='(2I10,2X,A8,2X,A8,3F20.10,2X,A8,2X,A8,F20.10)'
void simulation::read_born_radii(char * fname, double * per_atom_born_radii)
{
    FILE * input;
    char buffer[255];
    int junk, ires,iatom;
    char resname[4],aname[5];
    const char * fmt1 = "%d %d %3s %4s";
    bool flag; //reached ext
    double alpha;
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open file %s.\n",fname);
        die();
    }
    flag=false;
    while (!feof(input)) {
        fgets(buffer,sizeof(buffer),input);
        if (strstr(buffer,"EXT")!=NULL) break;
    }
    if (feof(input)) {
        printf("Error reading Born radii file %s.\n",fname);
        die();
    }
    while (!feof(input)) {
        fgets(buffer,sizeof(buffer),input);
        sscanf(buffer,fmt1,&junk,&ires,resname,aname);
        iatom=top->find_atom(' ',ires,aname);//won't work for multichain
        if (iatom<0) {
            printf("could not find atom %d %s from Born radii file %s.\n",ires,aname,fname);
            die();
        }
        sscanf(buffer+120,"%lg",&alpha);
        per_atom_born_radii[iatom]=alpha;
    }
}


//Solve the equation  - (1-1/eps_w) Kcoul * sum_{i<j} q_i q_j (r_ij^2 + alpha_k exp(-r_ij^2 /(4 alpha_k^2))^(-1/2) for alpha_k, given q_i and r_ij
//dgx = delta G / tau = delta G / [(1/epsp - 1/epsw]*COUL_CONST]
//to do: use the spline system here

double fragmenttype::get_average_born_radius(forcefield * ffield, double dgx)
{
    int i,j,iter;
    double alpha,alpha2,qq,dgest,ddgdalpha,r2,aux,aux2,aux3,fGB,dfGBdalpha,newalpha,conv,x,spl,dspl,lambda;
    alpha=10;
    iter=1;
    if (!is_charged) return 0;
    while (iter<=GB_MAX_ITER) {
        dgest=0; //est. of deltaG
        ddgdalpha=0; //d(deltaG)/d(alpha);
        alpha2=alpha*alpha;
        for (i=0; i<natom; i++) for (j=i; j<natom; j++) {
            qq=(ffield->chargeParams[types[i]])*(ffield->chargeParams[types[j]]);
            if (qq==0) continue;
            if (i==j) {
                dgest+=qq/(2.0*alpha);
                ddgdalpha+=(-qq/(2.0*alpha2));
            } else {
                r2=d2matrix[i][j];
                //spline function?
                /*aux=exp(-r2/(GB_KAPPA*alpha2));
                fGB=1/sqrt(r2+alpha2*aux);
                aux2=2*alpha+(2*r2)/(GB_KAPPA*alpha);
                dfGBdalpha=-0.5*fGB*fGB*fGB*aux*aux2;*/ //this 0.5 comes from the power rule
                //If using the spline function, need to take into account that f_GB = 1/(alpha_k spline(r^2/alpha_k^2)) and differentiate properly
                x=r2/alpha2;
                splintGB_diff(x,&spl,&dspl);
                fGB=spl/alpha;
                dfGBdalpha=-(spl+2*x*dspl)/alpha2;
                dgest+=qq*fGB;
                ddgdalpha+=qq*dfGBdalpha;
            }
        }
        //print newalpha,dgest,dgx,ddgdalpha;
        //newton's method, scale deviation if necessary to ensure alpha stays >0
        lambda=1.0;
        do {
            newalpha=alpha-lambda*((dgest-dgx)/(ddgdalpha));
            lambda=0.5*lambda;
        } while (((lambda<1) && (newalpha<0.001)) || ((lambda>=1) && (newalpha<0))); //numerical safety in case alpha is too small
        /*if (newalpha>0) conv=dgest-dgx; else {
            newalpha=0.5*alpha;*/
        conv=fabs(dgest-dgx);
        alpha=newalpha;
        //if (conv<0) conv=-conv;
        if (conv<GB_RADIUS_CONV) break;
        //print newalpha,dgest,dgx,ddgdalpha,conv;
        //print alpha,dgest,ddgdalpha;
        //alpha+=0.1;
        iter++;
    }
    if (iter==GB_MAX_ITER) {
        printf("Warning: Newton iteration for determining average born radius failed to converge\n");
    }
    return alpha;
}

//The per_atom_born_radii array must be arranged in order by fragment atom.
double fragmenttype::get_average_born_radius(forcefield * ffield, double * per_atom_born_radii)
{
    int iatom,jatom,k;
    double dgx,r2,fGB,aux,aux2,aux3,aiaj,qq,x;
    dgx=0.0;
    if (!is_charged) return 0;
    for (iatom=0; iatom<natom; iatom++) for (jatom=iatom; jatom<natom; jatom++) {
        qq=(ffield->chargeParams[types[iatom]])*(ffield->chargeParams[types[jatom]]);
	if (qq==0) continue;
        if (iatom==jatom) {
            dgx+=qq/(2.0*per_atom_born_radii[iatom]);
        } else {
            r2=d2matrix[iatom][jatom];
            aiaj=per_atom_born_radii[iatom]*per_atom_born_radii[jatom];
            /*aux=exp(-r2/(GB_KAPPA*aiaj));
            fGB=1/sqrt(r2+aiaj*aux);*/
            x=r2/aiaj;
            fGB=splintGB(x)/sqrt(aiaj);
            dgx+=qq*fGB;
        }
    }
#ifdef DEBUG
    printf(" %.10f\n",dgx);
#endif
    return get_average_born_radius(ffield,dgx);
}

//Compute per-fragment born radii from per-atom born radii.
void topology::convert_born_radii(forcefield * ffield, const double * per_atom_born_radii, double * per_fragment_born_radii)
{
    int ifrag,iatom,iactualatom,nknownatoms,k;
    fragmenttype * fragtype;
    double tempradii[MAX_ATOMS_PER_FRAGMENT];
    for (ifrag=0; ifrag<nfrag; ifrag++) {
        fragtype=fragtypes[frags[ifrag].type];
        //Rearrange born radii into fragment order.
        //Also figure out how many atoms have known coordinates.
        for (iatom=0; iatom<fragtype->natom; iatom++) {
            iactualatom=frags[ifrag].atoms[iatom];
            tempradii[iatom]=per_atom_born_radii[iactualatom];
            //mass[iactualatom]=atoms[iactualatom].mass;
        }
#ifdef DEBUG
        printf("convert_born_radii: %d",ifrag); //matches with another print statement in fragmenttype::get_average_born_radius
#endif
        per_fragment_born_radii[ifrag]=fragtype->get_average_born_radius(ffield,tempradii);
    }
}


//read facts parameter file
//bool per_atom = true if they are per atom parameters, false otherwise
//if doing per fragment parameters, need to have read the definitions file first (for nfragtypes)
void simulation::read_facts_params(bool per_atom, char * fname, forcefield * ffield, double tau)
{
    FILE * input;
    char fragname[MAX_FRAGMENT_NAME],buffer[255];
    char * p;
    int i,itype,iclass,ntypes;
    double radius,dG_solv_free,sasa_free;
    facts_param_info info;
    if (per_atom) ntypes=MAX_NUM_OF_ATOM_CLASSES; else ntypes=top->nfragtypes;
    facts_params=ALLOCATE(ntypes,facts_param_info);
    for (i=0; i<ntypes; i++) facts_params[i].volume=0;
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open FACTS parameter file %s.\n",fname);
        die();
    }
    while (!feof(input)) {
        fgets(buffer,sizeof(buffer),input);
        //if (buffer[0]=='#') continue;
        //remove any comment
        if (feof(input)) break;
        p=strchr(buffer,'#');
        if (p!=NULL) *p='\0';
        if (strlen(buffer)==0) continue;
        if (strncasecmp(buffer,"END",3)==0) break;

        if (per_atom) {
            //type number, b1, b2, a2, a3, Rsphere, d1, d2, c2, c3
            sscanf(buffer,"%d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",&itype,&info.b1,&info.b2,&info.a0,&info.a2,&info.a3,&info.rsphere2,
                &info.d1,&info.d2,&info.c0,&info.c2,&info.c3);
            radius=ffield->vdwParams[itype].sigma;
            if (itype<=3) radius=1.0; //FACTS parameters are derived with a vdw radius for hydrogens of 1.0
            info.volume=(4*M_PI/3)*radius*radius*radius;
            //info.a0=-(tau/(2*radius))*(1+exp(-info.a2*info.a3)); //see embedded equation in FACTS paper, between eqs. 7 and 8, check this formula!
            dG_solv_free=-tau/(2*radius); //this is missing a factor of q_i^2 which is provided in finish_born_radii
            //info.a0=dG_solv_free*(1-1/(1+exp(info.a2*info.a3)));
            sasa_free=4*M_PI*(radius+WATER_RADIUS)*(radius+WATER_RADIUS);
        } else {
            //fragment type name, volume, dG_solv(free fragment), b1, b2, a2, a3, Rsphere, d1, d2, c2, c3
            sscanf(buffer,"%32s %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",fragname,&info.volume,&sasa_free,&dG_solv_free,&info.b1,&info.b2,&info.a2,&info.a3,&info.rsphere2,
                &info.d1,&info.d2,&info.c2,&info.c3);
            itype=top->frag_type_by_name(fragname);
            info.a0=dG_solv_free*(1-1/(1+exp(info.a2*info.a3)));
            info.c0=sasa_free*(1+exp(-info.c2*info.c3));//embedded equation between eqs 11 and 12, check this formula!
        }
        info.rsphere2=info.rsphere2*info.rsphere2;
        info.a1=-info.a0;
        info.c1=-info.c0;
        facts_params[itype]=info;
    }
    //this assumes a system is in place, check to make sure all parameters are defined
    //maybe this should be a separate membor function of topology?
    if (per_atom) {
        for (i=0; i<top->natom; i++) if (facts_params[top->atoms[i].classx].volume<=0) {
            printf("Error: missing FACTS parameters for atom class %d\n",i);
            die();
        }
    } else {
        for (i=0; i<top->nfrag; i++) if (facts_params[top->frags[i].type].volume<=0) {
            printf("Error: missing FACTS parameters for fragment type %s\n",top->fragtypes[top->frags[i].type]->fragname);
            die();
        }
    }
}

//This calculates the A_i and B_i using eqs. 3-5 from FACTS paper.  It can be used either on a per-atom or per-fragment basis.
void topology::partial_get_born_radii(bool per_atom, facts_param_info * facts_params, int n, const double * coords,  double * __restrict__ numer, double * __restrict__ denom, double * __restrict__ a, double * __restrict__ b)
{
    //double rsphere2;
    int i,j,k,itype,jtype;
    double theta,ivolume,jvolume,factor,rij[3],r2,r,bb,aux,c;
    //numer=(double *) checkalloc(3*n,sizeof(double)); //should we do this allocation early?
    //denom=(double *) checkalloc(n,sizeof(double));
    for (i=0; i<n; i++) {
        denom[i]=1;
        a[i]=0;
        for (k=0; k<3; k++) numer[3*i+k]=0;
    }
    //rsphere2=rsphere*rsphere;
    for (i=0; i<n; i++) {
        //aa=0;
        //for (k=0; k<3; k++) numer[k]=0;
        //denom=1; //in denominator, eq. 4
        //it may be possible to skip half the trips through this loop by symmetry
        //ivolume=fragtypes[itypes[i]].volume;
        if (per_atom) itype=atoms[i].classx; else itype=frags[i].type;
        for (j=i+1; j<n; j++) {
            if (per_atom) jtype=atoms[j].classx; else jtype=frags[j].type;
            for (k=0; k<3; k++) rij[k]=coords[3*j+k]-coords[3*i+k]; //need to make pbc sensitive
            r2=0;
            for (k=0; k<3; k++) r2+=rij[k]*rij[k];
            //if (r2>rsphere2) continue; //theta_ij = 0, no contributions
            if ((r2<=facts_params[itype].rsphere2) || (r2<=facts_params[jtype].rsphere2)) {
                r=sqrt(r2);
                for (k=0; k<3; k++) rij[k]/=r; //changes x_ij to x^_ij
            }
            if (r2<=facts_params[itype].rsphere2) {
                jvolume=facts_params[jtype].volume;
                theta=1-r2/facts_params[itype].rsphere2;
                theta=theta*theta; //eq. 5 FACTS paper
                a[i]+=jvolume*theta; //eq. 3, both sides
                denom[i]+=(jvolume*theta/r);
                for (k=0; k<3; k++) numer[3*i+k]+=(jvolume*theta*rij[k]/r);
            }
            if (r2<=facts_params[jtype].rsphere2) {
                ivolume=facts_params[itype].volume;
                theta=1-r2/facts_params[jtype].rsphere2;
                theta=theta*theta; //eq. 5 FACTS paper
                a[j]+=ivolume*theta;
                denom[j]+=(ivolume*theta/r);
                for (k=0; k<3; k++) numer[3*j+k]-=(ivolume*theta*rij[k]/r);
            }
        }
    }
    for (i=0; i<n; i++) {
        bb=0;
        for (k=0; k<3; k++) {
            aux=numer[3*i+k]/denom[i];
            bb+=aux*aux;
        }
        bb=sqrt(bb);
        b[i]=bb;
#ifdef DEBUG
        if ((b[i]<=0) || (b[i]>=1)) {
             printf("Error in FACTS, particle %d %.10f %.10f %.10f %.10f %.10f\n",i,b[i],numer[3*i],numer[3*i+1],numer[3*i+2],denom[i]);
             die();
        }
#endif
        //peptide's fit params,
        /*c=a[i]+825.7*b[i]-1.20688*a[i]*b[i];
        aux=1/(1+exp(-0.00764646*c));
        dgx[i]=-9.564*(1-aux);*/
    }
}

//this finishes the born radii given the A_i's and B_i's, also does the SASA term, using eqs. 2, 6, 7, 10, and 11
//need to add GB self term to this!
void topology::finish_born_radii(bool per_atom, facts_param_info * facts_params, int n,  forcefield * ffield, const double tau, const double * a, const double * b, double * __restrict__ born_radii, double * __restrict__ dg_self, double * __restrict__ sasa)
{
    double c,d,aux,aux2,q;
    int i,itype;
    //sasatot=0.0;
    //gb_self=0.0;
    for (i=0; i<n; i++) {
        if (per_atom) itype=atoms[i].classx; else itype=frags[i].type;
        c=a[i]+facts_params[itype].b1*b[i]+facts_params[itype].b2*a[i]*b[i]; //eq. 6 FACTS paper
        aux=facts_params[itype].a2*(c-facts_params[itype].a3);
        aux2=1-1/(1+exp(-aux));
        dg_self[i]=facts_params[itype].a0*aux2; //eq. 7, in all cases a1=-a0, still missing a factor of q^2
        if ((dg_self[i]>-0.001) && (per_atom || fragtypes[itype]->is_charged)) {
            if (per_atom) printf("error: solvation energy %g too high for atom %d res %d name %s\n",dg_self[i],i,atoms[i].resNum+1,atoms[i].name);
            else printf("error: solvation energy %g too high for fragment %d type %s\n",dg_self[i],i,fragtypes[frags[i].type]->fragname);
        }
        if (per_atom) {
            q=ffield->chargeParams[atoms[i].type];
            born_radii[i]=-tau/(2*dg_self[i]); //eq. 2 FACTS paper, the factor of q^2 cancels out
            dg_self[i]=dg_self[i]*q*q; //now we supply it
#ifdef DEBUG
            printf("gb_self: %d %s %d %s %d %.4f %.4f\n",i,atoms[i].resName,atoms[i].resNum+1,atoms[i].name,atoms[i].fragment,atoms[i].radius/born_radii[i],dg_self[i]);
#endif
        } else {
            //we're finding a per fragment born radius after a per-fragment determination, need to solve with newton's method
            born_radii[i]=fragtypes[itype]->get_average_born_radius(ffield,-dg_self[i]/tau); //do we need this minus sign?
        }
        //now do the sasa term, eqs. 10 and 11
        d=a[i]+facts_params[itype].d1*b[i]+facts_params[itype].d2*a[i]*b[i]; //eq. 10 FACTS paper
        aux=facts_params[itype].c2*(d-facts_params[itype].c3);
        if (aux>100) sasa[i]=0; else sasa[i]=facts_params[itype].c0*(1-1/(1+exp(-aux))); //eq. 11 (using c1=-c0)
        if (sasa[i]<0.000001) {
            if (per_atom) printf("error: SASA %g too low for atom %d res %d name %s\n",sasa[i],i,atoms[i].resNum+1,atoms[i].name);
            else printf("error: SASA %g too low for fragment %d type %s\n",sasa[i],i,fragtypes[frags[i].type]->fragname);
        }
        //sasatot+=sasa;
        //gb_self+=dg;
    }
}

//control routine. checks to see all arrays are allocated, if not allocates them, then calls above two routines as appropriate
//also deals with timers (switching back to TIMER_OTHER)
//if we are working with old/new born radii, this may need additional parameters
//we could have some interaction terms from forcefield::non_tabulated_energy in energies[EN_GB_EXACT]
void simulation::calculate_born_radii(const double * center, const double * coords, double * __restrict__ per_atom_born_radii, double * __restrict__ per_fragment_born_radii, double * __restrict__ dg_self, double * __restrict__ sasa)
{
    //double sasa,egb_self;
    int i;
    switch(gb_mode) {
        case GB_MODE_NONE: return;
        case GB_MODE_PER_ATOM:
        case GB_MODE_PER_ATOM_CONVERTED:
#ifdef TIMERS
            switch_timer(TIMER_PER_ATOM_BORN_RADII);
#endif
            if (facts_a==NULL) facts_a=ALLOCATE(top->natom,double);
            if (facts_b==NULL) facts_b=ALLOCATE(top->natom,double);
            if (facts_numer==NULL) facts_numer=ALLOCATE(3*top->natom,double);
            if (facts_denom==NULL) facts_denom=ALLOCATE(top->natom,double);
            //if (per_atom_born_radii==NULL) per_atom_born_radii=ALLOCATE(top->natom,double);
            //if (dg_self==NULL) dg_self=ALLOCATE(top->natom,double);
            //if (sasa==NULL) sasa=ALLOCATE(top->natom,double);
            top->partial_get_born_radii(true,facts_params,top->natom,coords,facts_numer,facts_denom,facts_a,facts_b);
            top->finish_born_radii(true,facts_params,top->natom,ffield,gb_params.tau,facts_a,facts_b,per_atom_born_radii,dg_self,sasa);
            if (gb_mode==GB_MODE_PER_ATOM_CONVERTED) {
#ifdef TIMERS
                switch_timer(TIMER_CONVERT_BORN_RADII);
#endif
                if (per_fragment_born_radii==NULL) per_fragment_born_radii=ALLOCATE(top->nfrag,double);
                top->convert_born_radii(ffield,per_atom_born_radii,per_fragment_born_radii);
            }
            break;
        case GB_MODE_PER_FRAGMENT:
#ifdef TIMERS
            switch_timer(TIMER_PER_FRAGMENT_BORN_RADII);
#endif
            if (facts_a==NULL) facts_a=ALLOCATE(top->nfrag,double);
            if (facts_b==NULL) facts_b=ALLOCATE(top->nfrag,double);
            if (facts_numer==NULL) facts_numer=ALLOCATE(3*top->nfrag,double);
            if (facts_denom==NULL) facts_denom=ALLOCATE(top->nfrag,double);
            //if (per_fragment_born_radii==NULL) per_fragment_born_radii=ALLOCATE(top->nfrag,double);
            //if (dg_self==NULL) dg_self=ALLOCATE(top->nfrag,double);
            //if (sasa==NULL) sasa=ALLOCATE(top->nfrag,double);
            top->partial_get_born_radii(false,facts_params,top->nfrag,center,facts_numer,facts_denom,facts_a,facts_b);
            top->finish_born_radii(false,facts_params,top->nfrag,ffield,gb_params.tau,facts_a,facts_b,per_fragment_born_radii,dg_self,sasa);
            break;
    }
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
}





//GB solvation energy using atomic born radii
double topology::born_solvation_energy1(forcefield * ffield, int ifrag, int jfrag, double tau, double * coords, double * per_atom_born_radii)
{
    int iatom, jatom, ifragatom, jfragatom,k;
    double egb, egbtot,qq,r2,x,aiaj;
    fragmenttype * ifragtype = fragtypes[frags[ifrag].type];
    fragmenttype * jfragtype = fragtypes[frags[jfrag].type];
    egbtot=0.0;
    for (ifragatom=0; ifragatom<ifragtype->natom; ifragatom++)
        for (jfragatom=0; jfragatom<jfragtype->natom; jfragatom++){
            iatom=frags[ifrag].atoms[ifragatom];
            jatom=frags[jfrag].atoms[jfragatom];
            qq=(ffield->chargeParams[ifragtype->types[ifragatom]])*(ffield->chargeParams[jfragtype->types[jfragatom]]);
            if (iatom>jatom) continue;
            else if (iatom==jatom) {
                egb=-qq/(2*per_atom_born_radii[iatom]);
            } else {
                r2=pbc_distance2(false,0,0,&coords[3*iatom],&coords[3*jatom]);
                aiaj=per_atom_born_radii[iatom]*per_atom_born_radii[jatom];
                x=splintGB(r2/aiaj);//=1/sqrt(r^2/aiaj + exp(-GB_KAPPA*r^2/aiaj))
                egb=-qq*x/sqrt(aiaj);
                //egb=-qq/sqrt(r2+aiaj*exp(-r2/(GB_KAPPA*aiaj)));
            }
            egbtot+=egb;
    }
    egbtot*=tau;
    return egbtot;
}

//use per-fragment born radii
double topology::born_solvation_energy2(forcefield * ffield, int ifrag, int jfrag, double tau, double rkl2, double * coords, double xx)
{
    int iatom, jatom, ifragatom, jfragatom,k;
    double egb, egbtot,qq,r2,x,akal;
    fragmenttype * ifragtype = fragtypes[frags[ifrag].type];
    fragmenttype * jfragtype = fragtypes[frags[jfrag].type];
    egbtot=0.0;
    //printf("%.2f %.2f\n",rkl2,xx);
    akal=(xx*xx)*rkl2;
    //printf("%.2f %.2f %.2f\n",rkl2,xx,akal);
    for (ifragatom=0; ifragatom<ifragtype->natom; ifragatom++)
        for (jfragatom=0; jfragatom<jfragtype->natom; jfragatom++){
            iatom=frags[ifrag].atoms[ifragatom];
            jatom=frags[jfrag].atoms[jfragatom];
            qq=(ffield->chargeParams[ifragtype->types[ifragatom]])*(ffield->chargeParams[jfragtype->types[jfragatom]]);
            if (iatom>jatom) continue;
            else if (iatom==jatom) {
                egb=-qq/(2*sqrt(akal));
            } else {
                r2=pbc_distance2(false,0,0,&coords[3*iatom],&coords[3*jatom]);
                x=splintGB(r2/akal);//=1/sqrt(r^2/aiaj + exp(-kappa*r^2/aiaj))
                egb=-qq*x/sqrt(akal);
            }
            egbtot+=egb;
    }
    egbtot*=tau;
    return egbtot;
}
//discretize sqrt(ak al)
#define ALPHAPOINTS     20
#define ALPHAMAX        10
double discretize(double x, double maxx, int nx)
{
    int ix;
    double dx = maxx/nx;
    //printf("discretize1: x, maxx, nx = %.2f %.2f %d\n",x,maxx,nx);
    ix=((int) (x/dx));
    if (ix>=nx) ix=nx-1;
    //printf("discretize2: x = %.2f\n",((ix+0.5)*dx));
    return ((ix+0.5)*dx);
    //printf("discretize2: x = %.2f\n",x);
}


/*double topology::born_solvation_energy4(forcefield * ffield, int ifrag, int jfrag, double tau, double rkl2, double * coords, double akal)
{
    int iatom, jatom, ifragatom, jfragatom,k;
    double egb, evacuum,qq,r2,x,factor;
    fragmenttype * ifragtype = fragtypes[frags[ifrag].type];
    fragmenttype * jfragtype = fragtypes[frags[jfrag].type];
    //evacuum=0.0;
    evacuum=exact_interaction_energy(ffield,false,0,0,1.0,false,ifrag,jfrag,coords); //this has been modified to be electrostatics only
    //rkl2=pbc_distance2(false,0,0,&center[3*ifrag],&center[3*jfrag]);
    x=splintGB(rkl2/akal);
    factor=sqrt(rkl2/akal)*x;
    egb=-factor*evacuum*tau;
    return egb;
}*/

/*void simulation::born_test(char * born_radii_file, int nalpha, char * output_file)
{
    int ifrag,jfrag;
    double egb1,egb2,egb3, egb4,evacuum, epsp, epsw, tau,akal,rkl2,x,xx,xdisc,factor,epseff;
    double egbtot1, egbtot2, egbtot3, egbtot4;
    FILE * output;
    unsigned long long start,end;
    bool close;
    output=fopen(output_file,"w");
    if (output==NULL) output=stdout;
    setup_gb_spline(4.0);//maybe do this from somewhere else
    per_atom_born_radii=(double *) checkalloc(top->natom,sizeof(double));
    per_fragment_born_radii=(double *) checkalloc(top->nfrag,sizeof(double));
    read_born_radii(born_radii_file,per_atom_born_radii);
#ifdef TIMERS
    start=get_rdtsc();
#endif
    top->convert_born_radii(ffield,per_atom_born_radii,per_fragment_born_radii);
#ifdef TIMERS
    end=get_rdtsc();
    printf("Clock cycles for born radii conversion: %lld\n",end-start);
#endif
    epsp=1.0;
    epsw=80.00;
    tau=(1/epsp-1/epsw)*COUL_CONST;
    egbtot1=0;
    egbtot2=0;
    egbtot3=0;
    egbtot4=0;
    for (ifrag=0; ifrag<top->nfrag; ifrag++)
        for (jfrag=0; jfrag<top->nfrag; jfrag++) {
            //if (!((top->fragtypes[frags[ifrag].type]->is_charged) && (fragtypes[frags[ifrag].type]->is_charged))) continue;
            egb1=top->born_solvation_energy1(ffield,ifrag,jfrag,tau,oldcoords,per_atom_born_radii);
            if (egb1==0) continue;
            if (ifrag==jfrag) {
                //For ifrag==jfrag, egb1==egb2 by construction.  There would be no use of tables,
                //so we would do the self term separately using the most exact value available.
            egbtot1+=egb1;
                egbtot2+=egb1; //egb2==egb1 by construction
                egbtot3+=egb1;
                egbtot4+=egb1;
                continue; //can't do the interaction solvation energies
            }
            rkl2=pbc_distance2(false,0,0,&oldcenter[3*ifrag],&oldcenter[3*jfrag]);
            akal=per_fragment_born_radii[ifrag]*per_fragment_born_radii[jfrag];
            xx=sqrt(akal/rkl2);
            //printf("in born_test: xx= %.2f\n",xx);
            egb2=top->born_solvation_energy2(ffield,ifrag,jfrag,tau,rkl2,oldcoords,xx);
            xdisc=discretize(xx,1,nalpha);
            //printf("in born_test: xdisc = %.2f\n",xdisc);
            egb3=top->born_solvation_energy2(ffield,ifrag,jfrag,tau,rkl2,oldcoords,xdisc);
            //egb4=top->born_solvation_energy4(ffield,ifrag,jfrag,tau,rkl2,oldcoords,akal);
            evacuum=top->exact_interaction_energy(ffield,false,0,0,1.0,ifrag,jfrag,oldcoords); //this has been modified to be electrostatics only and to be missing the COUL_CONST
            x=splintGB(rkl2/akal);
            factor=sqrt(rkl2/akal)*x;
            egb4=-factor*evacuum*tau;
            epseff=1/epsp+egb1/(evacuum*COUL_CONST); //"exact" effective dielectric constant
            epseff=1/epseff;
            //if (ifrag!=jfrag)
            close=top->closefragments[ifrag*top->nfrag+jfrag];
            if (ifrag!=jfrag) fprintf(output,"%d %s %d %s %c %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n",
                ifrag,top->fragtypes[top->frags[ifrag].type]->fragname,jfrag,top->fragtypes[top->frags[jfrag].type]->fragname,
                yesno(close),sqrt(rkl2),1/xx,epseff,egb1,egb2,egb3,egb4);
            egbtot1+=egb1;
            egbtot2+=egb2;
            if (close) { //we could calculate close interactions using per-fragmetn born radii
                egbtot3+=egb2;
                egbtot4+=egb2;
            } else {
                egbtot3+=egb3;
                egbtot4+=egb4;
            }
    }
    printf("GB totals: %.4f %.4f %.4f %.4f\n",egbtot1,egbtot2,egbtot3,egbtot4);
}*/
