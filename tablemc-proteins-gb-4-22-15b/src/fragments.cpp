#include <cstdio>
#include <cctype>
#include "fragments.h"
#include "ffield.h"
#include "rotations.h"
#include "util.h"
#include "solvation.h"

fragmenttype::fragmenttype(const char * name, const char * fname, forcefield * ffield)
{
    FILE * f;
    double xx,yy,zz,q,r2;
    double inertia[3][3],axes[3][3],cross[3],det,temp[3];
    char junk[255];
    char aname[4];
    int i,j,k,junk2;
    //printf("*** %d\n",strlen(fname));
    n_used=0;
    strncpy(fragname,name,sizeof(fragname));
    has_covalent_tables=((strstr(fragname,"PEPTIDE")!=NULL) || (strcmp(fragname,"N-METHYL-AMIDE")==0));
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("Could not open fragment file %s\n",fname);
        die();
    }
    fscanf(f,"%d\n",&natom);
    ///printf("*** %d\n",*natom);
    fgets(junk,sizeof(junk),f);
    xx=0.0;
    yy=0.0;
    zz=0.0;
    totalmass=0.0;
    for (i=0; i<natom; i++) {
        fscanf(f,"%s %lg %lg %lg %d\n",&names[i],&refgeom[3*i],&refgeom[3*i+1],&refgeom[3*i+2],&types[i]);
        mass[i]=ffield->atomTypeLookUp[types[i]].mass;
        xx+=mass[i]*refgeom[3*i];
        yy+=mass[i]*refgeom[3*i+1];
        zz+=mass[i]*refgeom[3*i+2];
        totalmass+=mass[i];
    }
    /*place barycenter at origin*/
    xx/= totalmass;
    yy/= totalmass;
    zz/= totalmass;
    for (i=0; i<natom; i++) {
        refgeom[3*i]-=xx;
        refgeom[3*i+1]-=yy;
        refgeom[3*i+2]-=zz;
    }
    //compute moment of inertia and align with axes
    for (i=0; i<3; i++) for (j=0; j<3; j++) inertia[i][j]=0.0;
    for (i=0; i<natom; i++) {
        xx=refgeom[3*i];
        yy=refgeom[3*i+1];
        zz=refgeom[3*i+2];
        inertia[0][0]+=mass[i]*(yy*yy+zz*zz);
        inertia[0][1]+=-mass[i]*xx*yy;
        inertia[0][2]+=-mass[i]*xx*zz;
        inertia[1][1]+=mass[i]*(xx*xx+zz*zz);
        inertia[1][2]+=-mass[i]*yy*zz;
        inertia[2][2]+=mass[i]*(xx*xx+yy*yy);
    }
    inertia[1][0]=inertia[0][1];
    inertia[2][0]=inertia[0][2];
    inertia[2][1]=inertia[1][2];
    jacobi(3,&inertia[0][0],&axes[0][0]);
    //the axes are {axes[0][0], axes[1][0], axes[2][0]}, etc.
    //make sure it's a rotation matrix (det axes = +1)
    cross[0]=axes[1][1]*axes[2][2]-axes[2][1]*axes[1][2];
    cross[1]=axes[2][1]*axes[0][2]-axes[0][1]*axes[2][2];
    cross[2]=axes[0][1]*axes[1][2]-axes[1][1]*axes[0][2];
    det=axes[0][0]*cross[0]+axes[1][0]*cross[1]+axes[2][0]*cross[2];
    if (det<0) for (i=0; i<3; i++) axes[i][2]*=-1.0;
    //this actually multiplies by the transpose of axes,
    for (i=0; i<natom; i++) {
        matmul(&axes[0][0],&refgeom[3*i],&temp[0]);
        for (k=0; k<3; k++) refgeom[3*i+k]=temp[k];
    }
    //calculate dipole moment for possible use -- may take this out (but leave hte intitialization of is_charged!)
    is_charged=false;
    qtot=0.0;
    for (k=0; k<3; k++) dipole[k]=0.0;
    for (i=0; i<natom; i++) {
        q=ffield->chargeParams[types[i]];
        if (fabs(q)>1e-6) is_charged=true;
        qtot+=q;
        for (k=0; k<3; k++) dipole[k]+=q*(refgeom[3*i+k]);
    }
    dipolemag=0.0;
    for (k=0; k<3; k++) dipolemag+=dipole[k]*dipole[k];
    dipolemag=sqrt(dipolemag);
    //calculate distance matrix
    for (i=0; i<natom; i++) {
        d2matrix[i][i]=0;
        for (j=i+1; j<natom; j++) {
            r2=0.0;
            for (k=0; k<3; k++) r2+=(refgeom[3*i+k]-refgeom[3*j+k])*(refgeom[3*i+k]-refgeom[3*j+k]);
            d2matrix[i][j]=r2;
            d2matrix[j][i]=r2;
        }
    }
#ifdef DEBUG_NON_TABULATED
    printf("Total charge for fragment %s: %g\n",fname,qtot);
    printf("Dipole moment for fragment %s: %.2f %.2f %.2f\n",fname,dipole[0],dipole[1],dipole[2]);
    printf("Magnitude of dipole:           %.2f\n",dipolemag);
#else
    printf("Fragment file %s loaded.\n",fname);
#endif
    fclose(f);
}

void fragmenttype::get_coords(double * center, double * orient, double * coords)
{
    double x[3],xx[3],rotmatrix[3][3];
    int i;
    quat_to_matrix(&orient[0],&rotmatrix[0][0]);
    for (i=0; i<natom; i++) {
        /*copy reference coordinates from fragment atom i */
        x[0]=refgeom[3*i];
        x[1]=refgeom[3*i+1];
        x[2]=refgeom[3*i+2];
        matmul(&rotmatrix[0][0],x,xx);
        coords[3*i]=xx[0]+center[0];
        coords[3*i+1]=xx[1]+center[1];
        coords[3*i+2]=xx[2]+center[2];
    }
}

//returns the RMSD
void fragmenttype::fit_fragment(double * coords, double * center, double * orient, double * rmsd)
{
    fit_fragment(coords,center,orient,mass,rmsd);
}

void fragmenttype::fit_fragment(double * coords, double * center, double * orient, double * weights, double * rmsd)
{
    //I think it will be ok if only 1 or 2 atoms.  Will check this later.
    //Do the fit using only heavy atoms.
    int i;
    rmsd_fit(natom,weights,refgeom,coords,center,orient,rmsd);
}
int fragmenttype::findatombyname(const char * name)
{
    int iatom;
    for (iatom=0; iatom<natom; iatom++)
        if (strncasecmp(name,names[iatom],sizeof(names[iatom]))==0) return iatom;
    return -1; //not found
}

//These are used for model building
double fragmenttype::get_bond_length(int iatom, int jatom)
{
    double r;
    int k;
    r=0.0;
    for (k=0; k<3; k++) r+=(refgeom[3*iatom+k]-refgeom[3*jatom+k])*(refgeom[3*iatom+k]-refgeom[3*jatom+k]);
    r=sqrt(r);
    return r;
}

double fragmenttype::get_angle(int iatom, int jatom, int katom)
{
    double rji[3],rjk[3];
    int l;
    for (l=0; l<3; l++) {
        rjk[l]=refgeom[3*katom+l]-refgeom[3*jatom+l];
        rji[l]=refgeom[3*iatom+l]-refgeom[3*jatom+l];
    }
    return angle(rji,rjk);
}

double fragmenttype::get_dihedral(int iatom, int jatom, int katom, int latom)
{
    double dihed, rji[3],rjk[3],rkl[3],cross[3],dot;
    int m;
    for (m=0; m<3; m++) {
        rji[m]=refgeom[3*iatom+m]-refgeom[3*jatom+m];
        rjk[m]=refgeom[3*katom+m]-refgeom[3*jatom+m];
        rkl[m]=refgeom[3*latom+m]-refgeom[3*katom+m];
    }
    dihed=dihedral(rji,rjk,rkl);
    //We need to determine the sign of the dihedral.
    /*cross[0]=rji[1]*rjk[2]-rji[2]*rjk[1];
    cross[1]=rji[2]*rjk[0]-rji[0]*rjk[2];
    cross[2]=rji[0]*rjk[1]-rji[1]*rjk[0];
    dot=cross[0]*rkl[0]+cross[1]*rkl[1]+cross[2]*rkl[2];
    if (dot>0) dihed=-dihed;*/
    return dihed;
}




