#include <cstdio>
#include <cstdlib>
#include "fragments.h"
#include "ffield.h"
#include "mc.h"
#include "rotations.h"
#include "mt.h"
#include "util.h"
//DCD related stuff.
#define AKMA               4.88882129E-02
#define CHARMM_VERSION     34
//Writing the random number generator state vector
const char * charmm_signature = "CORD";
//Taken from the ICNTRL array in subroutine WRITCV of dynio.src in CHARMM.
struct dcd_header {
    char signature[4]; //="CORD"
    unsigned int nframes;
    unsigned int begin; //Number of steps first frame.
    unsigned int skip; //Number of steps interval between frames.
    unsigned int nstep; //Total number of steps.
    unsigned int icntrl6to8[3];
    unsigned int ndegf; //Number of degrees of freedom.
    unsigned int nfixed; //Number of fixed atoms.
    float timestep_akma; //Timestep in AKMA units.
    unsigned int crystal; //Whether or not crystal data is written (side lengths for PBC
    unsigned int dim4; //whether or not 4-dimensinal dynamics
    unsigned int qcg; //whether or not charges written (CHEQ)
    unsigned int icntrl14to19[6];
    unsigned int version; //CHARMM version.
};

struct dcd_titles {
    unsigned int ntitle;
    char title[80];
};

//The offsets are one less than indicated on the PDB standard, because
//C arrays are zero-based.
const char * origin = "Produced by Tabulated Monte Carlo code";

void topology::read_pdb_file(char * fname, double * coords)
{
    FILE * f;
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("Could not open pdb file %s for reading\n",fname);
        die();
    }
    read_pdb_stream(f,coords);
    fclose(f);
}
void topology::read_pdb_stream(FILE * input, double * coords)
{
    FILE * f;
    char buf[255],buf2[8],aname[6],chain;
    int ires,iatom;
    double x,y,z;
    while (!feof(input)) {
        fgets(buf,sizeof(buf),input);
        if (strncasecmp("END",buf,3)==0) break;
        if (strncasecmp("ATOM  ",buf,6)!=0) continue;
        //if (strncasecmp("ENDMDL",buf,6)==0)
        sscanf(buf+12,"%4s",aname);
            aname[4]='\0';
        sscanf(buf+22,"%3d",&ires);
        chain=buf[21];
        sscanf(buf+30,"%lf%lf%lf",&x,&y,&z);
        iatom=find_atom(chain,ires,aname);
        if (iatom<0) {
            printf("Could not find atom %c %d %s from PDB file\n",chain,ires,aname);
            die();
        }
        coords[3*iatom]=x;
        coords[3*iatom+1]=y;
        coords[3*iatom+2]=z;
    }
}

void topology::write_pdb_stream(FILE * output, char * fmt, double * coords, forcefield * ffield)
{
    const char * pdbfmt = "ATOM  %5d %4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
    const char * pqrfmt = "ATOM  %5d %4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f %.4f %.4f # %d\n";
    int iatom,ires,iseg,mainfragcount,scfragcount,iclass;
    char chain;
    double color,charge,radius;
    fprintf(output,"REMARK %s\n",origin);
    mainfragcount=0;
    scfragcount=0;
    for (iatom=0; iatom<natom; iatom++) {
        //Determine what segment we are in.
        ires=atoms[iatom].resNum;
        chain=' ';
        for (iseg=0; iseg<nseg; iseg++) if ((ires>segstart[iseg]) && (ires<segend[iseg])) {
            chain=chaincodes[iseg];
            ires=ires-segstart[iseg];
        }
        color=0.0;
        //color main-chain and side-chain fragments separately
        if ((iatom>0) && (atoms[iatom].fragment!=atoms[iatom-1].fragment)) {
            if (atoms[iatom].is_side_chain) {
                scfragcount++;
            } else {
                mainfragcount++;
            }
        }
        if (atoms[iatom].is_side_chain) {
            color=(double) (scfragcount%2) +4.0;
        } else {
            color=(double) (mainfragcount%2);
        }
        if (strcasecmp(fmt,"pdb")==0) {
        //indices are zero-based, add 1 to each
            fprintf(output,pdbfmt,iatom+1,atoms[iatom].name,atoms[iatom].resName,chain,ires+1,
                coords[3*iatom],coords[3*iatom+1],coords[3*iatom+2],1.0,color);
        } else if (strcasecmp(fmt,"pqr")==0) {
            charge=ffield->chargeParams[atoms[iatom].type];
            iclass=atoms[iatom].classx;
            radius=ffield->vdwParams[iclass].sigma;
            fprintf(output,pqrfmt,iatom+1,atoms[iatom].name,atoms[iatom].resName,chain,ires+1,
                coords[3*iatom],coords[3*iatom+1],coords[3*iatom+2],charge,radius,atoms[iatom].fragment);
        } else {
            printf("Unrecognized output format %s,\n",fmt);
            die();
        }
    }
    fprintf(output,"END\n");

}
void topology::write_pdb_file(char * fname, char * fmt, double * coords, forcefield * ffield)
{
    FILE * output;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("Could not open pdb file %s for writing\n",fname);
        die();
    }
    write_pdb_stream(output,fmt,coords,ffield);
    fclose(output);
}


void topology::write_pdb_file2(forcefield * ffield, char * fname, double * coords)
{
    FILE * f;
    const char * pdbfmt = "ATOM  %5d %4.4s %3.3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
    int iatom,ires,iseg;
    char chain;
    f=fopen(fname,"w");
    if (f==NULL) {
        printf("Could not open pdb file %s for writing\n",fname);
        die();
    }
    fprintf(f,"REMARK %s\n",origin);
    for (iatom=0; iatom<natom; iatom++) {
        //Determine what segment we are in.
        ires=atoms[iatom].resNum;
        chain=' ';
        for (iseg=0; iseg<nseg; iseg++) if ((ires>segstart[iseg]) && (ires<segend[iseg])) {
            chain=chaincodes[iseg];
            ires=ires-segstart[iseg];
        }
        //indices are zero-based, add 1 to each
        fprintf(f,pdbfmt,iatom+1,atoms[iatom].name,atoms[iatom].resName,chain,ires+1,
            coords[3*iatom],coords[3*iatom+1],coords[3*iatom+2],1.0,ffield->chargeParams[atoms[iatom].type]);
    }
    fprintf(f,"END\n");
    fclose(f);
}


/* FORTRAN formats copied from subroutine PSFWR2 in source/io/psfres.src :
        heading='PSF EXT'
        fmt00='(/I10,A)'
        fmt01='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8)'
        fmt01a='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)'
        fmt02='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8)'
        fmt02a='(I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6)'
        fmt03='(8I10)'
        fmt04='(9I10)'
        fmt05='(/2I10,A)'
        fmt06='(2I10,3X,L1,3G14.6)'
*/
//This subroutine is incomplete.
void topology::write_psf_file(char * fname, forcefield * ffield)
{
    FILE * f;
    int iatom,ires,iseg;
    char chain;
    double q;
    const char * atomfmt = "%10d        %c %-8d %-8s %-8s %4d %14.6f%14.6f%8d%14.6f%14.6f\n"; //fmt02a above

    f=fopen(fname,"w");
    if (f==NULL) {
        printf("Could not open pdb file %s for writing\n",fname);
        die();
    }
    fprintf(f,"PSF EXT CMAP CHEQ\n");
    fprintf(f,"\n");
    fprintf(f,"%10d !NTITLE\n",1);
    fprintf(f,"%s\n",origin);
    fprintf(f,"\n");
    for (iatom=0; iatom<natom; iatom++) {
        ires=atoms[iatom].resNum;
        chain=' ';
        for (iseg=0; iseg<nseg; iseg++) if ((ires>segstart[iseg]) && (ires<segend[iseg])) {
            chain=chaincodes[iseg];
            ires=ires-segstart[iseg];
        }
        q=ffield->chargeParams[atoms[iatom].type];
        fprintf(f,atomfmt,iatom+1,chain,ires+1,atoms[iatom].resName,atoms[iatom].name,atoms[iatom].type,q,
            atoms[iatom].mass,0,0.0,0.0);
    }
    fclose(f);
}

void simulation::read_restart(char * fname)
{
    int ifrag,i,natomx,nfragx;
    double x[3],q[4];
    FILE * f;
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("Could not open restart file %s\n",fname);
        die();
    }
    fscanf(f,"%d %d %ld\n",&natomx,&nfragx,&nprevstep);
    if ((natomx!=top->natom) || (nfragx!=top->nfrag)){
        printf("Error in restart file %s.\n",fname);
        die();
    }
    printf("Reading RNG state vector from file %s.\n",fname);
    read_rng_state(f);
    /*Read the fragment, center, and quaternion.*/
    while (!feof(f)) {
    //for (i=0;i<nfrags;i++){
        fscanf(f,"%d %lg %lg %lg %lg %lg %lg %lg\n",&ifrag,&x[0],&x[1],&x[2],&q[0],&q[1],&q[2],&q[3]);
        ifrag=ifrag-1;
        normalize_quat(q);
        oldcenter[3*ifrag]=x[0];
        oldcenter[3*ifrag+1]=x[1];
        oldcenter[3*ifrag+2]=x[2];
        oldorient[4*ifrag]=q[0];
        oldorient[4*ifrag+1]=q[1];
        oldorient[4*ifrag+2]=q[2];
        oldorient[4*ifrag+3]=q[3];
        normalize_quat(&oldorient[4*ifrag]);
        //update_coords(ifrag,oldcenter,oldorient,oldcoords);
        //copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
    }
    fclose(f);
}

void simulation::write_restart(long int istep, char * fname)
{
    FILE * output;
    int ifrag,i;
    output=fopen(fname,"w");
    if (output==NULL) {
        printf("Could not open restart file %s for writing\n",fname);
        die();
    }
    printf("Writing restart file %s.\n",restartfname);
    fprintf(output,"%d %d %ld\n",top->natom,top->nfrag,istep);
    write_rng_state(output);
    for (ifrag=0; ifrag<top->nfrag; ifrag++) fprintf(output,"%d %.20f %.20f %.20f %.20f %.20f %.20f %.20f\n",
        ifrag+1,newcenter[3*ifrag],newcenter[3*ifrag+1],newcenter[3*ifrag+2],neworient[4*ifrag],neworient[4*ifrag+1],neworient[4*ifrag+2],neworient[4*ifrag+3]);
    fclose(output);
}

void write_pair_pdb2(FILE * output, fragmenttype * frag1, fragmenttype * frag2, double * coords1, double * coords2)
{
    char buf[32];
    int i,iatom,ii;
    const char * pdbcrysfmt = "CRYST1%9.3%9.3%9.3%7.2%7.2%7.2\n";
    const char * pdbatomfmt = "ATOM  %6d%4s F%02d %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n";
    /*fprintf(output,"%d\n",natom);
    fprintf(output,"step %ld\n",istep);*/
    //for (ifrag=0; ifrag<nfrag; ifrag++) {
    ii=1;
        //frag=ufragtypes[frags[ifrag]];
        for (i=0; i<frag1->natom; i++) {
            //iatom=fragtypes_tart[ifrag]+i;
            //snprintf(buf,sizeof(buf),"%s_%d",frag->names[i],ifrag+1);
            fprintf(output,pdbatomfmt,ii,frag1->names[i],0,1,coords1[3*i],coords1[3*i+1],coords1[3*i+2],1.0,0.0);
            ii++;
        }
    //}
        //frag=fragtypes[frags[jfrag]];
        for (i=0; i<frag2->natom; i++) {
            //iatom=fragstart[jfrag]+i;
            //snprintf(buf,sizeof(buf),"%s_%d",frag->names[i],jfrag+1);
            fprintf(output,pdbatomfmt,ii,frag2->names[i],1,2,coords2[3*i],coords2[3*i+1],coords2[3*i+2],1.0,0.0);
            ii++;
        }

    fprintf(output,"END\n");
    fflush(output);
}

//Fortran writes the record length as a 32-bit integer before and after writing each section of data.
//For DCDs we must emulate this behavior exactly.
void fortran_fwrite(const void * ptr, size_t size, size_t count, FILE * stream)
{
    unsigned int nbytes;
    nbytes=size*count;
    fwrite(&nbytes,sizeof(nbytes),1,stream);
    fwrite(ptr,size,count,stream);
    fwrite(&nbytes,sizeof(nbytes),1,stream);
}


void simulation::write_dcd_header(FILE * dcdfile)
{
    int i;
    dcd_header hdr;
    dcd_titles titles;
    strncpy(hdr.signature,charmm_signature,sizeof(hdr.signature));
    hdr.nframes = (nmcstep/nsave);
    hdr.begin = nprevstep + nsave; //I think this is correct for multiple files, but not absolutely sure
    hdr.skip = nsave;
    hdr.nstep = nmcstep;
    for (i=0; i<3; i++) hdr.icntrl6to8[i]=0;
    hdr.ndegf = 6*top->nfrag-6; //not necessarily exact degrees of freedom
    hdr.nfixed = 0;
    hdr.timestep_akma = 0.001f/AKMA; //Fake time step of 1 fs per MC step.
    hdr.crystal = 0; //maybe should set equal to "pbc" and write crystal data?
    hdr.dim4 = 0;
    hdr.qcg = 0;
    for (i=0; i<6; i++) hdr.icntrl14to19[i]=0;
    hdr.version = CHARMM_VERSION;
    titles.ntitle = 1;
    for (i=0; i<sizeof(titles.title); i++) titles.title[i]=' ';
    strncpy(titles.title,origin,sizeof(origin));
    fortran_fwrite(&hdr,sizeof(hdr),1,dcdfile);
    fortran_fwrite(&titles,sizeof(titles),1,dcdfile);
    fortran_fwrite(&top->natom,sizeof(top->natom),1,dcdfile);
    fflush(dcdfile);
}

void simulation::write_dcd_frame(FILE * dcdfile, double * coords)
{
    int i;
    for (i=0; i<top->natom; i++) {
        xwrite[i]=coords[3*i];
        ywrite[i]=coords[3*i+1];
        zwrite[i]=coords[3*i+2];
    }
    fortran_fwrite(xwrite,sizeof(float),top->natom,dcdfile);
    fortran_fwrite(ywrite,sizeof(float),top->natom,dcdfile);
    fortran_fwrite(zwrite,sizeof(float),top->natom,dcdfile);
    //fflush(dcdfile);
}

//allows rereading "center-quaternion" files to obtain their energies.  More exact than reading a pdb file.
//this assumes each frame is written as a block of lines equal to the number of fragments
void topology::read_frame_quat(FILE * input, long int * istep, double * center, double * orient)
{
    int ifrag,ifragx,k;
    double x[3],q[4];
    char fragname[MAX_FRAGMENT_NAME];
    for (ifragx=0; ifragx<nfrag; ifragx++){
        fscanf(input,"%ld %d %s %lg %lg %lg %lg %lg %lg %lg\n",istep,&ifrag,fragname,&x[0],&x[1],&x[2],&q[0],&q[1],&q[2],&q[3]);
        ifrag=ifrag-1;
        if ((ifrag>=nfrag) || (strncasecmp(fragname,fragtypes[frags[ifrag].type]->fragname,MAX_FRAGMENT_NAME)!=0)) {
            printf("error reading trajectory file\n");
            die();
        }
        normalize_quat(q);
        for (k=0; k<3; k++) center[3*ifrag+k]=x[k];
        for (k=0; k<4; k++) orient[4*ifrag+k]=q[k];
        normalize_quat(&orient[4*ifrag]);
    }
}


void topology::write_frame_quat(FILE * output, long int istep, double * center, double * orient)
{
    int ifrag;
    for (ifrag=0; ifrag<nfrag; ifrag++) fprintf(output,"%ld %d %s %.10f %.10f %.10f %.10f %.10f %.10f %.10f\n",
        istep,ifrag+1,fragtypes[frags[ifrag].type]->fragname,center[3*ifrag],center[3*ifrag+1],center[3*ifrag+2],orient[4*ifrag],orient[4*ifrag+1],orient[4*ifrag+2],orient[4*ifrag+3]);
    fflush(output);
}
