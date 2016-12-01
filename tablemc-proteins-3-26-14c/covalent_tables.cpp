#include <stdio.h>
#include "covalent_tables.h"
#include "ffield.h"
#include "fragments.h"
#include "rotations.h"
#include "util.h"

#define TABLESIZE 131072

covalent_table::covalent_table(const char * fname, int new_table)//, int newtable)
{
    FILE * f;
    energy=NULL;
    hdr.totalpoints=0;
    if (new_table) {
        generate_table(fname);
    } else {
        f=fopen(fname,"rb");
        if (f==NULL) { //still couldn't open?
            printf("Could not open table file %s.\n",fname);
            exit(1);
        }
        alloc_read_table(f,fname);
        fclose(f);
    }
}

covalent_table::~covalent_table()
{
    if (energy!=NULL) free(energy);
}


void covalent_table::alloc_read_table(FILE * f, const char * fname)
{
    printf("Loading covalent table file %s:\n",fname);
    //printf("*** %d\n",sizeof(*hdr));
    fread(&hdr,sizeof(hdr),1,f);
    if (ferror(f)) {
        printf("Error reading table file %s\n",fname);
        exit(1);
    }
    if ((hdr.headersize!=sizeof(hdr))) {
        printf("File format error in table file %s\n",fname);
        exit(1);
    }
    if ((sizeof(float)==sizeof(cov_energy_t)) && !(hdr.options & OP_SINGLE_PRECISION)) {
        printf("This table uses double precision numbers, but program is compiled for single precision tables.\n");
        exit(1);
    }
    if ((sizeof(double)==sizeof(cov_energy_t)) && !(hdr.options & OP_SINGLE_PRECISION)) {
        printf("This table uses single precision numbers, but program is compiled for double precision tables.\n");
        exit(1);
    }
    print_header_info();
    /*clash_table = (double *) malloc(hdr.clashpoints*sizeof(double));
    if (clash_table == NULL) {
        printf("Could not allocate memory for clash table in file %s\n",fname);
        exit(1);
    }*/
    energy = (cov_energy_t *) malloc(hdr.totalpoints*sizeof(cov_energy_t)); //generate
    if (energy == NULL) {
        printf("Could not allocate memory for table in file %s\n",fname);
        exit(1);
    }
    if (ferror(f)) {
        printf("Error reading table file %s\n",fname);
        exit(1);
    }
    //fread(clash_table,sizeof(double),hdr.clashpoints,f);
    fread(energy,sizeof(cov_energy_t),hdr.totalpoints,f);
    //fclose(f);
}

void covalent_table::write_table(const char * fname)
{
    FILE * f;
    f=fopen(fname,"wb");
    if (f==NULL) {
        printf("Could not open table file %s for writing\n",fname);
        exit(1);
    }
    fwrite(&hdr,sizeof(hdr),1,f);
    //fwrite(clash_table,sizeof(double),hdr.clashpoints,f);
    fwrite(energy,sizeof(cov_energy_t),hdr.totalpoints,f);
    fclose(f);
    printf("Table written to file %s.\n",fname);
}

double covalent_table::getsize(void)
{
    return ((double) hdr.totalpoints*sizeof(cov_energy_t))/((double) 1024*1024);
}

double covalent_table::get_energy(double phi, double psi, double alpha)
{
    long int iphi,ipsi,ialpha;
    long long int index;
    double en;
    /*phi=-M_PI+(iphi+0.5)*hdr.dphi;
                psi=-M_PI+(ipsi+0.5)*hdr.dpsi;
                alpha=hdr.minalpha+(ialpha+0.5)*hdr.dalpha;*/
    if ((alpha>=hdr.maxalpha)||(alpha<=hdr.minalpha)) return DUMMY_ENERGY;
    iphi=(long int) ((phi+M_PI)/hdr.dphi);
    if (iphi>=hdr.nphi) iphi=hdr.nphi-1;
    ipsi=(long int) ((psi+M_PI)/hdr.dpsi);
    if (ipsi>=hdr.npsi) ipsi=hdr.npsi-1;
    ialpha=(long int) ((alpha-hdr.minalpha)/hdr.dalpha);
    if (ialpha>=hdr.nalpha) ialpha=hdr.nalpha-1;
    index=calculate_index(iphi,ipsi,ialpha);
    en=energy[index];
    return en;
}

double covalent_table::table_interaction_energy(double * coords, int ic1, int in1, int ica1, int ic2, int in2, int ica2)
{
    int k;
    double rnc1[3],rnca1[3],rcan1[3],rcac1[3],rcn2[3],rnca2[3];
    double phi, psi, alpha;
    for (k=0; k<3; k++) {
        rnc1[k]=coords[3*ic1+k]-coords[3*in1+k];
        rnca1[k]=coords[3*ica1+k]-coords[3*in1+k];
        rcan1[k]=-rnca1[k];//Because of the way dihedral is written, we need both ways for this vector.
        rcac1[k]=coords[3*ic2+k]-coords[3*ica1+k];
        rcn2[k]=coords[3*in2+k]-coords[3*ic2+k];
        //rnca2[k]=coords[3*ica2+k]-coords[3*in2+k];
    }
    phi=dihedral(rnc1,rnca1,rcac1);
    psi=dihedral(rcan1,rcac1,rcn2);
    alpha=angle(rcan1,rcac1);
    return get_energy(phi,psi,alpha);
}

void covalent_table::print_header_info(void)
{
    char buffer[255];
    double mem;
    int iregion;
    printf("Phi and psi resolution:                           %.2f %.2f degrees\n",hdr.dphi*RAD_TO_DEG,hdr.dpsi*RAD_TO_DEG);
    /*printf("Omega range:                                      (+/-) %.2f-180.00 degrees\n",hdr.minomega*RAD_TO_DEG);
    printf("Omega resolution:                                 %.2f degrees\n",hdr.domega*RAD_TO_DEG);*/
    printf("N-Ca-C angle range:                               %.2f-%.2f degrees\n",hdr.minalpha*RAD_TO_DEG,hdr.maxalpha*RAD_TO_DEG);
    printf("N-Ca-C angle resolution:                          %.2f degrees\n",hdr.dalpha*RAD_TO_DEG);
    printf("Number of N-Ca-C angle points:                    %ld\n",hdr.nalpha);
    printf("Ca-C bond length:                                 %.2f A\n",hdr.rcac);
    printf("Ca-C-N angle:                                     %.2f degrees\n",hdr.thetacacn*RAD_TO_DEG);
    printf("Omega angle:                                      %.2f degrees\n",hdr.omega*RAD_TO_DEG);
    /*printf("Ca-C-N angle range:                               %.2f-%.2f degrees\n",hdr.minthetac*RAD_TO_DEG,hdr.maxthetac*RAD_TO_DEG);
    printf("Ca-C-N angle resolution:                          %.2f degrees\n",hdr.dthetac*RAD_TO_DEG);
    printf("Ca-C bond length range:                           %.2f-%.2f A\n",hdr.minr,hdr.maxr);
    printf("Ca-C bond length resolution:                      %.2f A\n",hdr.dr);
    printf("Number of phi/psi points:                         %ld %ld\n",hdr.nphi,hdr.npsi);
    printf("Number of omega points:                           %ld\n",hdr.nomega);
    printf("Number of N-Ca-C angle, Ca-C-N angle points:      %ld %ld\n",hdr.nalpha,hdr.nthetac);
    printf("Number of Ca-C bond length points:                %ld\n",hdr.nr);*/
    printf("Total points:                                     %lld\n",hdr.totalpoints);
    printf("Table size:                                       %.2f MB\n",getsize());
    printf("Dielectric constant:                              %.2f\n",hdr.eps);
    printf("Definitions file used for creation:               %s\n",hdr.defs);
    printf("Sequence used for creation:                       %s\n",hdr.sequence);
    printf("Force field:                                      %s\n",hdr.forcefield);
    printf("Fragment names:                                   %s %s\n",hdr.fragtype1,hdr.fragtype2);
    snprintf(buffer,sizeof(buffer),hdr.fragfmt,hdr.fragtype1);
    strlower(buffer);
    printf("Fragment file 1:                                  %s\n",buffer);
    snprintf(buffer,sizeof(buffer),hdr.fragfmt,hdr.fragtype2);
    strlower(buffer);
    printf("Fragment file 2:                                  %s\n",buffer);
    strncpy(buffer,ctime(&hdr.created),sizeof(buffer));
    printf("This table file created on:                       %s\n",buffer);
    if ((hdr.options & OP_BOLTZMANN_AVERAGED)!=0) {
        printf("Boltzmann averaging was used in the construction of this table.\n");
        //printf("Number of cells averaged (r, anguler, Euler):     %d %d %d\n",hdr.nsmoothr,hdr.nsmoothsph,hdr.nsmootheuler);
        //printf("Values of e1 and alpha1:                          %.2f %.2f kcal/mol\n",hdr.e1,hdr.alpha1);
        //printf("Values of e2 and alpha2:                          %.2f %.2f kcal/mol\n",hdr.e2,hdr.alpha2);
        //printf("Energy limit for smoothing:                       %.2f kcal/mol\n",hdr.en_smooth_limit);
        printf("Temperature for Boltzmann averaging:              %.2f K\n",1/(KBOLTZ*hdr.smooth_beta));
        printf("Regions for smoothing:\n");
        printf("Name  phi range (degrees)  psi range (degrees)  Angular scale:\n");
        for (iregion=0; iregion<hdr.n_smooth_region; iregion++)
            printf("%8s %.2f %.2f %.2f %.2f %.2f\n",hdr.smooth_regions[iregion].name,
                hdr.smooth_regions[iregion].philo*RAD_TO_DEG,hdr.smooth_regions[iregion].phihi*RAD_TO_DEG,
                hdr.smooth_regions[iregion].psilo*RAD_TO_DEG,hdr.smooth_regions[iregion].psihi*RAD_TO_DEG,
                hdr.smooth_regions[iregion].smooth_scale*RAD_TO_DEG);
        printf("Energy limit for smoothing disabled this version.\n");
        //printf("Seed used:                                   %ld\n",hdr.seed);
    }
    if ((hdr.options & OP_DIST_DEP_DIELECTRIC)!=0) {
        printf("This table used a distance dependent dielectric.\n");
    }
    if ((hdr.options & OP_SINGLE_PRECISION)!=0) {
        printf("This table contains single precision numbers.\n");
    }
}


void covalent_table::read_table_header_info(const char * fname)
{
    /*char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != NULL)
        printf("Current working dir: %s\n", cwd);*/
    char buffer[255];
    int rdie,iregion;
    hdr.headersize=sizeof(hdr);
    hdr.options=0;
    time(&hdr.created); /*current date/time*/
    printf("Reading control file: %s\n",fname);
    FILE * f;
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("FATAL ERROR: file %s is not found\n",fname);
        exit(1);
    }
    fscanf(f,"%lg %lg %lg %lg %lg\n",&hdr.dphi,&hdr.dpsi,&hdr.minalpha,&hdr.maxalpha,&hdr.dalpha);
    hdr.dphi=hdr.dphi*DEG_TO_RAD;
    hdr.dpsi=hdr.dpsi*DEG_TO_RAD;
    hdr.nphi=floor(2*M_PI/hdr.dphi);
    hdr.npsi=floor(2*M_PI/hdr.dpsi);
    hdr.dphi=2*M_PI/hdr.nphi;
    hdr.dpsi=2*M_PI/hdr.npsi;

    //fscanf(f,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg\n",&hdr.minomega,&hdr.domega,&hdr.minalpha,&hdr.maxalpha,&hdr.dalpha,
    //    &hdr.minthetac,&hdr.maxthetac,&hdr.dthetac,&hdr.minr,&hdr.maxr,&hdr.dr);
    /*hdr.minomega=hdr.minomega*DEG_TO_RAD;
    hdr.domega=hdr.domega*DEG_TO_RAD;
    hdr.nomega=floor((M_PI-hdr.minomega)/hdr.domega);
    hdr.domega=(M_PI-hdr.minomega)/hdr.nomega;
    hdr.nomega=hdr.nomega*2; //since we have both positive and negative omega*/
    hdr.minalpha=hdr.minalpha*DEG_TO_RAD;
    hdr.maxalpha=hdr.maxalpha*DEG_TO_RAD;
    hdr.dalpha=hdr.dalpha*DEG_TO_RAD;
    hdr.nalpha=floor((hdr.maxalpha-hdr.minalpha)/hdr.dalpha)+1;
    hdr.dalpha=(hdr.maxalpha-hdr.minalpha)/hdr.nalpha;
    /*hdr.minthetac=hdr.minthetac*DEG_TO_RAD;
    hdr.maxthetac=hdr.maxthetac*DEG_TO_RAD;
    hdr.dthetac=hdr.dthetac*DEG_TO_RAD;
    hdr.nthetac=floor((hdr.maxthetac-hdr.minthetac)/hdr.dthetac);
    hdr.dthetac=(hdr.maxthetac-hdr.minthetac)/hdr.nthetac;
    hdr.minr=hdr.minr*DEG_TO_RAD;
    hdr.maxr=hdr.maxr*DEG_TO_RAD;
    hdr.dr=hdr.dr*DEG_TO_RAD;
    hdr.nr=floor((hdr.maxr-hdr.minr)/hdr.dr);
    hdr.dr=(hdr.maxr-hdr.minr)/hdr.nr;*/
    hdr.totalpoints=hdr.nphi*hdr.npsi*hdr.nalpha; //hdr.nomega*hdr.nalpha*hdr.nthetac*hdr.nr;
    fscanf(f,"%lg %lg %lg\n",&hdr.rcac,&hdr.thetacacn,&hdr.omega);
    hdr.thetacacn=hdr.thetacacn*DEG_TO_RAD;
    hdr.omega=hdr.omega*DEG_TO_RAD;
    //hdr.totalpoints=(long long int) hdr.ntrans * (long long int) hdr.norient;
    fscanf(f,"%lg %d %d %lg\n",&hdr.eps,&rdie,&hdr.n_smooth_region,&hdr.smooth_beta); //, &hdr.n_invalid_dev, &hdr.en_clash,&hdr.en_invalid_margin);
    hdr.smooth_beta=1/(KBOLTZ*hdr.smooth_beta);
    if (hdr.n_smooth_region>0) {
        hdr.options = hdr.options | OP_BOLTZMANN_AVERAGED;
        for (iregion=0; iregion<hdr.n_smooth_region; iregion++) read_region_info(f,&hdr.smooth_regions[iregion]);
    }
    if (rdie) hdr.options = hdr.options | OP_DIST_DEP_DIELECTRIC;
    if (sizeof(cov_energy_t)==sizeof(float)) hdr.options = hdr.options | OP_SINGLE_PRECISION;
    fgets(buffer,sizeof(buffer),f);
    sscanf(buffer,"%s %s %s\n",hdr.defs,hdr.forcefield,hdr.fragfmt);
    trim_string(hdr.defs);
    trim_string(hdr.forcefield);
    fgets(buffer,sizeof(buffer),f);
    strncpy(hdr.sequence,buffer,sizeof(hdr.sequence));
    fclose(f);
}
//We allow iphihi and ipsihi to be out of range, so that we can wrap around the edges of the ramachandran plot
//when smoothing.
void covalent_table::read_region_info(FILE * f, rama_region * rgn)
{
    fscanf(f,"%8s %lg %lg %lg %lg %lg\n",&rgn->name,&rgn->philo,&rgn->phihi,
        &rgn->psilo,&rgn->psihi,&rgn->smooth_scale);
    rgn->philo=rgn->philo*DEG_TO_RAD;
    rgn->phihi=rgn->phihi*DEG_TO_RAD;
    rgn->psilo=rgn->psilo*DEG_TO_RAD;
    rgn->psihi=rgn->psihi*DEG_TO_RAD;
    rgn->smooth_scale=rgn->smooth_scale*DEG_TO_RAD;
    rgn->iphilo=round((rgn->philo+M_PI)/hdr.dphi);
    rgn->iphihi=round((rgn->phihi+M_PI)/hdr.dphi)-1;
    rgn->philo=-M_PI+(rgn->iphilo)*hdr.dphi;
    rgn->phihi=-M_PI+(rgn->iphihi-1)*hdr.dphi;
    if (rgn->iphilo>rgn->iphihi) rgn->iphihi+=hdr.nphi;
    rgn->ipsilo=round((rgn->psilo+M_PI)/hdr.dpsi);
    rgn->ipsihi=round((rgn->psihi+M_PI)/hdr.dpsi)-1;
    rgn->psilo=-M_PI+(rgn->ipsilo)*hdr.dpsi;
    rgn->psihi=-M_PI+(rgn->ipsihi-1)*hdr.dpsi;
    if (rgn->ipsilo>rgn->ipsihi) rgn->ipsihi+=hdr.npsi;
}

#ifdef DEBUG
//identical to simulation::print_energies, but needed here because there is no simulation
void covalent_table::print_energies(int hdr, const char * title, long int istep, double * energies, double etot)
{
    const char * hdrfmt = "%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n";
    const char * enfmt = "%15s %15ld %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n";
    if (hdr) {
        printf(hdrfmt,title,"Step","Total","Interaction","Bonds","Angles","Dihedrals","Impropers","Non-tab VDW","Non-tab Elec");
    }
    printf(enfmt,title,istep,etot,energies[EN_INTERACTION],energies[EN_BOND],energies[EN_ANGLE],energies[EN_DIHEDRAL],energies[EN_IMPROPER],energies[EN_VDW_EXACT],energies[EN_ELEC_EXACT]);
}
#endif

void covalent_table::fill_table(topology * top, forcefield * ffield,fragmenttype * fragtype1,fragmenttype * fragtype2, int nb_atom_list_size, atom_nb_entry * nb_atom_list)
{
    double * coords;
#ifdef DEBUG
    double * coords2;
    double q[8],x[6],rmsd;
#else
    double q[4],x[3];
#endif
    int k;
    int rdie;
    int ic1,in1,ica1,ic2,io2,ih2,in2, ica2,icf,icaf,inf,iof,ihf;
    long int iphi,ipsi,ialpha;
    long long int points, index;
    double r, phi, psi, alpha, thetacacn,omega;
    double rcn,rnca,rco,rnh,thetanco,thetacnca,thetacnh,dihocnca,dihocnh;
    double energies[EN_TERMS],etot;
    bool * movedatoms;
    rdie=(hdr.options & OP_DIST_DEP_DIELECTRIC)!=0;
    //save the pointers to the fragment types
    //Find the atoms we need to use for building.
    ic1=top->find_atom(0,"C");
    in1=top->find_atom(0,"N");
    ica1=top->find_atom(0,"CA");
    ic2=top->find_atom(1,"C");
    in2=top->find_atom(1,"N");
    ica2=top->find_atom(1,"CA");
    io2=top->find_atom(1,"O");
    ih2=top->find_atom(1,"H"); //may not exist if we are doing peptide-pro
    //these indices are relative to fragment #2, rather than whole system
    icf=top->atoms[ic2].fragatom;
    inf=top->atoms[in2].fragatom;
    icaf=top->atoms[ica2].fragatom;
    iof=top->atoms[io2].fragatom;
    if (ih2>0) ihf=top->atoms[ih2].fragatom; else ihf=-1;
    coords=(double *) checkalloc(3*top->natom,sizeof(double));
    movedatoms=(bool *) checkalloc(top->natom,sizeof(bool));
    //We use the moved_non_tabulated_energy function to evaluate energies,
    //passing it a "moved atoms" array in which fragment 2 has moved and fragment 1 has not.
    //In this way we obtain the total interaction energy between the two fragments.
    //This doesn't include the bond energy (but this doesn't matter anyway since we use the equilibrium bond length)
    for (k=0; k<fragtype1->natom; k++) movedatoms[k]=false;
    for (k=fragtype1->natom; k<top->natom; k++) movedatoms[k]=true;
    for (k=0; k<3*top->natom; k++) coords[k]=0.0;
    //Set first fragment at the origin, with reference orientation.
    for (k=0; k<3; k++) x[k]=0.0;
    q[0]=1.0;
    q[1]=0.0;
    q[2]=0.0;
    q[3]=0.0;
    //This fills in the first fragment's coordinates.
    top->update_coords(0,x,q,coords);
    //measure some bond distances and angles in the second fragment
    //We should have the fragment measure this from its reference geometry
    rcn=fragtype2->get_bond_length(icf,inf);
    rnca=fragtype2->get_bond_length(inf,icaf);
    rco=fragtype2->get_bond_length(icf,iof);
    thetacnca=fragtype2->get_angle(icf,inf,icaf);
    thetanco=fragtype2->get_angle(iof,icf,inf);
    dihocnca=fragtype2->get_dihedral(iof,icf,inf,icaf);
    if (ihf>0) {
        rnh=fragtype2->get_bond_length(inf,ihf);
        thetacnh=fragtype2->get_angle(icf,inf,ihf);
        dihocnh=fragtype2->get_dihedral(iof,icf,inf,ihf);
    }
    //set up parameters (testing)
    //phi=-180.0*DEG_TO_RAD;
    //psi=-180.0*DEG_TO_RAD;
    /*phi=-79.0*DEG_TO_RAD;
    psi=139.0*DEG_TO_RAD;
    alpha=117.50*DEG_TO_RAD;*/
    points=0;
/*#ifdef DEBUG
    iphi=15; ipsi=105; ialpha=11;
#else*/
    for (iphi=0; iphi<hdr.nphi; iphi++)
        for (ipsi=0; ipsi<hdr.npsi; ipsi++)
            for (ialpha=0; ialpha<hdr.nalpha; ialpha++) {
//#endif
                phi=-M_PI+(iphi+0.5)*hdr.dphi;
                psi=-M_PI+(ipsi+0.5)*hdr.dpsi;
                alpha=hdr.minalpha+(ialpha+0.5)*hdr.dalpha;
                build_atom(coords,ic1,in1,ica1,ic2,hdr.rcac,alpha,phi);
                build_atom(coords,in1,ica1,ic2,in2,rcn,hdr.thetacacn,psi);
                build_atom(coords,ica1,ic2,in2,ica2,rnca,thetacnca,hdr.omega); //assume trans configuration (omega=180)
                build_atom(coords,ica2,in2,ic2,io2,rco,thetanco,dihocnca);
                if (ih2>0) build_atom(coords,io2,ic2,in2,ih2,rnh,thetacnh,dihocnh);
#ifdef DEBUG
                /*top->write_pdb_file("test.pdb",coords);
                coords2=(double *) checkalloc(3*top->natom,sizeof(double));
                top->assemble_fragments(coords,x,q,coords2);
                top->write_pdb_file("test2.pdb",coords2);
                free(coords2);*/
#endif
                ffield->moved_non_tabulated_energy(hdr.eps,rdie,top->natom,top->atoms,movedatoms,nb_atom_list_size,nb_atom_list,coords,energies);
                etot=0.0;
                for (k=EN_BOND;k<EN_TERMS;k++) etot+=energies[k];
                index=calculate_index(iphi,ipsi,ialpha);
                energy[index]=etot;
#ifdef DEBUG
                //print_energies(true,"Int:",points,energies,etot);
#endif
                points++;
                if ((points%10000)==0) printf("%lld points completed\n",points);
//#ifndef DEBUG
            }
//#endif
    free(coords);
    free(movedatoms);
}

double covalent_table::partition_function(double beta,rama_region rgn)
{
    double sum, dv, totalvolume;
    long int iphi, ipsi, ialpha;
    long long int index;
    //this is incorrect
    totalvolume=4.0*M_PI*M_PI*(hdr.maxalpha-hdr.minalpha);
    dv=hdr.dphi*hdr.dpsi*hdr.dalpha;
    sum=0.0;
    for (ialpha=0; ialpha<hdr.nalpha; ialpha++)
        for (iphi=rgn.iphilo; iphi<=rgn.iphihi; iphi++)
            for (ipsi=rgn.ipsilo; ipsi<=rgn.ipsihi; ipsi++) {
                index=calculate_index(iphi%hdr.nphi,ipsi%hdr.npsi,ialpha);
                sum+=dv*exp(-beta*energy[index]);
            }
    return -(1.0/beta)*log(sum/totalvolume);
}


void covalent_table::boltzmann_average(rama_region rgn)
{
    double * prob;
    long int iphi, ipsi, ialpha, ineighphi, ineighpsi;
    long long int index, neighborindex, points;
    double sum, weightsum, weight,coeff,betainv,maxdist2;
    double phi, psi, neighphi, neighpsi, diffphi, diffpsi,dist2,ddist;
    double kernel[TABLESIZE];
    int bin;
    coeff=1/(rgn.smooth_scale*rgn.smooth_scale);
    ddist=2.0*M_PI*M_PI/TABLESIZE;
    for (bin=0; bin<TABLESIZE; bin++) {
        dist2=(bin+0.5)*ddist;
        kernel[bin]=exp(-coeff*dist2);
    }
    prob=(double *) malloc(hdr.totalpoints*sizeof(double));
    if (prob==NULL) {
        printf("Could not allocate space for Boltzmann factors.\n");
        exit(1);
    }
    printf("Computing Boltzmann factors...\n");
    for (index=0; index<hdr.totalpoints; index++) {
        if (energy[index]>50.0) prob[index]=0.0; else prob[index]=exp(-hdr.smooth_beta*energy[index]);
    }
    betainv=1.0/hdr.smooth_beta;
    points=0;
    printf("%ld %ld %ld %ld\n",rgn.iphilo,rgn.iphihi%hdr.nphi,rgn.ipsilo,rgn.ipsihi%hdr.npsi);
    for (ialpha=0; ialpha<hdr.nalpha; ialpha++)
        for (iphi=rgn.iphilo; iphi<=rgn.iphihi; iphi++)
            for (ipsi=rgn.ipsilo; ipsi<=rgn.ipsihi; ipsi++) {
                sum=0.0;
                weightsum=0.0;
                //iphi and ipsi may be out of range (see read_region_info) so we must reduce modulo
                //in order to wrap around.
                phi=-M_PI+(iphi%hdr.nphi+0.5)*hdr.dphi;
                psi=-M_PI+(ipsi%hdr.npsi+0.5)*hdr.dpsi;
                index=calculate_index(iphi%hdr.nphi,ipsi%hdr.npsi,ialpha);
		//if (ialpha==0) printf("%d %d\n",iphi%hdr.nphi,ipsi%hdr.npsi);
                for (ineighphi=rgn.iphilo; ineighphi<=rgn.iphihi; ineighphi++)
                    for (ineighpsi=rgn.ipsilo; ineighpsi<=rgn.ipsihi; ineighpsi++) {
                        neighborindex=calculate_index(ineighphi%hdr.nphi,ineighpsi%hdr.npsi,ialpha);
			//printf("%d %d %d %d\n",iphi%hdr.nphi,ipsi%hdr.npsi,ineighphi%hdr.nphi,ineighpsi%hdr.npsi);
                        neighphi=-M_PI+(ineighphi%hdr.nphi+0.5)*hdr.dphi;
                        neighpsi=-M_PI+(ineighpsi%hdr.npsi+0.5)*hdr.dpsi;
                        diffphi=neighphi-phi;
                        diffpsi=neighpsi-psi;
                        if (diffphi<-M_PI) diffphi+=2.0*M_PI;
                        if (diffphi>M_PI) diffphi-=2.0*M_PI;
                        if (diffpsi<-M_PI) diffpsi+=2.0*M_PI;
                        if (diffpsi>M_PI) diffpsi-=2.0*M_PI;
                        dist2=diffphi*diffphi+diffpsi*diffpsi;
                        bin=(int)(dist2/ddist);
                        weight=kernel[bin];
                        /*if ((index==1)||(neighborindex==1)) {
                                printf("%lld %lld %.4f %d %.4g\n",index,neighborindex,dist2,bin,weight);
                        }*/
                        sum+=weight*prob[neighborindex];
                        weightsum+=weight;
                    }
                sum=sum/weightsum;
                if (sum<0) {
                    printf("Error in smoothing for indices %ld %ld %ld\n",iphi,ipsi,ialpha);
                    exit(1);
                } else if (sum>0) energy[index]=-betainv*log(sum);
                points++;
                if ((points%1000)==0) {
			printf("%lld points averaged\n",points);
			fflush(stdout);
		}
    }
    free(prob);
}


/*void covalent_table::smooth(void)
{
    long int iphi, ipsi, ialpha;
    double enmin,e1,e2,en,a;
    long long int points,index;
    points=0;
    for (ialpha=0; ialpha<hdr.nalpha; ialpha++) {
        enmin=DUMMY_ENERGY;
        for (iphi=0; iphi<hdr.nphi; iphi++)
            for (ipsi=0; ipsi<hdr.npsi; ipsi++) {
                index=calculate_index(iphi,ipsi,ialpha);
                en=energy[index];
                if (en<enmin) enmin=en;
            }
        e1=enmin+hdr.e1;
        e2=enmin+hdr.e2;
        for (iphi=0; iphi<hdr.nphi; iphi++)
            for (ipsi=0; ipsi<hdr.npsi; ipsi++) {
		if ((iphi==115) && (ipsi==234) && (ialpha==11)) {
			printf("debugging pad\n");
		}
                index=calculate_index(iphi,ipsi,ialpha);
                en=energy[index];
                if (en>e1) {
			a=(e2-en)/hdr.alpha2;
			//printf("%.10f %.10f %.10f %.10f\n",en,e1,e2,a);
			if (a>-50.0) energy[index]=en-((en-e1)*(en-e1))/((hdr.alpha1+en-e1)*(1.0+exp(-a)));
		}
                points++;
                if ((points%1000)==0) {
                    printf("%lld points averaged\n",points);
                    fflush(stdout);
                }
            }
    }
}*/


void covalent_table::generate_table(const char * control_file)
{
    topology * top;
    forcefield * ffield;
    std::vector<atom_nb_entry> non_tab_list;
    double beta;
    int iregion;
    printf("Generating covalent table.\n");
    read_table_header_info(control_file);
    ffield = new forcefield(hdr.forcefield);
    top = new topology(hdr.defs,ffield,hdr.fragfmt);
    top->add_segment(' ',hdr.sequence);
    top->link_fragments(); //probably don't need this, but just in case
    top->create_angle_dihedral_lists(false);
    top->create_non_tab_list(false,&non_tab_list);
    ffield->find_parameters(top->natom,top->atoms);
    top->create_improper_dihedral_lists(false,ffield);
    top->print_detailed_info();
    //We assume the structure we built contains exactly two fragments.
    if (top->nfrag!=2) {
        printf("Wrong number of fragments in covalent table generation.\n");
        exit(1);
    }
    fragmenttype * fragtype1;
    fragmenttype * fragtype2;
    fragtype1=top->fragtypes[top->frags[0].type];
    fragtype2=top->fragtypes[top->frags[1].type];
    strncpy(hdr.fragtype1,fragtype1->fragname,sizeof(hdr.fragtype1));
    strncpy(hdr.fragtype2,fragtype2->fragname,sizeof(hdr.fragtype2));
    print_header_info();
    energy=(cov_energy_t *) malloc(hdr.totalpoints*sizeof(cov_energy_t));
    if (energy == NULL) {
        printf("Could not allocate memory for table\n");
        exit(1);
    }
    fill_table(top,ffield,fragtype1,fragtype2,non_tab_list.size(),&non_tab_list[0]);
    beta=1/(KBOLTZ*300.0);
    //printf("Partition function at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*beta),partition_function(beta));
    if (hdr.options & OP_BOLTZMANN_AVERAGED){
        //smooth();
        for (iregion=0; iregion<hdr.n_smooth_region; iregion++) {
            printf("Smoothing region %8s\n",hdr.smooth_regions[iregion].name);
            printf("Partition function at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*beta),partition_function(beta,hdr.smooth_regions[iregion]));
            boltzmann_average(hdr.smooth_regions[iregion]);
            printf("Partition function at temperature %.2f K: %.6f kcal/mol\n",1/(KBOLTZ*beta),partition_function(beta,hdr.smooth_regions[iregion]));
        }
    }
    //free(coords2);
    //free(weights);
    delete top;
    delete ffield;
}

void covalent_table::write_energies(double alpha, const char * fname)
{
    FILE * output;
    long int iphi,ipsi,ialpha;
    long long int index;
    double phi,psi;
    output=fopen(fname,"w");
    alpha=alpha*DEG_TO_RAD;
    ialpha=(long int) ((alpha-hdr.minalpha)/hdr.dalpha);
    for (iphi=0; iphi<hdr.nphi; iphi++)
        for (ipsi=0; ipsi<hdr.npsi; ipsi++) {
            phi=-M_PI+(iphi+0.5)*hdr.dphi;
            psi=-M_PI+(ipsi+0.5)*hdr.dpsi;
            index=calculate_index(iphi,ipsi,ialpha);
            fprintf(output,"%.4f %.4f %.4f\n",phi*RAD_TO_DEG,psi*RAD_TO_DEG,energy[index]);
    }
    printf("Energies for plotting written to file %s.\n",fname);
    fclose(output);
}




