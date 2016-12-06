#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include "tables.h"
#include "mc.h"
#include "mt.h"
#include <math.h>
#include "rotations.h"
#include "ffield.h"
#include "util.h"
#include "fragments.h"
#if defined(__unix__) || defined(__unix) || \
        (defined(__APPLE__) && defined(__MACH__))
#define UNIX
#include <sys/types.h>
#include <unistd.h>
#endif

#if defined(PARALLEL) || defined(EXCHANGE)
simulation::simulation(const char * command, const char * fname, int _mynod, int _numnod)
#else
simulation::simulation(const char * command, const char * fname)
#endif
{
    FILE * f;
    unsigned long seed;
    double listcutoff;
    double ptot,mctemp,p,size;
    double total_table_size;
    char buffer[255],word[255],gbword[255],tablefmt[255],covtablefmt[255],fragfmt[255],struct3fname[255],factsfname[255];
    char * seq;
    char * pch;
    bool start,restart,do_mc,do_energy;
    int itype,ifrag,i,jtype,move;
    double energies[EN_TERMS],etot;
#if defined(PARALLEL) || defined(EXCHANGE)
    mynod=_mynod;
    numnod=_numnod;
#endif
    do_mc=(strcasecmp(command,"run")==0);
    do_energy=(strcasecmp(command,"energy")==0);
    tables=NULL;
    en_by_table=NULL;
    frag_nblist=NULL;
    old_per_fragment_born_radii=NULL;
    old_per_atom_born_radii=NULL;
    new_per_fragment_born_radii=NULL;
    new_per_atom_born_radii=NULL;
    facts_params=NULL;
    facts_a=NULL;
    facts_b=NULL;
    facts_numer=NULL;
    facts_denom=NULL;
    old_dg_self=NULL;
    old_sasa=NULL;
    new_dg_self=NULL;
    new_sasa=NULL;
    printf("Reading control file: %s\n",fname);
    f=fopen(fname,"r");
    if (f==NULL) {
        printf("FATAL ERROR: file %s is not found\n",fname);
        die();
    }
    //The force field
    fgets(buffer,sizeof(buffer),f);
#ifdef EXCHANGE
    sscanf(buffer,"%s %s %s\n",deffname,fragfmt,forcefieldfname);
#else
    sscanf(buffer,"%s %s %s %lg %d\n",deffname,fragfmt,forcefieldfname,&tables_lambda,&use_cov_tables);
#endif
    trim_string(forcefieldfname);
    printf("Loading parameter file %s.\n",forcefieldfname);
    ffield = new forcefield(forcefieldfname);
    trim_string(deffname);
    printf("Loading definitions file %s.\n",deffname);
    top = new topology(deffname,ffield,fragfmt);
#ifndef EXCHANGE //otherwise, this will be done in exchange_init
    if (tables_lambda<0.0) tables_lambda=0.0;
    if (tables_lambda>1.0) tables_lambda=1.0;
    //0 = fully exact, 1 = fully tabulated.
    use_std_tables=(tables_lambda>0.0);
    if ((tables_lambda>0.0) && (tables_lambda<1.0)) {
        printf("Will mix exact/tabulated energies in this simulation.  Fraction %.2f exact and %.2f tabulated.\n",(1.0-tables_lambda),tables_lambda);
    } else if (tables_lambda==0.0) {
        printf("Will calculate all nonbonded interactions exactly.\n");
    } else if (tables_lambda==1.0) {
        printf("Will use noncovalent tables in this simulation.\n");
    }
    if (use_cov_tables) printf("Will use covalent tables for peptide groups in this simulation.\n"); else printf("Will calculate all peptide covalent interactions exactly.\n");
#endif
    //Read the sequence and assemble all the topology info.
    //For now, read only one chain.  There is support for multiple chains, but need to work on chain "names", etc.
    //provide support for multiple lines
    seq=read_multiline(f);
    printf("Sequence: %s\n",seq);
    top->add_segment(' ',seq); //temporary
    free(seq);
    top->link_fragments();

    top->create_angle_dihedral_lists(use_cov_tables);
    top->create_non_tab_list(use_cov_tables,&non_tab_list);
    ffield->find_parameters(top->natom,top->atoms);
    top->create_improper_dihedral_lists(use_cov_tables,ffield);
#ifdef DEBUG_NON_TABULATED
    top->print_detailed_info();
#else
    top->print_summary_info();
#endif

    //Allocate all the coordinate arrays.
    oldcenter=ALLOCATE(3*top->nfrag,double);
    oldorient=ALLOCATE(4*top->nfrag,double);
    oldcoords=ALLOCATE(3*top->natom,double);
    newcenter=ALLOCATE(3*top->nfrag,double);
    neworient=ALLOCATE(4*top->nfrag,double);
    newcoords=ALLOCATE(3*top->natom,double);
    if (do_mc) {
        //Read start/restart and initial structure.
        fgets(buffer,sizeof(buffer),f);
        sscanf(buffer,"%s %s %s\n",word,structfname,struct2fname);
        initcoords=NULL;
        start=(strncasecmp("START",word,5)==0);
        restart=(strncasecmp("RESTART",word,7)==0);
        if (start) {
            //We're starting from a PDB file.
            printf("Will read PDB file %s and assemble fragments to fit it.\n",structfname);
            initcoords=(double *) checkalloc(3*top->natom,sizeof(double));
            top->read_pdb_file(structfname,initcoords);
            top->assemble_fragments(initcoords,oldcenter,oldorient,oldcoords);
            if (strlen(struct2fname)>0) {
                printf("Writing fitted structure to file %s.\n",struct2fname);
                top->write_pdb_file(struct2fname,oldcoords);
            }
            //PQR files: it is a format string with a "%d" that designates a fragment number for which charged fragment.
            /*if (strstr(struct2fname,"%d")!=NULL) {
                for (ifrag=0; ifrag<top->nfrag; ifrag++) {
                   snprintf(struct3fname,sizeof(struct3fname),struct2fname,ifrag);
                   printf("Writing PQR file with fragment %d charged to file %s\n",ifrag,struct3fname);
                   top->write_pqr_file(struct3fname,oldcoords,ifrag,ffield);
                }
            } else {
                printf("Writing PQR file with all fragments charged to file %s\n",struct2fname);
                top->write_pqr_file(struct2fname,oldcoords,-1,ffield);
            }*/
            nprevstep=0;
        } else if (restart) {
            initcoords=NULL;
#ifdef EXCHANGE
            snprintf(struct3fname,sizeof(struct3fname),structfname,mynod+1);//1-based replicas
            printf("Will read restart file %s.\n",struct3fname);
            read_restart(struct3fname);
#else
            printf("Will read restart file %s.\n",structfname);
            read_restart(structfname);
#endif
            printf("Number of previous steps: %ld\n",nprevstep);
            for (ifrag=0; ifrag<top->nfrag; ifrag++) top->update_coords(ifrag,oldcenter,oldorient,oldcoords);
            if (strlen(struct2fname)>0) {
                printf("Reading initial coordinates from PDB file %s.\n",struct2fname);
                initcoords=(double *) checkalloc(3*top->natom,sizeof(double));
                top->read_pdb_file(struct2fname,initcoords);
            }
        } else {
            printf("Must specify start or restart in input file.\n");
            die();
        }
        for (ifrag=0; ifrag<top->nfrag; ifrag++) top->copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
        //read the rest of MC related parameters
    	fgets(buffer,sizeof(buffer),f);
#ifdef EXCHANGE
        sscanf(buffer,"%ld %ld %ld %ld %ld %d\n",&nmcstep,&nsave,&nprint,&ncheck,&seed, &enwrite);
#else
    	sscanf(buffer,"%ld %ld %ld %ld %ld %d %lg\n",&nmcstep,&nsave,&nprint,&ncheck,&seed, &enwrite, &mctemp);
#endif
    	if (!restart) { //If restarting, we do not reinitialize the random number generator; read_restart put back the state.
    	   if (seed==0) {
       	   	seed=time(NULL);
#ifdef UNIX
                seed^=getpid(); //To ensure uniqueness even if many are started at same time.
#endif
           	printf("Initialized seed based on time.  Seed used is %ld\n",seed);
    	   } else {
           	printf("Seed specified in file.  Seed used is %ld\n",seed);
           }
    	   init_genrand(seed);
        }
    } //end of if (do_mc)
    fgets(buffer,sizeof(buffer),f);
    sscanf(buffer,"%d %lg %lg %lg\n",&pbc,&boxsize,&cutoff,&listcutoff);
    //Create nonbond list object.  Set listcutoff=cutoff or less to disable the nb list.
    use_nb_list=(listcutoff>cutoff);
    if (use_nb_list) frag_nblist=new fragment_nblist(top->nfrag,listcutoff);
    halfboxsize=boxsize*0.5;
    cutoff2=cutoff*cutoff;
    if (pbc) {
        printf("PBC is on. Box size =      %.2f A\n",boxsize);
    } else {
        printf("PBC is off.\n");
    }
    /*rdie=false;
    if (rdie) {
        printf("Distance dependent dielectric will be used.\n");
    }*/
    if (use_nb_list) {
        printf("Nonbond list will be used.\n");
        printf("Spherical/list cutoff         %.2f %.2f\n",cutoff,listcutoff);
    } else {
        printf("Nonbond list will not be used.\n");
        printf("Spherical cutoff              %.2f\n",cutoff);
    }
    fgets(buffer,sizeof(buffer),f);
    sscanf(buffer,"%s %s %lg %lg %lg %lg\n",gbword,factsfname,&gb_params.epsp,&gb_params.epsw,&gb_params.kappa,&sasa_coef);
    if (strcasecmp(gbword,"vacuum")==0) {
        gb_mode=GB_MODE_NONE;
    } else if (strcasecmp(gbword,"peratom")==0) {
        gb_mode=GB_MODE_PER_ATOM;
    } else if (strcasecmp(gbword,"peratomconv")==0) {
        gb_mode=GB_MODE_PER_ATOM_CONVERTED;
    } else if (strcasecmp(gbword,"perfrag")==0) {
        gb_mode=GB_MODE_PER_FRAGMENT;
    } else {
        printf("Invalid GB mode selection.\n"); 
        die();
    }
    switch(gb_mode) {
        case GB_MODE_NONE: printf("No solvation contribution.  Will be done in vacuum.\n"); break;
        case GB_MODE_PER_ATOM: printf("Born radii will be calculated on a per-atom basis.\n"); break;
        case GB_MODE_PER_ATOM_CONVERTED: printf("Born radii will be calculated on a per-atom basis, then converted to per-fragment Born radii.\n"); break;
        case GB_MODE_PER_FRAGMENT: printf("Born radii will be calculated on a per-fragment basis.\n"); break;
        default: printf("Invalid GB mode selection.\n"); die(); break;
    }
    if (gb_mode!=GB_MODE_NONE) {
        trim_string(factsfname);
        printf("FACTS parameter file name:               %s\n",factsfname);
        printf("Dielectric constants (solute, solvent):  %.2f %.2f\n",gb_params.epsp,gb_params.epsw);
        printf("Exponent in Still formula:               %.2f\n",gb_params.kappa);
        printf("Hydrophobic surface tension:             %.3f kcal/(mol-A^2)\n",sasa_coef);
        gb_params.tau=COUL_CONST*(1/gb_params.epsp-1/gb_params.epsw);
        read_facts_params(gb_mode!=GB_MODE_PER_FRAGMENT,factsfname,ffield,gb_params.tau);
        setup_gb_spline(gb_params.kappa);
    } else {
        printf("Dielectric constant:                     %.2f\n",gb_params.epsp);
    }
    if ((gb_mode==GB_MODE_PER_ATOM) || (gb_mode==GB_MODE_PER_ATOM_CONVERTED)) {
        old_per_atom_born_radii=ALLOCATE(top->natom,double);
        new_per_atom_born_radii=ALLOCATE(top->natom,double);
        old_dg_self=ALLOCATE(top->natom,double);
        old_sasa=ALLOCATE(top->natom,double);
        new_dg_self=ALLOCATE(top->natom,double);
        new_sasa=ALLOCATE(top->natom,double);
    }
    if ((gb_mode==GB_MODE_PER_FRAGMENT) || (gb_mode==GB_MODE_PER_ATOM_CONVERTED)) {
        old_per_fragment_born_radii=ALLOCATE(top->natom,double);
        new_per_fragment_born_radii=ALLOCATE(top->natom,double);
    }
    if (gb_mode==GB_MODE_PER_FRAGMENT) {
        old_dg_self=ALLOCATE(top->nfrag,double);
        old_sasa=ALLOCATE(top->nfrag,double);
        new_dg_self=ALLOCATE(top->nfrag,double);
        new_sasa=ALLOCATE(top->nfrag,double);
    }
    //keyword, (relative) probability, and size
    if (do_mc) {
    	for (i=1; i<=NUM_MOVES; i++) prob[i]=0.0;
    	while (TRUE) {
            fgets(buffer,sizeof(buffer),f);
            sscanf(buffer,"%s %lg %lg\n",word,&p,&size);
            move=-1;
            for (i=1; i<=NUM_MOVES; i++)
                if (strncasecmp(mc_move_names[i],word,strlen(mc_move_names[i]))==0) move=i;
            if (move<0) break;
            prob[move]=p;
            movesize[move]=size*DEG_TO_RAD;
       }
       ptot=0.0;
       for(i=1;i<=NUM_MOVES;i++)ptot+=prob[i];
       for(i=1;i<=NUM_MOVES;i++)prob[i]/=ptot;
       cumprob[1]=prob[1];
       for(i=2;i<=NUM_MOVES;i++)cumprob[i]=cumprob[i-1]+prob[i];
       for(i=1;i<=NUM_MOVES;i++)printf("%.10s moves:  Maximum size %.2f degrees  Fraction %.2f%%\n",
          mc_move_names[i],movesize[i]*RAD_TO_DEG,prob[i]*100.0);
       printf("Number of monte carlo steps: %ld\n",nmcstep);
       printf("Save and print frequencies:  %ld %ld\n",nsave,nprint);
#ifdef EXCHANGE
        fflush(stdout);
        exchange_init(f,fragfmt); //.will read the rest of the input file
    } //if (do_mc)
#else
       printf("Temperature:                 %.2f K\n",mctemp);
       beta=1/(KBOLTZ*mctemp);

       fgets(buffer,sizeof(buffer),f);
       sscanf(buffer,"%s %s %s\n",xyzfname,quatfname,restartfname);
       trim_string(xyzfname);
       trim_string(quatfname);
       trim_string(restartfname);
       if (start && (strlen(restartfname)>0)) write_restart(nprevstep,restartfname); //so that initial fragment centers/orientations are available. for analysis
    } 
    //And now... load all the tables.
    if (use_std_tables || enwrite) {
        fgets(tablefmt,sizeof(tablefmt),f); //A format string for file names.
        tables = (table * *) checkalloc(top->nfragtypes*top->nfragtypes,sizeof(table *));
        for (i=0; i<top->nfragtypes*top->nfragtypes; i++) tables[i]=NULL;
        top->load_tables(tablefmt,fragfmt,tables);
    } //else tablefmt[0]='\0';
    if (use_cov_tables || enwrite) {
        fgets(covtablefmt,sizeof(covtablefmt),f); //A format string for file names.
        covalent_tables = (covalent_table * *) checkalloc(top->nfragtypes*top->nfragtypes,sizeof(covalent_table *));
        for (i=0; i<top->nfragtypes*top->nfragtypes; i++) covalent_tables[i]=NULL;
        top->load_covalent_tables(covtablefmt,covalent_tables);
    } //else covtablefmt[0]='\0';
    total_table_size=0.0;
    if (use_std_tables || enwrite) for (i=0; i<top->nfragtypes*top->nfragtypes; i++) if (tables[i]!=NULL) total_table_size+=tables[i]->getsize();
    if (use_cov_tables || enwrite) for (i=0; i<top->nfragtypes*top->nfragtypes; i++) if (covalent_tables[i]!=NULL) total_table_size+=covalent_tables[i]->getsize();
    if (use_std_tables || use_cov_tables || enwrite) printf("Total table size:     %.2f MB\n",total_table_size);
    if (enwrite) {
        printf("Will calculate both exact and table-based energies and write them to energy.dat.\n");
        energy_output=fopen("energy.dat","w");
        pairs_output=fopen("pairs.pdb","w");
    }
#ifdef DEBUG_NON_TABULATED
    if (!do_energy) {
    	total_energy(oldcenter,oldorient,initcoords,energies,&etot);
    	print_energies(stdout,true,"init:",0,energies,etot);
    }
#endif
    //don't need this anymore, we will do the energies in do_energies
    /*if (do_energy) {
       total_energy(oldcenter,oldorient,oldcoords,energies,&etot);
       print_energies(stdout,TRUE,"Energy: ",0,energies,etot);
    }*/
#endif //EXCHANGE
    fclose(f); //Close control file. We're done with it.
}




simulation::~simulation()
{
    int i;
    if (tables!=NULL) {
        for (i=0; i<top->nfragtypes*top->nfragtypes; i++)
        	if (tables[i]!=NULL) delete tables[i];
    }
    if (en_by_table!=NULL) free(en_by_table);
    free(tables);
    free(initcoords);
    free(oldcoords);
    free(oldcenter);
    free(oldorient);
    free(newcoords);
    free(newcenter);
    free(neworient);
    if (frag_nblist!=NULL) delete frag_nblist;
    if (old_per_fragment_born_radii!=NULL) free(old_per_fragment_born_radii);
    if (old_per_atom_born_radii!=NULL) free(old_per_atom_born_radii);
    if (new_per_fragment_born_radii!=NULL) free(new_per_fragment_born_radii);
    if (new_per_atom_born_radii!=NULL) free(new_per_atom_born_radii);
    if (facts_params!=NULL) free(facts_params);
    delete ffield;
    delete top;
}


//Supervises the loading of tables. fmt is a format string with the names of fragments.
//void simulation::load_tables(char * fmt)
void topology::load_tables(const char * fmt, const char * fragfmt, table * * tables)
{
    int num_tables, count, frags_in_use, ifrag,jfrag,itype,jtype,i;
    double total_size;
    FILE * f;
    bool * need_table;
    char fname[255];
    need_table=(bool *) malloc(nfragtypes*nfragtypes*sizeof(bool));
    num_tables=0;
    for (i=0; i<nfragtypes*nfragtypes; i++) need_table[i]=false;
    for (ifrag=0; ifrag<nfrag; ifrag++)
        for (jfrag=0; jfrag<nfrag; jfrag++)
            if ((ifrag!=jfrag) && (!closefragments[ifrag*nfrag+jfrag])) {
                itype=frags[ifrag].type;
                jtype=frags[jfrag].type;
                need_table[itype*nfragtypes+jtype]=true;
            }
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=itype; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) num_tables++;
    //total_size=0.0;
    frags_in_use=0;
    for (ifrag=0; ifrag<nfragtypes; ifrag++) if (fragtypes[ifrag]->n_used>0) frags_in_use++;
    //num_tables=frags_in_use*(frags_in_use+1)/2;
    printf("Fragment types in use: %d\n",frags_in_use);
    printf("Need to load a total of %d tables.\n",num_tables);
    count=0;
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=itype; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) {
                //Maybe this code for finding the file name should be moved to table's constructor.
                count++;
                printf("\n");
                printf("--------------------------------------------------------------------------------\n");
                printf("Loading interaction table file %d of %d.\n",count,num_tables);
                //ifrag is the lesser fragment type.  Load the table! (Constructor will try both possible names.
                tables[itype*nfragtypes+jtype]=new table(fmt,fragfmt,fragtypes[itype]->fragname,fragtypes[jtype]->fragname,this);
                //tables[jtype*nfragtypes+itype]=tables[itype*nfragtypes+jtype];
                //total_size+=tables[itype*nfragtypes+jtype]->getsize();
		fflush(stdout);
            }
    free(need_table);
    printf("Total interaction tables loaded: %d\n",count);

    //printf("Total standard table size:     %.2f MB\n",total_size);
}

void topology::load_covalent_tables(const char * covtablefmt, covalent_table * * cov_tables)
{
    int num_cov_tables, count, frags_in_use, seg, ifrag,jfrag,itype,jtype,i;
    double total_size;
    FILE * f;
    bool * need_table;
    char fname[255];
    need_table=(bool *) malloc(nfragtypes*nfragtypes*sizeof(bool));
    num_cov_tables=0;
    for (i=0; i<nfragtypes*nfragtypes; i++) need_table[i]=false;
    for (seg=0; seg<nseg; seg++) {
        ifrag=first_main_chain_frag[seg];
        jfrag=frags[ifrag].main_chain_next;
        while (jfrag>0) {
            itype=frags[ifrag].type;
            jtype=frags[jfrag].type;
            //we rely on short-circuit evaluaton here
            if (use_covalent_table(itype,jtype)) need_table[itype*nfragtypes+jtype]=true;
            ifrag=jfrag;
            jfrag=frags[jfrag].main_chain_next;
        }
    }
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=0; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) num_cov_tables++;
    //total_size=0.0;
    printf("Need to load a total of %d covalent tables.\n",num_cov_tables);
    count=0;
    for (itype=0; itype<nfragtypes; itype++)
        for (jtype=0; jtype<nfragtypes; jtype++)
            if (need_table[itype*nfragtypes+jtype]) {
                count++;
                printf("\n");
                printf("--------------------------------------------------------------------------------\n");
                printf("Loading covalent table file %d of %d.\n",count,num_cov_tables);
                snprintf(fname,sizeof(fname),covtablefmt,fragtypes[itype]->fragname,fragtypes[jtype]->fragname);
                strlower(fname);
                trim_string(fname);
                //ifrag is the lesser fragment type.  Load the table! (Constructor will try both possible names.
                cov_tables[itype*nfragtypes+jtype]=new covalent_table(fname,false);
                //tables[jtype*nfragtypes+itype]=tables[itype*nfragtypes+jtype];
                //total_size+=cov_tables[itype*nfragtypes+jtype]->getsize();
		fflush(stdout);
            }
    free(need_table);
    printf("Total covalent tables loaded: %d\n",count);
    //printf("Total covalent table size:     %.2f MB\n",total_size);
}

void simulation::do_energies(char * type, char * fname,  char * enfname, char * trajfname)
{
    long int frame,istep;
    int ifrag;
    FILE * input;
    FILE * enoutput;
    FILE * trajoutput;
    double * rmsds;
    double energies[EN_TERMS],etot;
    top->chaincodes[0]='P'; //hack to make compatible with pdb files derived from catdcd
    rmsds=ALLOCATE(top->nfrag,double);
    initcoords=ALLOCATE(3*top->natom,double);
    en_by_table=ALLOCATE(top->nfragtypes*top->nfragtypes,double);
    if ((enfname==NULL) || (strlen(enfname)==0)) enoutput=stdout; else enoutput=fopen(enfname,"w");
    if (enoutput==NULL) {
        printf("Could not open energies file %s\n",enfname);
        die();
    }
    if (strcasecmp(type,"rest")==0) {  //Calculate energies from restart file. Needed for replica exchange.
        read_restart(fname);
        for (ifrag=0; ifrag<top->nfrag; ifrag++) top->update_coords(ifrag,oldcenter,oldorient,oldcoords);
        calculate_born_radii(oldcenter,oldcoords,old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa);
        total_energy(old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa,oldcenter,oldorient,oldcoords,energies,&etot);
        print_energies(stdout,false,"Energy: ",0,energies,etot);
        if ((trajfname!=NULL) && (strlen(trajfname)>0)) top->write_pqr_file(trajfname,oldcoords,-1,ffield); 
        return;
    }
    input=fopen(fname,"r");
    if (input==NULL) {
        printf("Could not open trajectory file %s\n",fname);
        die();
    }
    frame=1;
    while (!feof(input)) {
        printf("Reading frame %ld\n",frame);
        if (strcasecmp(type,"pdb")==0) {
            //No need for diagnostic messages.  Later we can implement statistics on rmsds.
            top->read_pdb_stream(input,initcoords);
            top->fit_all_fragments(initcoords,oldcenter,oldorient,oldcoords,rmsds);
        } else if (strcasecmp(type,"dat")==0) {
            read_frame_quat(input,&istep,oldcenter,oldorient);
            printf("Step number %ld\n",istep);
            for (ifrag=0; ifrag<top->nfrag; ifrag++) top->update_coords(ifrag,oldcenter,oldorient,oldcoords);
        }
        //if (feof(input)) break;
        calculate_born_radii(oldcenter,oldcoords,old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa);
        total_energy(old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa,oldcenter,oldorient,oldcoords,energies,&etot);
        print_energies(enoutput,(frame==1),"Energy: ",frame,energies,etot);
        //we could write an actual trajectory or a set of pqr files -- explore htis later
        if ((frame==1) && (trajfname!=NULL) && (strlen(trajfname)>0)) top->write_pqr_file(trajfname,oldcoords,-1,ffield);
        if (use_std_tables) print_energies_by_table();
        if (feof(input)) break;
        frame++;
    }
    fclose(input);
    fclose(enoutput);
    free(rmsds);
}

