#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#if defined(PARALLEL) || defined(EXCHANGE)
#include <mpi.h>
#endif
#include "mc.h"
//#include "tables.h"
#include "mt.h"
#include "rotations.h"
#include "util.h"
#include "covalent_tables.h"
#include "solvation.h"

double topology::exact_interaction_energy(forcefield * ffield, int pbc, double halfboxsize, double boxsize, int gb_mode, gb_param_info * gb_params,
    int ifrag, int jfrag, double * coords, double * born_radii)
{
    double en;
    double tempcoords1[3*MAX_ATOMS_PER_FRAGMENT],tempcoords2[3*MAX_ATOMS_PER_FRAGMENT];
    double born_radii1[MAX_ATOMS_PER_FRAGMENT],born_radii2[MAX_ATOMS_PER_FRAGMENT],a1a2;
    fragmenttype * ifragtype;
    fragmenttype * jfragtype;
    int iatom, iactualatom,k;
#ifdef TIMERS
    switch_timer(TIMER_INT_EXACT_PREP);
#endif
    for (iatom=0; iatom<3*MAX_ATOMS_PER_FRAGMENT; iatom++) {
         tempcoords1[iatom]=NAN; tempcoords2[iatom]=NAN;
    }
    //printf("Doing exact interaction, ifrag=%d jfrag=%d gb_mode=%d\n",ifrag,jfrag,gb_mode);
    ifragtype=fragtypes[frags[ifrag].type];
    jfragtype=fragtypes[frags[jfrag].type];
    for (iatom=0; iatom<ifragtype->natom; iatom++) {
        iactualatom=frags[ifrag].atoms[iatom];
        for (k=0; k<3; k++) tempcoords1[3*iatom+k]=coords[3*iactualatom+k];
        if (gb_mode==GB_MODE_PER_ATOM) born_radii1[iatom]=born_radii[iactualatom];
    }
    for (iatom=0; iatom<jfragtype->natom; iatom++) {
        iactualatom=frags[jfrag].atoms[iatom];
        for (k=0; k<3; k++) tempcoords2[3*iatom+k]=coords[3*iactualatom+k];
        if (gb_mode==GB_MODE_PER_ATOM) born_radii2[iatom]=born_radii[iactualatom];
    }
#ifdef TIMERS
    switch_timer(TIMER_INT_EXACT);
#endif
    switch(gb_mode) {
        case GB_MODE_NONE:
            en=ffield->exact_interaction_energy(pbc,halfboxsize,boxsize,gb_params->epsp,ifragtype->natom,&ifragtype->types[0],&tempcoords1[0],
                jfragtype->natom,&jfragtype->types[0],&tempcoords2[0]);
            break;
        case GB_MODE_PER_ATOM:
            en=ffield->exact_interaction_energy_gb_per_atom(pbc,halfboxsize,boxsize,gb_params,ifragtype->natom,&ifragtype->types[0],&tempcoords1[0],&born_radii1[0],
                jfragtype->natom,&jfragtype->types[0],&tempcoords2[0],&born_radii2[0]);
            break;
        case GB_MODE_PER_ATOM_CONVERTED:
        case GB_MODE_PER_FRAGMENT:
            a1a2=born_radii[ifrag]*born_radii[jfrag];
            //printf("a1a2=%.5f\n",a1a2);
            //fflush(stdout);
            if (a1a2>0) {
                en=ffield->exact_interaction_energy_gb_per_frag(pbc,halfboxsize,boxsize,gb_params,ifragtype->natom,&ifragtype->types[0],&tempcoords1[0],
                   jfragtype->natom,&jfragtype->types[0],&tempcoords2[0],a1a2);
            } else { //we do not need to do gb, both fragments are uncharged, this prevents floating point exceptions possibly due to overoptimization
                en=ffield->exact_interaction_energy(pbc,halfboxsize,boxsize,gb_params->epsp,ifragtype->natom,&ifragtype->types[0],&tempcoords1[0],
                   jfragtype->natom,&jfragtype->types[0],&tempcoords2[0]);
            }
            break;
    }
#ifdef TIMERS
    switch_timer(TIMER_INT_OTHER);
#endif
    return en;
}

double topology::covalent_table_energy(double * coords, bool * moved, covalent_table * * covalent_tables)
{
    double encov;
    int seg, ifrag, jfrag, itype, jtype;
    covalent_table * covtbl;
#ifdef TIMERS
    switch_timer(TIMER_COV_TABLES);
#endif
    encov=0.0;
    for (seg=0; seg<nseg; seg++) {
        ifrag=first_main_chain_frag[seg];
        jfrag=frags[ifrag].main_chain_next;
        while (jfrag>0) {
            itype=frags[ifrag].type;
            jtype=frags[jfrag].type;
            //we rely on short-circuit evaluaton here
            if (use_covalent_table(itype,jtype) && ((moved==NULL) || (moved[ifrag]^moved[jfrag]))) {
                covtbl=covalent_tables[itype*nfragtypes+jtype];
                if (covtbl==NULL) {
                    printf("Covalent table for fragment types %s and %s not loaded.\n",fragtypes[itype]->fragname, fragtypes[jtype]->fragname);
                    die();
                }
                encov+=covtbl->table_interaction_energy(coords,frags[ifrag].atc,frags[ifrag].atn,frags[ifrag].atca,
                    frags[jfrag].atc,frags[jfrag].atn,frags[jfrag].atca);
            }
            ifrag=jfrag;
            jfrag=frags[jfrag].main_chain_next;
        }
    }
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
    return encov;
}
//This subroutine does tabulated interaction energies and exact interaction energies for fragments not having a 1,4-relationship
//If the two fragments are too close for tabulation,all the vdw/electrostatic terms are handled through the bottom section of "internal_energy" or "moved_internal_energy"
double simulation::interaction_energy(int ifrag, int jfrag, double * center, double * orient,double * coords, double * born_radii)
{
    //Switches between exact_interaction_energy and table_interaction_energy.
    int reffrag,otherfrag,tableindex,itype,jtype,iatomstart,jatomstart,inatom,jnatom,k;
    bool close;
    table * tbl;
    double en,enexact,entable,rij[3],r2,rdiff;    //Master spherical cutoff ensures same for both table-based and exact simulations.
#ifdef TIMERS
    switch_timer(TIMER_CHECK_CUTOFF);
#endif
    for (k=0; k<3; k++) {
        rij[k]=center[3*jfrag+k]-center[3*ifrag+k];
        if (pbc) {
            if (rij[k]>halfboxsize) rij[k]-=boxsize;
            if (rij[k]<-halfboxsize) rij[k]+=boxsize;
        }
    }
    r2=rij[0]*rij[0]+rij[1]*rij[1]+rij[2]*rij[2];
    if (r2>cutoff2) {
#ifdef TIMERS
        switch_timer(TIMER_INT_OTHER);
#endif
        return 0.0;
    }
    if (r2<minr2) minr2=r2;
    enevalcount++;
    enexact=0.0;
    entable=0.0;
    //tables_lambda: 0 = fully exact, 1 = fully tabulated.
    if ((tables_lambda>0.0) || enwrite) {
#ifdef TIMERS
        switch_timer(TIMER_INT_PREP);
#endif
        //figure out which fragment will be the reference.
        //ensure interactions are always taken in a consistent order
        itype=top->frags[ifrag].type;
        jtype=top->frags[jfrag].type;
        if (itype==jtype) { //it really doesn't matter which is the reference fragment
            tbl=tables[itype*top->nfragtypes+jtype];
            if (tbl==NULL) {
                printf("Table for fragment types %s and %s not loaded.\n",top->fragtypes[itype]->fragname, top->fragtypes[jtype]->fragname);
                die();
            }
            if (ifrag<jfrag) {
                reffrag=ifrag;
                otherfrag=jfrag;
            } else if (ifrag>jfrag) {
                reffrag=jfrag;
                otherfrag=ifrag;
                for (k=0; k<3; k++) rij[k]=-rij[k];
            } else { //shouldn't happen
                printf("Interaction energy error.\n");
                die();
            }
        } else {
            //We need to consult the table to see which fragment should be the reference.
            //This guarantees consistency in the choice of which fragment is the reference between table generation and table use.
            tbl=tables[itype*top->nfragtypes+jtype];
            if (tbl==NULL) tbl=tables[jtype*top->nfragtypes+itype];
            if (tbl==NULL) {
                printf("Table for fragment types %s and %s not loaded.\n",top->fragtypes[itype]->fragname, top->fragtypes[jtype]->fragname);
                die();
            }
            if (itype==tbl->reffragtype) {
                reffrag=ifrag;
                otherfrag=jfrag;
            } else if (jtype==tbl->reffragtype) {
                reffrag=jfrag;
                otherfrag=ifrag;
                for (k=0; k<3; k++) rij[k]=-rij[k];
            } else { //shouldn't happen
                printf("Interaction energy error.\n");
                die();
            }
        }
        //double table_interaction_energy(double r, double * rij, int reftype, double * qref, int othertype, double * qother)
        entable=tbl->table_interaction_energy(enwrite,interp,r2,&rij[0],&orient[4*reffrag],&orient[4*otherfrag],&rdiff);
        if (en_by_table!=NULL) {
            if (itype<=jtype) tableindex=itype*top->nfragtypes+jtype; else tableindex=jtype*top->nfragtypes+itype;
            en_by_table[tableindex]+=entable;
        }
        entablecount++;
    }
    if ((tables_lambda<1.0) || enwrite || (entable<=INVALID_ENERGY)) {
        //double exact_interaction_energy(double eps, int natom1, int * types1, double * coords1, int natom2, int * types2, double * coords2)
        //if (ifrag<jfrag) {
        //the appropriate set of born_radii has been selected by moved_energy or total_energy
        enexact=top->exact_interaction_energy(ffield,pbc,halfboxsize,boxsize,gb_mode,&gb_params,ifrag,jfrag,coords,born_radii);
        enexactcount++;
    }
    /*if (enwrite) {
        fprintf(energy_output,"%d %d %.4f %.10f %.10f %.10f\n",ifrag,jfrag,rdiff,enexact,entable,en);
        if (fabs(enexact-entable)>5.0) write_pair_pdb(pairs_output,ifrag,jfrag,coords);
    }*/
    en=(1.0-tables_lambda)*enexact+tables_lambda*entable;
#ifdef DEBUG_NON_TABULATED
        printf("Interaction energy: %d %d %.5f\n",ifrag,jfrag,en);
#endif
    //if (exact) return enexact; else return entable;
    return en;
}

//Change in energy upon moving fragment.
//this assumes only one fragment has been moved.
//check nonbond list.  We only need to count interactions between moved and non-moved fragments,
//since each move is a rigid transformation of the moved fragments.
//Note on the GB term.  It is sum(ij) tau*q_i*q_j/fGB over all atoms i and j. It is divided into three parts.
//The case where both atoms i and j belong to the same fragment is covered with dG_k in calculate_born_radii and is placed in energies[EN_GB_EXACT].
//If atoms i and j are in different fragments, but
/*void simulation::moved_energy(bool * moved, bool * movedatoms, bool new_born_radii, double * center, double * orient, double * coords, double * energies, double * etot)
{
    double energy;
    int i,j,itype,jfrag,ifrag,n;
    double * born_radii;
    //eold=0.0;
    //for (j=0; j<nfrags; j++)
    //    if (j!=imoved){
    for (i=0; i<EN_TERMS; i++) energies[i]=0.0;
    if (en_by_table!=NULL) for (itype=0; itype<top->nfragtypes*top->nfragtypes; itype++) en_by_table[itype]=0.0;
    if (gb_mode!=GB_MODE_NONE) {
        if (new_born_radii) calculate_born_radii(center,coords); //this results in doing this twice but is simple for now.
        if (gb_mode==GB_MODE_PER_ATOM) n=top->natom; else n=top->nfrag;
        for (i=0; i<n; i++) {
            energies[EN_GB_SELF]+=dg_self[i];
            energies[EN_SASA]+=sasa[i];
        }
    }
    energies[EN_SASA]*=sasa_coef;    //Figure out the correct born radii to use.
    switch(gb_mode) {
        case GB_MODE_NONE: born_radii=NULL; break;
        case GB_MODE_PER_ATOM:
            born_radii=per_atom_born_radii; break;
        case GB_MODE_PER_ATOM_CONVERTED:
        case GB_MODE_PER_FRAGMENT:
            born_radii=per_fragment_born_radii; break;
    }
    ffield->moved_non_tabulated_energy(top->natom,top->atoms,movedatoms,non_tab_list.size(),&non_tab_list[0],gb_mode,gb_params,coords,born_radii,energies);
    if (use_cov_tables) energies[EN_COVALENT_TABLE]=top->covalent_table_energy(coords,moved,covalent_tables);
#ifdef TIMERS
    switch_timer(TIMER_INT_OTHER);
#endif
    if (use_nb_list) {
        for (ifrag=0; ifrag<top->nfrag; ifrag++) if (moved[ifrag]) {
            for (j=0; j<frag_nblist->nb_list_count[ifrag]; j++) {
                jfrag=frag_nblist->nonbond_list[frag_nblist->nb_list_per_frag*ifrag+j];
                if (top->closefragments[ifrag*top->nfrag+jfrag]) continue; //the fragments are too close for tabulation
                if (!moved[jfrag]) energies[EN_INTERACTION]+=interaction_energy(ifrag,jfrag,center,orient,coords,born_radii);
            }
        }
    } else { //no nb list
        for (ifrag=0; ifrag<top->nfrag; ifrag++)  {
            for (jfrag=(ifrag+1); jfrag<top->nfrag; jfrag++) {
                //if (moved[ifrag]==moved[jfrag]) continue;  //The interaction needs to be computed only if exactly one of iatom/jatom  has been moved (since all moves are rigid on subsets of atoms)
                if (top->closefragments[ifrag*top->nfrag+jfrag]) continue; //the fragments are too close for tabulation
                energies[EN_INTERACTION]+=interaction_energy(ifrag,jfrag,center,orient,coords,born_radii);
            }
        }
    }
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
    *etot=0.0;
    for (i=0; i<EN_TERMS; i++) *etot += energies[i];
}*/

//this only works on "new" coordinates.
void simulation::total_energy(double * per_atom_born_radii, double * per_fragment_born_radii, double * dg_self, double * sasa, double * center, double * orient, double * coords, double * energies, double * etot)
{
    int i,j,itype,jtype,n;
    double energy, inte, evdwint, eelecint;
    double * born_radii;
    for (i=0; i<EN_TERMS; i++) energies[i]=0.0;
    if (en_by_table!=NULL) for (itype=0; itype<top->nfragtypes*top->nfragtypes; itype++) en_by_table[itype]=0.0;
    if (gb_mode!=GB_MODE_NONE) {
        //if (new_born_radii) calculate_born_radii(center,coords); //this results in doing this twice but is simple for now.
        if ((gb_mode==GB_MODE_PER_ATOM) || (gb_mode==GB_MODE_PER_ATOM_CONVERTED)) n=top->natom; else n=top->nfrag;
        for (i=0; i<n; i++) {
            energies[EN_GB_SELF]+=dg_self[i];
            energies[EN_SASA]+=sasa[i];
        }
    }
    energies[EN_SASA]*=sasa_coef;
    switch(gb_mode) {
        case GB_MODE_NONE: born_radii=NULL; break;
        case GB_MODE_PER_ATOM:
	case GB_MODE_PER_ATOM_CONVERTED:
            born_radii=per_atom_born_radii; break;
        //case GB_MODE_PER_ATOM_CONVERTED:
        case GB_MODE_PER_FRAGMENT:
            born_radii=per_fragment_born_radii; break;
    }
    ffield->non_tabulated_energy(top->natom,top->atoms,non_tab_list.size(),&non_tab_list[0],gb_mode,gb_params,coords,born_radii,energies);
#ifdef TIMERS
    switch_timer(TIMER_INT_OTHER);
#endif
    if (use_cov_tables) energies[EN_COVALENT_TABLE]=top->covalent_table_energy(coords,NULL,covalent_tables);
    /*energies[EN_VDW_INT]=evdw_internal;
    energies[EN_ELEC_INT]=eelec_internal;*/
#ifdef TIMERS
    switch_timer(TIMER_INT_OTHER);
#endif
    for (i=0; i<top->nfrag; i++)
        for (j=(i+1); j<top->nfrag; j++)
            if (!top->closefragments[i*top->nfrag+j]) energies[EN_INTERACTION]+=interaction_energy(i,j,center,orient,coords,born_radii);
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
    *etot=0.0;
    for (i=0; i<EN_TERMS; i++) *etot += energies[i];
        /*{
            inte=interaction_energy(i,j,newcenter,neworient,newcoords);
            fprintf(energy_output,"total_energy: i, j, inte, energy = %d %d %.4f %.4f\n",i,j,inte,energy);
            energy+=inte;
        }*/
#ifdef TIMERS
    switch_timer(TIMER_OTHER);
#endif
}


void simulation::print_energies_by_table(void)
{
    int itype,jtype,idx;
    printf("Total interaction energies by table.\n");
    for (itype=0; itype<top->nfragtypes; itype++)
    for (jtype=itype; jtype<top->nfragtypes; jtype++) {
        idx=itype*top->nfragtypes+jtype;
        if (tables[idx]!=NULL) printf("Types %s and %s: %.6f\n",top->fragtypes[itype]->fragname,top->fragtypes[jtype]->fragname,en_by_table[idx]);
    }
}


//Recenter fragments within the box.
/*void simulation::recenter(void)
{
    int ifrag,k,flag;
    for (ifrag=0; ifrag<nfrags; ifrag++)
    {
        flag=FALSE;
        for (k=0; k<3; k++) {
            if (oldcenter[3*ifrag+k]>halfboxsize) {
                oldcenter[3*ifrag+k]-=boxsize;
                flag=TRUE;
            }
            if (oldcenter[3*ifrag+k]<-halfboxsize) {
                oldcenter[3*ifrag+k]+=boxsize;
                flag=TRUE;
            }
            if (flag) {
                update_coords(ifrag,oldcenter,oldorient,oldcoords);
                copy_frag(ifrag,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
            }
        }
    }
}*/
//Fake loop in which monte carlo moves are made and then accepted.  Test monte carlo moves/i/o routines.
void simulation::fakeloop(void)
{
    long int istep;
    bool * moved;
    bool * movedatoms;
    int movetype,ifrag;
    double fresh_energy, fresh_energies[EN_TERMS];
    xyzoutput = fopen(xyzfname,"wb"); //DCD file
    if (xyzoutput==NULL) {
        printf("Could not open DCD trajectory file %s\n",xyzfname);
        die();
    }
    write_dcd_header(xyzoutput);
    quatoutput = fopen(quatfname,"w");
    if (quatoutput==NULL) {
        printf("Could not open center/orientation file %s\n",quatfname);
        die();
    }
    printf("Will write trajectory output to file %s.\n",xyzfname);
    printf("Will write centers/quaternions to file %s.\n",quatfname);

    moved = (bool *) checkalloc(top->nfrag,sizeof(bool));
    movedatoms = (bool *) checkalloc(top->natom,sizeof(bool));
    for (istep=1; istep<=nmcstep; istep++) {
        mcmove(&movetype,moved,movedatoms,newcenter,neworient,newcoords);
	//total_energy(newcenter,neworient,newcoords,fresh_energies,&fresh_energy);
        //print_energies(FALSE,"Energy:",istep,fresh_energies,fresh_energy);
        printf("step %ld:  %s move.\n",istep,mc_move_names[movetype]);
        //accept all moves
        for (ifrag=0; ifrag<top->nfrag; ifrag++) if (moved[ifrag]) top->copy_frag(ifrag,newcenter,neworient,newcoords,oldcenter,oldorient,oldcoords);
        if ((istep%nsave)==0) {
            printf("step %ld\n",istep);
            write_dcd_frame(xyzoutput,newcoords);
            write_frame_quat(quatoutput,istep,newcenter,neworient);
        }
    }
    fclose(xyzoutput);
    fclose(quatoutput);
    free(moved);
    free(movedatoms);
}

void simulation::print_energies(FILE * output, int hdr, const char * title, long int istep, double * energies, double etot)
{
    const char * hdrfmt = "%15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s %15s\n";
    const char * enfmt = "%15s %15ld %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f %15.6f\n";
    if (hdr) {
        fprintf(output,hdrfmt,title,"Step","Total","Interaction","Cov. intxn","Bonds","Angles","Dihedrals","Impropers","Non-tab VDW","Non-tab Elec","GB Self","Non-tab GB","SASA");
    }
    fprintf(output,enfmt,title,istep,etot,energies[EN_INTERACTION],energies[EN_COVALENT_TABLE],energies[EN_BOND],energies[EN_ANGLE],energies[EN_DIHEDRAL],energies[EN_IMPROPER],
        energies[EN_VDW_EXACT],energies[EN_ELEC_EXACT],energies[EN_GB_SELF],energies[EN_GB_EXACT],energies[EN_SASA]);
    fflush(output);
}


/*The Monte Carlo Loop.*/
//void mcloop(FILE * xyzoutput,FILE * quatoutput)
void simulation::mcloop(void)
{
    long int istep,nacc[NUM_MOVES+1],natt[NUM_MOVES+1];
    int movetype,i;
    bool new_nb_list;
    double cum_energy,fresh_energy,enew,eold,de,de2,r,p,accrate,deviation;
    double cputime,elapsedtime;
    double oldenergies[EN_TERMS],newenergies[EN_TERMS],cum_energies[EN_TERMS],fresh_energies[EN_TERMS];
    bool * moved;
    bool * movedatoms;
    double * backbone;
    double q[4],c[3],rmsd;
    clock_t starttime;
    time_t start,end;
    char buffer[255];
    xyzoutput = fopen(xyzfname,"wb"); //DCD file
    if (xyzoutput==NULL) {
        printf("Could not open DCD trajectory file %s\n",xyzfname);
        die();
    }
    write_dcd_header(xyzoutput);
    quatoutput = fopen(quatfname,"w");
    if (quatoutput==NULL) {
        printf("Could not open center/orientation file %s\n",quatfname);
        die();
    }
    printf("Will write trajectory output to file %s.\n",xyzfname);
    printf("Will write centers/quaternions to file %s.\n",quatfname);
    moved = (bool *) checkalloc(top->nfrag,sizeof(bool));
    movedatoms = (bool *) checkalloc(top->natom,sizeof(bool));
    backbone = (double *) checkalloc(top->natom,sizeof(double));
    for (i=0; i<top->natom; i++) if (strstr("N CA C",top->atoms[i].name)!=NULL) backbone[i]=1.0; else backbone[i]=0.0;
    //nb_atom_list=new std::vector<atom_nb_entry>;
    if (use_nb_list) frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,top->nfrag,newcenter);
    //top->create_overlap_list(nb_list_per_frag,nb_list_count,nonbond_list,oldcoords,&overlap_list);
    //top->create_nb_atom_exact_list(exact,nb_list_per_frag,nb_list_count,nonbond_list,&nb_atom_list);
    /*fresh_energy=total_energy();
    cum_energy=fresh_energy;*/
    for (i=1; i<=NUM_MOVES; i++) {
        nacc[i]=0;
        natt[i]=0;
    }
    enexactcount=0;
    entablecount=0;
    enevalcount=0;
    minr2=1000.0;
    //printf("Step %ld: Cumulative energy = %.4f   Fresh energy = %.4f\n",0,cum_energy,fresh_energy);
    time(&start);
#ifdef TIMERS
    init_timers();
#endif
    printf("Starting Monte Carlo at %s\n",ctime(&start));
    starttime=clock();
    calculate_born_radii(oldcenter,oldcoords,old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa);
    total_energy(old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa,oldcenter,oldorient,oldcoords,fresh_energies,&fresh_energy);
    print_energies(stdout,TRUE,"Energy:",0,fresh_energies,fresh_energy);
    for (i=0; i<EN_TERMS; i++) cum_energies[i]=fresh_energies[i];
    for (i=0; i<EN_TERMS; i++) oldenergies[i]=fresh_energies[i];
    cum_energy=fresh_energy;
    //die();

    for (istep=1; istep<=nmcstep; istep++) {
#ifdef TIMERS
         switch_timer(TIMER_MC_MOVE);
#endif
         mcmove(&movetype,moved,movedatoms,newcenter,neworient,newcoords);
         //calculate_born_radii(center,coords,
         //eold=moved_energy(movedfrag,oldcenter,oldorient,oldcoords);
         natt[movetype]++;
         //eold=0.0;
         //for (i=0; i<EN_TERMS; i++) eold+=oldenergies[i];
         //moved_energy(moved,movedatoms,oldcenter,oldorient,oldcoords,oldenergies,&eold);
         //calculate_born_radii(oldcenter,oldcoords,old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa);
         total_energy(old_per_atom_born_radii,old_per_fragment_born_radii,old_dg_self,old_sasa,oldcenter,oldorient,oldcoords,oldenergies,&eold);
         if (use_nb_list) {
#ifdef TIMERS
             switch_timer(TIMER_NB_LIST);
#endif

             new_nb_list=frag_nblist->check_nb_list(pbc,halfboxsize,boxsize,cutoff,top->nfrag,moved,newcenter);
            if (new_nb_list) {
                //printf("New nonbond list at step %ld\n",istep);
                //if (pbc) recenter();
                frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,top->nfrag,newcenter);
                //top->create_nb_atom_exact_list(exact,nb_list_per_frag,nb_list_count,nonbond_list,&nb_atom_list);
            }
#ifdef TIMERS
             switch_timer(TIMER_OTHER);
#endif
         }
         //moved_energy(moved,movedatoms,newcenter,neworient,newcoords,newenergies,&enew);
         calculate_born_radii(newcenter,newcoords,new_per_atom_born_radii,new_per_fragment_born_radii,new_dg_self,new_sasa);
         total_energy(new_per_atom_born_radii,new_per_fragment_born_radii,new_dg_self,new_sasa,newcenter,neworient,newcoords,newenergies,&enew);
         de=enew-eold;
         //fresh_energy=total_energy(newcenter,neworient);
         //de2=fresh_energy-cum_energy;
         //printf("mcloop: step, de, de2 = %d %.4f %.4f\n",istep,de,de2);
         //This guards against floating point errors.
         if (de<0) p=1.0; else p=exp(-beta*de);
         r=genrand_real3();
         if (r<p) {
             //ACCEPTED
             cum_energy+=de;
             for (i=0; i<EN_TERMS; i++) cum_energies[i]+=newenergies[i]-oldenergies[i];
             nacc[movetype]++;
             for (i=0; i<top->nfrag; i++) if (moved[i]) top->copy_frag(i,newcenter,neworient,newcoords,oldcenter,oldorient,oldcoords);
             //for (i=0; i<EN_TERMS; i++) oldenergies[i]=newenergies[i];
             if ((gb_mode==GB_MODE_PER_ATOM) || (gb_mode==GB_MODE_PER_ATOM_CONVERTED)) for (i=0; i<top->natom; i++) {
                old_per_atom_born_radii[i]=new_per_atom_born_radii[i];
                old_dg_self[i]=new_dg_self[i];
                old_sasa[i]=new_sasa[i];
             }
             if ((gb_mode==GB_MODE_PER_FRAGMENT) || (gb_mode==GB_MODE_PER_ATOM_CONVERTED)) for (i=0; i<top->nfrag; i++) {
                old_per_fragment_born_radii[i]=new_per_fragment_born_radii[i];
             }
             if (gb_mode==GB_MODE_PER_FRAGMENT) for (i=0; i<top->nfrag; i++) {
                old_dg_self[i]=new_dg_self[i];
                old_sasa[i]=new_sasa[i];
             }
             //if (de<=-0.5*DUMMY_ENERGY) cum_energy=total_energy(); //This guards against numerical errors related to "declashing."
         } else {
             //REJECTED
             for (i=0; i<top->nfrag; i++) if (moved[i]) top->copy_frag(i,oldcenter,oldorient,oldcoords,newcenter,neworient,newcoords);
             //If we made a new nonbond list, we now need to get rid of it.  This is duplicative, but avoids the need to manage two nonbond lists.
             if (use_nb_list && new_nb_list) {
#ifdef TIMERS
                 switch_timer(TIMER_NB_LIST);
#endif
                 frag_nblist->create_nonbond_list(pbc,halfboxsize,boxsize,top->nfrag,newcenter);
#ifdef TIMERS
                 switch_timer(TIMER_OTHER);
#endif
                 //top->create_nb_atom_exact_list(exact,nb_list_per_frag,nb_list_count,nonbond_list,&nb_atom_list);
             }
             /*if ((gb_mode==GB_MODE_PER_ATOM) || (gb_mode==GB_MODE_PER_ATOM_CONVERTED)) for (i=0; i<top->natom; i++) new_per_atom_born_radii[i]=old_per_atom_born_radii[i];
             if (gb_mode==GB_MODE_PER_ATOM) for (i=0; i<top->natom; i++) {
                new_dg_self[i]=old_dg_self[i];
                new_sasa[i]=old_dg_self[i];
             }
             if ((gb_mode==GB_MODE_PER_FRAGMENT) || (gb_mode==GB_MODE_PER_ATOM_CONVERTED)) for (i=0; i<top->nfrag; i++) {
                new_per_fragment_born_radii[i]=old_per_fragment_born_radii[i];
                new_dg_self[i]=old_dg_self[i];
                new_sasa[i]=old_dg_self[i];
             }*/
         }
         //check_nb_list(istep);
         //Now, "new" and "old" should be the same again, and on the MC trajectory.
         //Handle saving and printing.
         if ((istep%nsave)==0) {
             write_dcd_frame(xyzoutput,newcoords);
             write_frame_quat(quatoutput,istep,newcenter,neworient);
         }
         if ((istep%nprint)==0) {
             fflush(xyzoutput);
             fflush(quatoutput);
             //fresh_energy=total_energy();
             //print_energies(FALSE,"Cum: ",istep,cum_energies,cum_energy);
             calculate_born_radii(newcenter,newcoords,new_per_atom_born_radii,new_per_fragment_born_radii,new_dg_self,new_sasa);
             total_energy(new_per_atom_born_radii,new_per_fragment_born_radii,new_dg_self,new_sasa,newcenter,neworient,newcoords,fresh_energies,&fresh_energy);
             print_energies(stdout,FALSE,"Energy:",istep,fresh_energies,fresh_energy);
             deviation=cum_energy-fresh_energy;
             printf("Energy deviation: %.6f kcal/mol\n",deviation);
             if (fabs(deviation)>ENERGY_TOLERANCE) {
                 print_energies(stdout,FALSE,"Cum. energy:",istep,cum_energies,cum_energy);
                 printf("Too much deviation between cumulative and fresh energies.\n");
                 if (strlen(restartfname)>0) write_restart(nprevstep+istep,restartfname);
                 die();
             }
             //printf("Step %ld: Cumulative energy = %.4f   Fresh energy = %.4f Deviation %.4f\n",istep,cum_energy,fresh_energy,cum_energy-fresh_energy);
             for (i=1; i<=NUM_MOVES; i++) {
                if (natt[i]>0) accrate=((double) nacc[i]/(double) natt[i])*100.0; else accrate=0.0;
                printf("%s moves   Attempted %ld   Accepted %ld   Acceptance rate=%.2f%%\n",mc_move_names[i],natt[i],nacc[i],accrate);
             }
             //if (nacc[MOVE_BACKRUB]>0) die();
             //printf("Minimum r: %.3f\n",sqrt(minr2));
             printf("Exact/table/overall energy evaluations: %ld %ld %ld\n",enexactcount,entablecount,enevalcount);
             printf("Fraction exact/table energy evaluations: %.2f%% %.2f%%\n",100*((double)enexactcount/(double)enevalcount),100*((double)entablecount/(double)enevalcount));
             if (initcoords!=NULL) {
                 rmsd_fit(top->natom,backbone,initcoords,newcoords,c,q,&rmsd);
                 printf("Backbone RMSD: %.3f A\n",rmsd);
             }
         }
#ifndef NO_CRC
         if (use_std_tables && (ncheck>0) && ((istep%ncheck)==0)) {
             printf("Verifying table checksums.\n");
             //check table checksums to see if memory corrupted
             for (i=0; i<top->nfragtypes*top->nfragtypes; i++)
                if ((tables[i]!=NULL) && !tables[i]->verify_checksum()) {
                    printf("Bad checksum in table.\n");
                    die();
                }
         }
#endif
#ifdef EXCHANGE
         if ((istep%exchfreq)==0) {
		exchange(istep/exchfreq,cum_energies,&cum_energy,newcenter,neworient,newcoords);
		//print_energies(stdout,FALSE,"Cum: ",istep,cum_energies,cum_energy);
                //total_energy(newcenter,neworient,newcoords,fresh_energies,&fresh_energy);
		//The cum_energy and cum_energies are now up to date with the newly exchanged configuration (or not, as the case may be)
 		//therefore, we need to keep the "fresh" energies up to date as well.
                for (i=0; i<EN_TERMS; i++) fresh_energies[i]=cum_energies[i];
                fresh_energy=cum_energy;
             	//print_energies(stdout,FALSE,"Energy:",istep,fresh_energies,fresh_energy);
	}

#endif
    } //end of main loop "for(istep=1; istep<=nmcstep; istep++)"
#if defined(PARALLEL) || defined(EXCHANGE)
    //printf("node %d before final barrier\n",mynod);
    //fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
    //printf("node %d after final barrier\n",mynod);
    //fflush(stdout);
#endif
#ifdef TIMERS
    switch_timer(TIMER_NONE);
#endif
    cputime=((double)(clock()-starttime))/CLOCKS_PER_SEC;
    time(&end);
    elapsedtime=difftime(end,start);
    printf("Total CPU time %.2f seconds.\n",cputime);
    printf("%.2f MC steps per second.\n",nmcstep/cputime);
    printf("Total elapsed time %.2f seconds.\n",elapsedtime);
#ifdef TIMERS
    print_timers();
#endif
    if (strlen(restartfname)>0) write_restart(nprevstep+nmcstep,restartfname);
    fclose(xyzoutput);
    fclose(quatoutput);
    if (enwrite) {
        fclose(energy_output);
        fclose(pairs_output);
    }
#ifdef EXCHANGE
    //printf("node %d closing replica log\n",mynod);
    //fflush(stdout);
    MPI_File_close(&rexlog);
    //printf("node %d closed replica log.\n",mynod);
    //fflush(stdout);
#endif
    free(moved);
    free(movedatoms);
    free(backbone);
}

