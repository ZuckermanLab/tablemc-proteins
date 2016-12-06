#include <cstdio>
#include "fragments.h"
#include "rotations.h"
#include "mc.h"
#include "ffield.h"
#include "mt.h"
#include "util.h"
//Rotate a selection of fragments about an axis through two atoms.  Does not support PBC.
void simulation::rotate_fragments_by_axis(const bool * moved, const int atom1, const int atom2, const double angle, double * center, double * orient)
{
  int ifrag,k;
  double axis[3],point[3],axism,quat[4],newquat[4];
  double rotmatrix[3][3];
  double disp[3],newdisp[3];
  //Compute the axis
  axism=0.0;
  for (k=0; k<3; k++) {
      axis[k]=oldcoords[3*atom2+k]-oldcoords[3*atom1+k];
      axism+=axis[k]*axis[k];
      point[k]=oldcoords[3*atom1+k];
  }
  axism=sqrt(axism);
  for (k=0; k<3; k++) axis[k]/=axism;
  axisangle_to_quat(angle,axis,quat);
  quat_to_matrix(quat,&rotmatrix[0][0]);
  for (ifrag=0; ifrag<top->nfrag; ifrag++) if (moved[ifrag]) {
      //Rotate center about axis.
      for (k=0; k<3; k++) disp[k]=center[3*ifrag+k]-point[k];
      matmul(&rotmatrix[0][0],disp,newdisp);
      for (k=0; k<3; k++) center[3*ifrag+k]=point[k]+newdisp[k];
      //Rotate orientation.
      multiply_quat(&orient[4*ifrag],quat,newquat);
      normalize_quat(newquat);
      for (k=0; k<4; k++) orient[4*ifrag+k]=newquat[k];
  }
}


void topology::create_backbone_move(bool * moved, int * atom1, int * atom2)
{
    int res,bond,iatom,jatom,ok,ifrag,jfrag,mainfrag,sidefrag;
    //Choose residue at random.  Make sure the residue has at least one rotatable bond.
      do {
        res=(int) (genrand_real3()*nres);
        ok=resinfo[res].nbbrot>0;
    } while (!ok);
    //Choose a rotatable bond at random. At least one of the atoms has to belong to side chain.
    bond=(int) (genrand_real3()*resinfo[res].nbbrot);
    //ensure that the "j" atom is on the side chain
    iatom=resinfo[res].ibbrot[bond];
    jatom=resinfo[res].jbbrot[bond];
    ifrag=atoms[iatom].fragment;
    jfrag=atoms[jatom].fragment;
    //Either ifrag and jfrag are same, or they aren't.  If they are, that fragment is not rotated (phi bond), but any side chain fragments attached to it are.
    //If ifrag is not the same as jfrag,  jfrag is the first to be rotated.
    mainfrag=jfrag;
    moved[mainfrag]=(ifrag!=jfrag);
    while (mainfrag>0) {
        sidefrag=frags[mainfrag].side_chain_next;
        while (sidefrag>0) {
            moved[sidefrag]=true;
            sidefrag=frags[sidefrag].side_chain_next;
        }
        mainfrag=frags[mainfrag].main_chain_next;
        if (mainfrag>0) moved[mainfrag]=true;
    }
    *atom1 = iatom;
    *atom2 = jatom;
}

void topology::create_sidechain_move(bool * moved, int * atom1, int * atom2)
{
    int res,bond,iatom,jatom,ok,ifrag,jfrag,frag,swap;
    bool special_case,is_proline;
    //Choose residue at random.  Make sure the residue has at least one rotatable bond.
    /*do {
        res=(int) (genrand_real3()*nres);
        is_proline=(strcasecmp(resdef[resinfo[res].restype].name,"PRO")==0);
        ok=(resinfo[res].nscrot>0)||(is_proline);
    } while (!ok);
    //Deal with proline specially.
    if (is_proline) {
        iatom=find_atom(res,"CB");
        jatom=find_atom(res,"CD");
    } else {*/
        //Choose a rotatable bond at random. At least one of the atoms has to belong to side chain.
    bond=(int) (genrand_real3()*nscrot);
        //ensure that the "j" atom is on the side chain
    iatom=iscrot[bond];
    jatom=jscrot[bond];
    //}
    if (!atoms[jatom].is_side_chain) {
        swap=iatom;
        iatom=jatom;
        jatom=swap;
    }
//#ifdef DEBUG
//    printf("Side chain axis involving atoms %d %d\n",iatom,jatom);
//#endif
    ifrag=atoms[iatom].fragment;
    jfrag=atoms[jatom].fragment;
    //Either ifrag and jfrag are same, or they aren't.  If they are, that fragment is to be rotated; if they aren't, jfrag is the first to be rotated.
    frag=jfrag;
    //special case: LYS CG-CD bond does not rotate the fragment, even though both atoms are within it.
    special_case=((strcasecmp(atoms[iatom].resName,"LYS")==0) && (strcasecmp(atoms[iatom].name,"CG")==0));
    if (special_case) frag=frags[frag].side_chain_next;
    while (frag>0) {
        moved[frag]=true;
        frag=frags[frag].side_chain_next;
    }
    *atom1 = iatom;
    *atom2 = jatom;
}

void topology::create_backrub_move(bool * moved, int * atom1, int * atom2)
{
    int res1,seg,res2,bond,iatom,jatom,ok,ifrag,jfrag,mainfrag,sidefrag,swap;
    bool is_proline;
    //Choose two residues at random from the same segment.
    do {
        res1=(int) (genrand_real3()*(nres-2));
        seg=resinfo[res1].whichseg;
        //Backrub moves are not allowed to end on proline because this changes the N-CD bond with the new frag scheme.
        //Must choose a new res1 if this happens, otherwise we could get stuck if the only residues after res1 are prolines.
    	res2=res1+2+((int) (genrand_real3()*(segend[seg]-res1-2)));
        is_proline=(strcasecmp(resdef[resinfo[res2].restype].name,"PRO")==0);
    } while (is_proline);
    //printf("Backrub move between res. %d %d\n",res1,res2);
    iatom=resinfo[res1].branchatom;
    jatom=resinfo[res2].branchatom;
    ifrag=atoms[iatom].fragment;
    jfrag=atoms[jatom].fragment;
    //We assume ifrag and jfrag are main chain fragments.
    if (frags[ifrag].is_side_chain || frags[jfrag].is_side_chain) {
        printf("Branch atoms not on backbone.\n");
        die();
    }
    mainfrag=ifrag;
    //All side chains between ifrag and jfrag (but not those attached to ifrag and jfrag themselves)
    //and all main chain fragments between ifrag and jfrag (including jfrag but not ifrag)
    do {
        if (mainfrag!=ifrag) {
            sidefrag=frags[mainfrag].side_chain_next;
            while (sidefrag>0) {
                moved[sidefrag]=true;
                sidefrag=frags[sidefrag].side_chain_next;
            }
        }
        mainfrag=frags[mainfrag].main_chain_next;
        if (mainfrag>0) moved[mainfrag]=true;
    } while (mainfrag!=jfrag);
    *atom1 = iatom;
    *atom2 = jatom;
}

void simulation::mcmove(int * movetype, bool * moved, bool * movedatoms, double * center, double * orient, double * coords)
{
    int ifrag,k,move,atom1,atom2,iatom;
    double angle,r;
    //initialize moved array to all false
    for (ifrag=0; ifrag<top->nfrag; ifrag++) moved[ifrag]=false;
    //pick a move type
    r=genrand_real3();
    for (move=1; move<=NUM_MOVES; move++) {
        if (r<=cumprob[move]) break;
    }
    *movetype=move;
    //have the topology file create the move, choosing two atoms for the axis and
    switch (move) {
        case MOVE_BACKBONE:
            top->create_backbone_move(moved,&atom1,&atom2);
            break;
        case MOVE_SIDECHAIN:
            top->create_sidechain_move(moved,&atom1,&atom2);
            break;
        case MOVE_BACKRUB:
            top->create_backrub_move(moved,&atom1,&atom2);
            break;
        default:
            printf("Error in switch statement.\n");
            die();
    }
    angle=(2.0*genrand_real3()-1.0)*movesize[move];
    rotate_fragments_by_axis(moved,atom1,atom2,angle,center,orient);
    for (ifrag=0; ifrag<top->nfrag; ifrag++) if (moved[ifrag]) top->update_coords(ifrag,center,orient,coords);
    top->get_moved_atoms(moved,movedatoms);
}

void topology::get_moved_atoms(bool * movedfrag, bool * movedatoms)
{
    int ifrag, iatom;
    for (iatom=0; iatom<natom; iatom++) movedatoms[iatom]=movedfrag[atoms[iatom].fragment];
}

