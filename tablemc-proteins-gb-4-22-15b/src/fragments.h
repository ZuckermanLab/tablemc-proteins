#ifndef FRAGMENTS_H_INCLUDED
#define FRAGMENTS_H_INCLUDED

#include "ffield.h"
#include <vector>
//#include "tables.h"

#define MAX_ATOMS_PER_FRAGMENT 20
#define MAX_FRAGMENT_NAME      32
//#define MAX_FRAGMENT_TYPES 10
//#define MAX_FRAGMENTS 1000
//#define MAX_ATOMS MAX_FRAGMENTS*MAX_ATOMS_PER_FRAGMENT
//#define MAX_MAINCHAIN_FRAGMENTS 4
//#define MAX_SIDECHAIN_FRAGMENTS 4


class table;
class covalent_table;

//This actually represents a fragment type.
class fragmenttype {
    friend class topology; //Allows topology to access private members.
private:
    double refgeom[3*MAX_ATOMS_PER_FRAGMENT]; /*Reference geometry*/
    double dipolemag; /*total charge, magnitude of dipole moment*/
    double dipole[3]; /*dipole moment in reference geometry*/
    double totalmass;
    double mass[MAX_ATOMS_PER_FRAGMENT]; //Masses of atoms, for RMSD fitting.  Also used by topology::assemble_fragments for overall RMSD fit, so must be public.
    //distance matrix. used for finding born radii
    double d2matrix[MAX_ATOMS_PER_FRAGMENT][MAX_ATOMS_PER_FRAGMENT];
public:
    int natom; /*Number of atoms in the fragment*/
//    int in_use; //At least one of these fragments is used in the system. For determining if tables need to be loaded.
    int n_used;
    double qtot;
    bool has_covalent_tables; /*True if either a peptide or n-methyl-acetamide*/
    bool is_charged;
    char fragname[MAX_FRAGMENT_NAME];
    char names[MAX_ATOMS_PER_FRAGMENT][4]; /*names of all the atoms */
    int types[MAX_ATOMS_PER_FRAGMENT]; /*Types of atoms.  Some table stuff needs this.*/
    fragmenttype(const char * name, const char * fname, forcefield * ffield);
    void get_coords(double * center, double * orient, double * coords);
    int findatombyname(const char * name);
    void fit_fragment(double * coords, double * center, double * orient,double * rmsd);
    void fit_fragment(double * coords, double * center, double * orient, double * weights, double * rmsd);
    double get_bond_length(int iatom, int jatom);
    double get_angle(int iatom, int jatom, int katom);
    double get_dihedral(int iatom, int jatom, int katom, int latom);
    double get_average_born_radius(forcefield * ffield, double dgx);
    double get_average_born_radius(forcefield * ffield, double * per_atom_born_radii);
    //void internal_interaction_energy(forcefield * ffield, double eps, int rdie, double * evdwint, double * eelecint);
};




void write_pair_pdb2(FILE * output, fragmenttype * frag1, fragmenttype * frag2, double * coords1, double * coords2);
#endif // FRAGMENTS_H_INCLUDED
