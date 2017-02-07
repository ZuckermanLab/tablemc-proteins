#ifndef FRAGMENTS_H_INCLUDED
#define FRAGMENTS_H_INCLUDED

#include <vector>
#include "ffield.h"
#include "solvation.h"

//#include "tables.h"

#define MAX_ATOMS_PER_FRAGMENT 20
//#define MAX_FRAGMENT_TYPES 10
//#define MAX_FRAGMENTS 1000
//#define MAX_ATOMS MAX_FRAGMENTS*MAX_ATOMS_PER_FRAGMENT
//#define MAX_MAINCHAIN_FRAGMENTS 4
//#define MAX_SIDECHAIN_FRAGMENTS 4
#define MAX_FRAGMENTS_PER_RESIDUE 4
#define MAX_SEGMENTS              10
#define MAX_ATOMS_PER_RESIDUE     100
#define MAX_BONDS_PER_RESIDUE     MAX_ATOMS_PER_RESIDUE*4
#define MAX_FRAGMENT_NAME         32
#define UNKNOWN_COORD             9999.0
#define RT_NOT_ROTATABLE          0
#define RT_BACKBONE               1
#define RT_SIDECHAIN              2

class table;
class covalent_table;

//This actually represents a fragment type.
class fragmenttype {
    friend class topology; //Allows topology to access private members.
private:
    double refgeom[3*MAX_ATOMS_PER_FRAGMENT]; /*Reference geometry*/

    double totalmass;
    double mass[MAX_ATOMS_PER_FRAGMENT]; //Masses of atoms, for RMSD fitting.  Also used by topology::assemble_fragments for overall RMSD fit, so must be public.
public:
    int natom; /*Number of atoms in the fragment*/
//    int in_use; //At least one of these fragments is used in the system. For determining if tables need to be loaded.
    int n_used;
    double qtot;
    double dipolemag; /*total charge, magnitude of dipole moment*/
    double dipole[3]; /*dipole moment in reference geometry*/
    bool is_charged; /* does this fragment have a net charge or dipole? */
    //solvation_radius is the radius such that the GK energy of the charge/dipole equals the PB solvation energy of the fragment
    double solvation_radius, volume_radius, volume;
    bool has_covalent_tables; /*True if either a peptide or n-methyl-acetamide*/
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
    //void internal_interaction_energy(forcefield * ffield, double eps, int rdie, double * evdwint, double * eelecint);
};


struct fragment {
    int type; //Fragment type index.
    int is_side_chain; //Whether or not the fragment is part of the side chain.
    //Each of the following residues is -1 if
    int main_chain_prev; //Previous residue along main chain.
    int main_chain_next; //Next residue along main chain.
    int side_chain_prev; //Previous residue along side chain (towards main chain)
    int side_chain_next; //Next residue along side chain (away from main chain)
    int atoms[MAX_ATOMS_PER_FRAGMENT]; //Atoms in the fragment.
    int atn, atca, atc; //N, CA, and C atoms, for covalent table lookup
};
//Information on fragments within residues.  The arrays in this structure are indexed by the index of an atom within fragment.
struct resfraginfo {
    int fragtype;
    int is_side_chain; //is this fragment part of the side chain or the backbone?
    //Names of atoms within residue that are part of the fragment.
    char atomnames[MAX_ATOMS_PER_FRAGMENT][6];
    //int overlap; //Does this fragment overlap with a fragment in the previous residue?
    //int atomsprevious[MAX_ATOMS_PER_FRAGMENT]; //Atoms that overlap in previous residue.
    int offset[MAX_ATOMS_PER_FRAGMENT]; //Residue offset.  If the residue overlaps, this designates the atoms that are actually in the previous residue.
};

/*struct covtablelookup {
    int ifrag,jfrag;
    int itype,jtype;
    int */
//Definition of a residue.
struct residuedef {
    char name[4];
    int natom; //Number of atoms per residue.
    char atomnames[MAX_ATOMS_PER_RESIDUE][6];
    int nbond; //Number of bonds.
    char iname[MAX_BONDS_PER_RESIDUE][6];//Atoms pair for each bond.
    char jname[MAX_BONDS_PER_RESIDUE][6];
    int rottype[MAX_BONDS_PER_RESIDUE]; //Is this bond rotatable? If so, MC moves are generated based on rotation about this bond.
    int joffset[MAX_BONDS_PER_RESIDUE]; //Residue offset on "j" atom.
    int nfrag; //Number of fragments intersecting this residue.
    resfraginfo frags[MAX_FRAGMENTS_PER_RESIDUE];
    int branchatom; //Name of the branch atom.
};

struct reslookup {
    int restype; //Index of the definition.
    int whichseg; //Which segment is it a part of?
    int branchatom; //Alpha carbon, or other branching atom.
    int atomstart;
    int atomend;
    //int nscrot; //Number of backbone and sidechain rotatable bonds.
    int nbbrot;
    /*int iscrot[MAX_BONDS_PER_RESIDUE];
    int jscrot[MAX_BONDS_PER_RESIDUE];*/
    int ibbrot[MAX_BONDS_PER_RESIDUE];
    int jbbrot[MAX_BONDS_PER_RESIDUE]; //Is the bond part of the side chain?
    //int peptidebond; //Number of the actual peptide bond fragment.
    //int sidechainstart; //Starting fragment of the side chain.
    //int sidechainend; //Ending fragment of the side chain.
};

struct topology {
    int nfragtypes;
    fragmenttype * * fragtypes; //[MAX_FRAGMENT_TYPES];
    /*number of actual fragments, total number of atoms*/
    int nfrag,natom,nseg,nres,nscrot;
    /*which type of fragments each fragment is*/
    fragment * frags; //[MAX_FRAGMENTS];
    //virtual_site * virtual_sites;
    /*starting atom number for each fragment*/
    //int * fragstart;  //[MAX_FRAGMENTS];
    /*Atom information*/
    ATOMS * atoms;  //[MAX_ATOMS];
    int * iscrot;
    int * jscrot;
    int nresdef;
    residuedef * resdef;
    reslookup * resinfo;
    //int * sequence; //Sequence of residue index numbers;
    //Start and end of each segment, segment index for each residue.
    int * segstart;
    int * segend;
    int * first_main_chain_frag;
    //Chain code for each segment. To match PDB file.
    char * chaincodes;
    bool * closefragments;
    double qsystem;
    //int * whichseg;
    topology(const char * commandfile, forcefield * ffield, const char * fragfmt);
    ~topology();
    int fragtypebyname(const char * name);
    int fragtypebyfile(const char * fname);
    int resdefbyname(const char * name);

    void read_residue_definition(FILE * f, residuedef * def);
    int find_atom(int actual_res, const char * aname);
    int find_atom(char chain, int res, const char * aname);
    //int is_bonded(int ifrag, int jfrag);
    //void addseg(int nsegres, char * seq);
    void print_detailed_info(void);
    void print_summary_info(void);
    void insert_residue(const char * res);
    void link_fragments(void);
    void create_angle_dihedral_lists(bool using_cov_tables);
    void create_improper_dihedral_lists(bool using_cov_tables, forcefield * ffield);
    void create_non_tab_list(bool using_cov_tables,std::vector<atom_nb_entry> * atom_nb_list);
    bool use_covalent_table(int itype, int jtype);
    bool term_in_covalent_tables(int iatom, int jatom);
    bool term_in_covalent_tables(int iatom, int jatom, int katom);
    bool term_in_covalent_tables(int iatom, int jatom, int katom, int latom);
    //void create_nb_atom_exact_list(int exact, int nb_list_per_frag, int * nb_list_count, int * nonbond_list, std::vector<atom_nb_entry> * atom_nb_list);
    void add_segment(char chain, const char * sequence);
    void assemble_fragments(double * orig_coords, double * center, double * orient, double * new_coords);
    void fit_all_fragments(double * orig_coords, double * center, double * orient, double * new_coords, double * rmsds);
    void load_tables(const char * fmt, const char * fragfmt, table * * tables);
    void load_covalent_tables(const char * covtablefmt, covalent_table * * cov_tables);
    void update_coords(int ifrag, double * center, double * orient, double * coords);
    void copy_frag(int ifrag, double * center1, double * orient1, double * coords1, double * center2, double * orient2, double * coords2);
    double exact_interaction_energy(forcefield * ffield, int pbc, double halfboxsize, double boxsize, double eps,  int rdie, int frag1, int frag2, double * coords);
    double covalent_table_energy(double * coords, bool * moved, covalent_table * * covalent_tables);
    double dipolar_gb_energy(gb_param_info * gb_params, double * per_frag_born_radii, double * center, double * orient);
    void calculate_born_radii(gb_param_info * params, double * center, double * per_frag_born_radii);
    void read_gb_params(char * line, gb_param_info * params, double cutoff);
    void print_gb_info(gb_param_info * params);
    //void total_internal_interaction_energy(forcefield * ffield, double eps, int rdie, double * evdw_internal, double * eelec_internal);
    //I/O routines.
    void read_pdb_stream(FILE * input, double * coords);
    void read_pdb_file(char * fname, double * coords);
    void write_pdb_stream(FILE * output, char * fmt, double * coords, forcefield * ffield);
    void write_pdb_file(char * fname, char * fmt, double * coords, forcefield * ffield);
    void write_pdb_file2(forcefield * ffield, char * fname, double * coords);
    void write_psf_file(char * fname, forcefield * ffield);
    void read_frame_quat(FILE * input, long int * istep, double * center, double * orient);
    void write_frame_quat(FILE * output, long int istep, double * center, double * orient);
    //Monte Carlo move generation.
    void create_backbone_move(bool * moved, int * atom1, int * atom2);
    void create_sidechain_move(bool * moved, int * atom1, int * atom2);
    void create_backrub_move(bool * moved, int * atom1, int * atom2);
    void get_moved_atoms(bool * movedfrag, bool * movedatoms);
    //SASA stuff.
    void add_to_overlap_list(int ifrag, int jfrag, double * coords, std::vector<atom_nb_entry> * overlap_list);
    void create_overlap_list(int nb_list_per_frag, int * nb_list_count, int * nonbond_list, double * coords, std::vector<atom_nb_entry> * overlap_list);
};

void write_pair_pdb2(FILE * output, fragmenttype * frag1, fragmenttype * frag2, double * coords1, double * coords2);
#endif // FRAGMENTS_H_INCLUDED
