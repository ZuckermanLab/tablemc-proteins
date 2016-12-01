#ifndef COVALENT_TABLES_H_INCLUDED
#define COVALENT_TABLES_H_INCLUDED

#include "time.h"
#include "ffield.h"
//#include "mc.h"
#include "fragments.h"

#define RAD_TO_DEG                  57.2957795//converts radians to degrees
#define DEG_TO_RAD                  1.74532925e-2//converts degrees to radians
#define KBOLTZ                      1.987E-3
#define OP_BOLTZMANN_AVERAGED       0x1
#define OP_DIST_DEP_DIELECTRIC      0x4
#define OP_SINGLE_PRECISION         0x8
#define MAX_SMOOTH_REGIONS          10

typedef float cov_energy_t;

/*A rectangular region of the Ramachandran plane.*/
struct rama_region {
    char name[8];
    double philo, phihi, psilo, psihi;
    long int iphilo, iphihi, ipsilo, ipsihi;
    double smooth_scale;
};
//The degrees of freedom are:
//phi = C-N-CA-C dihedral
//psi = N-CA-C-N dihedral
//alpha = N-CA-C angle.
//omega = CA-C-N-CA dihedral (need to keep track of the sign change)
//thetac = CA-C-N angle.
//r = CA-C bond length.

struct covalent_table_header {
    size_t headersize; /*size of the header*/
    long int options; /*for bit flags*/
    time_t created; /*Date/time of creation*/
    long long int totalpoints;
    long int nphi, npsi,nalpha; // nr, nalpha, nthetac, nomega;
    double dphi, dpsi,  minalpha, maxalpha, dalpha; //minr, maxr, dr, minthetac, maxthetac, dthetac, minomega, domega;
    double rcac, thetacacn, omega; //Fixed values of other parameters.
    double eps; /*dielectric constant used for electrostatic interactions here*/
    int n_smooth_region;
    double smooth_beta;
    rama_region smooth_regions[MAX_SMOOTH_REGIONS];
    //double e1,e2,alpha1,alpha2; //parameters for eq. 3 in the Sinko paper, e1 and e2 are relative to the minimum for a given alpha
    char defs[255]; /*definitions*/
    char fragfmt[255];
    char sequence[32]; /*sequence*/
    char fragtype1[MAX_FRAGMENT_NAME]; /*Fragment 1 file name*/
    char fragtype2[MAX_FRAGMENT_NAME]; /*Fragment 2 file name*/
    char forcefield[255]; /*Force field file name*/
};

class covalent_table {
private:
    covalent_table_header hdr;
    cov_energy_t * energy;
    void read_table_header_info(const char * fname);
    void read_region_info(FILE * f, rama_region * rgn);
    void print_header_info(void);
    //ordering from highest to lowest order indices: r, thetac, alpha, phi, psi
    /*inline long long int calculate_index(const int iphi, const int ipsi, const int iomega, const int ialpha, const int ithetac, const int ir)
    {
        return ((((ir*hdr.nthetac+ithetac)*hdr.nalpha+ialpha)*hdr.nomega+iomega)*hdr.nphi+iphi)*hdr.npsi+ipsi;
    }*/
    inline long long int calculate_index(const int iphi, const int ipsi, const int ialpha)
    {
        return (ialpha*hdr.nphi+iphi)*hdr.npsi+ipsi;
    }
    void alloc_read_table(FILE * f, const char * fname);
    double get_energy(double phi, double psi, double alpha);
    void fill_table(topology * top, forcefield * ffield,fragmenttype * fragtype1,fragmenttype * fragtype2, int nb_atom_list_size, atom_nb_entry * nb_atom_list);
    double partition_function(double beta,rama_region rgn);
    //void smooth(void);
    void boltzmann_average(rama_region rgn);
#ifdef DEBUG
    void print_energies(int hdr, const char * title, long int istep, double * energies, double etot);
#endif
public:
    covalent_table(const char * fname, int newtable);
    ~covalent_table();
    double table_interaction_energy(double * coords, int ic1, int in1, int ica1, int ic2, int in2, int ica2);
    void generate_table(const char * control_file);
    void write_table(const char * fname);
    void write_energies(double alpha, const char * fname);
    double getsize(void);
};

#endif // COVALENT_TABLES_H_INCLUDED
