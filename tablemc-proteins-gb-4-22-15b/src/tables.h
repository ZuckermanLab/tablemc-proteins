#ifndef TABLES_H_INCLUDED
#define TABLES_H_INCLUDED

#include "time.h"
#include "ffield.h"
//#include "mc.h"
#include "fragments.h"
#include "topology.h"

#define RAD_TO_DEG                  57.2957795//converts radians to degrees
#define DEG_TO_RAD                  1.74532925e-2//converts degrees to radians
#define SPHERICAL_DIFFUSION_CONST   2.000
#define SO3_DIFFUSION_CONST         2.000
#define KBOLTZ                      1.987E-3
#define VERSION                     3.0
#define OP_BOLTZMANN_AVERAGED       0x1
#define OP_DIST_DEP_DIELECTRIC      0x4
#define OP_SINGLE_PRECISION         0x8
#define OP_DECIMATED                0x10

#ifndef __cplusplus
#error Not C++
#endif
/*This structure defines the format of the header of table files.  If it is changed old table files will no longer be readable.*/

typedef float energy_t;

struct table_header {
    size_t headersize; /*size of the header*/
    long int options; /*for bit flags*/
    double version;/*Version that created this table*/
    time_t created; /*Date/time of creation*/
    unsigned int crc32_data; /*CRC-32 of table entries.*/
    double minr,maxr,dr,dsph,deuler; /*Maximum r, min. resolution in r, spherical theta/phi, euler angles*/
    //double min_clash, max_clash;
    double rfactor, logrfactor;
    long int nr,nsphtheta, nsphphi, ntheta, nphipsi; /*Number of r's, spherical theta,phi, Euler theta, phipsi*/
    long int ntrans,norient; /*Number of translational/orientational points*/
    long long int totalpoints, clashpoints,invalidpoints; /*Total number of points */
    int decr, decsph, deceuler; /*Decimation factors*/
    /*The following information is not used directly, but serves for recordkeeping*/
    /*int nboltzmannaverage;
    long int seed;
    double boltzmannbeta;*/
    double eps; /*Dielectric constant used*/
    //int nsmoothr, nsmoothsph, nsmootheuler; /*Smoothing range for radius, angular, orientational*/
    double beta,sph_alpha,euler_alpha;
    double en_smooth_limit,en_invalid_margin;
    int n_invalid_dev;
    char ref_frag_fname[255]; /*Fragment 1 file name*/
    char other_frag_fname[255]; /*Fragment 2 file name*/
    char forcefield[255]; /*Force field file name*/
};

class table {
private:
    table_header hdr;
    energy_t * energy;
#ifndef NO_MMAP_TABLES
    void * map; //points to the memory map of the table file
    size_t mapsize;
#endif
    double * radial_table;
    double one_minr, one_logrfactor, one_dsph, one_deuler;
    //double * clash_table;

    void read_table_header_info(const char * fname);
    void print_header_info(void);
    inline long long int calculate_index(const int ir, const int isphtheta, const int isphphi, const int iphi, const int itheta, const int ipsi)
    {
        long int itrans,iorient;
        //itrans=hdr.nx*hdr.nx*ix+hdr.nx*iy+iz;
        itrans=((hdr.nsphphi*isphtheta)+isphphi)*hdr.nr+ir;
        iorient=(hdr.ntheta*iphi+itheta)*hdr.nphipsi+ipsi;
        //printf("*** %ld %ld\n",itrans,iorient);
        return hdr.ntrans*((long long int) iorient)+itrans;
    }
    //inline long long int clash_table_index(long int isphtheta, long int isphphi, long int iphi, long int itheta, long int ipsi);
    void alloc_read_table(FILE * f, const char * fname);

    /*double get_clash_radius(double sphtheta, double sphphi, double phi, double theta, double psi);
    double interpolate_clash_radius(double sphtheta, double sphphi, double phi, double theta, double psi);*/
    //double get_energy(int enwrite, double r, double sphtheta, double sphphi, double phi, double theta, double psi);
    //void fill_clash_table(void);
    void fill_table(forcefield * ffield, fragmenttype * frag1, fragmenttype * frag2);
    void fake_fill_table();
    void boltzmann_average_trans(void);
    void boltzmann_average_orient(void);
    double get_pair_energy(forcefield * ffield, fragmenttype * frag1, fragmenttype * frag2, double * coords1, double * coords2, double r, double sphtheta, double sphphi, double phi, double theta, double psi);
    int find_neighbor(long int * indices, long int * disp, long int * neighbor);
    double assoc_free_energy(double beta);

    inline double get_energy(int enwrite, double r, double sphtheta, double sphphi, double phi, double theta, double psi)
    {
        long int ir,isphtheta,isphphi,iphi,itheta,ipsi;
        long int irmin,irmax,irmid;
        long long int index;
        int i,boundscheck;
        double u;
        energy_t * p;
        //More efficient to multiply than divide.
        //if (r<hdr.min_clash) return DUMMY_ENERGY;
        isphtheta=(long int) (sphtheta*one_dsph);
        if (isphtheta>=hdr.nsphtheta) isphtheta=hdr.nsphtheta-1;
        isphphi=((long int) (sphphi*one_dsph + 0.5));
        if (isphphi>=hdr.nsphphi) isphphi-=hdr.nsphphi;
        itheta=(long int) (theta*one_deuler);
        if (itheta>=hdr.ntheta) itheta=hdr.ntheta-1;
        iphi=((long int) (phi*one_deuler+0.5));
        if (iphi>=hdr.nphipsi) iphi-=hdr.nphipsi;
        ipsi=((long int) (psi*one_deuler+0.5));
        if (ipsi>=hdr.nphipsi) ipsi-=hdr.nphipsi;
        //printf("indices: %d %d %d %d %d %d\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
        //printf(energy_output,"indices: %d %d %d %d %d %d\n",ir,isphtheta,isphphi,iphi,itheta,ipsi);
        u=log(r*one_minr)*one_logrfactor;
        ir=(long int) (u+0.5);
        //if (ir>=hdr.nr) ir=hdr.nr-1;
        index=calculate_index(ir,isphtheta,isphphi,iphi,itheta,ipsi);
        //if (table[index]<=INVALID_ENERGY) printf("hit invalid energy cell.\n");
#ifdef TIMERS       
        switch_timer(TIMER_INT_LOOKUP);
#endif
        return energy[index];
        //return 0.0;
    }

public:
    //These are indices into topoloy's array of fragment types, which indicate which fragment type was the reference.
    int reffragtype, otherfragtype;
    table(const char * fname, int newtable);
    table(const char * fmt, const char * fragfmt, const char * fragname1, const char * fragname2, topology * top);
    table(table * orig_table, int decr, int decsph, int deceuler);
    ~table();
    void get_random_cell(double * r, double * sphtheta, double * sphphi, double * phi, double * theta, double * psi);
    double table_interaction_energy(int enwrite, int interp, double r2, double * rij,  double * qref,  double * qother, double * rdiff);
    void generate_table(const char * control_file);
    void write_table(const char * fname);
    void write_dx(char * fname, double phi, double theta, double psi, double xmax, double dx);
    void write_dx_exact(char * fname, double phi, double theta, double psi, double xmax, double dx);
    void write_dx_orient(char * fname, double r, double sphtheta, double sphphi, double dx);
    double getsize(void);
    bool verify_checksum(void);
    void volume_test(int ntest);
};


#endif // TABLES_H_INCLUDED
