#ifndef SOLVATION_H_INCLUDED
#define SOLVATION_H_INCLUDED

//possible modes for GB calculation
#define GB_MODE_NONE                0 //vacuum calculation
#define GB_MODE_PER_ATOM            1 //calculate and use per atom born radii
#define GB_MODE_PER_ATOM_CONVERTED  2 //calculate per atom born radii, convert them to per fragment radii
#define GB_MODE_PER_FRAGMENT        3 //calculate per fragment born radii
//#define GB_KAPPA                    4 //exponenet in Still equation
#define WATER_RADIUS          1.4 //in A
//convergence criterion for delta G in fragment-based born radius solution
#define GB_RADIUS_CONV        0.001
//maximum number of Newton iterations
#define GB_MAX_ITER           1000
//These paraemters are for the Still formula and SASA term.
//kappa -- in the exponent of the still formula (desginated c_f in Schneiders and Ponder)
//c0, c1 -- coefficients in kirkwood expansion (eq. 19 of the paper)
//alpha, beta, gamma -- for the tanh formula
//radii_offset and radii_scale_factor -- used in calculating Born radii (radii_scale_factor corresponds to S_kk' in the GB-HCT paper
struct gb_param_info {
    double epsp, epsw, kappa, c0, c1, sasa_coef, rgbmax, rgbmaxpsmax2, alpha, beta, gamma, radii_offset, radii_scale_factor;
};

//These parameters are for the determination of a Born radius using FACTS
/*struct facts_param_info {
    //int type; //could be an atom class (see ffield.h) or a fragment type
    double volume, b1, b2, a0, a1, a2, a3, rsphere2, d1, d2, c0, c1, c2, c3;
};*/



void setup_gb_spline(void);
double splintGB(double x);
void splintGB_diff(const double x, double * y, double * dy);
#endif // SOLVATION_H_INCLUDED
