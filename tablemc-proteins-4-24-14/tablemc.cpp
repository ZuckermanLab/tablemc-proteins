#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "string.h"
#include "mc.h"
#include "tables.h"
#include "covalent_tables.h"
#include "rotations.h"
#include "fragments.h"
#include "util.h"
#include "ffield.h"

#ifdef __unix__
#define _GNU_SOURCE   // gives us feenableexcept on older gcc's
#define __USE_GNU     // gives us feenableexcept on newer gcc's
#include <fenv.h>
#else
#include <float.h>
//#ifndef _EM_OVERFLOW
//#define _EM_OVERFLOW EM_OVERFLOW
//#endif
#endif

//Command line syntax: tablemc input_file  OR tablemc generate input_file table_file]
int main(int argc, char * argv[])
{
    int imoved,i,mynod,numnod;
    //double en,en2,de;
    char fname[255],buffer[255];
    time_t now;

    table * newtable;
    table * newdectable;
    covalent_table * newcovtable;
    simulation * sim;
#ifdef __unix__
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
#else
//    _controlfp(0, _EM_ZERODIVIDE | _EM_INVALID | _EM_OVERFLOW);
#endif
#ifndef NO_TRIG_TABLES
    fill_trig_tables();
#endif
    if (argc<2) {
        printf("Syntax: tablemc input_file\n");
        die();
    }
#if defined(PARALLEL) || defined(EXCHANGE)
    parallel_init(&mynod,&numnod,argv[2]);
#endif
    time(&now);
    strncpy(buffer,ctime(&now),sizeof(buffer));
    printf("Tabulated Monte Carlo starting at %s\n",buffer);
    printf("Executable name: %s\n",argv[0]);
    fflush(stdout);
#if defined(PARALLEL) || defined(EXCHANGE)
//Run simulations only in parallel mode (as of right now)
        sim=new simulation("run",argv[1],mynod,numnod);
        fflush(stdout);
        sim->mcloop();
        delete sim;
#else
    if (strcasecmp(argv[1],"generate")==0) {
        if (argc<4) {
           printf("Syntax: tablemc generate input_file table_file\n");
           die();
        }
        newtable=new table(argv[2],TRUE);
        newtable->write_table(argv[3]);
        delete newtable;
    } else if (strcasecmp(argv[1],"decimate")==0) {
        newtable=new table(argv[2],FALSE);
        newdectable=new table(newtable,atoi(argv[3]),atoi(argv[4]),atoi(argv[5]));
        newdectable->write_table(argv[6]);
        delete newtable;
        delete newdectable;
    } else if (strcasecmp(argv[1],"run")==0) {
        sim=new simulation(argv[1],argv[2]);
        //sim->fakeloop();
        sim->mcloop();
        delete sim;
    } else if (strcasecmp(argv[1],"energy")==0) {
        sim=new simulation(argv[1],argv[2]); //Constructor sets up energies only, does not read MC-related stuff.
        sim->do_energies(argv[3],argv[4],argv[5]);
        delete sim;
    } else if (strcasecmp(argv[1],"dx-orient")==0) {
        newtable=new table(argv[2],FALSE);
        newtable->write_dx_orient(argv[7],atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]));
        delete newtable;
    } else if (strcasecmp(argv[1],"dx-exact")==0) {
        newtable=new table(argv[2],FALSE);
        newtable->write_dx_exact(argv[8],atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]),atof(argv[7]));
        delete newtable;
    } else if (strcasecmp(argv[1],"dx")==0) {
        newtable=new table(argv[2],FALSE);
        newtable->write_dx(argv[8],atof(argv[3]),atof(argv[4]),atof(argv[5]),atof(argv[6]),atof(argv[7]));
        delete newtable;
    /*} else if (strcasecmp(argv[1],"volume")==0) {
        newtable=new table(argv[2],FALSE);
        newtable->volume_test(atoi(argv[3]));
        delete newtable;*/
    } else if (strcasecmp(argv[1],"gencov")==0) {
        newcovtable=new covalent_table(argv[2],TRUE);
        newcovtable->write_table(argv[3]);
        delete newcovtable;
    } else if (strcasecmp(argv[1],"plotcov")==0) {
        newcovtable=new covalent_table(argv[2],FALSE);
        newcovtable->write_energies(atof(argv[3]),argv[4]);
        delete newcovtable;
    } else {
        //Default to running a simulation.
        sim=new simulation(argv[1],argv[2]);
        sim->mcloop();
        delete sim;
    }
#endif
    time(&now);
    strncpy(buffer,ctime(&now),sizeof(buffer));
    printf("Finished at %s\n",buffer);
#if defined(PARALLEL) || defined(EXCHANGE)
    parallel_finish();
#endif
    return 0;
}
