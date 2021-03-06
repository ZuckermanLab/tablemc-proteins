#include <cstdio>
#include "ffield.h"
#include "util.h"
#include "rotations.h"

//Outside this limit on dot product, dihedrals will be 1
#define LIMIT    0.999999
forcefield::forcefield(char * fname)
{
  char buffer[3000];
  char name[300];
  //int devNumber;
  char *str=NULL;
  char delim[5]=":";
  int i,j,k,m;
  int ij;
  int index;
  int count;
  int countBonds,countAngles,countImprops,countDiheds;
  int frame;
  int numOfAtomsFrag;
  int startAtom,endAtom;
  float eng[3];
  double r;
  double num[12];
  double eps, sigma, sigma3, sigma6, sigma12;
  int mult;

  FILE* fHand;
    if((fHand = fopen(fname,"r"))==NULL){
    printf("FATAL ERROR: file %s is not found\n",fname);
    die();
  }

  countBonds=0;
  countAngles=0;
  countImprops=0;
  countDiheds=0;
  while(fgets(buffer,sizeof(buffer),fHand)){
    name[0] = '\0';
    sscanf(buffer,"%s",name);

    if(strcmp(name,"atom")==0){
      sscanf(buffer,"%s %d %d",name,&index,&i);
      atomTypeLookUp[index].classx = i;

      str=strtok(buffer,"\"");
      str=strtok(NULL,"\"");
      str=strtok(NULL,"\"");
      sscanf(str,"%d %lf %d",&i,&num[0],&j);
      atomTypeLookUp[index].atomicNum = i;
      atomTypeLookUp[index].mass = num[0];
      //printf("index: %d, atom class: %d, atomic number: %d\n",index,atomTypeLookUp[index].classx,atomTypeLookUp[index].atomicNum);
    }else if(strcmp(name,"vdw")==0){
      sscanf(buffer,"%s %d %lf %lf",name,&i,&num[0],&num[1]);
      vdwParams[i].sigma=num[0];
      vdwParams[i].epsilon=num[1];
      vdwParams[i].sigma14=num[0];//assume the 1-4 parameters are the same, unless specified differently
      vdwParams[i].epsilon14=num[1];
    }else if(strcmp(name,"vdw14")==0){
      sscanf(buffer,"%s %d %lf %lf",name,&i,&num[0],&num[1]);
      vdwParams[i].sigma14=num[0];
      vdwParams[i].epsilon14=num[1];
    }else if(strcmp(name,"bond")==0){
      sscanf(buffer,"%s %d %d %lf %lf",name,&i,&j,&num[0],&num[1]);
      bondParams[countBonds].atomClass[0] = i;
      bondParams[countBonds].atomClass[1] = j;
      bondParams[countBonds].K  = num[0];
      bondParams[countBonds].r0 = num[1];
      countBonds++;
    }else if(strcmp(name,"angle")==0){
      sscanf(buffer,"%s %d %d %d %lf %lf",name,&i,&j,&k,&num[0],&num[1]);
      angleParams[countAngles].atomClass[0] = i;
      angleParams[countAngles].atomClass[1] = j;
      angleParams[countAngles].atomClass[2] = k;
      angleParams[countAngles].K      = num[0];//it seems that tinker's K are in units of kcal/(mol*rad^2)
      angleParams[countAngles].theta0 = num[1]*DEG_TO_RAD;//convert degrees to radians
      countAngles++;
    }else if((strcmp(name,"improper")==0) || (strcmp(name,"imptors")==0)){
      sscanf(buffer,"%s %d %d %d %d %lf %lf",name,&i,&j,&k,&m,&num[0],&num[1]);
      impropParams[countImprops].atomClass[0] = i;
      impropParams[countImprops].atomClass[1] = j;
      impropParams[countImprops].atomClass[2] = k;
      impropParams[countImprops].atomClass[3] = m;
       //if((i==23 && j==1 && k==21 && m==22) || (i==1 && j==21 && k==23 && m==24)){
		impropParams[countImprops].K          = num[0];
		impropParams[countImprops].chi0     = num[1]*DEG_TO_RAD;
       /*}else{
	   impropParams[countImprops].V          = 0.5*num[0];
      }*/
      //it seems that tinker does not account for 1/2 in the improper dihedral potential energy function.
      //there may be some issues of symmetry. I will just keep the energy the same for now as in tinker but tinker may be wrong!
      //printf("%d %d %d %d  %f\n",impropParams[countImprops].atomClass[0],impropParams[countImprops].atomClass[1],impropParams[countImprops].atomClass[2],impropParams[countImprops].atomClass[3],impropParams[countImprops].V);
      countImprops++;
    }else if(strcmp(name,"torsion")==0){
      dihedParams[countDiheds].V[0] = 0;
      dihedParams[countDiheds].V[1] = 0;
      dihedParams[countDiheds].V[2] = 0;
      dihedParams[countDiheds].V[3] = 0;
      for (i=0; i<=3; i++) dihedParams[countDiheds].phase[i]=false;
      for(i=0;i<12;i++){
	num[i]=0;

      }

      sscanf(buffer,"%s %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
         name,&i,&j,&k,&m,&num[0],&num[1],&num[2],&num[3],&num[4],&num[5],&num[6],&num[7],&num[8],&num[9],&num[10],&num[11]);
      dihedParams[countDiheds].atomClass[0] = i;
      dihedParams[countDiheds].atomClass[1] = j;
      dihedParams[countDiheds].atomClass[2] = k;
      dihedParams[countDiheds].atomClass[3] = m;
      //this assumes there is only ONE dihedral per line, consistent with the CHARMM19 force field.
      //num[0]=magnitude, num[1]=phase, num[2]=multiplicity.
      mult=(int) num[2];
      if (mult>0) {
	dihedParams[countDiheds].phase[mult-1]=(num[1]>0.0); //assumes either 0 or 180
        dihedParams[countDiheds].V[mult-1]=num[0];
      }


      /*if(num[2]==1){
	dihedParams[countDiheds].V[0] = 0.5*num[0];//account for 1/2 in the dihedral potential energy function
      }else if(num[2]==2){
	dihedParams[countDiheds].V[1] = 0.5*num[0];
      }else if(num[2]==3){
        dihedParams[countDiheds].V[2] = 0.5*num[0];
      }else if(num[2]==4){
	dihedParams[countDiheds].V[3] = 0.5*num[0];
      }

      if(num[5]==2){
	dihedParams[countDiheds].V[1] = 0.5*num[3];
      }else if(num[5]==3){
	dihedParams[countDiheds].V[2] = 0.5*num[3];
	}else if(num[5]==4){
	dihedParams[countDiheds].V[3] = 0.5*num[3];
      }

      if(num[8]==3){
        dihedParams[countDiheds].V[2] = 0.5*num[6];
      }else if(num[8]==4){
	dihedParams[countDiheds].V[3] = 0.5*num[6];
      }

      if(num[11]==4){
        dihedParams[countDiheds].V[3] = 0.5*num[9];
      }
      //printf("%f %f %f\n",dihedParams[countDiheds].V[0],dihedParams[countDiheds].V[1],dihedParams[countDiheds].V[2]);*/

      countDiheds++;
    }else if(strcmp(name,"charge")==0){
      sscanf(buffer,"%s %d %lf",name,&i,&num[0]);
      chargeParams[i]=num[0];
    }
  }
  numOfBondParams=countBonds;
  printf("the number of bond params: %d\n",numOfBondParams);
  numOfAngleParams=countAngles;
  printf("the number of angle params: %d\n",numOfAngleParams);
  numOfImpropParams=countImprops;
  printf("the number of improper dihedral params: %d\n",numOfImpropParams);
  numOfDihedParams=countDiheds;
  printf("the number of dihedral params: %d\n",numOfDihedParams);

  fclose(fHand);
  //fabs() to protect against floating point exceptions on uninitialized data.
    for(i=1;i<MAX_NUM_OF_ATOM_CLASSES;i++){
        for(j=1;j<MAX_NUM_OF_ATOM_CLASSES;j++){
            eps = sqrt(fabs(vdwParams[i].epsilon*vdwParams[j].epsilon));
            //sigma = sqrt(fabs(vdwParams[i].sigma*vdwParams[j].sigma));
            //The CHARMM force field uses an additive rule for the vdW interactions.
            sigma = vdwParams[i].sigma+vdwParams[j].sigma;
            sigma3 = sigma*sigma*sigma;
            sigma6 = sigma3*sigma3;
            sigma12= sigma6*sigma6;
            //this change for the CHARMM vdw function
            vdwAFact[i][j] = eps*sigma12; //4*eps*sigma12;
            vdwBFact[i][j] = 2*eps*sigma6; //4*eps*sigma6;
            eps = sqrt(fabs(vdwParams[i].epsilon14*vdwParams[j].epsilon14));
            //sigma = sqrt(fabs(vdwParams[i].sigma*vdwParams[j].sigma));
            //The CHARMM force field uses an additive rule for the vdW interactions.
            sigma = vdwParams[i].sigma14+vdwParams[j].sigma14;
            sigma3 = sigma*sigma*sigma;
            sigma6 = sigma3*sigma3;
            sigma12= sigma6*sigma6;
            //this change for the CHARMM vdw function
            vdwAFact14[i][j] = eps*sigma12; //4*eps*sigma12;
            vdwBFact14[i][j] = 2*eps*sigma6;
        }
    }
}


double MLJ(double coefA,double coefB, double r6)
{
    if (r6==0) return DUMMY_ENERGY;
    double ir6 = 1./r6;
    double ir12 = ir6*ir6;
    double ene =(coefA*ir12 - coefB*ir6);
    //double ene = coefA*ir12;
    //if(ene > 100000.){return 100000.;}
    return ene;

}

double MCE(double chara,double charb, double rad)
{
    if (rad==0) return DUMMY_ENERGY;
    //double irad = 1./rad;
    double ene =(chara*charb/rad);
    /*if(ene > 100000.){return 100000.;}
    else if(ene <-100000.){return -100000.;}
    else{return ene;}*/
        return ene;

}

//changed this to use trig-free approach
//if c = cos(phi) then 1+cos(2phi) = 2.0*c^2, 1+cos(3phi) = 1-3*c+4c^3
//the CHARMM 19 force field doesn't have 4th order dihedrals, so dropped that term
//rewrote to accommodate phases more generally, but assumes they are either 0 or 180
double forcefield::MDE(double cosphi, int type){
	double cosphi1,cosphi2,cosphi3;
	if (dihedParams[type].phase[0]) cosphi1=-cosphi; else cosphi1=cosphi;
	if (dihedParams[type].phase[1]) cosphi2=2.0*(1-cosphi*cosphi); else cosphi2=2.0*cosphi*cosphi;
	if (dihedParams[type].phase[2]) cosphi3=1+cosphi*(3+4*cosphi*cosphi); else cosphi3=(1+cosphi*(-3+4*cosphi*cosphi));
	return dihedParams[type].V[0]*cosphi1 + dihedParams[type].V[1]*cosphi2+dihedParams[type].V[2]*cosphi3;
}

double angle(double a[], double b[]){
  double theta;

  theta = acos((a[0]*b[0] + a[1]*b[1] + a[2]*b[2])/sqrt((a[0]*a[0]+a[1]*a[1]+a[2]*a[2])*(b[0]*b[0]+b[1]*b[1]+b[2]*b[2])));

  return theta;
}

//justin altered this to return the cosine of the dihedral
double cos_dihedral(double a[],double b[],double c[]){
  double e[3];
  double f[3];
  double sign;
  double dot;
  double dihed;

  //vector product: e = a x b
  e[0] = a[1]*b[2] - a[2]*b[1];
  e[1] = a[2]*b[0] - a[0]*b[2];
  e[2] = a[0]*b[1] - a[1]*b[0];

  //vector product: f = c x b
  f[0] =-b[1]*c[2] + b[2]*c[1];
  f[1] =-b[2]*c[0] + b[0]*c[2];
  f[2] =-b[0]*c[1] + b[1]*c[0];

  //dot product: e dot f.
  dot = e[0]*f[0] + e[1]*f[1] + e[2]*f[2];
  dot = dot/sqrt((e[0]*e[0] + e[1]*e[1] + e[2]*e[2])*(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]));
  //asm("rsqrtss %1, %0":"=x"(r2):"x"(r2));

  //to avoid acos(dot)=Nan check if dot>1.0 and dot<-1.0
  if(dot>LIMIT){
    dot=1.0;
  }else if(dot<-LIMIT){
    dot=-1.0;
  }


  return dot;
  //return dihed;
}

double dihedral(double a[],double b[],double c[]){
  double e[3];
  double f[3];
  double g[3];
  double sign;
  double ee,ef,ff,gg,ec,bb,dot;
  double dihed;

  //vector product: e = a x b
  e[0] = a[1]*b[2] - a[2]*b[1];
  e[1] = a[2]*b[0] - a[0]*b[2];
  e[2] = a[0]*b[1] - a[1]*b[0];

  //vector product: f = c x b
  f[0] =-b[1]*c[2] + b[2]*c[1];
  f[1] =-b[2]*c[0] + b[0]*c[2];
  f[2] =-b[0]*c[1] + b[1]*c[0];

  //dot product: e dot f.
  ef = e[0]*f[0] + e[1]*f[1] + e[2]*f[2];
  ee = e[0]*e[0] + e[1]*e[1] + e[2]*e[2];
  ff = f[0]*f[0] + f[1]*f[1] + f[2]*f[2];
  //dot = dot/sqrt((e[0]*e[0] + e[1]*e[1] + e[2]*e[2])*(f[0]*f[0] + f[1]*f[1] + f[2]*f[2]));
  dot = ef/sqrt(ee*ff);
  //asm("rsqrtss %1, %0":"=x"(r2):"x"(r2));

  //to avoid acos(dot)=Nan check if dot>1.0 and dot<-1.0
  //justin: change this to use LIMIT (above) to avoid numerical problems with changing table entries near end of domain
  /*if(dot>LIMIT){
    dot=1.0;
    dihed=0.0;
  }else if(dot<-LIMIT){
    dot=-1.0;
    dihed=M_PI;
  } else dihed=acos(dot);*/
  //avoid using acos near end of domain
  /*if (fabs(dot)>0.99999) {
     dihed=sqrt(2*(1-dot));
     if (dot<0) dihed=-dihed;
     printf("dihedral: %.20f %.20f %.20f %.20f %.20f\n",dot,dihed,ef,ee,ff);
     printf("dihedral: e = %.20f %.20f %.20f\n",e[0],e[1],e[2]);
     printf("dihedral: f = %.20f %.20f %.20f\n",f[0],f[1],f[2]);
  } else dihed=acos(dot);*/
  //vector product: a = e x f. reuse variable a
  //e x f = (a x b) x (c x b) = -[(a x b) . c] b = (e.c) b
  //so (e.c) {b} = {e} |f| sin(phi) with opposite sign
  //gg = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]); //this is equal to |e| |f| sin(phi); ef=|e| |f| cos(phi)
  ec = e[0]*c[0] + e[1]*c[1] + e[2]*c[2];
  bb = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
  dihed = -atan2(ec*bb,ef); //this minus sign compensates for using c x b above
  //dihed = acos(dot); //this is quadrant-sensitive
  //if (ec>0) dihed=-dihed; // (a x b) . c is negative
  //check sign by dot product a dot b. reuse variable cx
  //sign = g[0]*b[0] + g[1]*b[1] + g[2]*b[2];

  /*if(sign<0){
    dihed = -dihed;
  }*/
  return dihed;
}

//Nonbonded energy for an individual atom pair.  The electrostatic energy returned by this subroutine does not have the Coulomb constant or dielectric constant for efficiency reasons.
void forcefield::nonbond_energy( int rdie, int type1,  int type2, int is14, double dx, double dy, double dz, double * evdw, double * eelec)
{
    int class1,class2;
    double q1,q2,r2,r6;

    r2=dx*dx+dy*dy+dz*dz;
    r6=r2*r2*r2;
    class1=atomTypeLookUp[type1].classx;
    class2=atomTypeLookUp[type2].classx;
    q1=chargeParams[type1];
    q2=chargeParams[type2];
    if (rdie) {
        *eelec=MCE(q1,q2,r2);
    } else {
        *eelec=MCE(q1,q2,sqrt(r2));
    }
    if (is14) {
        *evdw=MLJ(vdwAFact14[class1][class2],vdwBFact14[class1][class2],r6);
        *eelec*=BONDED14_SCALE;
    } else {
        *evdw=MLJ(vdwAFact[class1][class2],vdwBFact[class1][class2],r6);
    }
}

//This is only for pairs of fragments that do not share a 1-4 bond.
double forcefield::exact_interaction_energy(int pbc, double halfboxsize, double boxsize, double eps,  int rdie, int natom1, int * types1, double * coords1, int natom2, int * types2, double * coords2)
{
    int i,j,class1,class2;
    double dx,dy,dz,eelec, evdw,eelec1,evdw1;
    eelec=0.0;
    evdw=0.0;
    for (i=0; i<natom1; i++) {
        for (j=0; j<natom2; j++) {
            dx=coords1[3*i]-coords2[3*j];
            dy=coords1[3*i+1]-coords2[3*j+1];
            dz=coords1[3*i+2]-coords2[3*j+2];
            if (pbc) {
                if (dx>halfboxsize) dx-=boxsize;
                if (dx<-halfboxsize) dx+=boxsize;
                if (dy>halfboxsize) dy-=boxsize;
                if (dy<-halfboxsize) dy+=boxsize;
                if (dz>halfboxsize) dz-=boxsize;
                if (dz<-halfboxsize) dz+=boxsize;
            }
            nonbond_energy(rdie,types1[i],types2[j],FALSE,dx,dy,dz,&evdw1,&eelec1);
            eelec+=eelec1;
            evdw+=evdw1;
        }
    }
    //printf("%g %g\n",eelec*eps*COUL_CONST,evdw);
    return eelec*COUL_CONST/eps+evdw;
    //return evdw;
}


/*int term_needed(bool * moved, ATOMS * atoms, int a, int b, int c)
{
    int one_moved = (moved[atoms[a].fragment] || moved[atoms[b].fragment] || moved[atoms[c].fragment]);
    int one_not_moved = (!moved[atoms[a].fragment] || !moved[atoms[b].fragment] || !moved[atoms[c].fragment]);
    return one_moved && one_not_moved;
}

int term_needed(bool * moved, ATOMS * atoms, int a, int b, int c, int d)
{
    int one_moved = (moved[atoms[a].fragment] || moved[atoms[b].fragment] || moved[atoms[c].fragment] || moved[atoms[d].fragment]);
    int one_not_moved = (!moved[atoms[a].fragment] || !moved[atoms[b].fragment] || !moved[atoms[c].fragment] || !moved[atoms[d].fragment]);
    return one_moved && one_not_moved;
}*/
//Optimized bit twiddling.  A term is needed only if at least one atom is moved, and at least one atom isn't.
//The bit is set if the atom is moved.  Consequently the term is needed if the mask is not 000 and not 111 (or 1111 in the case of 4 atoms).
//Overloaded twice -- once for 3 atoms, once for 4.
inline bool term_needed(const bool * movedatoms, const int a, const int b, const int c)
{
    int mask;
    mask=0;
    if (movedatoms[a]) mask|=1;
    if (movedatoms[b]) mask|=2;
    if (movedatoms[c]) mask|=4;
    return (mask!=0) && (mask!=7);
}

inline bool term_needed(const bool * movedatoms, const int a, const int b, const int c, const int d)
{
    int mask;
    mask=0;
    if (movedatoms[a]) mask|=1;
    if (movedatoms[b]) mask|=2;
    if (movedatoms[c]) mask|=4;
    if (movedatoms[d]) mask|=8;
    return (mask!=0) && (mask!=15);
}

//Internal energy function for backbone, sidechain, and backrub moves.  Each is a rigid transformation of a subset of fragments.
//Therefore, we don't need bonds at all, and we only need terms that cross the boundary between moved/unmoved fragments.
//Also includes all non-tabulated interaction terms.
void forcefield::moved_non_tabulated_energy(double eps, int rdie, int numOfAtoms, ATOMS * atoms, bool * movedatoms, int nb_atom_list_size, atom_nb_entry * nb_atom_list,  double * coords, double * energies)
{
  int i,j,k;
  int ij;
  int ires;
  int iatom,jatom;
  int a,b,c,d;
  int type;
  int iclass,jclass;//M
  int ilibsc;
  int icsc;
  int count;
  int needed;
  double nbflag;
  double fi,temp;
  double ba[3],bc[3],cd[3],ad[3];
  double dx,dy,dz;
  double r2,r6;
  double K;
  double r,r0;
  double theta,theta0,chi0;
  double dihed;
  double en;
  double evdw,evdw2,eelec; /*for corrections*/
  bool is14;
  //zero out all energy terms excapt fragment interactions, which may have been previously calculated

  for (i=EN_BOND; i<EN_TERMS; i++) energies[i]=0.0;

  //calculate bond energies.
  //Since all moves implemented so far are rigid transformations of a selected subset of fragments, none of them changed
  //calculate bond angle energies
#ifdef TIMERS
   switch_timer(TIMER_NT_ANGLES);
#endif
  count=0;

  for(iatom=0;iatom<numOfAtoms;iatom++){
    for(j=0;j<atoms[iatom].numOfAngles;j++){
      a = atoms[iatom].angleAtomList[3*j];
      b = atoms[iatom].angleAtomList[3*j+1];
      c = atoms[iatom].angleAtomList[3*j+2];

      if(term_needed(movedatoms,a,b,c) && ((iatom==a && iatom<c) || (iatom==c && iatom<a))){//avoid double counting bond angle energies
        ba[0] = coords[3*a]   - coords[3*b];
        ba[1] = coords[3*a+1] - coords[3*b+1];
        ba[2] = coords[3*a+2] - coords[3*b+2];

        bc[0] = coords[3*c]   - coords[3*b];
        bc[1] = coords[3*c+1] - coords[3*b+1];
        bc[2] = coords[3*c+2] - coords[3*b+2];

        theta = angle(ba,bc);

        type   = atoms[iatom].angleParamType[j];
        K      = angleParams[type].K;
        theta0 = angleParams[type].theta0;

        energies[EN_ANGLE] += K*(theta - theta0)*(theta - theta0);
        count++;
#ifdef DEBUG_NON_TABULATED
        printf("angle: %d %d %d\n",a,b,c);
#endif
        //We assume here that if the term is needed, then fragments "a" and "c" are not the same.
        //if (atoms[a].fragment!=atoms[c].fragment) {
        /*&    dx = coords[3*c]   - coords[3*a];
            dy = coords[3*c+1] - coords[3*a+1];
            dz = coords[3*c+2] - coords[3*a+2];
            nonbond_energy(rdie,atoms[a].type,atoms[c].type,FALSE,dx,dy,dz,&evdw,&eelec);
            //printf("1-3 interaction %d %d %d %d %.4f\n",a,c,atoms[a].fragment,atoms[c].fragment,evdw);
            energies[EN_VDW1213]+=evdw;
            energies[EN_ELEC1213]+=eelec;*/
        //}
      }
    }
  }
  //printf("angle bending: %f kcal/mol, the number of interactions: %d\n",engAngle,count);

  //calculate improper dihedral energy
#ifdef TIMERS
  switch_timer(TIMER_NT_IMPROPERS);
#endif
  count=0;
  for(iatom=0;iatom<numOfAtoms;iatom++){
    for(j=0;j<atoms[iatom].numOfImprops;j++){
      a = atoms[iatom].impropAtomList[4*j];
      b = atoms[iatom].impropAtomList[4*j+1];
      c = atoms[iatom].impropAtomList[4*j+2];
      d = atoms[iatom].impropAtomList[4*j+3];
 
      //if(term_needed(movedatoms,a,b,c,d) && (iatom==a)&&(a<d)){//avoid multi-counting improper dihedral energies
      if(term_needed(movedatoms,a,b,c,d) && (iatom==c)){
          //The CHARMM 19 force field needs to be defined like this, replacing the "BC" axis with the "AD" axis.
        ba[0] = coords[3*a]   - coords[3*b];
        ba[1] = coords[3*a+1] - coords[3*b+1];
        ba[2] = coords[3*a+2] - coords[3*b+2];
        ad[0] = coords[3*a]   - coords[3*d];
        ad[1] = coords[3*a+1] - coords[3*d+1];
        ad[2] = coords[3*a+2] - coords[3*d+2];

        /*bc[0] = coords[3*c]   - coords[3*b];
        bc[1] = coords[3*c+1] - coords[3*b+1];
        bc[2] = coords[3*c+2] - coords[3*b+2];*/

        cd[0] = coords[3*d]   - coords[3*c];
        cd[1] = coords[3*d+1] - coords[3*c+1];
        cd[2] = coords[3*d+2] - coords[3*c+2];

        //dihed = acos(cos_dihedral(ba,ad,cd));
        dihed=dihedral(ba,ad,cd);
        //dihed = M_PI-fabs(dihed);

        type = atoms[iatom].impropParamType[j];
        K = impropParams[type].K;
        chi0 = impropParams[type].chi0;
#ifdef DEBUG_NON_TABULATED
        printf("improper: %d %d %d %d \n",a,b,c,d);//,dihed*DEG_TO_RAD,
#endif
        //energies[EN_IMPROPER] += impropParams[type].K*(1.0-cos(2.0*(dihed-chi0)));
        energies[EN_IMPROPER] += K*(dihed-chi0)*(dihed-chi0);
        //printf("%d %d %d %d  %d  %f  %f  %.4f\n",a+1,b+1,c,d,type,impropParams[type].V,dihed*RAD_TO_DEG,eng);

        count++;
      }
    }
  }

  //printf("improper dihedral: %f kcal/mol, the number of interactions: %d\n",engImprop,count);

  //calculate dihedral energy
#ifdef TIMERS
  switch_timer(TIMER_NT_DIHEDRALS);
#endif
  count=0;
  for(iatom=0;iatom<numOfAtoms;iatom++){
    for(j=0;j<atoms[iatom].numOfBonded14Atoms;j++){
      a = atoms[iatom].bonded14AtomList[4*j];
      b = atoms[iatom].bonded14AtomList[4*j+1];
      c = atoms[iatom].bonded14AtomList[4*j+2];
      d = atoms[iatom].bonded14AtomList[4*j+3];

      if(term_needed(movedatoms,a,b,c,d) && (d>a)){//avoid double counting dihedral energies
        ba[0] = coords[3*a]   - coords[3*b];
        ba[1] = coords[3*a+1] - coords[3*b+1];
        ba[2] = coords[3*a+2] - coords[3*b+2];

        bc[0] = coords[3*c]   - coords[3*b];
        bc[1] = coords[3*c+1] - coords[3*b+1];
        bc[2] = coords[3*c+2] - coords[3*b+2];

        cd[0] = coords[3*d]   - coords[3*c];
        cd[1] = coords[3*d+1] - coords[3*c+1];
        cd[2] = coords[3*d+2] - coords[3*c+2];

        dihed = cos_dihedral(ba,bc,cd);
        type  = atoms[iatom].bonded14DihedParamType[j];
        /*en = MDE(dihed,type);
        printf("%d %s %s %s %s %d %d %d %d %.4f %.4f %.4f %.4f %.2f %.4f\n",count,atoms[a].name,atoms[b].name,atoms[c].name,atoms[d].name,
            atoms[a].classx,atoms[b].classx,atoms[c].classx,atoms[d].classx,
            dihedParams[type].V[0],dihedParams[type].V[1],dihedParams[type].V[2],dihedParams[type].V[3],dihed*RAD_TO_DEG,en);*/
        energies[EN_DIHEDRAL] += MDE(dihed,type);
        count++;
#ifdef DEBUG_NON_TABULATED
        printf("dihedral: %d %d %d %d\n",a,b,c,d);
#endif
        //Include 1-4 nonbonded interactions between atoms belonging to the same fragment.  Make use of special 1-4 parameters when necessary.
        //Subtract any possibly erroneous VDW interactions involving atoms not belonging to the same fragment.
        //if (atoms[a].fragment!=atoms[d].fragment) {
        /*    dx = coords[3*d]   - coords[3*a];
            dy = coords[3*d+1] - coords[3*a+1];
            dz = coords[3*d+2] - coords[3*a+2];
            nonbond_energy(rdie,atoms[a].type,atoms[d].type,TRUE,dx,dy,dz,&evdw2,&eelec);
            //energies[EN_VDW14]+=BONDED14_SCALE*evdw;
            nonbond_energy(rdie,atoms[a].type,atoms[d].type,FALSE,dx,dy,dz,&evdw,&eelec);
            energies[EN_VDW14]+=(evdw2-evdw);
            //if (atoms[a].fragment!=atoms[d].fragment)
            //energies[EN_VDW14]-=evdw;
            //Already been counted if both atoms are in different fragments.

            //if (atoms[a].fragment==atoms[d].fragment) {
            //energies[EN_ELEC14]+=BONDED14_SCALE*eelec;
            //if (atoms[a].fragment!=atoms[d].fragment)
            //energies[EN_ELEC14]-=eelec;
        //}*/
      }
    }
  }
 //Calculate all interaction terms that need to be calculated exactly.
 #ifdef TIMERS
  switch_timer(TIMER_NT_VDW_ELEC);
#endif
  for (i=0; i<nb_atom_list_size; i++) {
      iatom=nb_atom_list[i].iatom;
      jatom=nb_atom_list[i].jatom;
      //The interaction needs to be computed if exactly one of iatom/jatom  has been moved (since all moves are rigid on subsets of atoms)
      //I think this is ok because the C++ standard requires that "bool" values be casted to int, with false=0 and true=1.
      if (!(movedatoms[iatom]^movedatoms[jatom])) continue;
      is14=nb_atom_list[i].is14;
      dx=coords[3*jatom]-coords[3*iatom];
      dy=coords[3*jatom+1]-coords[3*iatom+1];
      dz=coords[3*jatom+2]-coords[3*iatom+2];
      nonbond_energy(rdie,atoms[iatom].type,atoms[jatom].type,is14,dx,dy,dz,&evdw,&eelec);
      energies[EN_VDW_EXACT]+=evdw;
      energies[EN_ELEC_EXACT]+=eelec;
#ifdef DEBUG_NON_TABULATED
      printf("nonbond: %d %d\n",iatom,jatom);
#endif
  }
  energies[EN_ELEC_EXACT]=energies[EN_ELEC_EXACT]*COUL_CONST/eps;
#ifdef TIMERS
  switch_timer(TIMER_OTHER);
#endif
  //printf("dihedral: %f kcal/mol, the number of interactions: %d\n",engDihed,count);

}

void forcefield::non_tabulated_energy(double eps, int rdie, int numOfAtoms, ATOMS * atoms, int full_atom_list_count, atom_nb_entry * non_tab_list, double * coords, double * energies)
{
  int i,j,k;
  int ij;
  int ires;
  int iatom,jatom;
  int a,b,c,d;
  int type;
  int iclass,jclass;
  int ilibsc;
  int icsc;
  int count;
  double nbflag;
  double fi,temp;
  double ba[3],bc[3],cd[3],ad[3];
  double dx,dy,dz;
  double r2,r6;
  double K;
  double r,r0;
  double theta,theta0,chi0;
  double dihed;
  double en;
  double evdw,evdw2,eelec; /*for corrections*/
  bool is14;
  //zero out all energy terms excapt interactions, which may have been previously calculated

  //for (i=1; i<EN_TERMS; i++) energies[i]=0.0;

  //calculate bond energies
  count=0;
  for(iatom=0;iatom<numOfAtoms;iatom++){
    for(j=0;j<atoms[iatom].numOfBondedAtoms;j++){
      jatom = atoms[iatom].bondedAtomList[j];
      if(jatom>iatom){//avoid double counting bond energies b/c if iatom contains jatom, jatom contains iatom
        dx = coords[3*jatom]   - coords[3*iatom];
        dy = coords[3*jatom+1] - coords[3*iatom+1];
        dz = coords[3*jatom+2] - coords[3*iatom+2];
        r  = sqrt(dx*dx + dy*dy + dz*dz);

        type = atoms[iatom].bondedParamType[j];
        K    = bondParams[type].K;
        r0   = bondParams[type].r0;
	en = K*(r - r0)*(r - r0);
        energies[EN_BOND] += en;
        //only interfragment interactions need to be counted here!
        //if (atoms[iatom].fragment!=atoms[jatom].fragment) {
        /*    nonbond_energy(rdie,atoms[iatom].type,atoms[jatom].type,FALSE,dx,dy,dz,&evdw,&eelec);
            energies[EN_VDW1213]+=evdw;
            energies[EN_ELEC1213]+=eelec;*/
            //printf("1-2 interaction %d %d %d %d %.4f\n",iatom,jatom,atoms[iatom].fragment,atoms[jatom].fragment,evdw);
        //}
#ifdef DEBUG_NON_TABULATED
        printf("Bond: %d %d  %.10f %.10f\n",iatom,jatom,r,en);
#endif
        count++;
      }
    }
  }

  //printf("bond stretching: %f kcal/mol, the number of interactions: %d\n",engBond,count);

  //calculate bond angle energies

  count=0;
  for(iatom=0;iatom<numOfAtoms;iatom++){
    for(j=0;j<atoms[iatom].numOfAngles;j++){
      a = atoms[iatom].angleAtomList[3*j];
      b = atoms[iatom].angleAtomList[3*j+1];
      c = atoms[iatom].angleAtomList[3*j+2];

      if((iatom==a && iatom<c) || (iatom==c && iatom<a)){//avoid double counting bond angle energies
        ba[0] = coords[3*a]   - coords[3*b];
        ba[1] = coords[3*a+1] - coords[3*b+1];
        ba[2] = coords[3*a+2] - coords[3*b+2];

        bc[0] = coords[3*c]   - coords[3*b];
        bc[1] = coords[3*c+1] - coords[3*b+1];
        bc[2] = coords[3*c+2] - coords[3*b+2];

        theta = angle(ba,bc);

        type   = atoms[iatom].angleParamType[j];
        K      = angleParams[type].K;
        theta0 = angleParams[type].theta0;

        energies[EN_ANGLE] += K*(theta - theta0)*(theta - theta0);
        count++;
#ifdef DEBUG_NON_TABULATED
        printf("Angle: %d %d %d %d\n",count,a,b,c);
#endif
        //if (atoms[a].fragment!=atoms[c].fragment) {
        /*    dx = coords[3*c]   - coords[3*a];
            dy = coords[3*c+1] - coords[3*a+1];
            dz = coords[3*c+2] - coords[3*a+2];
            nonbond_energy(rdie,atoms[a].type,atoms[c].type,FALSE,dx,dy,dz,&evdw,&eelec);
            //printf("1-3 interaction %d %d %d %d %.4f\n",a,c,atoms[a].fragment,atoms[c].fragment,evdw);
            energies[EN_VDW1213]+=evdw;
            energies[EN_ELEC1213]+=eelec;*/
        //}
      }
    }
  }
  //printf("angle bending: %f kcal/mol, the number of interactions: %d\n",engAngle,count);

  //calculate improper dihedral energy

  count=0;
  for(iatom=0;iatom<numOfAtoms;iatom++){
    for(j=0;j<atoms[iatom].numOfImprops;j++){
      a = atoms[iatom].impropAtomList[4*j];
      b = atoms[iatom].impropAtomList[4*j+1];
      c = atoms[iatom].impropAtomList[4*j+2];
      d = atoms[iatom].impropAtomList[4*j+3];

      if((iatom==c)){//avoid multi-counting improper dihedral energies
          //The CHARMM 19 force field needs to be defined like this, replacing the "BC" axis with the "AD" axis.
        ba[0] = coords[3*a]   - coords[3*b];
        ba[1] = coords[3*a+1] - coords[3*b+1];
        ba[2] = coords[3*a+2] - coords[3*b+2];
        ad[0] = coords[3*a]   - coords[3*d];
        ad[1] = coords[3*a+1] - coords[3*d+1];
        ad[2] = coords[3*a+2] - coords[3*d+2];

        /*bc[0] = coords[3*c]   - coords[3*b];
        bc[1] = coords[3*c+1] - coords[3*b+1];
        bc[2] = coords[3*c+2] - coords[3*b+2];*/

        cd[0] = coords[3*d]   - coords[3*c];
        cd[1] = coords[3*d+1] - coords[3*c+1];
        cd[2] = coords[3*d+2] - coords[3*c+2];

        //dihed = acos(cos_dihedral(ba,ad,cd));
        dihed = dihedral(ba,ad,cd);
        //dihed = M_PI-fabs(dihed);

        type = atoms[iatom].impropParamType[j];
        K = impropParams[type].K;
        chi0 = impropParams[type].chi0;
        //printf("%d %d %d %d %d \n",count,a+1,b+1,c+1,d+1);//,dihed*DEG_TO_RAD,
        //energies[EN_IMPROPER] += impropParams[type].K*(1.0-cos(2.0*(dihed-chi0)));
        en= K*(dihed-chi0)*(dihed-chi0);
        energies[EN_IMPROPER] += en;
#ifdef DEBUG_NON_TABULATED
        printf("Improper: %d %d %d %d  %d %.10f %.10f \n",a,b,c,d,type,dihed*RAD_TO_DEG,en);
#endif

        count++;
      }
    }
  }

  //printf("improper dihedral: %f kcal/mol, the number of interactions: %d\n",engImprop,count);

  //calculate dihedral energy

  count=0;
  for(iatom=0;iatom<numOfAtoms;iatom++){
    for(j=0;j<atoms[iatom].numOfBonded14Atoms;j++){
      a = atoms[iatom].bonded14AtomList[4*j];
      b = atoms[iatom].bonded14AtomList[4*j+1];
      c = atoms[iatom].bonded14AtomList[4*j+2];
      d = atoms[iatom].bonded14AtomList[4*j+3];

      if (d>a){//avoid double counting dihedral energies
        ba[0] = coords[3*a]   - coords[3*b];
        ba[1] = coords[3*a+1] - coords[3*b+1];
        ba[2] = coords[3*a+2] - coords[3*b+2];

        bc[0] = coords[3*c]   - coords[3*b];
        bc[1] = coords[3*c+1] - coords[3*b+1];
        bc[2] = coords[3*c+2] - coords[3*b+2];

        cd[0] = coords[3*d]   - coords[3*c];
        cd[1] = coords[3*d+1] - coords[3*c+1];
        cd[2] = coords[3*d+2] - coords[3*c+2];

        dihed = cos_dihedral(ba,bc,cd);
        type  = atoms[iatom].bonded14DihedParamType[j];
        en = MDE(dihed,type);
        /*printf("%d %s %s %s %s %d %d %d %d %.4f %.4f %.4f %.4f %.2f %.4f\n",count,atoms[a].name,atoms[b].name,atoms[c].name,atoms[d].name,
            atoms[a].classx,atoms[b].classx,atoms[c].classx,atoms[d].classx,
            dihedParams[type].V[0],dihedParams[type].V[1],dihedParams[type].V[2],dihedParams[type].V[3],acos(dihed)*RAD_TO_DEG,en);*/
        energies[EN_DIHEDRAL] += MDE(dihed,type);
        count++;
        //printf("%d %d %d %d %d\n",count,a+1,b+1,c+1,d+1);
        //Include 1-4 nonbonded interactions between atoms belonging to the same fragment.  Make use of special 1-4 parameters when necessary.
        //Subtract any possibly erroneous VDW interactions involving atoms not belonging to the same fragment.
        //if (atoms[a].fragment!=atoms[d].fragment) {
        /*    dx = coords[3*d]   - coords[3*a];
            dy = coords[3*d+1] - coords[3*a+1];
            dz = coords[3*d+2] - coords[3*a+2];
            nonbond_energy(rdie,atoms[a].type,atoms[d].type,TRUE,dx,dy,dz,&evdw2,&eelec);
            //energies[EN_VDW14]+=BONDED14_SCALE*evdw;
            nonbond_energy(rdie,atoms[a].type,atoms[d].type,FALSE,dx,dy,dz,&evdw,&eelec);
            energies[EN_VDW14]+=evdw2-evdw;
            //if (atoms[a].fragment!=atoms[d].fragment) energies[EN_VDW14]-=evdw;
            //Already been counted if both atoms are in different fragments.

            //if (atoms[a].fragment==atoms[d].fragment) {
            energies[EN_ELEC14]+=(BONDED14_SCALE-1)*eelec;*/
            //if (atoms[a].fragment!=atoms[d].fragment) energies[EN_ELEC14]-=eelec;
        //}
      }
    }
  }

  //printf("dihedral: %f kcal/mol, the number of interactions: %d\n",engDihed,count);
  //Calculate all interaction terms that need to be calculated exactly.
  for (i=0; i<full_atom_list_count; i++) {
      iatom=non_tab_list[i].iatom;
      jatom=non_tab_list[i].jatom;
      is14=non_tab_list[i].is14;
      dx=coords[3*jatom]-coords[3*iatom];
      dy=coords[3*jatom+1]-coords[3*iatom+1];
      dz=coords[3*jatom+2]-coords[3*iatom+2];
      nonbond_energy(rdie,atoms[iatom].type,atoms[jatom].type,is14,dx,dy,dz,&evdw,&eelec);
      energies[EN_VDW_EXACT]+=evdw;
      energies[EN_ELEC_EXACT]+=eelec;
#ifdef DEBUG_NON_TABULATED
      printf("Nonbonded interaction: %d %d %c %d %d %.10f %.10f\n",iatom,jatom,yesno(is14),atoms[iatom].type,atoms[jatom].type,evdw,eelec);
#endif
  }
  energies[EN_ELEC_EXACT]=energies[EN_ELEC_EXACT]*COUL_CONST/eps;
}


void forcefield::find_parameters(int numOfAtoms, ATOMS * atoms)
{

  int i,j,k,m,n;
  int a,b,c,d;
  int ca,cc,cb;
  int ii,ij;
  int iclass,jclass,kclass,mclass;
  int iatom,jatom,katom,matom,natom,iiatom;
  int ires,jres,res;
  int resOld;
  int fragTypeOld;
  int num,num2;
  int count;
  int iFlag,nFlag;
  int flag,flag2;
  int chain;
  double sum,sum2;
  double r;
  double scale;
  double eps,sigma,sigma3,sigma6,sigma12;

  /*for(i=1;i<numOfAtoms;i++){
    for(j=0;j<5;j++){
      if(atoms[i].bondedAtomList[j]==0){
        atoms[i].numOfBondedAtoms=j;
        break;
      }
    }
    //printf("%d %d\n",i,atomComb[i].numOfBondedAtoms);
  }*/

  //fetch tinker atom class and atomic number for each atom
  for(i=0;i<numOfAtoms;i++){
    atoms[i].classx=atomTypeLookUp[atoms[i].type].classx;
    atoms[i].atomicNum=atomTypeLookUp[atoms[i].type].atomicNum;
    //printf("atom: %d, type: %d, class: %d\n",i,atoms[i].type,atoms[i].classx);
  }

  //fetch non-bonded params for each atom
  for(i=0;i<numOfAtoms;i++){
    atoms[i].radius   = vdwParams[atoms[i].classx].sigma;

    //printf("atom: %d, charge: %f, sigma: %f, epsilon: %f\n",i,atoms[i].q,atoms[i].sigma,atoms[i].epsilon);
  }

  //calculate the number of atoms within each residue
  /*res=1;
  count=0;
  for(iatom=0;iatom<numOfAtoms;iatom++){
    if(atoms[iatom].resNum!=res){
      resLookUp[res].numOfAtoms = count;
      res++;
      count=0;
    }else if(iatom==numOfAtoms){
      resLookUp[res].numOfAtoms = count+1;
      count=0;
    }
    count++;
  }*/
  for(iatom=0;iatom<numOfAtoms;iatom++)
  {
    //regular bonded1-2
    for(j=0;j<atoms[iatom].numOfBondedAtoms;j++)
    {
      jatom=atoms[iatom].bondedAtomList[j];
      for(k=0;k<numOfBondParams;k++)
      {
        if(atoms[iatom].classx==bondParams[k].atomClass[0] && atoms[jatom].classx==bondParams[k].atomClass[1] || atoms[iatom].classx==bondParams[k].atomClass[1] && atoms[jatom].classx==bondParams[k].atomClass[0])
        {
          atoms[iatom].bondedParamType[j] = k;
          //printf("%d %d  %d\n",iatom,jatom,k);
          break;
        }

      }
    }
  }
//prepare bonded 1-2, 1-3, 1-4 arrays

 //take care of 1-3 bonded atoms


  for (iatom=0; iatom<numOfAtoms; iatom++) {
    for (num=0; num<atoms[iatom].numOfAngles; num++) {
        jatom=atoms[iatom].angleAtomList[3*num+1];
        katom=atoms[iatom].angleAtomList[3*num+2];
        for(m=0;m<numOfAngleParams;m++)
        {
            if(atoms[iatom].classx==angleParams[m].atomClass[0] && atoms[jatom].classx==angleParams[m].atomClass[1] && atoms[katom].classx==angleParams[m].atomClass[2] || atoms[iatom].classx==angleParams[m].atomClass[2] && atoms[jatom].classx==angleParams[m].atomClass[1] && atoms[katom].classx==angleParams[m].atomClass[0])
            {
                atoms[iatom].angleParamType[num] = m;
                //printf("%d %d %d %d  %d\n",i,iatom,jatom,katom,k);
            break;
            }
        }
    }
    }

  printf("setup bond angle interactions: passed\n");


  /*
  //print out improper dihedrals
  for(iatom=0;iatom<numOfAtoms;iatom++){
    printf("atom: %d, involved in the following improper dihedrals: ",iatom);
    for(j=0;j<atoms[iatom].numOfImprops;j++){
      printf("%d %d %d %d  %d   ",atoms[iatom].impropAtomList[4*j],atoms[iatom].impropAtomList[4*j+1],atoms[iatom].impropAtomList[4*j+2],atoms[iatom].impropAtomList[4*j+3],atoms[iatom].impropParamType[j]);
    }
    printf("\n");
  }
  */

  printf("setup improper torsion interactions: passed\n");


		  //fetch forcefield parameters for bonded 1-4 interactions
  for(i=0;i<numOfAtoms;i++){
    for(j=0;j<atoms[i].numOfBonded14Atoms;j++){
      iatom = atoms[i].bonded14AtomList[4*j];
      jatom = atoms[i].bonded14AtomList[4*j+1];
      katom = atoms[i].bonded14AtomList[4*j+2];
      matom = atoms[i].bonded14AtomList[4*j+3];
      for(k=0;k<numOfDihedParams;k++){
        if(atoms[iatom].classx==dihedParams[k].atomClass[0] && atoms[jatom].classx==dihedParams[k].atomClass[1] && atoms[katom].classx==dihedParams[k].atomClass[2] && atoms[matom].classx==dihedParams[k].atomClass[3] || atoms[iatom].classx==dihedParams[k].atomClass[3] && atoms[jatom].classx==dihedParams[k].atomClass[2] && atoms[katom].classx==dihedParams[k].atomClass[1] && atoms[matom].classx==dihedParams[k].atomClass[0]){
          atoms[i].bonded14DihedParamType[j] = k;
          //printf("%d %d %d %d %d  %d\n",i,iatom,jatom,katom,matom,k);
          break;
        }
      }
    }
  }
    /*for(j=0;j<atoms[i].numOfBonded14AddAtoms;j++){
      iatom = atoms[i].bonded14AddAtomList[4*j];
      jatom = atoms[i].bonded14AddAtomList[4*j+1];
      katom = atoms[i].bonded14AddAtomList[4*j+2];
      matom = atoms[i].bonded14AddAtomList[4*j+3];
      for(k=0;k<numOfDihedParams;k++){
        if(atoms[iatom].classx==dihedParams[k].atomClass[0] && atoms[jatom].classx==dihedParams[k].atomClass[1] && atoms[katom].classx==dihedParams[k].atomClass[2] && atoms[matom].classx==dihedParams[k].atomClass[3] || atoms[iatom].classx==dihedParams[k].atomClass[3] && atoms[jatom].classx==dihedParams[k].atomClass[2] && atoms[katom].classx==dihedParams[k].atomClass[1] && atoms[matom].classx==dihedParams[k].atomClass[0]){
          atoms[i].bonded14AddDihedParamType[j] = k;
          break;
        }
      }
      //printf("iatom: %d, %d %d %d %d\n",i,iatom,jatom,katom,matom);
    }*/

  /*
  //print out bonded 1-4 atoms
  for(iatom=0;iatom<numOfAtoms;iatom++){
    printf("atom: %d, has bonded 1-4 atoms: ",iatom);
    for(j=0;j<atoms[iatom].numOfBonded14Atoms;j++){
      printf("%d %d %d %d   ",atoms[iatom].bonded14AtomList[4*j],atoms[iatom].bonded14AtomList[4*j+1],atoms[iatom].bonded14AtomList[4*j+2],atoms[iatom].bonded14AtomList[4*j+3]);
    }
    printf("\n");
  }
  */

  printf("setup bonded1-4 interactions: passed\n");

}

//void forcefield::virtual_hydrogen_reconstruct(ATOMS * atoms,
