/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_ttemp.h"
#include "atom.h"
#include "atom_masks.h"
#include "accelerator_kokkos.h"
#include "update.h"
#include "modify.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */
/*Originally taken from fix addforce*/
FixTwoTemp::FixTwoTemp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), dTi_dt(NULL)
{

  if (narg != 19) error->all(FLERR,"Illegal fix ttemp command");

  if (strcmp(arg[3],"K0") == 0) K0 = force->numeric(FLERR,arg[4]);

  if (strcmp(arg[5],"dt") == 0) HeatTransportTimeStep = force->numeric(FLERR,arg[6]);

  if (strcmp(arg[7],"Tmax") == 0) Tmax = force->numeric(FLERR,arg[8]);

  if (strcmp(arg[9],"loop") == 0) loop = force->numeric(FLERR,arg[10]);

  if (strcmp(arg[11],"rdcut") == 0) rdcut = force->numeric(FLERR,arg[12]);

  if (strcmp(arg[13],"G") == 0) G = force->numeric(FLERR,arg[14]);

  if (strcmp(arg[15],"Vatom") == 0) Vatom = force->numeric(FLERR,arg[16]);
  
  if (strcmp(arg[17],"lambda") == 0) lambda = force->numeric(FLERR,arg[18]);

  if (narg != 19) error->all(FLERR,"Illegal fix ttemp command");

  comm_forward = 1;
  if (force->newton_pair) comm_reverse = 1; //size from ghost atom to own atoms

  dynamic_group_allow = 1;
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;

  xstr = ystr = zstr = NULL;

  xvalue = 0.0;
  xstyle = CONSTANT;

  yvalue = 0.0;
  ystyle = CONSTANT;

  zvalue = 0.0;
  zstyle = CONSTANT;

  // optional args

  nevery = 1;
  iregion = -1;
  idregion = NULL;
  estr = NULL;

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;

  nmax = 0;
}

/* ---------------------------------------------------------------------- */

FixTwoTemp::~FixTwoTemp()
{
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  delete [] estr;
  delete [] idregion;
  memory->destroy(dTi_dt);
}

/* ---------------------------------------------------------------------- */

int FixTwoTemp::setmask()
{
  datamask_read = datamask_modify = 0;

  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTwoTemp::init()
{

  // need a full neighbor list, built whenever re-neighboring occurs

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->fix = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  estyle = NONE;

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (varflag == CONSTANT && estyle != NONE)
    error->all(FLERR,"Cannot use variable energy with "
               "constant force in fix addforce");
  if ((varflag == EQUAL || varflag == ATOM) &&
      update->whichflag == 2 && estyle == NONE)
    error->all(FLERR,"Must use variable energy with fix addforce");

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTwoTemp::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTwoTemp::min_setup(int vflag)
{
  post_force(vflag);
}

void FixTwoTemp::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */


// Function call for exponential calculation


double FixTwoTemp::exp16(double exp_temp)
{
	
	exp_temp = 1.0 + exp_temp / 16.0;
	exp_temp *= exp_temp; exp_temp *= exp_temp; 
	exp_temp *= exp_temp; exp_temp *= exp_temp;
	
	return exp_temp;	
	
}

/* ---------------------------------------------------------------------- */


void FixTwoTemp::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  double mvv2e = force->mvv2e;
  double **v = atom->v;
  double *temp = atom->mxe_temperature; //Needed for max-ent
  double *mass = atom->mass; //Needed for max-ent
  int *type = atom->type;
  int *mask = atom->mask;
  int newton_flag = force->newton_pair;
  imageint *image = atom->image;
  int nlocal = atom->nlocal;
  bigint natoms = atom->natoms;

  int i,j,ii,jj,inum,jnum,iloop,itype,count_ij;   // count_ij = count the number of atoms within the cutoff for temp_latt and ken
  double dx,dy,dz,rsq; // This value should be given in the script!!!
  double ken, Ke=0.0, xi = 0.0, e_dt = 0.0, dT = 0.0;  //dTi_dt = 0.0,
  double Ce = 0.0, Kt = 0.0, kB = force->boltz; //eV/sec

  int 	 temp_flag[nlocal];
  double temp_latt[nlocal];
  double com[4], COM[4], lKe = 0.0; //Center of motion velocity (vx, vy, vz, Mt. com: local, COM: Global). lKe: Local Kinetic energy of the node
  com[0] = 0.0, com[1] = 0.0, com[2] = 0.0, com[3] = 0.0;
  COM[0] = 0.0, COM[1] = 0.0, COM[2] = 0.0, COM[3] = 0.0;

  // grow dTi_dt array if necessary  
  if (atom->nmax > nmax) {
    memory->destroy(dTi_dt);
    nmax = atom->nmax;
    memory->create(dTi_dt,nmax,"ttemp");
  }

  // set to zeros the arrays 
  for(i=0; i<nmax; i++) dTi_dt[i]=0.0;

  // set to zeros the arrays 
  for(i=0; i<nlocal; i++) {
	temp_latt[i]=0.0;
	temp_flag[i]=1;
  }

  double tfactor;

  int *ilist, *jlist,*numneigh,**firstneigh; 

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  if (update->ntimestep % nevery) return;

  ///////////////////// MAIN MODIFICATION //////////////////////////////// 
  //variable initialization
  double t = 0.0, t_temp = 0.0, scalar = 0.0;
  ken = 0.0;
  double e_temp = 0.0, E_Temp = 0.0;
  double c_atom = 0;
  double C_Atom = 0;

  // calculate lattice temp, electronic temp, atom number and kinetic energy  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      t_temp = (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * mass[type[i]];
      ken += 0.5 * mvv2e * t_temp;
      t += t_temp;
      e_temp += temp[i];
    }
  }

  double data[3], Data[3];
  data[0] = t;
  data[1] = e_temp;
  data[2] = ken;
 
  MPI_Allreduce(&data[0],&Data[0],3,MPI_DOUBLE,MPI_SUM,world);
 
  scalar = Data[0];
  E_Temp = Data[1];
  Ke     = Data[2];

  tfactor = force->mvv2e / (natoms*3*force->boltz);
  scalar *= tfactor;

  // Calculate temperature dependent Ce and Kt
  // Ce = lambda*Te;
  // Ke = K0*Te/Tl; Kt = (2*Ke*d)/(Ce*Z*b*b) = (2*K0*Te*d)/(Ce*Z*Tl*b*b)
  // User will provide in the input script K0 = K = (2*K0*d)/(Z*b*b) [These are all constants]
  // Hence Kt = (K0*Te)/(Ce*Tl) 
  
  Ce = lambda*E_Temp/natoms;             
  Kt = (K0*E_Temp)/(Ce*scalar*natoms);   

   /////////////////////\\\\\\\\\\\\\\\\\\\\\\//////////////////////////////// 

  // variable timestep for electronic heat conduction 
  e_dt = (update->dt)/loop;

  //Scan atoms that are going to exchange energy
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
 
    if (mask[i] & groupbit) {    	
      com[0] += mass[type[i]]*v[i][0];
      com[1] += mass[type[i]]*v[i][1];
      com[2] += mass[type[i]]*v[i][2];
      com[3] += mass[type[i]];
      ken = 0.5 * mvv2e * mass[type[i]] * (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
      if(ken < 1.0e-4) temp_flag[i] = 0;
    }
  }

  MPI_Allreduce(&com[0],&COM[0],4,MPI_DOUBLE,MPI_SUM,world);

  //V_j = sum_i^N v_j^i * m^i / Mt
  COM[0] /= COM[3];
  COM[1] /= COM[3];
  COM[2] /= COM[3];

  double KeCOM = 0.5 * COM[3] * (COM[0]*COM[0] + COM[1]*COM[1] + COM[2]*COM[2]); //1/2 mass v_com^2
  double TCOM = KeCOM * force->mvv2e / (natoms*3.0*force->boltz);

  //Compute electronic heat conduction
  for (iloop = 0; iloop < loop; iloop++) {

  comm->forward_comm_fix(this);  

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    dTi_dt[i] = 0.0;
    temp_latt[i] = 0.0;
    count_ij = 0;
    t = 0.0;
    
    if (mask[i] & groupbit) {
    	
      t += (v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]) * mass[type[i]];
      count_ij++;
    	
    	for (jj = 0; jj < jnum; jj++) {
    		
    	  j = jlist[jj];
    	  j &= NEIGHMASK;

          dx = x[i][0] - x[j][0];
          dy = x[i][1] - x[j][1];
          dz = x[i][2] - x[j][2];
          rsq = dx*dx + dy*dy + dz*dz;

          t += (v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]) * mass[type[j]];
          count_ij++;

          if (rsq < rdcut*rdcut) {
      	      
      	      // Electronic heat conduction //
      	      double theta = -2.0*(temp[j] - temp[i])/(Tmax);

      	      dTi_dt[i] += Kt*( temp[j]/Tmax*(1.0-temp[i]/Tmax)*exp16( theta) 
		              - temp[i]/Tmax*(1.0-temp[j]/Tmax)*exp16(-theta) );
      	     }
          } // end jj
                    
          // Modifiation of equation of motion for TTM
          tfactor = force->mvv2e / (count_ij*3.0*force->boltz);
          t *= tfactor;	
          temp_latt[i] = t - TCOM; //Kelvin
     }
  } //end ii
	
  // reverse comm of dx_dt -> from ghost atoms to own atoms
  if (force->newton_pair){
   	     comm->reverse_comm_fix(this);
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (mask[i] & groupbit) {
    	
	if(temp_flag[i] == 1) {

	    	if (HeatTransportTimeStep > 0.0) {
	    		temp[i] += (Tmax*dTi_dt[i]*HeatTransportTimeStep) - ((G/Ce)*(temp[i] - temp_latt[i]))*HeatTransportTimeStep;
	    	} else {        	
			// variable timestep from lammps
			temp[i] += (Tmax*dTi_dt[i]*e_dt) - ((G/Ce)*(temp[i] - temp_latt[i]))*e_dt;    
		}
	}
	else {
	    	if (HeatTransportTimeStep > 0.0) {
	    		temp[i] += (Tmax*dTi_dt[i]*HeatTransportTimeStep)*HeatTransportTimeStep;
	    	} else {        	
			// variable timestep from lammps
			temp[i] += (Tmax*dTi_dt[i]*e_dt)*e_dt;    
		}
	}
    }  // end if
  }  // end for
 } // end loop

  //// Change/Addition Here -------------- Mohammad

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    if (mask[i] & groupbit) {

      if(temp_flag[i] == 1) {

       ken = 0.5 * mvv2e * mass[type[i]] * (v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);

       //dimensionless Note: G*Vatom should	 be in the same units (Angstrom^3) 	
       xi = G*Vatom*(temp[i] - temp_latt[i])/(2.0*ken);

//       if ( fabs(xi) > 1.0e5) printf("Warning: xi is %8.4e, temp_e %8.4e, temp_l %8.4e, ken %8.4e \n", xi, temp[i], temp_latt[i], ken);

       f[i][0] += xi*mvv2e*mass[type[i]]*v[i][0];
       f[i][1] += xi*mvv2e*mass[type[i]]*v[i][1];
       f[i][2] += xi*mvv2e*mass[type[i]]*v[i][2];
      }
    }
  }

}

/* ---------------------------------------------------------------------- */

void FixTwoTemp::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTwoTemp::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixTwoTemp::compute_scalar()
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[0];
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixTwoTemp::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal,foriginal_all,4,MPI_DOUBLE,MPI_SUM,world);
    force_flag = 1;
  }
  return foriginal_all[n+1];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixTwoTemp::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

int FixTwoTemp::pack_forward_comm(int n, int *list, double *buf,
                                  int pbc_flag, int *pbc)
{
  int m;
    for(m = 0; m < n; m++) buf[m] = atom->mxe_temperature[list[m]];
  return n;
}

/* ---------------------------------------------------------------------- */

void FixTwoTemp::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m;

    for(m = 0, i = first; m < n; m++, i++) atom->mxe_temperature[i] = buf[m];
}

/* ---------------------------------------------------------------------- */

int FixTwoTemp::pack_reverse_comm(int n, int first, double *buf)
{    
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = dTi_dt[i];

  return m;
  
}

/* ---------------------------------------------------------------------- */

void FixTwoTemp::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
	j = list[i];
	dTi_dt[j] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */
