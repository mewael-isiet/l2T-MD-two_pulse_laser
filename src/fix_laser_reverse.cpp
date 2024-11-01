/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   MXE implementation: M. Ponga (UBC)
------------------------------------------------------------------------- */

#include <string.h>
#include <algorithm>
#include <vector>
#include <stdlib.h>
#include "fix_laser_reverse.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "domain.h"
#include "comm.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLaserReverse::FixLaserReverse(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  // fix 1 all laser_reverse I0 1.0 R 0.1 Lp 50.0 t0 0.50 sigma0 0.1 x 0.0
  if (narg != 15) error->all(FLERR,"Illegal fix laser reverse command");

  dynamic_group_allow = 1;
#if 0
	S(z,t) = I_0 (1-R)L_p^{-1} \exp[-z/L_p] \exp[-(t-t_0)^2/(2*\sigma_0^2)]
	This fix implements a heat flow in the atoms given by the following expression:
	$q_{i} = i_o \exp(-x_i/l_p)\exp(-(t-t_0)/(2*\sigma^2))$ [K $\cdot$ ps] 
	and the rate of change og the energy associated to each atom is $\dot{e}_i = k_B q_i$ 
	where $k_B$ is the Boltzmann constant. The amount of energy absorbed by the sample is 
	e_i(x_i) = \int_0^t \dot{e}_i dt = \int_0^t  i_o \exp(-x_i/l_p)\exp(-(t-t_0)/(2*\sigma^2)) dt.
#endif
  if (strcmp(arg[3],"F") == 0) F = force->numeric(FLERR,arg[4]);
  if (strcmp(arg[5],"R" ) == 0) R = force->numeric(FLERR,arg[6]);
  if (strcmp(arg[7],"Lp") == 0) Lp = force->numeric(FLERR,arg[8]);
  if (strcmp(arg[9],"t0") == 0) t0 = force->numeric(FLERR,arg[10]);
  if (strcmp(arg[11],"lambda") == 0) lambda = force->numeric(FLERR,arg[12]);
  if (strcmp(arg[13],"x") && strcmp(arg[13],"y") && strcmp(arg[13],"z"))
    error->all(FLERR,"Illegal laser reverse command");
  axis = arg[13][0];

  if (axis == 'x') {
    dim_flag = 0;
    c1 = force->numeric(FLERR,arg[14]);
  } else if (axis == 'y') {
    dim_flag = 1;
    c1 = force->numeric(FLERR,arg[14]);
  } else if (axis == 'z') {
    dim_flag = 2;
    c1 = force->numeric(FLERR,arg[14]);
  }

  time  = 0.0;
  time0 = 0.0;
}

/* ---------------------------------------------------------------------- */

FixLaserReverse::~FixLaserReverse()
{

}

/* ---------------------------------------------------------------------- */

int FixLaserReverse::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLaserReverse::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixLaserReverse::init()
{
    // Some checks
    
    if (force->pair == NULL)
        error->all(FLERR,"Fix optfreq requires a pair style be defined");
    
    int count = 0;
    for (int i = 0; i < modify->nfix; i++) {
      if (strcmp(modify->fix[i]->style,"optfreq") == 0) count++;
      if (count > 1 && comm->me == 0)
         error->warning(FLERR,"More than one fix optfreq");
    }
    time0 = update->dt*update->ntimestep;
}

/* ---------------------------------------------------------------------- */

void FixLaserReverse::initial_integrate(int vflag)
{
  double **x = atom->x;
  double *mxe_temperature = atom->mxe_temperature;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt = 0.0;
  double value;
  time = update->dt*update->ntimestep-time0;
  double local_min_x, local_max_x, local_min_y, local_max_y, local_min_z, local_max_z;
  double global_min_x, global_max_x, global_min_y, global_max_y, global_min_z, global_max_z;
  std::vector<double> vec_x, vec_y, vec_z;
  double Ce[nlocal];

  for(int l = 0; l < nlocal; l++){
    if (mask[l] & groupbit) {
          vec_x.push_back(x[l][0]);
          vec_y.push_back(x[l][1]);
          vec_z.push_back(x[l][2]);
    }
  }

  local_min_x = *std::min_element(vec_x.begin(), vec_x.end());
  local_max_x = *std::max_element(vec_x.begin(), vec_x.end());
  local_min_y = *std::min_element(vec_y.begin(), vec_y.end());
  local_max_y = *std::max_element(vec_y.begin(), vec_y.end());
  local_min_z = *std::min_element(vec_z.begin(), vec_z.end());
  local_max_z = *std::max_element(vec_z.begin(), vec_z.end());
  vec_x.clear();
  vec_y.clear();
  vec_z.clear();

  MPI_Allreduce(&local_min_x, &global_min_x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max_x, &global_max_x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&local_min_y, &global_min_y, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max_y, &global_max_y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&local_min_z, &global_min_z, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&local_max_z, &global_max_z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
	    Ce[i] = lambda*mxe_temperature[i];
	    if(dim_flag == 0) mxe_temperature[i] += (0.9394372787*F*(1-R)/(Lp*t0*Ce[i]))*exp(-(c1 - (x[i][0] - global_min_x)/(global_max_x - global_min_x))/(Lp/(global_max_x - global_min_x)))*exp(-2.7725887222*(time-2*t0)*(time-2*t0)/(t0*t0))*update->dt;
	    if(dim_flag == 1) mxe_temperature[i] += (0.9394372787*F*(1-R)/(Lp*t0*Ce[i]))*exp(-(c1 - (x[i][1] - global_min_x)/(global_max_x - global_min_x))/(Lp/(global_max_x - global_min_x)))*exp(-2.7725887222*(time-2*t0)*(time-2*t0)/(t0*t0))*update->dt;
	    if(dim_flag == 2) mxe_temperature[i] += (0.9394372787*F*(1-R)/(Lp*t0*Ce[i]))*exp(-(c1 - (x[i][2] - global_min_x)/(global_max_x - global_min_x))/(Lp/(global_max_x - global_min_x)))*exp(-2.7725887222*(time-2*t0)*(time-2*t0)/(t0*t0))*update->dt;    
    }
  } 
}
/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixLaserReverse::memory_usage()
{
  double bytes = 0.0;
  return bytes;
}
