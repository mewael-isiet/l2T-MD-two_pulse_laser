/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   MXE implementation: M. Ponga (UBC) and J.P. Mendez (Caltech)
------------------------------------------------------------------------- */

#include <string.h>
#include <stdlib.h>
#include "fix_heatflow.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
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

FixHeatFlow::FixHeatFlow(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  // fix fix_id group_id optfreq xvalue yvalue zvalue 
  // [0]   [1]    [2]     [3]     [4]    [5]     [6] 
  if (narg < 6 || narg > 6) error->all(FLERR,"Illegal fix optfreq command");

  dynamic_group_allow = 1;

  //Q_in is the heat flow in/out at the i-th atom. The fix actually does $T_i [K] = q_in [K/ps] * \Delta t [ps]$
  qin =    force->numeric(FLERR,arg[3]);
  //yvalue = force->numeric(FLERR,arg[4]);
  //zvalue = force->numeric(FLERR,arg[5]);
  
}

/* ---------------------------------------------------------------------- */

FixHeatFlow::~FixHeatFlow()
{

}

/* ---------------------------------------------------------------------- */

int FixHeatFlow::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::init()
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

}

/* ---------------------------------------------------------------------- */

void FixHeatFlow::initial_integrate(int vflag)
{

  double *freq_force = atom->freq_force;
  double *frequency = atom->frequency;
  double *mxe_temperature = atom->mxe_temperature;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double ptime = xvalue; //pseudo-time for dynamic relaxation
  double pdamping = yvalue; //pseudo-damping for dynamic relaxation
  double dfreq = 0.0; // Change in frequency
  double dt = 0.0;   
  double value;

  dt = (update->dt);

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      mxe_temperature[i] += qin*dt;

    }
  }
  
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixHeatFlow::memory_usage()
{
  double bytes = 0.0;
  return bytes;
}
