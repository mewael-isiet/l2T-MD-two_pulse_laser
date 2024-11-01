/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
   l2T-MD implementation: M. Ponga (UBC) 
------------------------------------------------------------------------- */

#include <math.h>
#include <string.h>
#include "compute_ee.h"
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "modify.h"
#include "fix.h"
#include "memory.h"
#include "error.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeEE::ComputeEE(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg),
  energy(NULL)
{
  //compute ee all ee Ge 1.5e7 V 100.0
  if (narg != 7) error->all(FLERR,"Illegal compute ee command");

  if (strcmp(arg[3],"lambda") == 0) lambda = force->numeric(FLERR,arg[4]);

  if (strcmp(arg[5],"V") == 0) V = force->numeric(FLERR,arg[6]);

  scalar_flag = 1;
  size_peratom_cols = 0;
  peatomflag = 1;
  timeflag = 1;
  comm_reverse = 1;
  atom_style = atom->atom_style;
  
  char *atom_style = atom->atom_style;
  if ( strcmp(atom_style,"atomic_mxe")!= 0 ) 
      error->all(FLERR,"Illegal atom style for fe command");

  nmax = 0;  
  
}

/* ---------------------------------------------------------------------- */

ComputeEE::~ComputeEE()
{
  memory->destroy(energy);
}

/* ---------------------------------------------------------------------- */

double ComputeEE::compute_scalar()
{
  int i;

  invoked_peratom = update->ntimestep;
  if (update->eflag_atom != invoked_peratom)
    error->all(FLERR,"Electronic energy was not tallied on needed timestep");

  // grow local energy array if necessary
  // needs to be atom->nmax in length

  if (atom->nmax > nmax) {
    memory->destroy(energy);
    nmax = atom->nmax;
    memory->create(energy,nmax,"ee:energy");
    vector_atom = energy;
  }

  int nlocal = atom->nlocal;

  double elec_energy=0.0, EE = 0.0;
  double *temp = atom->mxe_temperature;
  int *mask = atom->mask;
  bigint natoms = atom->natoms;
 
  //Compute electronic energy due to electronic temperature ee[i]=C_eT_i^e*V_{atom}
  //if (!(mask[i] & groupbit)) 
  	for (i = 0; i < nlocal; i++) elec_energy += (temp[i]*temp[i]);

  // Reduce
  MPI_Allreduce(&elec_energy,&scalar,1,MPI_DOUBLE,MPI_SUM,world);

  scalar *= lambda*V;

  return scalar;
}

/* ---------------------------------------------------------------------- */

int ComputeEE::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) buf[m++] = energy[i];
  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeEE::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    energy[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeEE::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
