/* -*- c++ -*- ----------------------------------------------------------
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

#ifdef COMPUTE_CLASS

ComputeStyle(ee,ComputeEE)

#else

#ifndef LMP_COMPUTE_EE_H
#define LMP_COMPUTE_EE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeEE : public Compute {
 public:
  ComputeEE(class LAMMPS *, int, char **);
  ~ComputeEE();
  void init() {}
  double compute_scalar();
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int pairflag,bondflag,angleflag,dihedralflag,improperflag;
  int kspaceflag,fixflag;
  int nmax;
  double *energy, lambda, V;
  char *atom_style;
  int mt_flag;
  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Per-atom energy was not tallied on needed timestep

You are using a thermo keyword that requires potentials to
have tallied energy, but they didn't on this timestep.  See the
variable doc page for ideas on how to make this work.

*/
