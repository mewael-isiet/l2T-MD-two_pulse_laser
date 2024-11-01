/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(heatflow,FixHeatFlow)

#else

#ifndef LMP_FIX_HEAT_FLOW_H
#define LMP_FIX_HEAT_FLOW_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHeatFlow : public Fix {
 public:
  FixHeatFlow(class LAMMPS *, int, char **);
  virtual ~FixHeatFlow();
  int setmask();
  virtual void init();
  virtual void initial_integrate(int);
  //virtual void post_force(int);
  double memory_usage();

 protected:
  double qin, xvalue,yvalue,zvalue;

};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix setforce does not exist

Self-explanatory.

E: Variable name for fix setforce does not exist

Self-explanatory.

E: Variable for fix setforce is invalid style

Only equal-style variables can be used.

E: Cannot use non-zero forces in an energy minimization

Fix setforce cannot be used in this manner.  Use fix addforce
instead.

*/
