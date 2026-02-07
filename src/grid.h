/*
	Siam Quantum [SQ] is a program package that performs quantum
	mechanical calculations on atomic systems.

	Copyright (C) 2015 Teepanis Chachiyo <teepanisc@nu.ac.th>,
	The Institute for Fundamental Study (IF),
	Naresuan University, Thailand.

	This file is a part of SQ.
	                                                       
	This program is free software; you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation; either version 2, or (at your option) 
	any later version.                                                  
*/
#ifndef GRID_H
#define GRID_H

#include "util.h"
#include "basis.h"

// quadrature grid for molecule
struct MolGrid_t{
	int nPoint;   // the number of grid points
	double *w;    // weight
	double *x;    // sampling points
	double *y;    // sampling points
	double *z;    // sampling points
};

// basis function grid point
struct BGrid_t{
	int nBasis;    // number of basis
	long nPoint;   // number of points
	double *maxr2; // maximum r2
	double *val;   // basis function value
	double *gx;    // x-direction gradient
	double *gy;    // y-direction gradient
	double *gz;    // z-direction gradient
	double *x;     // x-coordinate center
	double *y;     // y-coordinate center
	double *z;     // z-coordinate center
};

// radial/angular 21/50
#define GRIDSIZE_S_NRADIAL  21  // enough for milli-hartree or chmical accuracy
#define GRIDSIZE_S_NLABEDEV 4   // but unstable for geometry optimization

// radial/angular 47/146
#define GRIDSIZE_M_NRADIAL  47  // default for energy: enough for a few micro-hartree,
#define GRIDSIZE_M_NLABEDEV 8   // may not be good enough for geometry optmization

// radial/angular 75/302
#define GRIDSIZE_L_NRADIAL  75  // default for optimization: sub micro-hartree range,
#define GRIDSIZE_L_NLABEDEV 13  // good for testing numerical accuracies 

// radial/angular 99/590
// GAMESS uses 99/590 as default for meta-GGA 
#define GRIDSIZE_XL_NRADIAL 99  // probably too large, a good way to test cpu speed
#define GRIDSIZE_XL_NLABEDEV 19 // or parallelization performance

// radial/angular 155/974
#define GRIDSIZE_XXL_NRADIAL 155  // for ultimate testing only, just to be sure
#define GRIDSIZE_XXL_NLABEDEV 25  // this is the same as JANS keywords for GAMESS

#define BALL_RADIUS   2.5584    // parameter for opt->whichGrid = GRID_BALLS
// 2.5584
//   Fitted to CO geometry from correlation paper 
//   6-31G* /BeckeGrid/Chachiyo-XC
//   NRADIAL 47 NLABEDEV 8
//   on the condition that the forces 
//   are equal in magnitude
//
//   Since BeckeGrid depends on R[A]. This means we
//   we treat all atoms equally.
//
//   With NRADIAL 47 NLABEDEV 8, error is about 1e-5 
//   compared to the fine grid 75/13
//
//   With NRADIAL 47 NLABEDEV 8 + RA[A] = get_bragg...,
//   the forces were not equal and geometry
//   optimization failed. The H2O optmization
//   failed also.
//
// 2.42
//   Fitted to LiH instead of CO
//
// Teepanis Chachiyo 15 Aug, 2019


//
// NRADIAL 90, NLABEDEV 8 + Becke + Bragg radius 
//   - caused the forces on CO not
//     being equal in magnitude
//
//   - the force on oxygen in H2O is
//     unstable and OPT do not converge
//
//   - the problem goes away if we use
//     constant atom radius
//

struct MolGrid_t * genMolGrid(
const int whichGrid,          // grid type in opt->whichGrid
const int nradial,            // number of radial grid point
const int nrule,              // order of angular grid
const int nBasis,             // number of basis
const struct GTOBasis_t *gto, // basis info
const struct Molecule_t *mol, // molecule info
const double cutoff);         // cut-off

struct MolGrid_t * cleanMolGrid(struct MolGrid_t *grid);

int LabedevP(int rule);
void LabedevQ(int npoint, double *x, double *y, double *z, double *w);

#endif
