/*
	Siam Quantum [SQ] is a program pacakage that performs quantum
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
#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include "basis.h"
#include "mol.h"
#include "option.h"

// convergence criterion
//#define CONV_FORCEMAX   0.000450
//#define CONV_FORCERMS   0.000300
//#define CONV_DISPMAX    0.001800
//#define CONV_DISPRMS    0.001200
#define OPT_CONV_FORCEMAX   0.000100
#define OPT_CONV_FORCERMS   0.000080
#define OPT_CONV_DISPMAX    0.000800
#define OPT_CONV_DISPRMS    0.000400
// Setting FORCEMAX 0.00450 and RMSMAX=0.00300 causes 
// OPT C2H6/6-31G* giving incorrect lowest vibation frequency.
// (319 versus the correct value of 326).
// Setting the criterian to .0002 and 0.0001 fixes this issue.
//
// Teepanis Chachiyo - 28 March, 2020

#define MAXSTEPSIZE     0.5

void invHessian_BFGS(
	int nDim,            // matrix dimension
	const double *dR,    // displacement vector
	const double *dGrad, // change of gradient vector,
	double *invHessian); // original inverse of the Hessian

void stepVector_Newton(
	int nDim,            // dimension of the matrix
	const double *invH,  // pointer to inverse of the Hessian matrix
	const double *grad,  // pointer to gradient vector
	double maxStepSize,  // maximum step size in each direction
	double *dR);         // returned step vector

void deleteTranslation(
	int nDim,     // vector dimension
	double *dR); // displacement 

void optimize(
	int dbSize,                         // number of record in basis set db
	const struct GTOBasisSet_t *basisDB,// pointer to basis set database
	struct Molecule_t *mol,             // returned molecular coordinate
	struct option_t *opt);              // options

#endif
