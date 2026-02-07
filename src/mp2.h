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

#ifndef MP2_H
#define MP2_H

#include <stdlib.h>
#include <stdio.h>

#include "basis.h"
#include "option.h"
#include "mol.h"

double mp2_rhf_aqbj(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_rhf_aqij(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_rhf_aqij_Parallel_Contribute(
	int a,                         // starting orbital index
	int maxCorr,                   // number of correlated orbitals
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const double *Schwarz,         // schwarz inequality
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_rhf_aqij_Parallel(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_rhf_direct(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_rhf_semi_direct_aqbj(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_rhf_semi_direct_aqij(
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_uhf_direct(
	int nBasis,                    // number of basis function
	int nEA,                       // number of occupied spin up orbitals
	int nEB,                       // number of occupied spin dn orbitals
	const double *eA,              // eigen values
	const double *eB,              // eigen values
	const double *CA,              // molecular orbitals
	const double *CB,              // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_uhf_semi_direct_aqij(
	int nBasis,                    // number of basis function
	int nEA,                       // number of occupied spin up orbitals
	int nEB,                       // number of occupied spin dn orbitals
	const double *eA,              // eigen values
	const double *eB,              // eigen values
	const double *CA,              // molecular orbitals
	const double *CB,              // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

double mp2_uhf_aqij(
	int nBasis,                    // number of basis function
	int nEA,                       // number of occupied spin up orbitals
	int nEB,                       // number of occupied spin dn orbitals
	const double *eA,              // eigen values
	const double *eB,              // eigen values
	const double *CA,              // molecular orbitals
	const double *CB,              // molecular orbitals
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt);   // global options

#endif
