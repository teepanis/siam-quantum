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

#ifndef UHF_H
#define UHF_H
#include <math.h>
#include "mol.h"
#include "basis.h"
#include "option.h"

void uhf_getDMatrix(int nBasis, int nOcc, double *C, double *P);

void uhf_getGMatrix(
	int nBasis,
	const struct GTOBasis_t *gto,
	const double *Schwarz,
	double cutoff,
	const double *PA, const double *PB,
	double *JT, double *KA, double *KB,
	struct option_t *opt);

double uhf_rho(int nBasis,              // number of basis function
               struct GTOBasis_t *gto,  // function structure
               double *PA, double *PB,  // density matrix
               double x, double y, double z);

void uhf_rho_XYPlane(
	int nBasis,                         // number of basis function
    const struct GTOBasis_t *gto,       // basis function structure
    const double *PA, const double *PB, // density matrix
	double cutoff,                      // ignore prefactor
	double x0, double y0, double z0,    // origin of the grid
	double dx, double dy,               // step size
	int nx, int ny,                     // number of points
	double *rhoXY);                     // returned matrix of size ny*nx

void uhf_potential_XYPlane(
	int nBasis,                         // number of basis function
	const struct GTOBasis_t *gto,       // basis function structure
	const struct Molecule_t *mol,       // molecule structure
	const double *PA, const double *PB, // density matrix
	double cutoff,                      // ignore prefactor
	double x0, double y0, double z0,    // origin of the grid
	double dx, double dy,               // step size
	int nx, int ny,                     // number of points
	int inc_nuclei,                     // flag to include nuclei contribution
	double *phiXY);                     // returned matrix of size ny*nx

double uhf_mo(int nBasis,              // number of basis function
              struct GTOBasis_t *gto,  // function structure
              double *C,               // molecular orbital
              int n,                   // orbital index (1st orbital is zero)
              double x, double y, double z);

double uhf_potential(int nBasis,                    // number of basis function
                     const struct GTOBasis_t *gto,  // function structure
                     const struct Molecule_t *mol,  // molecular structure info
                     const double *PA,              // alpha density matrix
                     const double *PB,              // beta density matrix
                     int inc_nuclei,                // flag to include nuclei contribution
                     double x, double y, double z);

void uhf_efield(int nBasis,                          // number of basis function
                const struct GTOBasis_t *gto,        // function structure
                const struct Molecule_t *mol,        // molecular structure info
                const double *PA,                    // alpha density matrix
                const double *PB,                    // beta density matrix
                double x, double y, double z,        // specified point
                double *ex, double *ey, double *ez); // return field

double uhf_getValue(
const int type,                                  // type of value to return
const int nBasis,                                // number of basis function
const struct option_t *opt,                      // global options
const struct GTOBasis_t *gto,                    // function structure
const struct Molecule_t *mol,                    // molecule structure
const double *PA, const double *PB,              // density matrix
const double x, const double y, const double z); // coordinate to get the value

double uhf(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nEA,                 // total number of spin up electrons
	int nEB,                 // total number of spin down electrons
	double *CA,              // returned molecular alpha spin orbital
	double *CB,              // returned molecular beta spin orbital 
	double *eA,              // returned eigen values
	double *eB,              // returned eigen values
	struct option_t *opt);   // global option

#endif
