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

#ifndef GRAD_H
#define GRAD_H

#define GRAD_SCHWARZ_CUTOFF 1.0E-10
void printForce(const struct Molecule_t * mol,       // molecule info
                const double *Fx,                    // force in x direction 
                const double *Fy,                    // force in y direction
                const double *Fz);                   // force in z direction

void uhf_force(
	int nBasis,              // number of basis functions
	struct GTOBasis_t * gto, // pointer to function structure
	struct Molecule_t * mol, // pointer to molecule structure
	int nEA,                 // total number of spin up electrons
	int nEB,                 // total number of spin down electrons
	double *CA,              // molecular alpha spin orbital
	double *CB,              // molecular beta spin orbital 
	double *eA,              // eigen values
	double *eB,              // eigen values
	struct option_t *opt,    // global option
	double *Fx,              // returned force in x direction
	double *Fy,              // returned force in y direction
	double *Fz);             // returned force in z direction

void GradEri_ShellSet_Parallel(
	int childID,                    // this child id
	int nCPU,                       // number of CPUs
	int nBasis,                     // number of basis function
	const struct GTOBasis_t *gto,   // basis function database
	const int *basis2Atom,          // basis to atom mapping
	const double *PT,               // total density matrix
	const double *PS,               // spin density matrix
	const double Kfactor,           // exchange weight [0,1]
	double *Gx,                     // returned gradient in x direction
	double *Gy,                     // returned gradient in y direction
	double *Gz);                    // returned gradient in z direction
#endif
