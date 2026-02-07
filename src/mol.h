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

#ifndef MOL_H
#define MOL_H
struct Molecule_t{
	int   nAtom; // number of nuclei
	int   Q;     // total molecular charge
	int   *Z;    // atomic number array
	double *x;   // Cartesian x-coordinate
	double *y;   // Cartesian y-coordinate
	double *z;   // Cartesian z-coordinate
};

void printMolecule_XYZ(
	const struct Molecule_t * mol,// molecule info
	FILE *fd);                    // output file pointer
struct Molecule_t *readMolecule_XYZ(FILE *inFile);
struct Molecule_t *cleanMolecule(struct Molecule_t *mol);
int get_nElectron(struct Molecule_t *mol);
double nuclei_coulomb(const struct Molecule_t *mol);
int get_nEA(struct Molecule_t *mol, int M);
int get_nEB(struct Molecule_t *mol, int M);
#endif

