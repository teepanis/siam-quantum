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

#ifndef BASIS_H
#define BASIS_H
#include "mol.h"

struct GTOBasis_t{
	int           nContract; // number of contracted functions
	int           l;         // x-coordinate angular index
	int           m;         // y-coordinate angular index
	int           n;         // z-coordinate angular index
	double        x0;        // x-coordinate center
	double        y0;        // y-coordinate center
	double        z0;        // z-coordinate center
	double        *coef;     // array for coeffients
	double        *exp;      // array for exponent
	double        *norm;     // array for normalization factor
};


#define MAX_BASIS 10000
#define BASIS_NAME 256
// GTOBasisSet_t is the structure that hold a primitive
// gaussian basis set, meaning it does not recognize the
// concept of shell. Therefore, building this structure
// from, says, a GAMESS basis set format, you need to 
// construct the px,py,pz yourself.
struct GTOBasisSet_t{
	char          name[BASIS_NAME];
	int           Z;               // atomic number it belongs to
	int           nContract;       // number of contracted functions
	int           l;               // x-coordinate angular index
	int           m;               // y-coordinate angular index
	int           n;               // z-coordinate angular index
	double        *exp;            // array for exponent
	double        *coef;           // array for coeffients
};

int printGTOBasis(FILE *outFile,                // pointer to file
                  int nBasis,                   // number of basis
                  struct GTOBasis_t *gtoBasis); // pointer to array

struct GTOBasis_t *readGTOInput(FILE *inFile,    // pointer to file
                                int  *nBasis);   // number of basis

struct GTOBasis_t *cleanGTOBasis(struct GTOBasis_t *gtoBasis, int nBasis);

struct GTOBasisSet_t *read_GAMESS_BasisSet(
	FILE *inFile,     // input file pointer
	char *name,       // name to assign
	int *nBasisSet);  // return also number of basis set read

void print_BasisSet(int nSet, struct GTOBasisSet_t *dbSet);
void print_BasisSetCode(int nSet, struct GTOBasisSet_t *dbSet);

struct GTOBasis_t * genBasis(
	struct Molecule_t *mol,               // pointer to molecular structure
	int *nBasis,                          // return the number of basis created
	int dbSize,                           // number of record in basis set database
	const struct GTOBasisSet_t *basisDB); // pointer to basis set database

int get_nPrim(int nBasis, const struct GTOBasis_t *gto);

#endif
