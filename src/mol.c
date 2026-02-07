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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mol.h"
#include "util.h"

// get_nElectron : return total number of electron in the 
// molecule. This already takes into account of the total
// molecular charge Q.
//
// Mar 1, 2008 - Teepanis Chachiyo
//     Initial implementation
//
int get_nElectron(struct Molecule_t *mol){
	int i;
	int Q=0;
	for(i=0; i < mol->nAtom; i++)
		Q += mol->Z[i];
	return Q - mol->Q;
}

// get_nEA and get_nEB : return total number of spin up electron in the 
// molecule. This already takes into account of the total
// molecular charge Q.
//
// July 12, 2010 - Teepanis Chachiyo
//     Initial implementation
//
int get_nEA(struct Molecule_t *mol, int M){
	int N;
	N = get_nElectron(mol);
	if( (N+M-1)%2 == 1){
		printf("get_nEA - error number of electron does not match multiplicity\n");
		exit(-1);
	}
	return (N+M-1)/2;
}
int get_nEB(struct Molecule_t *mol, int M){
	int N;
	N = get_nElectron(mol);
	if((N-M+1)%2 == 1){
		printf("get_nEB - error number of electron does not match multiplicity\n");
		exit(-1);
	}
	if((N-M+1) < 0){
		printf("get_nEB - error multiplicity cannot larger than number of electron\n");
		exit(-1);
	}
	return (N-M+1)/2;
}


// readMolecule_XYZ : read molecule in XYZ format and 
// return Molecule_t structure. In this case, we are
// using Cartesian basis a=x b=y and c=z
//
// Feb 29, 2008 - Teepanis Chachiyo
//     Initial implementation
//
// Mar  1, 2008 - Teepanis Chachiyo
//     Setting default total molecular charge mol->Q
//
// Dec 13, 2008 - Teepanis Chachiyo
//	   Do not use reduced coordinate
//
// Dec 31, 2009 - Teepanis Chachiyo
//     Potential bug fix: reading comment lines after the number of atoms
//
// Jan 08, 2010 - Teepanis Chachiyo
//     Fix bug - xyz format now does not have to have comments (blank line)
//
// Oct 5, 2012 - Teepanis Chachiyo
//     Adjust the print so that it is identical to printMolecule_XYZ
//
struct Molecule_t *readMolecule_XYZ(FILE *inFile){
	int i;
	int nItem;
	struct Molecule_t *mol;
	char str[256];

	// allocate memory
	mol = (struct Molecule_t *) calloc(1, sizeof(struct Molecule_t));
	if(mol == NULL){
		printf("readMolecule_XYZ: Error - Cannot allocate memory\n");
		exit(-1);
	}

	// reading number of atoms in unit cell
	if(fgets(str, 256, inFile)!=str){
		printf("readMolecule_XYZ: Error - Reading number of atoms\n");
		exit(EXIT_FAILURE);
	}
	nItem = sscanf(str, "%d", &(mol->nAtom));
	if(nItem != 1){
		printf("readMolecule_XYZ: Error - Reading number of atoms\n");
		exit(-1);
	}

	printf("Reading molecule in XYZ file format\n");
	printf("Detected %d atoms\n", mol->nAtom);

	// allocate atomic number and coordinates
	mol->Z = calloc(mol->nAtom, sizeof(int));
	if(mol->Z == NULL){
		printf("readMolecule_XYZ: Error - Cannot allocate atomic number\n");
		exit(EXIT_FAILURE);
	}
	mol->x = calloc(mol->nAtom, sizeof(double));
	mol->y = calloc(mol->nAtom, sizeof(double));
	mol->z = calloc(mol->nAtom, sizeof(double));
	if(mol->x == NULL || mol->y == NULL || mol->z == NULL){
		printf("readMolecule_XYZ: Error - Cannot allocate coordinates\n");
		exit(EXIT_FAILURE);
	}

	// reading comments to buffer
	if(fgets(str, 256, inFile)!=str){
		printf("readMolecule_XYZ: Error - Reading comment from XYZ file\n");
		exit(EXIT_FAILURE);
	}

	// reading atomic number and coordinates
	for(i=0; i < mol->nAtom; i++){
		nItem = fscanf(inFile, "%s %lf %lf %lf", str,
		                                         mol->x+i,
		                                         mol->y+i,
		                                         mol->z+i);
		if(nItem != 4){
			printf("readMolecule_XYZ: Error - Reading coordinate\n");
			exit(EXIT_FAILURE);
		}
		// parse atom name in short format
		mol->Z[i] = sym2Z(str,SYMB_SHORTNAME);

		// convert angstrom to bohr
		mol->x[i] = (mol->x[i]) * ANGSTROM2BOHR;
		mol->y[i] = (mol->y[i]) * ANGSTROM2BOHR;
		mol->z[i] = (mol->z[i]) * ANGSTROM2BOHR;
	}

	// print molecular structure
	printMolecule_XYZ(mol, stdout);

	// set the total molecular charge to zero by default
	mol->Q = 0;

	return mol;

}


// cleanMolecule : cleans up memory allocated by
// readMolecule function. It returns null on sucess.
//
// Feb 10, 2008 - Teepanis Chachiyo
//     Initial implementation
//
struct Molecule_t *cleanMolecule(
	struct Molecule_t *mol){  // pointer to molecule structure

	// clear atom info
	if(mol->Z)    free(mol->Z);
	if(mol->x)    free(mol->x);
	if(mol->y)    free(mol->y);
	if(mol->z)    free(mol->z);

	// clear the entire structure
	if(mol)       free(mol);
	
	return NULL;
}


// nuclei_coulomb : compute nuclei coulomb repulsion energy
//
// Dec 13, 2008 - Teepanis Chachiyo
// 		Initial implementation
//
double nuclei_coulomb(const struct Molecule_t *mol){
	double U;
	int i,j;

	U = 0.0;
	for(i=0; i < mol->nAtom; i++)
	for(j=0; j < i; j++)
		U += mol->Z[i]*mol->Z[j] / 
		     sqrt((mol->x[i] - mol->x[j])*(mol->x[i] - mol->x[j]) +
		          (mol->y[i] - mol->y[j])*(mol->y[i] - mol->y[j]) +
		          (mol->z[i] - mol->z[j])*(mol->z[i] - mol->z[j]));
	return U;
}


// printMolecule_XYZ print molecule in XYZ format
// 
// Oct 22, 2010 - Teepanis Chachiyo
//   Initial implementation
//
void printMolecule_XYZ(
	const struct Molecule_t * mol,// molecule info
	FILE *fd){                    // output file pointer

	int A;        // atom index;
	char str[64]; // atom name in short format

	fprintf(fd,
"-------------------------------------------------------------\n"
"                     Atom Coordinates (Angstroms)            \n"
"     Atom         X              Y             Z             \n"
"-------------------------------------------------------------\n");
#define CAPVALUE(a) (fabs(a)<1.0E-6)?0.0:a
	for(A=0; A < mol->nAtom; A++){
		Z2SymShort(mol->Z[A],str);
		fprintf(fd, "   %5s %15.6f%15.6f%15.6f\n",
		        str,
		        CAPVALUE(mol->x[A])*BOHR2ANGSTROM,
		        CAPVALUE(mol->y[A])*BOHR2ANGSTROM,
		        CAPVALUE(mol->z[A])*BOHR2ANGSTROM);
	}
	fprintf(fd,
"-------------------------------------------------------------\n");
#undef CAPVALUE
}
