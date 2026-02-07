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
#include <string.h>

#include "int.h"
#include "util.h"
#include "basis.h"


// printGTOBasis : prints out basis function information
// to file which can be standard output. The printing 
// format is the same as in the file sqmc.py 
// class: GTOBasis.printInfo()
//
// It can also read this format thru function readGTOInput()
//
// Feb 22, 2008 - Teepanis Chachiyo
//     Initial implementation.
//
int printGTOBasis(FILE *outFile,                // pointer to file
                  int nBasis,                   // number of basis
                  struct GTOBasis_t *gtoBasis){ // pointer to array

	int i;     // basis function loop
	int j;     // contracted function loop
	int nItem; // number of item pushed to file

	fprintf(outFile, "nBasis: %6d\n", nBasis);
	// loop over basis functions
	for(i=0; i < nBasis; i++){

		nItem = fprintf(outFile,
		                "%2d%2d%2d %15.8E %15.8E %15.8E %6d\n",
		                gtoBasis[i].l,  gtoBasis[i].m,  gtoBasis[i].n,
		                gtoBasis[i].x0, gtoBasis[i].y0, gtoBasis[i].z0,
		                gtoBasis[i].nContract);

		if(nItem < 0){
			printf("printGTOBasis: Error - Printing basis function\n");
			exit(-1);
		}

		for(j=0; j < gtoBasis[i].nContract; j++){

			nItem = fprintf(outFile,
			                "                       %15.8E %15.8E\n",
			                gtoBasis[i].exp[j], gtoBasis[i].coef[j]);

			if(nItem < 0){
				printf("printGTOBasis: Error - "
				       "Print contracted function\n");
				exit(-1);
			}
		}
	}

	// return 1 on success
	return 1;
}


// readGTOInput : retrives GTO basis information
// thru file descriptor inFile. It then constructs
// GTODatabase for further usage.
//
// Feb 22, 2008 - Teepanis Chachiyo
//     Initial implementation
//
struct GTOBasis_t *readGTOInput(FILE *inFile,    // pointer to file
                                int  *nBasis){   // number of basis

	int i;      // basis function loop
	int j;      // contracted function loop
	int nItem;  // number of item read from file

	struct GTOBasis_t *gtoBasis; // pointer to array
	double norm;                 // normalization constant

	// read number of basis function
	nItem = fscanf(inFile, "%*s %d", nBasis);
	if(nItem != 1){
		printf("readGTOInput: Error - Reading number of basis\n");
		exit(-1);
	}
	if(*nBasis < 1){
		printf("readGTOInput: Error - nBasis < 1.\n");
		exit(-1);
	}

	// allocate memory
	gtoBasis = (struct GTOBasis_t *) 
	           calloc( *nBasis, sizeof(struct GTOBasis_t));

	if(inFile == NULL){
		printf("readGTOInput: Error - inFile == NULL\n");
		exit(-1);
	}
	if(gtoBasis == NULL){
		printf("readGTOInput: Error - gtoBasis == NULL\n");
		exit(-1);
	}

	// loop to read basis information
	for(i=0; i < *nBasis; i++){

		// reading angular index
		nItem = fscanf(inFile, "%d %d %d", &(gtoBasis[i].l),
		                                   &(gtoBasis[i].m),
		                                   &(gtoBasis[i].n));
		if(nItem != 3){
			printf("readGTOInput: Error - Reading GTO angular\n");
			exit(-1);
		}

		// reading center in reduced coordinate
		nItem = fscanf(inFile, "%lf %lf %lf", &(gtoBasis[i].x0),
		                                      &(gtoBasis[i].y0),
		                                      &(gtoBasis[i].z0));
		if(nItem != 3){
			printf("readGTOInput: Error - Reading GTO center");
			exit(-1);
		}


		// reading number of contracted functions
		nItem = fscanf(inFile, "%d", &(gtoBasis[i].nContract));
		if(nItem != 1){
			printf("readGTOInput: Error - Reading GTO nContract\n");
			exit(-1);
		}

		// allocate memory for contracted functions
		gtoBasis[i].coef  = calloc(gtoBasis[i].nContract,
		                           sizeof(double));
		gtoBasis[i].exp   = calloc(gtoBasis[i].nContract,
		                           sizeof(double));
		gtoBasis[i].norm  = calloc(gtoBasis[i].nContract,
		                           sizeof(double));

		if(gtoBasis[i].coef  == NULL || 
		   gtoBasis[i].exp   == NULL || 
		   gtoBasis[i].norm  == NULL){
			printf("readGTOInput: Error - Cannot allocate memory\n");
			exit(-1);
		}

		// loop to read contracted functions
		for(j=0; j < gtoBasis[i].nContract; j++){
			nItem = fscanf(inFile, "%lf %lf", &(gtoBasis[i].exp[j]),
			                                  &(gtoBasis[i].coef[j]));
			if(nItem != 2){
				printf("readGTOInput: Error - "
				       "Reading contracted functions");
				exit(-1);
			}
		}
	}

	//
	// normalize contracted coefficient
	//
	// loop over basis set
	for(i=0; i < *nBasis; i++){
		// loop over contracted functions
		for(j=0; j < gtoBasis[i].nContract; j++){
			norm = overlap(gtoBasis[i].exp[j],
			               gtoBasis[i].l,
			               gtoBasis[i].m,
			               gtoBasis[i].n,
			               gtoBasis[i].x0,
			               gtoBasis[i].y0,
			               gtoBasis[i].z0,

			               gtoBasis[i].exp[j],
			               gtoBasis[i].l,
			               gtoBasis[i].m,
			               gtoBasis[i].n,
			               gtoBasis[i].x0,
			               gtoBasis[i].y0,
			               gtoBasis[i].z0);

			gtoBasis[i].norm[j] = 1.0/sqrt(norm);
		}
	}

	return gtoBasis;
}


//
// cleanGTOBasis : cleans up memory allocation by
// readGTOBasis function. It returns null on success.
//
// Feb 22, 2008 - Teepanis Chachiyo
//     Initial implementation
//
struct GTOBasis_t *cleanGTOBasis(struct GTOBasis_t *gtoBasis, int nBasis){

	int i;
	
	// loop over the number of basis function
	for(i=0; i < nBasis; i++){
		if(gtoBasis[i].exp)   free(gtoBasis[i].exp);
		if(gtoBasis[i].coef)  free(gtoBasis[i].coef);
		if(gtoBasis[i].norm)  free(gtoBasis[i].norm);
	}
	
	// free the entire structure
	if(gtoBasis) free(gtoBasis);

	return NULL;
}


// print_BasisSet : print out basis set informatin
// We are not using L shell here, only S and P separatedly.
//
// Feb 29, 2008 - Teepanis Chachiyo
//     Initial implementaion
//
void print_BasisSet(int nSet, struct GTOBasisSet_t *dbSet){
	int i,j;
	char str[256];

	for(i=0; i < nSet; i++){

		Z2Sym(dbSet[i].Z, str);
		printf("%s %d %d %d\n", str, dbSet[i].l, dbSet[i].m, dbSet[i].n);
		for(j=0; j < dbSet[i].nContract; j++){
			printf("    %20.8lf %20.8lf\n", dbSet[i].exp[j],
			                                dbSet[i].coef[j]);
		}
	}
}

// print_BasisSetCode : print out C source code for including this
// basis set into the program.
//
// Feb 29, 2008 - Teepais Chachiyo
//     Initial implementation
//
void print_BasisSetCode(int nSet, struct GTOBasisSet_t *dbSet){
	int i,j;

	for(i=0; i < nSet; i++){

		// exponent
		printf("static double BASIS_%s_EXP_%d[] = {",dbSet[i].name,i);
		for(j=0; j < dbSet[i].nContract; j++){
			printf("%.8lf", dbSet[i].exp[j]);
			if(j != dbSet[i].nContract-1)
				printf(",");
		}
		printf("};\n");

		// coefficient
		printf("static double BASIS_%s_COEF_%d[] = {",dbSet[i].name,i);
		for(j=0; j < dbSet[i].nContract; j++){
			printf("%.8lf", dbSet[i].coef[j]);
			if(j != dbSet[i].nContract-1)
				printf(",");
		}
		printf("};\n");
	}

	// structure
	printf("struct GTOBasisSet_t  BASIS_DB[%d] = {\n",nSet);
	for(i=0; i < nSet; i++){
		printf("{\"%s\",%d,%d,%d,%d,%d,BASIS_%s_EXP_%d,BASIS_%s_COEF_%d}",
		    dbSet[i].name,
		    dbSet[i].Z,
		    dbSet[i].nContract,
		    dbSet[i].l,
		    dbSet[i].m,
		    dbSet[i].n,
		    dbSet[i].name,i,
		    dbSet[i].name,i);
		if(i!=nSet-1)
			printf(",\n");
	}
	printf("};\n");
}


// read_GAMESS_BasisSet : read basis set information
// in GAMESS US format and build a structure containing
// the basis set information.
//
// Feb 25, 2008 - Teepanis Chachiyo
//     Initial implementation.
//
// Feb 28, 2008 - Teepanis Chachiyo
//     Finalize
//
// Dec 24, 2014 - Teepanis Chachiyo
//		Handle character 'D' in scientific format such as 1.0D-2
//
#define MAX_BASIS_PER_FILE 20000
struct GTOBasisSet_t *read_GAMESS_BasisSet(
	FILE *inFile,     // input file pointer
	char *name,       // name to assign
	int *nBasisSet){  // return also number of basis set read

	int i;
	int Z=0;
	int nContract;

	// problem is that when reading file we do not know how
	// many basis set is in there before hand.
	struct GTOBasisSet_t * basisDB;     // at returning stage
	int nDBSet=0;                       // number of db items
	char keyword[1024];                 // current keyword
	double bas_exp[1024];               // buffer for exponent
	double bas_coef_a[1024];            // buffer for coefficient
	double bas_coef_b[1024];            // 2nd buffer if needed

	// allocate memeory at maximum
	basisDB = calloc(MAX_BASIS_PER_FILE,sizeof(struct GTOBasisSet_t));

	// search for word "$DATA"
	if(findf(inFile,1,"$DATA")==EOF){
		printf("read_GAMESS_BasisSet - "
		       "Error: Cannot find $DATA keyword\n");
		       exit(EXIT_FAILURE);
	}

	do{
		// read keyword
		if(fscanf(inFile,"%s",keyword)!=1){
			printf("read_GAMESS_BasisSet - "
			       "Error: File corrupted\n");
			exit(EXIT_FAILURE);
		}


		// keyword could be atom name
		if(sym2Z(keyword, SYMB_LONGNAME) > 0){
			Z = sym2Z(keyword, SYMB_LONGNAME);
			continue;
		}

		// keyword could be orbital index
		if(strcmp(keyword,"S")==0 ||
		   strcmp(keyword,"L")==0 ||
		   strcmp(keyword,"P")==0 ||
		   strcmp(keyword,"D")==0 ||
		   strcmp(keyword,"F")==0){
		
			// read number of contracted function
			if(fscanf(inFile,"%d",&nContract)!=1){
				printf("read_GAMESS_BasisSet - "
				"Error: No contracted function number\n");
				exit(EXIT_FAILURE);
			}


			// read in
			for(i=0; i < nContract; i++){

				// special for L: both coefficients at the same time
				if(strcmp(keyword,"L")==0){
					// read string into buffer for pre-processing
					char str1[1024], str2[1024], str3[1024];
					if(fscanf(inFile,"%*d %s %s %s",str1,str2,str3)!=3){
						printf("read_GAMESS_BasisSet - "
					               "Error: Reading coefficient\n");
						exit(EXIT_FAILURE);
					}
					// convert character D --> E if exists
					int n;
					for(n=0; n<strlen(str1); n++) if(str1[n]=='D') str1[n]='E';
					for(n=0; n<strlen(str2); n++) if(str2[n]=='D') str2[n]='E';
					for(n=0; n<strlen(str3); n++) if(str3[n]=='D') str3[n]='E';
					// convert string to double
					bas_exp   [i] = strtod(str1,NULL);
					bas_coef_a[i] = strtod(str2,NULL);
					bas_coef_b[i] = strtod(str3,NULL);
					// validate range
					if(bas_exp[i] <= 0.0 || bas_coef_a[i] == 0.0 || bas_coef_b[i] == 0.0){
						printf("read_GAMESS_BasisSet - "
					               "Error: Reading parsing coefficient\n");
						exit(EXIT_FAILURE);
					}

					//if(fscanf(inFile,"%*d %lf %lf %lf",bas_exp+i,
					//                                   bas_coef_a+i,
					//                                   bas_coef_b+i)!=3){
					//	printf("read_GAMESS_BasisSet - "
					//               "Error: Reading coefficient\n");
					//	exit(EXIT_FAILURE);
					//}
				}else{
					// process normal orbitals

					// read string into buffer for pre-processing
					char str1[1024], str2[1024];
					if(fscanf(inFile,"%*d %s %s",str1,str2)!=2){
						printf("read_GAMESS_BasisSet - "
					               "Error: Reading coefficient\n");
						exit(EXIT_FAILURE);
					}
					// convert character D --> E if exists
					int n;
					for(n=0; n<strlen(str1); n++) if(str1[n]=='D') str1[n]='E';
					for(n=0; n<strlen(str2); n++) if(str2[n]=='D') str2[n]='E';
					// convert string to double
					bas_exp   [i] = strtod(str1,NULL);
					bas_coef_a[i] = strtod(str2,NULL);
					// validate range
					if(bas_exp[i] <= 0.0 || bas_coef_a[i] == 0.0){
						printf("read_GAMESS_BasisSet - "
					               "Error: Reading parsing coefficient\n");
						exit(EXIT_FAILURE);
					}

					//if(fscanf(inFile,"%*d %lf %lf",bas_exp+i,
					//                               bas_coef_a+i)!=2){
					//	printf("read_GAMESS_BasisSet - "
					//               "Error: Reading coefficient\n");
					//	exit(EXIT_FAILURE);
					//}
				}
			}

#define ALLOCATE                                                 \
basisDB[nDBSet].nContract = nContract;                           \
basisDB[nDBSet].coef      = calloc(nContract, sizeof(double));   \
basisDB[nDBSet].exp       = calloc(nContract, sizeof(double))

#define SETINFO(li,lj,lk,ee,cc)      \
strcpy(basisDB[nDBSet].name,name);   \
basisDB[nDBSet].Z = Z;               \
basisDB[nDBSet].l = li;              \
basisDB[nDBSet].m = lj;              \
basisDB[nDBSet].n = lk;              \
for(i=0; i < nContract; i++){        \
    basisDB[nDBSet].coef[i] = cc[i]; \
    basisDB[nDBSet].exp[i]  = ee[i]; \
}

			// process orbital
			if(strcmp(keyword,"S")==0){
				ALLOCATE; SETINFO(0,0,0,bas_exp, bas_coef_a); nDBSet++;
			}

			// process p orbital
			if(strcmp(keyword,"L")==0){
				ALLOCATE; SETINFO(0,0,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1,0,0,bas_exp, bas_coef_b); nDBSet++;
				ALLOCATE; SETINFO(0,1,0,bas_exp, bas_coef_b); nDBSet++;
				ALLOCATE; SETINFO(0,0,1,bas_exp, bas_coef_b); nDBSet++;
			}

			// process p orbital
			if(strcmp(keyword,"P")==0){
				ALLOCATE; SETINFO(1,0,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,1,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,0,1,bas_exp, bas_coef_a); nDBSet++;
			}

			// process d orbital
			if(strcmp(keyword,"D")==0){
				ALLOCATE; SETINFO(2,0,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,2,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,0,2,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1,1,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1,0,1,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,1,1,bas_exp, bas_coef_a); nDBSet++;
			}

			if(strcmp(keyword,"F")==0){
				ALLOCATE; SETINFO(3,0,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,3,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,0,3,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(2,1,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1,2,0,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(2,0,1,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1,0,2,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,2,1,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(0,1,2,bas_exp, bas_coef_a); nDBSet++;
				ALLOCATE; SETINFO(1,1,1,bas_exp, bas_coef_a); nDBSet++;
			}

#undef SETINFO
#undef ALLOCATE
			continue;
		}

	}while( strcmp(keyword,"$END")!=0 );

	// truckate memory
	basisDB = realloc(basisDB, sizeof(struct GTOBasisSet_t) * nDBSet);

	// output and return value
	*nBasisSet = nDBSet;
	return basisDB;
}


// genBasis : generate basis set structure from the specified molecular
// structure.
//
// Feb 29, 2008 - Teepanis Chachiyo
//     Initial implementation
//
// Nov 04, 2009 - Theerapon Kumla
//     Check error if atomic number is not found in basis set database
//
// Oct 22, 2010 - Teepanis Chachiyo
//     Change basisDB to constant specifier
//
struct GTOBasis_t * genBasis(
	struct Molecule_t *mol,               // pointer to molecular structure
	int *nBasis,                          // return the number of basis created
	int dbSize,                           // number of record in basis set database
	const struct GTOBasisSet_t *basisDB){ // pointer to basis set database

	int i, b, j, n=0;
	struct GTOBasis_t *gto;
	double norm;

	// allocate memory
	gto = calloc(MAX_BASIS, sizeof(struct GTOBasis_t));
	if(gto==NULL){
		printf("genBasis: Error - Cannot allocate memory\n");
		exit(EXIT_FAILURE);
	}

	// loop thru all atoms
	for(i=0; i < mol->nAtom; i++){
	// make sure we have at least one matching
	for(b=0; b < dbSize; b++) if(basisDB[b].Z == mol->Z[i]) break;
	if(b==dbSize){
		printf("genBasis: Error - Cannot find atom %d in basis file.\n",
		       mol->Z[i]);
		exit(EXIT_FAILURE);
	}
	// find all matched atomic number
	for(b=0; b < dbSize; b++){
	if(basisDB[b].Z == mol->Z[i]){
		// copy values
		gto[n].nContract = basisDB[b].nContract;
		gto[n].l         = basisDB[b].l;
		gto[n].m         = basisDB[b].m;
		gto[n].n         = basisDB[b].n;
		gto[n].x0        = mol->x[i];
		gto[n].y0        = mol->y[i];
		gto[n].z0        = mol->z[i];
		// allocate memory
		gto[n].exp       = calloc(gto[n].nContract, sizeof(double));
		gto[n].coef      = calloc(gto[n].nContract, sizeof(double));
		gto[n].norm      = calloc(gto[n].nContract, sizeof(double));
		if(gto[n].exp ==NULL ||
		   gto[n].coef==NULL ||
		   gto[n].norm==NULL){
			printf("genBasis: Error - Cannot allocate memory\n");
			exit(EXIT_FAILURE);
		}
		// copy values
		for(j=0; j < gto[n].nContract; j++){
			gto[n].exp[j]  = basisDB[b].exp[j];
			gto[n].coef[j] = basisDB[b].coef[j];
		}
		n++;
	}
	}
	}

	// shrink memory
	gto     = realloc(gto, sizeof(struct GTOBasis_t) * n);
	*nBasis = n;

	//
	// normalize contracted coefficient
	//
	// loop over basis set
	for(i=0; i < *nBasis; i++){
		// loop over contracted functions
		for(j=0; j < gto[i].nContract; j++){
			norm = overlap(gto[i].exp[j],
			               gto[i].l,
			               gto[i].m,
			               gto[i].n,
			               gto[i].x0,
			               gto[i].y0,
			               gto[i].z0,

			               gto[i].exp[j],
			               gto[i].l,
			               gto[i].m,
			               gto[i].n,
			               gto[i].x0,
			               gto[i].y0,
			               gto[i].z0);

			gto[i].norm[j] = 1.0/sqrt(norm);
		}
	}

	// report stat
	printf("There are %d Cartesian-Gaussian contracted basis functions\n", n);
	i = get_nPrim(*nBasis, gto);
	printf("There are %d Primitive GTO functions\n", i);

	return gto;
}

// get_nPrim : return the number of primitive function in the
// basis set.
//
// Feb 28, 2008 - Teepanis Chachiyo
//
// May 15, 2011 - Teepanis Chachiyo
//                Change *gto to const modifier
//
int get_nPrim(int nBasis, const struct GTOBasis_t *gto){
	int i;
	int nPrim=0;

	for(i=0; i < nBasis; i++)
		nPrim += gto[i].nContract;
	return nPrim;
}

