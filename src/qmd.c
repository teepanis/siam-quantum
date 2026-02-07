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

#include "uhf.h"
#include "grad.h"
#include "util.h"
#include "mol.h"
#include "pop.h"
#include "qmd.h"


// change atomic number to mass in Dalton (or amu)
double Z2mass(int Z){
	double mass;
	switch(Z){
		case 1 : 	mass=1.008; 	break;
		case 2 :	mass=4.003; 	break;
		case 3 : 	mass=6.941; 	break;
		case 4 : 	mass=9.012; 	break;
		case 5 : 	mass=10.81; 	break;
		case 6 : 	mass=12.01; 	break;
		case 7 : 	mass=14.01; 	break;
		case 8 : 	mass=16.00; 	break;
		case 9 : 	mass=19.00; 	break;
		case 10 : 	mass=20.18; 	break;
		case 11 : 	mass=22.99; 	break;
		case 12 : 	mass=24.31; 	break;
		case 13 : 	mass=26.98; 	break;
		case 14 : 	mass=28.09; 	break;
		case 15 : 	mass=30.97; 	break;
		case 16 : 	mass=32.07; 	break;
		case 17 : 	mass=35.45; 	break;
		case 18 : 	mass=39.95; 	break;
		case 19 : 	mass=39.10; 	break;
		case 20 : 	mass=40.08; 	break;
		case 21 : 	mass=44.96; 	break;
		case 22 : 	mass=47.87; 	break;
		case 23 : 	mass=50.94; 	break;
		case 24 : 	mass=52.00; 	break;
		case 25 : 	mass=54.94; 	break;
		case 26 : 	mass=55.84; 	break;
		case 27 : 	mass=58.93; 	break;
		case 28 : 	mass=58.69; 	break;
		case 29 : 	mass=63.55; 	break;
		case 30 : 	mass=65.39; 	break;
		case 31 : 	mass=69.72; 	break;
		case 32 : 	mass=72.61; 	break;
		case 33 : 	mass=74.92; 	break;
		case 34 : 	mass=78.96; 	break;
		case 35 : 	mass=79.90; 	break;
		case 36 : 	mass=83.80; 	break;
		case 37 : 	mass=85.47; 	break;
		case 38 : 	mass=87.62; 	break;
		case 39 : 	mass=88.91; 	break;
		case 40 : 	mass=91.22; 	break;
		case 41 : 	mass=92.91; 	break;
		case 42 : 	mass=95.94; 	break;
		case 43 : 	mass=99.00; 	break;
		case 44 : 	mass=101.07;	break;
		case 45 : 	mass=102.91; 	break;
		case 46 : 	mass=106.42; 	break;
		case 47 : 	mass=107.87; 	break;
		case 48 : 	mass=112.41; 	break;
		case 49 : 	mass=114.82; 	break;
		case 50 : 	mass=118.71; 	break;
		case 51 : 	mass=121.76; 	break;
		case 52 : 	mass=127.60; 	break;
		case 53 : 	mass=126.90; 	break;
		case 54 : 	mass=131.29; 	break;
		case 55 : 	mass=132.91; 	break;
		case 56 : 	mass=137.33; 	break;
		case 57 : 	mass=138.91; 	break;
		case 58 : 	mass=140.12; 	break;
		case 59 : 	mass=140.91; 	break;
		case 60 : 	mass=144.24; 	break;
		case 61 : 	mass=145.00; 	break;
		case 62 : 	mass=150.36; 	break;
		case 63 : 	mass=151.96; 	break;
		case 64 : 	mass=157.25; 	break;
		case 65 : 	mass=158.93; 	break;
		case 66 : 	mass=162.50; 	break;
		case 67 : 	mass=164.93; 	break;
		case 68 : 	mass=167.26; 	break;
		case 69 : 	mass=168.93; 	break;
		case 70 : 	mass=173.04; 	break;
		case 71 : 	mass=174.97; 	break;
		case 72 : 	mass=178.49; 	break;
		case 73 : 	mass=180.95; 	break;
		case 74 : 	mass=183.84; 	break;
		case 75 : 	mass=186.21; 	break;
		case 76 : 	mass=190.23; 	break;
		case 77 : 	mass=192.22; 	break;
		case 78 : 	mass=195.08; 	break;
		case 79 : 	mass=196.97; 	break;
		case 80 : 	mass=200.59; 	break;
		case 81 : 	mass=204.38; 	break;
		case 82 : 	mass=207.2; 	break;
		case 83 : 	mass=208.98; 	break;
		case 84 : 	mass=209; 	 	break;
		case 85 : 	mass=210; 	 	break;
		case 86 : 	mass=222; 	 	break;
		case 87 : 	mass=223; 	 	break;
		case 88 : 	mass=226; 	 	break;
		case 89 : 	mass=227; 	 	break;
		case 90 : 	mass=232.04; 	break;
		case 91 : 	mass=231.04; 	break;
		case 92 : 	mass=238.03; 	break;
		case 93 : 	mass=237; 		break;
		case 94 : 	mass=244; 		break;
		case 95 : 	mass=243; 		break;
		case 96 : 	mass=247; 		break;
		case 97 : 	mass=247; 		break;
		case 98 : 	mass=251; 		break;
		case 99 : 	mass=252; 		break;
		case 100 : 	mass=257; 		break;
		case 101 : 	mass=258; 		break;
		case 102 : 	mass=259; 		break;
		case 103 : 	mass=262; 		break;
		// data from http://en.wikipedia.org/wiki/Atomic_weight/Table
		default: printf("Z2kg not supported atom type\n");
		         exit(-1);
	}
	return mass;
}


double KineticEnergy(
int nAtom, 
double *mass, 
double *vx, double *vy, double *vz){
	int n;
	double K=0.0;

	for(n=0; n < nAtom; n++){
		K += 0.5 * mass[n] * vx[n]*vx[n];
		K += 0.5 * mass[n] * vy[n]*vy[n];
		K += 0.5 * mass[n] * vz[n]*vz[n];
	}
	return K;
}


void EulerCromerUpdate(
	int nAtom,                                // number of atoms
	double dt,                                // time step in AU.
	double *mass,                             // atomic mass in AU.
	double *x,  double *y,  double *z,        // atomic positions
	double *vx, double *vy, double *vz,       // atomic velocities
	double *fx, double *fy, double *fz){      // atomic forces

	int n;

	// update velocity
	for(n=0; n < nAtom; n++){
		vx[n] += dt*fx[n]/mass[n];
		vy[n] += dt*fy[n]/mass[n];
		vz[n] += dt*fz[n]/mass[n];
	}

	// update position
	for(n=0; n < nAtom; n++){
		x[n] += dt*vx[n];
		y[n] += dt*vy[n];
		z[n] += dt*vz[n];
	}
}

//
// July 11, 2017 - Teepanis Chachiyo & Hathaithip Chachiyo
//   Step by step temperature control instead of one-shot 
//   rescaling. It takes away/pump in maximum of 1% of velocities
//   per time step. 
//
void VerletUpdate(
	int nAtom,                                // number of atoms
	double dt,                                // time step in AU.
	double *mass,                             // atomic mass in AU.
	double temp,                              // temperature in Kelvin
	double *x,  double *y,  double *z,        // atomic positions
	double *px, double *py, double *pz,       // previous positions
	double *vx, double *vy, double *vz,       // atomic velocities
	double *fx, double *fy, double *fz){      // atomic forces

	int n;
	double *cx, *cy, *cz; // current position
	
	// allocation
	cx    = calloc(nAtom,sizeof(double));
	cy    = calloc(nAtom,sizeof(double));
	cz    = calloc(nAtom,sizeof(double));

	/*
	for(n=0; n < nAtom; n++){

		// temporary store current position
		cx[n] = x[n];
		cy[n] = y[n];
		cz[n] = z[n];

		// update position
		x[n] = 2.0*cx[n] - px[n] + dt*dt*fx[n]/mass[n];
		y[n] = 2.0*cy[n] - py[n] + dt*dt*fy[n]/mass[n];
		z[n] = 2.0*cz[n] - pz[n] + dt*dt*fz[n]/mass[n];

		// compute velocity
		vx[n] = (x[n]-px[n])/(2.0*dt);
		vy[n] = (y[n]-py[n])/(2.0*dt);
		vz[n] = (z[n]-pz[n])/(2.0*dt);

		// shift previous position
		px[n] = cx[n];
		py[n] = cy[n];
		pz[n] = cz[n];
	}
	*/

	// temporary store current position
	for(n=0; n < nAtom; n++){
		cx[n] = x[n];
		cy[n] = y[n];
		cz[n] = z[n];
	}

	// compute velocity
	for(n=0; n < nAtom; n++){
		vx[n] = (cx[n]-px[n])/dt + 0.5*dt*fx[n]/mass[n];
		vy[n] = (cy[n]-py[n])/dt + 0.5*dt*fy[n]/mass[n];
		vz[n] = (cz[n]-pz[n])/dt + 0.5*dt*fz[n]/mass[n];
	}

	// rescale velocity to keep temperature if requested
	if(temp>0.0){
		double alpha = sqrt(1.5*MD_kB*temp*nAtom/KineticEnergy(nAtom,mass,vx,vy,vz));

		// cap rescaling factor to 1% difference
		//
		// Tested by H.Chachiyo on the Carbonic Acid reaction
		// to be presented by the Center of Excellence in
		// Quantum Modelling in Naresuan Research Fair 2017
		//
		// Teepanis Chachiyo - 11/7/2017
		if(fabs(alpha) > 1.01) alpha = 1.01*alpha/fabs(alpha);
		if(fabs(alpha) < 0.99) alpha = 0.99*alpha/fabs(alpha);

		for(n=0; n < nAtom; n++){
			vx[n] *= alpha;
			vy[n] *= alpha;
			vz[n] *= alpha;
		}
	}

	// update position
	for(n=0; n < nAtom; n++){
		x[n] = cx[n] + dt*vx[n] + 0.5*dt*dt*fx[n]/mass[n];
		y[n] = cy[n] + dt*vy[n] + 0.5*dt*dt*fy[n]/mass[n];
		z[n] = cz[n] + dt*vz[n] + 0.5*dt*dt*fz[n]/mass[n];
	}

	// shift previous position
	for(n=0; n < nAtom; n++){
		px[n] = cx[n];
		py[n] = cy[n];
		pz[n] = cz[n];
	}

	free(cx);
	free(cy);
	free(cz);
}

void printVelocity(
	int nAtom,                         // number of atom
	int *Z,                            // atomic number
	double *vx, double *vy, double *vz,// velocity
	FILE *fd){                         // output file pointer

	int A;        // atom index;
	char str[64]; // atom name in short format

	fprintf(fd,
"-------------------------------------------------------------\n"
"                     Atom Velocities (nm/ps)                 \n"
"     Atom         Vx             Vy             Vz           \n"
"-------------------------------------------------------------\n");
#define CAPVALUE(a) (fabs(a)<1.0E-6)?0.0:a
	for(A=0; A < nAtom; A++){
		Z2SymShort(Z[A],str);
		fprintf(fd, "   %5s %15.6f%15.6f%15.6f\n",
		        str,
		        CAPVALUE(vx[A]),
		        CAPVALUE(vy[A]),
		        CAPVALUE(vz[A]));
	}
	fprintf(fd,
"-------------------------------------------------------------\n");
#undef CAPVALUE
}


struct Molecule_t *readVelocity_XYZ(FILE *inFile){
	int i;
	int nItem;
	struct Molecule_t *mol;
	char str[256];

	// allocate memory
	mol = (struct Molecule_t *) calloc(1, sizeof(struct Molecule_t));
	if(mol == NULL){
		printf("readVelocity_XYZ: Error - Cannot allocate memory\n");
		exit(-1);
	}

	// reading number of atoms in unit cell
	if(fgets(str, 256, inFile)!=str){
		printf("readVelocity_XYZ: Error - Reading number of atoms\n");
		exit(EXIT_FAILURE);
	}
	nItem = sscanf(str, "%d", &(mol->nAtom));
	if(nItem != 1){
		printf("readVelocity_XYZ: Error - Reading number of atoms\n");
		exit(-1);
	}

	printf("Reading velocity in XYZ file format\n");
	printf("Detected %d atoms\n", mol->nAtom);

	// allocate atomic number and coordinates
	mol->Z = calloc(mol->nAtom, sizeof(int));
	if(mol->Z == NULL){
		printf("readVelocity_XYZ: Error - Cannot allocate atomic number\n");
		exit(EXIT_FAILURE);
	}
	mol->x = calloc(mol->nAtom, sizeof(double));
	mol->y = calloc(mol->nAtom, sizeof(double));
	mol->z = calloc(mol->nAtom, sizeof(double));
	if(mol->x == NULL || mol->y == NULL || mol->z == NULL){
		printf("readVelocity_XYZ: Error - Cannot allocate velocity\n");
		exit(EXIT_FAILURE);
	}

	// reading comments to buffer
	if(fgets(str, 256, inFile)!=str){
		printf("readVelocity_XYZ: Error - Reading comment from XYZ file\n");
		exit(EXIT_FAILURE);
	}

	// reading atomic number and velocities
	for(i=0; i < mol->nAtom; i++){
		nItem = fscanf(inFile, "%s %lf %lf %lf", str,
		                                         mol->x+i,
		                                         mol->y+i,
		                                         mol->z+i);
		if(nItem != 4){
			printf("readVelocity_XYZ: Error - Reading velocity\n");
			exit(EXIT_FAILURE);
		}
		// parse atom name in short format
		mol->Z[i] = sym2Z(str,SYMB_SHORTNAME);

	}

	// print output
	printVelocity(mol->nAtom, mol->Z, mol->x, mol->y, mol->z, stdout);

	// set the total molecular charge to zero by default
	mol->Q = 0;

	return mol;

}


//
// qmd: perform quantum molecular dynamics simulation
//
// July 11, 2017 - Teepanis Chachiyo & Hathaithip Chachiyo
//   Report Mullikan charge after UHF to delineate chemical species 
//
// March 10, 2016 - Teepanis Chachiyo & Chutchawan Jaisuk
//   Initial implementation and testing
//
/////////////////////////////////////////////////////////
////////INTERNAL VARIABLES MUST BE IN MOLECULAR DYNMAICS UNIT
//
// DISANCE  NANO-METER
// TIME     PICO-SECOND
// VELOCITY NM/PS
// MASS     DALTON
// FORCE    KJ/MOL/NANO-METER
// TEMPERATURE KELVIN
// BOLTZMAN 0.0083144621 KJ/MOL/TEMP
/////////////////////////////////////////////////////////
void qmd(
	int dbSize,                         // number of record in basis set db
	const struct GTOBasisSet_t *basisDB,// pointer to basis set database
	struct Molecule_t *mol,             // returned molecular coordinate
	struct option_t *opt){              // options

	int nBasis;                 // number of basis function
	struct GTOBasis_t *gto;     // pointer to basis function storage
	double *x, *y, *z;          // position in nm
	double *px, *py, *pz;       // previous position in nm
	double *fx,*fy,*fz;         // forces in kJ/mol/nm
	double *vx,*vy,*vz;         // velocity in nm/ps
	double *mass;               // mass of atoms in Dalton.
	FILE *fd;                   // file pointer

	double *CA, *eA, *CB, *eB;// molecular orbitals and eigen values
	double Etot;              // total electronic energy
	double Kinetic;           // total kinetic energy
	int savedSCFGuess;        // global guess option

	int nIter;                // iteration index
	int i;                    // loop index

	///////////////////
	// report options
	///////////////////
	printf("Time Step:               %10.4f ps\n",opt->tStep);
	printf("Number of Iterations:    %10d steps\n",opt->maxQMDStep);
	if(opt->Ex != 0.0 || opt->Ey != 0.0 || opt->Ez != 0.0){
	printf("Uniform Electric Field: (%9.4f,%9.4f,%9.4f) AU\n"
	       "                        (%9.4f,%9.4f,%9.4f) V/nm\n",
		opt->Ex, opt->Ey, opt->Ez,
		opt->Ex*AU2VOLTPERNM, opt->Ey*AU2VOLTPERNM, opt->Ez*AU2VOLTPERNM);
	if(opt->qmdFreq != 0.0){
	printf("Oscillate Electric Field with Frequency: %12.6f THz\n"
	       "          Equivalent to the Wavelength:  %12.6f nm\n",
		opt->qmdFreq, MD_LIGHTSPEED/opt->qmdFreq);
	}
	}
	if(opt->keepTemp > 0.0)
	printf("Maintain Temperature:    %10.2f K\n", opt->keepTemp);
	if(opt->initTemp > 0.0)
	printf("Initial Random Velocity: %10.4f K\n",opt->initTemp);
	printf("Output Trajectory File:  %10s\n",opt->traj);
	if(opt->initVel[0] != '\0')
	printf("Initial Velocity File:   %10s\n",opt->initVel);
	
	printf("\n");fflush(stdout);

	// save global guess option
	savedSCFGuess = opt->SCFGuess;

	// memory allocation
	x    = calloc(mol->nAtom,sizeof(double));
	y    = calloc(mol->nAtom,sizeof(double));
	z    = calloc(mol->nAtom,sizeof(double));
	px   = calloc(mol->nAtom,sizeof(double));
	py   = calloc(mol->nAtom,sizeof(double));
	pz   = calloc(mol->nAtom,sizeof(double));
	fx   = calloc(mol->nAtom,sizeof(double));
	fy   = calloc(mol->nAtom,sizeof(double));
	fz   = calloc(mol->nAtom,sizeof(double));
	vx   = calloc(mol->nAtom,sizeof(double));
	vy   = calloc(mol->nAtom,sizeof(double));
	vz   = calloc(mol->nAtom,sizeof(double));
	mass = calloc(mol->nAtom,sizeof(double));

	// determine atomic mass in Dalton
	for(i=0; i < mol->nAtom; i++)
		mass[i] = Z2mass(mol->Z[i]);

	/////////////////////////////
	// handle initial velocities
	/////////////////////////////

	// read initial velocity in nm/ps
	if(opt->initVel[0] != '\0'){
		struct Molecule_t *vel;

		// read file
		fd = fopen(opt->initVel,"r");
		if(fd == NULL){
			printf("qmd: Cannot open file %s\n", opt->initVel);
			exit(-1);
		}
		vel = readVelocity_XYZ(fd);

		// validation
		if(vel->nAtom != mol->nAtom){
			printf("qmd: error - different number of atoms in velocity files");
			exit(-1);
		}
		for(i=0; i < mol->nAtom; i++)
		if(vel->Z[i] != mol->Z[i]){
			printf("qmd: error - different atom types at atom %d\n",i+1);
			exit(-1);
		}

		// load value to vx,vy,vz
		for(i=0; i < vel->nAtom; i++){
			vx[i] = vel->x[i];
			vy[i] = vel->y[i];
			vz[i] = vel->z[i];
		}
		
		// clean memory
		fclose(fd);
		vel = cleanMolecule(vel);
	}

	// generate velocity
	if(opt->initTemp>0.0){
		for(i=0; i < mol->nAtom; i++){
			double speed;
			double theta, phi;

			// compute speed
			speed = sqrt(3.0*MD_kB*opt->initTemp/mass[i]);
			// random orientation
			theta = (double)rand()/RAND_MAX*M_PI;
			phi   = (double)rand()/RAND_MAX*M_2_PI;
			vx[i] = speed*sin(theta)*cos(phi);
			vy[i] = speed*sin(theta)*sin(phi);
			vz[i] = speed*cos(theta);
		}
	}

	// load mol (in Bohr) to x,y,z in nm
	for(i=0; i < mol->nAtom; i++){
		x[i] = mol->x[i]*AU2MD_DISTANCE;
		y[i] = mol->y[i]*AU2MD_DISTANCE;
		z[i] = mol->z[i]*AU2MD_DISTANCE;
	}

	// allocate memory for molecular orbitals and their eigen values
	gto   = genBasis(mol, &nBasis, dbSize, basisDB);
	cleanGTOBasis(gto,nBasis);
	CA = calloc(nBasis*nBasis,sizeof(double));
	eA = calloc(nBasis,sizeof(double));
	CB = calloc(nBasis*nBasis,sizeof(double));
	eB = calloc(nBasis,sizeof(double));
	if(CA==NULL || eA==NULL || CB==NULL || eB==NULL){
		printf("qmd: error - cannot allocate memory\n");
		exit(-1);
	}

	// open trajectory file
	fd = fopen(opt->traj,"w");
	if(fd == NULL){
		printf("qmd: Cannot open file %s\n", opt->traj);
		exit(-1);
	}

	nIter=0;
	do{
		///////////////////////////////////////////////
		// perform scf calculation and compute forces
		///////////////////////////////////////////////
		printf(
		"                                                             \n"
		"                                                             \n"
		"-------------------------------------------------------------\n"
		"-----       QUANTUM MOLECULAR DYNAMICS Step %5d       -----\n"
		"-------------------------------------------------------------\n",
		nIter+1);
		fflush(stdout);

		// generate basis function
		gto   = genBasis(mol, &nBasis, dbSize, basisDB);

		//
		// scf calculation
		//

		// load orbital from previous cycle to form initial guess
		if(nIter>0) opt->SCFGuess = SCFGUESS_CACB;

		Etot = uhf(nBasis, gto, mol, 
		    get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		    CA, CB, eA, eB, opt);
		if(Etot==0.0){
			printf("qmd: error SCF calculation did not converge\n");
			exit(-1);
		}

		mulliken(nBasis, gto, mol, 
		         get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity),
		         CA, CB, eA, eB, opt);

		// compute force
		uhf_force(nBasis, gto, mol, 
		      get_nEA(mol,opt->multiplicity), get_nEB(mol,opt->multiplicity), 
		      CA, CB, eA, eB, opt, fx, fy, fz);

		// convert force to kJol/mol/nm
		for(i=0; i < mol->nAtom; i++){
			fx[i] = fx[i]*AU2MD_FORCE;
			fy[i] = fy[i]*AU2MD_FORCE;
			fz[i] = fz[i]*AU2MD_FORCE;
		}

		////////////////////////////////////
		// update new molecular coordinate 
		////////////////////////////////////
		if(nIter==0){
			for(i=0; i < mol->nAtom; i++){
				px[i] = x[i];
				py[i] = y[i];
				pz[i] = z[i];
			}
			// EulerCromer advance velocity step so 
			// we need to print it before
			printVelocity(mol->nAtom,mol->Z,vx,vy,vz,stdout);
			Kinetic = KineticEnergy(mol->nAtom, mass, vx, vy, vz);
			EulerCromerUpdate(mol->nAtom, opt->tStep, mass, 
				x, y, z, vx, vy, vz, fx, fy, fz);
		}else{
			VerletUpdate(mol->nAtom, opt->tStep, mass, opt->keepTemp,
				x, y, z, px, py, pz, vx, vy, vz, fx, fy, fz);
			printVelocity(mol->nAtom,mol->Z,vx,vy,vz,stdout);
			Kinetic = KineticEnergy(mol->nAtom, mass, vx, vy, vz);
		}

		// write to trajectory
		fprintf(fd, "%d\nETOTAL=%.2f KINETIC=%.2f\n",mol->nAtom,
			Etot*HARTREE2KJMOL+Kinetic,Kinetic);
		for(i=0; i < mol->nAtom; i++){
			char str[256];
			Z2SymShort(mol->Z[i],str);
			fprintf(fd,"%s %12.5f %12.5f %12.5f\n", str, 
			mol->x[i]/ANGSTROM2BOHR, 
			mol->y[i]/ANGSTROM2BOHR, 
			mol->z[i]/ANGSTROM2BOHR);
		}
		fflush(fd);

		// report convergence status
		printf("\n");
		printf("ELECTRONIC     ENERGY        %20.2f kJ/mol\n", Etot*HARTREE2KJMOL);
		printf("NUCLEI KINETIC ENERGY        %20.2f kJ/mol\n", Kinetic);
		printf("TOTAL          ENERGY        %20.2f kJ/mol\n", Etot*HARTREE2KJMOL+Kinetic);
		fflush(stdout);

		// copy coordinate from nm to mol (Bohr)
		for(i=0; i < mol->nAtom; i++){
			mol->x[i] = x[i]/AU2MD_DISTANCE ;
			mol->y[i] = y[i]/AU2MD_DISTANCE;
			mol->z[i] = z[i]/AU2MD_DISTANCE;
		}

		// free memory
		cleanGTOBasis(gto,nBasis);

		fflush(stdout);

		nIter++;

	}while(nIter < opt->maxQMDStep);

	// close trajectory file
	fclose(fd);

	// free memory
	free(CA);
	free(eA);
	free(CB);
	free(eB);

	free(x);
	free(y);
	free(z);
	free(px);
	free(py);
	free(pz);
	free(fx);
	free(fy);
	free(fz);
	free(vx);
	free(vy);
	free(vz);
	free(mass);

	// set the guess option to original value
	opt->SCFGuess = savedSCFGuess;
}
