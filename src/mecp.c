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
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "mecp.h"
#include "lin.h"
#include "grad.h"
#include "uhf.h"
#include "mol.h"
#include "util.h"
#include "optimize.h"


// runGaussian: call Gaussian program and compute force and energy
//
// March 19, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
double runGaussian(
	struct Molecule_t *mol, // pointer to molecular structure
	char *exe,              // pointer to execution string
	char *input,            // pointer to input file name
	double *fx,             // returned force in x-direction
	double *fy,             // returned force in y-direction
	double *fz){            // returned force in z-direction

	double E;          // energy
	int i;             // generic integer
	char str[256];     // short buffer
	char buf[1024];    // buffer string
	char inFile[256];  // real input file name
	FILE *fdIn;        // file pointer of the original input
	FILE *fdOut;       // file pointer of the real input 

	// open input file
	sprintf(str,"%s.com",input);
	if((fdIn=fopen(str,"r"))==NULL){
		printf("runGaussian - error cannot open file %s.com\n",input);
		exit(-1);
	}

	// open real input file
	sprintf(inFile,"SQ_%s.com",input);
	if((fdOut=fopen(inFile,"w+"))==NULL){
		printf("runGaussian - error cannot open file %s\n",inFile);
		exit(-1);
	}

	// validate input file for NOSYMM, FORCE, and SQ_GEOMETRY keyword
	rewind(fdIn);
	if(findf(fdIn,1,"NOSYMM")==EOF){
		printf("runGaussian - error cannot find keyword NOSYMM in file %s.com\n",input);
		exit(-1);
	}
	rewind(fdIn);
	if(findf(fdIn,1,"FORCE")==EOF){
		printf("runGaussian - error cannot find keyword FORCE in file %s.com\n",input);
		exit(-1);
	}
	rewind(fdIn);
	if(findf(fdIn,1,"SQ_GEOMETRY")==EOF){
		printf("runGaussian - error cannot find keyword FORCE in file %s.com\n",input);
		exit(-1);
	}

	// begin copy line by line until reaching SQ_GEOMETRY keyword
	rewind(fdIn);
	for(;fgets(buf,1024,fdIn)!=NULL;){
		
		// check SQ_GEOMETRY matching
		if(strstr(buf,"SQ_GEOMETRY")!=NULL) break;

		// write to the real input file
		fprintf(fdOut,"%s",buf);
	}

	// print molecular structure
#define CAPVALUE(a) (fabs(a)<1.0E-8)?0.0:a
	for(i=0; i < mol->nAtom; i++){
		Z2SymShort(mol->Z[i],str);
		fprintf(fdOut, "   %5s %15.8f%15.8f%15.8f\n",
		        str,
		        CAPVALUE(mol->x[i])*BOHR2ANGSTROM,
		        CAPVALUE(mol->y[i])*BOHR2ANGSTROM,
		        CAPVALUE(mol->z[i])*BOHR2ANGSTROM);
	}

	// copy the rest
	for(;fgets(buf,1024,fdIn)!=NULL;){

		// write to the real input file
		fprintf(fdOut,"%s",buf);
	}

	// close file
	fclose(fdIn);
	fclose(fdOut);

	// execute gaussian
	sprintf(str,"%s %s",exe,inFile);
	system(str);

	// open output from Gaussian
	sprintf(str,"SQ_%s.log",input);
	if((fdOut=fopen(str,"r"))==NULL){
		printf("runGaussian - error cannot open Gaussian output %s\n",str);
		exit(-1);
	}

	// read energy
	if(findf(fdOut,1,"Done:")==EOF){
		printf("runGaussian - error cannot find total energy in %s.log\n",input);
		exit(-1);
	}
	if(fscanf(fdOut,"%*s %*s %lf",&E)!=1){
		printf("runGaussian - error cannot find total energy in %s.log\n",input);
		exit(-1);
	}
	//printf("Total Energy is %20.8lf Hartrees\n",E);

	// read forces
	if(findf(fdOut,1,"(Hartrees/Bohr)")==EOF){
		printf("runGaussian - error cannot find total energy in %s.log\n",input);
		exit(-1);
	}
	fscanf(fdOut,"%*s %*s %*s %*s %*s %*s");
	for(i=0; i< mol->nAtom; i++){
		// scan force
		if(fscanf(fdOut,"%*s %*s %lf %lf %lf",fx+i,fy+i,fz+i)!=3){
			printf("runGaussian - error while reading forces\n");
			exit(-1);
		}
	}

	// print out
	//printForce(mol,fx,fy,fz);

	// erase the files
	sprintf(str,"SQ_%s.com",input);
	if(unlink(str) != 0){
		printf("runGaussian - error cannot delete file %s\n",str);
		exit(-1);
	}
	sprintf(str,"SQ_%s.log",input);
	if(unlink(str) != 0){
		printf("runGaussian - error cannot delete file %s\n",str);
		exit(-1);
	}

	// close
	fclose(fdOut);

	return E;
}


// Hessian_BFGS: compute the Hessian matrix using BFBS schemme
// as explained in the equation (11) in
//
// Chachiyo, Teepanis, and Jorge H. Rodriguez. "A direct method for locating 
// minimum-energy crossing points (MECPs) in spin-forbidden transitions and 
// nonadiabatic reactions." The Journal of chemical physics 123 (2005): 094711.
//
// ***Note the typo in the paper:the last term should have been a minus sign
//
//
//                      T               T
// H = H + dGrad * dGrad       HdR * HdR
//         -------------   -  -----------
//             alpha             beta
//
//
// where constant alpha = dGrad DOT dR
//       constant  beta = Hdr   DOT dR
//       vector     HdR = H * dR  
//
// March 17, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void Hessian_BFGS(
	int nDim,            // matrix dimension
	const double *dR,    // displacement vector
	const double *dGrad, // change of gradient vector,
	double *H){          // original the Hessian and the returned updated values

	int i,j;           // loop index
	double alpha;      // dGrad DOT dR
	double beta;       // Hdr   DOT dR
	double *HdR;       // auxiliary vector

	// memory allocation
	HdR=calloc(nDim, sizeof(double));
	if(HdR==NULL){
		printf("Hessian_BFGS: error cannnot allocate memory\n");
		exit(-1);
	}

	// compute alpha
	alpha = 0.0;
	for(i=0; i < nDim; i++)
		alpha += dGrad[i] * dR[i];

	// compute HdR
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++)
		HdR[i] += H[i*nDim+j] * dR[j];

	// compute beta
	beta = 0.0;
	for(i=0; i < nDim; i++)
		beta += HdR[i] * dR[i];

	// check if alpha is zero
	if(alpha==0.0){
		printf("Hessian_BFGS: error dR perpendicular to dGrad\n");
		exit(-1);
	}

	// check if beta is zero
	if(beta==0.0){
		printf("Hessian_BFGS: error dR perpendicular to HdR\n");
		exit(-1);
	}

	// update the Hessian
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++)
		H[i*nDim+j] += dGrad[i]*dGrad[j]/alpha - HdR[i]*HdR[j]/beta;

	// free memory
	free(HdR);
}


// stepVector: calculate stepping vector using eigen vector method as 
// described in Equ 17-21 in
//
// Chachiyo, Teepanis, and Jorge H. Rodriguez. "A direct method for locating 
// minimum-energy crossing points (MECPs) in spin-forbidden transitions and 
// nonadiabatic reactions." The Journal of chemical physics 123 (2005): 094711.
//
// March 17, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void stepVector_Eigen(
	int nDim,            // dimension of the matrix
	double C,            // constraint C
	const double *H,     // pointer to the Hessian matrix
	const double *g,     // pointer to gradient vector
	const double *J,     // pointer to Jacobian vector
	double maxStepSize,  // maximum step size in each direction
	double *gL,          // returned lagrange gradient
	double *dR){         // returned step vector

	int i,j;        // loop index
	double r;       // generic double precision
	double *I;      // identity matrix
	double *h;      // eigen vector of hessian
	double *e;      // eigen value of hessian
	double *alpha;  // constant in equation (20)
	double *beta;   // constant in equation (21)
	double gamma;   // constant in equation (18)

	// allocate memory
	I     = calloc(nDim*nDim, sizeof(double));
	h     = calloc(nDim*nDim, sizeof(double));
	e     = calloc(nDim,      sizeof(double));
	alpha = calloc(nDim,      sizeof(double));
	beta  = calloc(nDim,      sizeof(double));
	if(I==NULL || h==NULL || e==NULL || alpha==NULL || beta==NULL){
		printf("stepVector_Eigen - error cannot allocate memory\n");
		exit(-1);
	}

	// compute eigen values and eigen vector of hessian
	for(i=0; i < nDim; i++) I[i*nDim+i] = 1.0;
	gen_sym_eigen(nDim, H, I, e, h);

	// compute alpha and beta
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++){
		alpha[i] += g[j]*h[i*nDim+j];
		beta[i]  += J[j]*h[i*nDim+j];
	}

	// compute gamma
	gamma=0.0;
	r=0.0;
	for(i=0; i < nDim; i++)
		if(fabs(e[i]) > 0.0){
			gamma += alpha[i]*beta[i]/e[i];
			r     +=  beta[i]*beta[i]/e[i];	
		}
	if(r==0.0) gamma=0.0; else gamma = (C - gamma)/r;

	// reset step vector
	for(i=0; i < nDim; i++) dR[i] = 0.0;

	// compute step vector
	for(i=0; i < nDim; i++)
		if(fabs(e[i]) > 0.0){
			r = ( alpha[i] + gamma * beta[i] ) / e[i];
			for(j=0; j < nDim; j++)
				dR[j] -= r * h[i*nDim+j];
		}

	// validate step size
	if(maxStepSize <= 0.0){
		printf("stepVector_Eigen: invalid range of maxStepSize");
		exit(-1);
	}

	// truncate the step size to maxStepSize
	r = 0;
	for(i=0; i < nDim; i++) r+= dR[i]*dR[i];
	r = sqrt(r);

	if(r > maxStepSize)
	for(i=0; i < nDim; i++)
		dR[i] = dR[i]*maxStepSize/r;

	// compute lagrange gradient for checking convergence
	for(i=0; i < nDim; i++)
		gL[i] = g[i] + gamma * J[i];

	// free memory
	free(I);
	free(h);
	free(e);
	free(alpha);
	free(beta);
}


//
// mecp: find the minimum energy crossing point as described in
//
// Chachiyo, Teepanis, and Jorge H. Rodriguez. "A direct method for locating 
// minimum-energy crossing points (MECPs) in spin-forbidden transitions and 
// nonadiabatic reactions." The Journal of chemical physics 123 (2005): 094711.
//
//
// March 17, 2013 - Teepanis Chachiyo
//   Initial implementation and testing
//
// Jan 1, 2016 - Teepanis Chachiyo
//   Delete translation from dR
//
void mecp(
	int dbSize,                         // number of record in basis set db
	const struct GTOBasisSet_t *basisDB,// pointer to basis set database
	struct Molecule_t *mol,             // returned molecular coordinate
	struct option_t *opt){              // options

	int nBasis;               // number of basis function
	struct GTOBasis_t *gto;   // pointer to basis function storage
	int nDim;                 // degree of freedoms
	double *HA,*HB,*H;        // Hessian of state A, B, and both
	double *gA,*gB,*g;        // gradient vector of state A, B, and both
	double *J;                // Jacobi vector
	double *dGA,*dGB;         // change of gradient vector
	double *dR;               // steping vector
	double *fx,*fy,*fz;       // forces acting on nuclei
	double *CA, *eA, *CB, *eB;// molecular orbitals and eigen values
	double EtotA,EtotB;       // total energy
	double C;                 // constrant

	int nIter;                // iteration index
	double force_max;         // maximum force on nuclei
	double force_rms;         // root mean square of force
	double step_max;          // maximum step size
	double step_rms;          // root mean square of step size

	int i,j;                  // loop index
	char str[256];            // string buffer

	// set degree of freedom
	nDim = mol->nAtom * 3;

	// memory allocation
	HA   = calloc(nDim*nDim, sizeof(double));
	HB   = calloc(nDim*nDim, sizeof(double));
	H    = calloc(nDim*nDim, sizeof(double));
	gA   = calloc(nDim,      sizeof(double));
	gB   = calloc(nDim,      sizeof(double));
	g    = calloc(nDim,      sizeof(double));
	J    = calloc(nDim,      sizeof(double));
	dR   = calloc(nDim,      sizeof(double));
	dGA  = calloc(nDim,      sizeof(double));
	dGB  = calloc(nDim,      sizeof(double));
	fx   = calloc(mol->nAtom,sizeof(double));
	fy   = calloc(mol->nAtom,sizeof(double));
	fz   = calloc(mol->nAtom,sizeof(double));
	if(HA==NULL || HB==NULL || H==NULL ||
	   gA==NULL || gB==NULL || g==NULL || J==NULL ||
	   dR==NULL || dGA==NULL || dGB==NULL ||
	   fx==NULL || fy==NULL || fz==NULL){
		printf("mecp - error cannot allocate memory\n");
		exit(-1);
	}

	// set initial Hessian to identity matrix
	for(i=0; i < nDim; i++)
	for(j=0; j < nDim; j++){
		if(i==j) HA[i*nDim+j] = 1.0; else HA[i*nDim+j] = 0.0;
		if(i==j) HB[i*nDim+j] = 1.0; else HB[i*nDim+j] = 0.0;
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
		"-----       MECP CALCULATIONS Step %5d                -----\n"
		"-------------------------------------------------------------\n",
		nIter+1);
		fflush(stdout);

		//////////////////////////////////////////////////////
		///// using external program: Gaussian if requested
		//////////////////////////////////////////////////////
		if(strlen(opt->gaussEXE)>0){
			
			// compute energy and force of state A
			printf("\n[Begin] executing Gaussian program for state A\n"); fflush(stdout);
			EtotA = runGaussian(mol, opt->gaussEXE, opt->gaussINA, fx, fy, fz);
			printf("[Done]  executing Gaussian program for state A\n"); fflush(stdout);

			// construct gradient vector for state A
			for(i=0; i < mol->nAtom; i++){
				// change of gradient vector
				dGA[i*3+0] = -fx[i] - gA[i*3+0];
				dGA[i*3+1] = -fy[i] - gA[i*3+1];
				dGA[i*3+2] = -fz[i] - gA[i*3+2];

				// current gradient vector
				gA[i*3+0] = -fx[i];
				gA[i*3+1] = -fy[i];
				gA[i*3+2] = -fz[i];
			}

			// update Hessian for state A if all info are available
			if(nIter > 0) Hessian_BFGS(nDim, dR, dGA, HA);

			// compute energy and force of state B
			printf("\n[Begin] executing Gaussian program for state B\n"); fflush(stdout);
			EtotB = runGaussian(mol, opt->gaussEXE, opt->gaussINB, fx, fy, fz);
			printf("[Done]  executing Gaussian program for state B\n"); fflush(stdout);

			// construct gradient vector for state B
			for(i=0; i < mol->nAtom; i++){
				// change of gradient vector
				dGB[i*3+0] = -fx[i] - gB[i*3+0];
				dGB[i*3+1] = -fy[i] - gB[i*3+1];
				dGB[i*3+2] = -fz[i] - gB[i*3+2];

				// current gradient vector
				gB[i*3+0] = -fx[i];
				gB[i*3+1] = -fy[i];
				gB[i*3+2] = -fz[i];
			}

			// update Hessian for state B if all info are available
			if(nIter > 0) Hessian_BFGS(nDim, dR, dGB, HB);

		}else{
		/////////////////////////////////////////////////////////
		// using Siam Quantum to compute energies and forces
		/////////////////////////////////////////////////////////

			//
			// scf preparations
			//

			// generate basis function
			gto   = genBasis(mol, &nBasis, dbSize, basisDB);

			// allocate meory for molecular orbitals and their eigen values
			CA = calloc(nBasis*nBasis,sizeof(double));
			eA = calloc(nBasis,sizeof(double));
			CB = calloc(nBasis*nBasis,sizeof(double));
			eB = calloc(nBasis,sizeof(double));
			if(CA==NULL || eA==NULL || CB==NULL || eB==NULL){
				printf("mecp: error - cannot allocate memory\n");
				exit(-1);
			}

			//
			// compute energy and force of state A
			//
#define SWAPSTR(a,b) strcpy(str,a); strcpy(a,b); strcpy(b,str);
			// replace density matrix and checkpoint file name
			SWAPSTR(opt->DMatrixFile,opt->DMatrixFileA);
			SWAPSTR(opt->CheckFile,opt->CheckFileA);

			EtotA = uhf(nBasis, gto, mol, 
			    get_nEA(mol,opt->mecpMA), get_nEB(mol,opt->mecpMA), 
			    CA, CB, eA, eB, opt);
			if(EtotA==0.0){
				printf("mecp: error SCF calculation did not converge\n");
				exit(-1);
			}

			// put back the density matrix and checkpoint file name
			SWAPSTR(opt->DMatrixFile,opt->DMatrixFileA);
			SWAPSTR(opt->CheckFile,opt->CheckFileA);

			// compute force for state A
			uhf_force(nBasis, gto, mol, 
			      get_nEA(mol,opt->mecpMA), get_nEB(mol,opt->mecpMA), 
			      CA, CB, eA, eB, opt, fx, fy, fz);

			//
			// construct gradient vector for state A
			//
			for(i=0; i < mol->nAtom; i++){
				// change of gradient vector
				dGA[i*3+0] = -fx[i] - gA[i*3+0];
				dGA[i*3+1] = -fy[i] - gA[i*3+1];
				dGA[i*3+2] = -fz[i] - gA[i*3+2];

				// current gradient vector
				gA[i*3+0] = -fx[i];
				gA[i*3+1] = -fy[i];
				gA[i*3+2] = -fz[i];
			}

			//
			// update Hessian for state A if all info are available
			//
			if(nIter > 0) Hessian_BFGS(nDim, dR, dGA, HA);


			//
			// compute energy and force of state A
			//

			// replace density matrix and checkpoint file name
			SWAPSTR(opt->DMatrixFile,opt->DMatrixFileB);
			SWAPSTR(opt->CheckFile,opt->CheckFileB);

			EtotB = uhf(nBasis, gto, mol, 
			    get_nEA(mol,opt->mecpMB), get_nEB(mol,opt->mecpMB), 
			    CA, CB, eA, eB, opt);
			if(EtotB==0.0){
				printf("mecp: error SCF calculation did not converge\n");
				exit(-1);
			}

			// put back the density matrix and checkpoint file name
			SWAPSTR(opt->DMatrixFile,opt->DMatrixFileB);
			SWAPSTR(opt->CheckFile,opt->CheckFileB);

			// compute force for state B
			uhf_force(nBasis, gto, mol, 
			      get_nEA(mol,opt->mecpMB), get_nEB(mol,opt->mecpMB), 
			      CA, CB, eA, eB, opt, fx, fy, fz);


			//
			// construct gradient vector for state B
			//
			for(i=0; i < mol->nAtom; i++){
				// change of gradient vector
				dGB[i*3+0] = -fx[i] - gB[i*3+0];
				dGB[i*3+1] = -fy[i] - gB[i*3+1];
				dGB[i*3+2] = -fz[i] - gB[i*3+2];

				// current gradient vector
				gB[i*3+0] = -fx[i];
				gB[i*3+1] = -fy[i];
				gB[i*3+2] = -fz[i];
			}

			//
			// update Hessian for state B if all info are available
			//
			if(nIter > 0) Hessian_BFGS(nDim, dR, dGB, HB);

			//
			// free memory
			//
			cleanGTOBasis(gto,nBasis);
			free(CA);
			free(eA);
			free(CB);
			free(eB);	
		}


		////////////////////////////////////
		// update new molecular coordinate 
		////////////////////////////////////

		// compute total Hessian
		for(i=0; i < nDim; i++)
		for(j=0; j < nDim; j++)
			H[i*nDim+j] = (HA[i*nDim+j] + HB[i*nDim+j])/2.0;

		// compute total gradient
		for(i=0; i < nDim; i++)
			g[i] = (gA[i] + gB[i])/2.0;

		// compute Jacobi vector for n=1
		for(i=0; i < nDim; i++)
			J[i] = gA[i]-gB[i];

		// compute constraint residual
		C = EtotA - EtotB;

		// compute steping vector
		stepVector_Eigen(nDim, C, H, g, J, MAXSTEPSIZE, g, dR);

		// delete simple translation
		deleteTranslation(nDim, dR);

		// update molecular gemetry
		for(i=0; i < mol->nAtom; i++){
			mol->x[i] = mol->x[i] + dR[i*3+0];
			mol->y[i] = mol->y[i] + dR[i*3+1];
			mol->z[i] = mol->z[i] + dR[i*3+2];
		}

		// compute convergence criteria
		force_max = 0.0; force_rms = 0.0;
		step_max = 0.0; step_rms = 0.0;
		for(i=0; i < nDim; i++){
			if(fabs(g[i])  > force_max) force_max = fabs(g[i]);
			if(fabs(dR[i]) > step_max)  step_max  = fabs(dR[i]);
			force_rms += g[i]*g[i];
			step_rms  += dR[i]*dR[i];
		}
		force_rms = sqrt(force_rms/nDim);
		step_rms  = sqrt(step_rms/nDim);

		// report total energy
		printf("\n");
		printf("MECP Information\n");
		printf("----------------\n");
		printf("State A Total Energy is %20.8lf Hartrees\n",EtotA);
		printf("State B Total Energy is %20.8lf Hartrees\n",EtotB);
		printf("\n");
		printMolecule_XYZ(mol, stdout);

		// report convergence status
		printf("Convergence Criterion         Value        Threshold\n");
		printf("  Max Lagrange Force        %10.6f   %10.6f  ", force_max, MECP_CONV_FORCEMAX);
		if(force_max > MECP_CONV_FORCEMAX) printf("NO\n"); else printf("YES\n");
		printf("  RMS Lagrange Force        %10.6f   %10.6f  ", force_rms, MECP_CONV_FORCERMS);
		if(force_rms > MECP_CONV_FORCERMS) printf("NO\n"); else printf("YES\n");
		printf("  Max Lagrange Displacement %10.6f   %10.6f  ", step_max, MECP_CONV_DISPMAX);
		if(step_max > MECP_CONV_DISPMAX) printf("NO\n"); else printf("YES\n");
		printf("  RMS Lagrange Displacement %10.6f   %10.6f  ", step_rms, MECP_CONV_DISPRMS);
		if(step_rms > MECP_CONV_DISPRMS) printf("NO\n"); else printf("YES\n");
		fflush(stdout);

		nIter++;

		if(nIter > opt->mecpMax){
			printf("MECP has not converged, nIter >= %d\n",opt->mecpMax);
			exit(-1);
		}

	}while(force_max > MECP_CONV_FORCEMAX || force_rms > MECP_CONV_FORCERMS ||
	       step_max  > MECP_CONV_DISPMAX  || step_rms  > MECP_CONV_DISPRMS);

	// free memory
	free(HA);
	free(HB);
	free(H);
	free(gA);
	free(gB);
	free(g);
	free(J);
	free(dR);
	free(dGA);
	free(dGB);
	free(fx);
	free(fy);
	free(fz);
}
