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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "uhf.h"
#include "mol.h"
#include "basis.h"
#include "util.h"
#include "option.h"
#include "matrix.h"

// print out density infomation in xsf format
// XSF is XCrySDEN format used in molecular 
// and crystal structure. The file format
// is reversed engineered from AbInit cut3d
// output file.
//
// Note! xmin,xmax,ymin,ymax,zmin,zmax in Bohr
//       set them to zero for automatic adjustment
//
//       set nx,ny,nz to zero for auto 
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation
//
// Dec 15, 2008 - Teepanis Chachiyo
//  	Automatically determine min-max if the 
//  	arguments are all zero
//
//	Fix Bohr > Angstrom bug
//
// Nov 4, 2009 - Theerapon Khamla
//      Fix another Angstrom bug when calling rhf_rho
//
// Dec 31, 2009 - Theerapon Khamla
//      1) Computing electron density using density matrix
//      which helps increase speed of the calculations.
//      2) Add plotting molecular orbitals selected using option_t
//
// July 11, 2010 - Teepanis Chachiyo
//      Migrate to unrestricted calculation
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//      Add potential volume information
//
// Oct 8, 2012 - Teepanis Chachiyo
//      - Change file name to 'volume.xsf'
//      - Adjust the code to coincide with the new struct option_t
//      - Compute the entire xy-plane at a time
//
// Jan 16, 2014 - Teepanis Chachiyo
//		- Handle the whicMO index in the option.c instead of here
//
//
// Nov, 18, 2017 - Teepanis Chachiyo
//      - Add GRADOVER43 and GRADRS options 
//      - Can now adjust the number of points per Angstrom from main option
//
void xsfden(int nBasis,               // number of basis function
            struct GTOBasis_t * gto,  // pointer to basis set information
            struct Molecule_t * mol,  // pointer to molecular structure
            double *CA,               // calculated molecular orbital coeff.
            double *CB,               // same but for beta spin
            double xmin, double xmax, // minimum and maximum in x-direction
            double ymin, double ymax, // minimum and maximum in y-direction
            double zmin, double zmax, // minimum and maximum in z-direction
            int nx, int ny, int nz,   // number of grid point
            struct option_t *option){ // print according to options

	int i,j,k;           // loop indexes in x,y,z direction
	double dx, dy, dz;   // step size
	double value;        // value to plot
	FILE *fd;            // file pointer
	double *PA;          // density matix
	double *PB;          // density matrix
	int nEA;             // number of occupied alpha electron
	int nEB;             // number of occupied beta electron
	double *xyPlane;     // grid points in xy-plane

	printf("\n[Begin] printing volume in XSF format to 'volume.xsf'\n");
	fflush(stdout);

	//
	// compute density matrix
	//

	// memory allocation
	PA = calloc(nBasis*nBasis, sizeof(double));
	PB = calloc(nBasis*nBasis, sizeof(double));               
	if(PA==NULL || PB==NULL){                                             
		printf("xsfden: Error - Cannot allocate memory\n");
		exit(EXIT_FAILURE);                                  
	}

	// check number of electrons
	//if(get_nElectron(mol) % 2 != 0){
	//	printf("xsfden: error RHF requires even number of electrons\n");
	//	exit(EXIT_FAILURE);
	//}

	// number of occupy orbital
	nEA = get_nEA(mol,option->multiplicity);
	nEB = get_nEB(mol,option->multiplicity);

	// compute densit matrix
	uhf_getDMatrix(nBasis, nEA, CA, PA);
	uhf_getDMatrix(nBasis, nEB, CB, PB);

	// determine min-max if all argument are zero
	if(xmin==0.0 && xmax==0.0 && 
	   ymin==0.0 && ymax==0.0 && 
	   zmin==0.0 && zmax==0.0){
		xmin = mol->x[0]; xmax = mol->x[0];
		ymin = mol->y[0]; ymax = mol->y[0];
		zmin = mol->z[0]; zmax = mol->z[0];

		// loop to get exact min-max
		for(i=1; i < mol->nAtom; i++){
			if(xmin > mol->x[i]) xmin = mol->x[i];
			if(xmax < mol->x[i]) xmax = mol->x[i];
			if(ymin > mol->y[i]) ymin = mol->y[i];
			if(ymax < mol->y[i]) ymax = mol->y[i];
			if(zmin > mol->z[i]) zmin = mol->z[i];
			if(zmax < mol->z[i]) zmax = mol->z[i];
		}

		// push boundary 4.0 Bohr
		xmin -= 4.0; xmax += 4.0;
		ymin -= 4.0; ymax += 4.0;
		zmin -= 4.0; zmax += 4.0;
	}

	// convert bohr to angstroms
	xmin = xmin*BOHR2ANGSTROM; xmax = xmax*BOHR2ANGSTROM;
	ymin = ymin*BOHR2ANGSTROM; ymax = ymax*BOHR2ANGSTROM;
	zmin = zmin*BOHR2ANGSTROM; zmax = zmax*BOHR2ANGSTROM;

	// automatic adjustment of nx,ny,nz
	if(nx==0 &&  ny==0 && nz==0){
		nx = ceil((xmax-xmin)*(double)option->outVolumeGrid);
		ny = ceil((ymax-ymin)*(double)option->outVolumeGrid);
		nz = ceil((zmax-zmin)*(double)option->outVolumeGrid);
	}

	// open file for writing
	fd = fopen("volume.xsf","w");
	if(fd==NULL){
		printf("xsfden : Cannot open volume file\n");
		exit(EXIT_FAILURE);
	}

	// sanitize inputs
	if(nx < 2 || ny < 2 || nz < 2){
		printf("Invalid range of grid points\n");
		exit(-1);
	}

	// compute translation displacement
	dx = (xmax-xmin)/(nx-1.);
	dy = (ymax-ymin)/(ny-1.);
	dz = (zmax-zmin)/(nz-1.);

	// place atom at the center of grid
	fprintf(fd, "ATOMS\n");
	for(i=0; i < mol->nAtom; i++)
		fprintf(fd, "%d %20.8E %20.8E %20.8E\n", mol->Z[i], 
		                                         mol->x[i] * BOHR2ANGSTROM, 
		                                         mol->y[i] * BOHR2ANGSTROM, 
		                                         mol->z[i] * BOHR2ANGSTROM);

	// header
	fprintf(fd, "BEGIN_BLOCK_DATAGRID_3D\n");
	fprintf(fd, "datagrids\n");

	fprintf(fd, "BEGIN_DATAGRID_3D_VOLUME\n");

	// number of grid points
	fprintf(fd, "%d %d %d\n", nx, ny, nz);

	// origin of grid
	fprintf(fd, "%20.8E%20.8E%20.8E\n", xmin, ymin, zmin);

	// displacement vector in cartesian
	fprintf(fd, "%20.8E%20.8E%20.8E\n", xmax-xmin, 0.0, 0.0);
	fprintf(fd, "%20.8E%20.8E%20.8E\n", 0.0, ymax-ymin, 0.0);
	fprintf(fd, "%20.8E%20.8E%20.8E\n", 0.0, 0.0, zmax-zmin);

	// allocate memory for points in xy-plane
	if((xyPlane=calloc(ny*nx,sizeof(double)))==NULL){
		printf("xsfden : error cannot allocate memeory\n");
		exit(-1);
	}

	switch(option->outVolumeType){

	case VOLUME_DENSITY_TOTAL:
		for(k=0; k<nz; k++){
			// compute the entire xy-plane at a time
			uhf_rho_XYPlane(nBasis, gto, PA, PB, option->outVolumeCut,
			                       xmin*ANGSTROM2BOHR,
			                       ymin*ANGSTROM2BOHR,
			                (zmin+dz*k)*ANGSTROM2BOHR,
			                         dx*ANGSTROM2BOHR,
			                         dy*ANGSTROM2BOHR, nx, ny, xyPlane);
			for(j=0; j<ny; j++)
			for(i=0; i<nx; i++)
				fprintf(fd, "%20.8E", xyPlane[j*nx+i]);

			// compute each point at a time
			//for(j=0; j<ny; j++)
			//for(i=0; i<nx; i++){
			//	value = uhf_rho(nBasis,gto,PA,PB,(xmin+dx*i)*ANGSTROM2BOHR,
			//	                                 (ymin+dy*j)*ANGSTROM2BOHR,
			//	                                 (zmin+dz*k)*ANGSTROM2BOHR);
			//	fprintf(fd, "%20.8E", value);
			//}
			fprintf(fd, "\n");
		}
	break;

	case VOLUME_POTENTIAL_TOTAL:
		for(k=0; k<nz; k++){
			// compute the entire xy-plane at a time
			uhf_potential_XYPlane(nBasis, gto, mol, PA, PB, option->outVolumeCut,
			                       xmin*ANGSTROM2BOHR,
			                       ymin*ANGSTROM2BOHR,
			                (zmin+dz*k)*ANGSTROM2BOHR,
			                         dx*ANGSTROM2BOHR,
			                         dy*ANGSTROM2BOHR, nx, ny, 0, xyPlane);
			for(j=0; j<ny; j++)
			for(i=0; i<nx; i++)
				fprintf(fd, "%20.8E", xyPlane[j*nx+i]);

			// compute each point at a time
			//for(j=0; j<ny; j++)
			//for(i=0; i<nx; i++){
			//	value = uhf_potential(nBasis,gto,mol,PA,PB,1,(xmin+dx*i)*ANGSTROM2BOHR,
			//				                                 (ymin+dy*j)*ANGSTROM2BOHR,
			//				                                 (zmin+dz*k)*ANGSTROM2BOHR);
			//	fprintf(fd, "%20.8E", value);
			//}
			fprintf(fd, "\n");
		}
	break;

	case VOLUME_MO_ALPHA:
		for(k=0; k<nz; k++){
			for(j=0; j<ny; j++)
			for(i=0; i<nx; i++){
				value = uhf_mo(nBasis,gto,CA,option->outWhichMO,(xmin+dx*i)*ANGSTROM2BOHR,
							                                    (ymin+dy*j)*ANGSTROM2BOHR,
							                                    (zmin+dz*k)*ANGSTROM2BOHR);
				fprintf(fd, "%20.8E", value);
			}
			fprintf(fd, "\n");
		}
	break;

	case VOLUME_MO_BETA:
		for(k=0; k<nz; k++){
			for(j=0; j<ny; j++)
			for(i=0; i<nx; i++){
				value = uhf_mo(nBasis,gto,CB,option->outWhichMO,(xmin+dx*i)*ANGSTROM2BOHR,
							                                    (ymin+dy*j)*ANGSTROM2BOHR,
							                                    (zmin+dz*k)*ANGSTROM2BOHR);
				fprintf(fd, "%20.8E", value);
			}
			fprintf(fd, "\n");
		}
	break;

	case VOLUME_GRADRHO_OVER_RHO43:
	case VOLUME_GRAD_RS:
		for(k=0; k<nz; k++){
			// compute each point at a time
			for(j=0; j<ny; j++)
			for(i=0; i<nx; i++){
				value = uhf_getValue(option->outVolumeType, nBasis,option,gto,mol,PA,PB,(xmin+dx*i)*ANGSTROM2BOHR,
				                                                                        (ymin+dy*j)*ANGSTROM2BOHR,
				                                                                        (zmin+dz*k)*ANGSTROM2BOHR);
				fprintf(fd, "%20.8E", value);
			}
			fprintf(fd, "\n");
		}
	break;

	default:
		printf("xsfden: error volume type not specified\n");
		exit(EXIT_FAILURE);
	}

	// footer
	fprintf(fd, "END_DATAGRID_3D\n");
	fprintf(fd, "END_BLOCK_DATAGRID_3D\n");

	// close file
	fclose(fd);

	// free memory
	free(xyPlane);
	free(PA);
	free(PB);

	printf("[Done]  printing volume in XSF format to 'volume.xsf'\n");
	fflush(stdout);
}


// print out volume infomation in xsf format
// XSF is XCrySDEN format used in molecular 
// and crystal structure. The file format
// is reversed engineered from AbInit cut3d
// output file.
//
// 2008 - Teepanis Chachiyo
// 		Initial implementation
//
// Oct 16, 2025 - Teepanis Chachiyo
//		Debug the value(double,double,double)
//
void xsfvol(double value(double, double, double),
            double xmin, double xmax, 
            double ymin, double ymax, 
            double zmin, double zmax,
            int nx, int ny, int nz){

	int i,j,k;
	double dx, dy, dz;

	// sanitize inputs
	if(nx < 2 || ny < 2 || nz < 2){
		printf("Invalid range of grid points\n");
		exit(-1);
	}

	// compute translation displacement
	dx = (xmax-xmin)/(nx-1.);
	dy = (ymax-ymin)/(ny-1.);
	dz = (zmax-zmin)/(nz-1.);

	// place dummy atom at the center of grid
	printf("ATOMS\n");
	printf("0 %20.8E %20.8E %20.8E\n", (xmax+xmin)/2.0, 
	                                   (ymax+ymin)/2.0, 
	                                   (zmax+zmin)/2.0);

	// header
	printf("BEGIN_BLOCK_DATAGRID_3D\n");
	printf("datagrids\n");
	printf("BEGIN_DATAGRID_3D_DENSITY\n");

	// number of grid points
	printf("%d %d %d\n", nx, ny, nz);

	// origin of grid
	printf("%20.8E%20.8E%20.8E\n", xmin, ymin, zmin);

	// displacement vector in cartesian
	printf("%20.8E%20.8E%20.8E\n", xmax-xmin, 0.0, 0.0);
	printf("%20.8E%20.8E%20.8E\n", 0.0, ymax-ymin, 0.0);
	printf("%20.8E%20.8E%20.8E\n", 0.0, 0.0, zmax-zmin);

	// generate data points	
	for (k=0; k<nz; k++){
		for (j=0; j<ny; j++)
		for (i=0; i<nx; i++)
			printf("%20.8E", value(xmin+dx*i, 
			                       ymin+dy*j, 
			                       zmin+dz*k));
		printf("\n");
	}

	// footer
	printf("END_DATAGRID_3D\n");
	printf("END_BLOCK_DATAGRID_3D\n");
}


// print out density infomation in cube format
// CUBE format is used when Gaussian compute electron
// density and molecular orbitals
//
// Note! xmin,xmax,ymin,ymax,zmin,zmax in *** Bohr ***
//       set them to zero for automatic adjustment
//
//       set nx,ny,nz to zero for auto 
//
// Nov 18, 2017 - Teepanis Chachiyo
//      - Adjust the number of points per Angstrom from the option
//
// Oct 8, 2012 - Teepanis Chachiyo
//      - Change the file name to 'volume.cube'
//      - Adjust the code to coincide with the new struct option_t
//      - Inner loop is now on the XY-plane
//      - Compute the entire XY-plane at a time
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//      Add potential volume information
//
// August 15, 2010 - Nanta Sophonrat
//      Fix bug: lines was not properly terminated during data printing
//      so that it can be sucessfully read using GabEdit
//
// July 11, 2010 - Teepanis Chachiyo
//      Migrate to unrestricted hartree fock
//
// Dec 31, 2009 - Theerapon Khamla
//      Initial implementation
//
void cubeden(int nBasis,              // number of basis function
            struct GTOBasis_t * gto,  // pointer to basis set information
            struct Molecule_t * mol,  // pointer to molecular structure
            double *CA,               // calculated molecular orbital coeff.
            double *CB,               // save but for beta spin
            double xmin, double xmax, // minimum and maximum in x-direction
            double ymin, double ymax, // minimum and maximum in y-direction
            double zmin, double zmax, // minimum and maximum in z-direction
            int nx, int ny, int nz,   // number of grid point
            struct option_t *option){ // print according to options

	int i,j,k;           // loop indexes in x,y,z direction
	double dx, dy, dz;   // step size
	double value;        // value to plot
	FILE *fd;            // file pointer
	double *PA, *PB;     // density matix
	int nEA, nEB;        // number of spin up and spin down electron
	double *xyPlane;     // volume information on xy-plane

	// cube specific variables
	double dx1, dy1, dz1; // first  step vector
	double dx2, dy2, dz2; // second step vector
	double dx3, dy3, dz3; // third  step vector
	double vox, voy, voz; // volume origin

	printf("\n[Begin] printing volume information in CUBE format to volume.cube\n");
	fflush(stdout);

	//
	// compute density matrix
	//

	// memory allocation
	PA = calloc(nBasis*nBasis, sizeof(double));            
	PB = calloc(nBasis*nBasis, sizeof(double));               
	if(PA==NULL || PB==NULL){                                             
		printf("cubeden: error - cannot allocate memory\n");
		exit(EXIT_FAILURE);                                  
	}

	// number of occupy orbital
	nEA = get_nEA(mol,option->multiplicity);
	nEB = get_nEB(mol,option->multiplicity);

	// compute densit matrix
	uhf_getDMatrix(nBasis, nEA, CA, PA);
	uhf_getDMatrix(nBasis, nEB, CB, PB);

	// determine min-max if all argument are zero
	if(xmin==0.0 && xmax==0.0 && 
	   ymin==0.0 && ymax==0.0 && 
	   zmin==0.0 && zmax==0.0){
		xmin = mol->x[0]; xmax = mol->x[0];
		ymin = mol->y[0]; ymax = mol->y[0];
		zmin = mol->z[0]; zmax = mol->z[0];

		// loop to get exact min-max
		for(i=1; i < mol->nAtom; i++){
			if(xmin > mol->x[i]) xmin = mol->x[i];
			if(xmax < mol->x[i]) xmax = mol->x[i];
			if(ymin > mol->y[i]) ymin = mol->y[i];
			if(ymax < mol->y[i]) ymax = mol->y[i];
			if(zmin > mol->z[i]) zmin = mol->z[i];
			if(zmax < mol->z[i]) zmax = mol->z[i];
		}

		// push boundary 4.0 Bohr
		xmin -= 4.0; xmax += 4.0;
		ymin -= 4.0; ymax += 4.0;
		zmin -= 4.0; zmax += 4.0;
	}

	// automatic adjustment of nx,ny,nz
	// 10 points / angstrom
	if(nx==0 &&  ny==0 && nz==0){
		nx = ceil((xmax-xmin)/ANGSTROM2BOHR*(double)option->outVolumeGrid);
		ny = ceil((ymax-ymin)/ANGSTROM2BOHR*(double)option->outVolumeGrid);
		nz = ceil((zmax-zmin)/ANGSTROM2BOHR*(double)option->outVolumeGrid);
	}

	// open file for writing
	fd = fopen("volume.cube","w");
	if(fd==NULL){
		printf("cubeden : error cannot open density file\n");
		exit(EXIT_FAILURE);
	}

	// sanitize inputs
	if(nx < 2 || ny < 2 || nz < 2){
		printf("Invalid range of grid points\n");
		exit(-1);
	}

	// compute translation displacement
	dx = (xmax-xmin)/(nx-1.);
	dy = (ymax-ymin)/(ny-1.);
	dz = (zmax-zmin)/(nz-1.);

	// compute step vector
	dx1 = 0.0; dy1 = 0.0; dz1 =  dz;
	dx2 = 0.0; dy2 =  dy; dz2 = 0.0;
	dx3 =  dx; dy3 = 0.0; dz3 = 0.0;

	// set volume origin
	vox = xmin; voy = ymin; voz = zmin;

	// cube file header
	fprintf(fd,"Siam Quantum (SQ) generated cube file\n");
	fprintf(fd,"data of molecular orbital\n");

	// number of atom and origin of volume in atomic unit
	// {nAtom} {vox} {voy} {voz}
	fprintf(fd," %d %20.8E %20.8E %20.8E\n",-1*mol->nAtom,vox,voy,voz);

	// number of voxel and axis vector
	// {nGrid} {vector x} {voctor y} {vector z}
	fprintf(fd," %d %20.8E %20.8E %20.8E\n",nz,dx1,dy1,dz1);    
	fprintf(fd," %d %20.8E %20.8E %20.8E\n",ny,dx2,dy2,dz2);
	fprintf(fd," %d %20.8E %20.8E %20.8E\n",nx,dx3,dy3,dz3);

	// molecule structure 
	// {Z} {charge} {x} {y} {z}
	for(i=0; i < mol->nAtom ;i++)
		fprintf(fd," %d %20.8E %20.8E %20.8E %20.8E\n",mol->Z[i], 
		                                               mol->Z[i]*1.0,
		                                               mol->x[i], 
		                                               mol->y[i], 
		                                               mol->z[i]);

	// number of volume data and molecular orbital index
	fprintf(fd,"1 %d\n",option->outWhichMO);

	// allocate memory for points in xy-plane
	if((xyPlane=calloc(ny*nx,sizeof(double)))==NULL){
		printf("cubeden : error cannot allocate memeory\n");
		exit(-1);
	}

	// generate data points
	switch(option->outVolumeType){

	case VOLUME_DENSITY_TOTAL:
		for (i=0; i<nz; i++){
			// compute the entire xy-plane at a time
			uhf_rho_XYPlane(nBasis, gto, PA, PB, option->outVolumeCut,
			                vox, voy, voz+dz1*i,
			                dx, dy, nx, ny, xyPlane);
			for (j=0; j<ny; j++){
				for (k=0; k<nx; k++){
					fprintf(fd, "%20.8E", xyPlane[j*nx+k]);
					if (k % 6 == 5) fprintf(fd,"\n");
				}
				fprintf(fd, "\n");
			}

			// compute one point at a time
			//for (j=0; j<ny; j++){
			//	for (k=0; k<nx; k++){
			//		value = uhf_rho(nBasis,gto,PA,PB,
			//		                (vox+dx1*i+dx2*j+dx3*k),
			//		                (voy+dy1*i+dy2*j+dy3*k),
			//		                (voz+dz1*i+dz2*j+dz3*k));
			//		fprintf(fd, "%20.8E", value);
			//		if (k % 6 == 5) fprintf(fd,"\n");
			//	}
			//	fprintf(fd, "\n");
			//}
		}
	break;

	case VOLUME_POTENTIAL_TOTAL:
		for (i=0; i<nz; i++){
			// compute the entire xy-plane at a time
			uhf_potential_XYPlane(nBasis, gto, mol, PA, PB, option->outVolumeCut,
			                      vox, voy, voz+dz1*i,
			                      dx, dy, nx, ny, 0, xyPlane);
			for (j=0; j<ny; j++){
				for (k=0; k<nx; k++){
					fprintf(fd, "%20.8E", xyPlane[j*nx+k]);
					if (k % 6 == 5) fprintf(fd,"\n");
				}
				fprintf(fd, "\n");
			}

			// compute one point at a time
			//for (j=0; j<ny; j++){
			//	for (k=0; k<nx; k++){
			//		value = uhf_potential(nBasis,gto,mol,PA,PB,0,
			//			                (vox+dx1*i+dx2*j+dx3*k),
			//			                (voy+dy1*i+dy2*j+dy3*k),
			//			                (voz+dz1*i+dz2*j+dz3*k));
			//		fprintf(fd, "%20.8E", value);
			//		if (k % 6 == 5) fprintf(fd,"\n");
			//	}
			//	fprintf(fd, "\n");
			//}
		}
	break;

	case VOLUME_MO_ALPHA:
		for (i=0; i<nz; i++)
			for (j=0; j<ny; j++){
				for (k=0; k<nx; k++){
					value = uhf_mo(nBasis,gto,CA,option->outWhichMO-1,
					                (vox+dx1*i+dx2*j+dx3*k),
					                (voy+dy1*i+dy2*j+dy3*k),
					                (vox+dz1*i+dz2*j+dz3*k));
					fprintf(fd, "%20.8E", value);
					if (k % 6 == 5) fprintf(fd,"\n");
				}
				fprintf(fd, "\n");
			}
	break;

	case VOLUME_MO_BETA:
		for (i=0; i<nz; i++)
			for (j=0; j<ny; j++){
				for (k=0; k<nx; k++){
					value = uhf_mo(nBasis,gto,CB,option->outWhichMO-1,
					                (vox+dx1*i+dx2*j+dx3*k),
					                (voy+dy1*i+dy2*j+dy3*k),
					                (vox+dz1*i+dz2*j+dz3*k));
					fprintf(fd, "%20.8E", value);
					if (k % 6 == 5) fprintf(fd,"\n");
				}
				fprintf(fd, "\n");
			}
	break;

	default:
		printf("cubeden - error do not recognize volume type\n");
		exit(-1);
	}

	// close file
	fclose(fd);

	// free memory
	free(xyPlane);
	free(PA);
	free(PB);

	printf("[Done]  printing volume information in CUBE format to volume.cube\n");
	fflush(stdout);
}

//
// export_gaussian: try to mimic output similar to that of Gaussian program
// so that the output can be interpreted by visualization software such as
// GabeEdit.
//
// Sep 7, 2010 - Teepanis Chaciyo
//   Initial implementation and testing
// 
void export_gaussian(int nBasis,      // number of basis functions
            struct GTOBasis_t * gto,  // pointer to basis set information
            struct Molecule_t * mol,  // pointer to molecular structure
            double *CA,               // calculated molecular orbital coeff.
            double *CB,               // save but for beta spin
            double *eA,               // spin up eigen value
            double *eB,               // spin down eigen value
            struct option_t *option){ // print according to options

	int i,n,c;     // generic index
	FILE *fd;      // file pointer
	int nEA, nEB;  // number of spin up and spin down electron

	printf("[Begin] exporting molecular orbital to Gaussian log format\n");
	fflush(stdout);

	// number of occupy orbital
	nEA = get_nEA(mol,option->multiplicity);
	nEB = get_nEB(mol,option->multiplicity);

	// open file for writing
	fd = fopen("gaussian.txt","w");
	if(fd==NULL){
		printf("export_gaussian : error cannot open gaussian log file\n");
		exit(EXIT_FAILURE);
	}

	// output signature of gaussian output format
	fprintf(fd," Entering Link 1\n");

	// output atomic positions
	fprintf(fd,
	"\n"
	"                         Standard orientation:                         \n"
	" ---------------------------------------------------------------------\n"
	" Center     Atomic     Atomic              Coordinates (Angstroms)\n"
	" Number     Number      Type              X           Y           Z\n"
	" ---------------------------------------------------------------------\n");

	for(i=0; i < mol->nAtom ;i++)
		fprintf(fd,"%5d%11d%14d    %12.6f%12.6f%12.6f\n",
		        i+1,mol->Z[i],0,mol->x[i]*BOHR2ANGSTROM
		                       ,mol->y[i]*BOHR2ANGSTROM
		                       ,mol->z[i]*BOHR2ANGSTROM);

	fprintf(fd,
	" ---------------------------------------------------------------------\n");

	// output basis function information
	fprintf(fd, "\n"
	            " AO basis set in the form of general basis input:\n");
	for(i=0; i < mol->nAtom; i++){
		// one entry for each atom
		fprintf(fd, "%3d  0\n", i+1);

		// loop for basis function at this center
		for(n=0; n<nBasis; n++){
			if(     gto[n].x0==mol->x[i] 
			     && gto[n].y0==mol->y[i] 
			     && gto[n].z0==mol->z[i]){

				switch(gto[n].l+gto[n].m+gto[n].n){
				case 0: // handle S type
					fprintf(fd, " S%4d 1.00       0.000000000000\n",
					        gto[n].nContract);
					for(c=0; c<gto[n].nContract; c++)
						fprintf(fd, "    %20.10E%20.10E\n",
						        gto[n].exp[c], gto[n].coef[c]);
				break;

				case 1: // handle P type

					// skip if printed this p-shell
					if( n>0 && isSameShell(nBasis,n,n-1,gto) &&
					   (gto[n-1].l+gto[n-1].m+gto[n-1].n)==1) continue;
					fprintf(fd, " P%4d 1.00       0.000000000000\n",
					        gto[n].nContract);
					for(c=0; c<gto[n].nContract; c++)
						fprintf(fd, "    %20.10E%20.10E\n",
						        gto[n].exp[c], gto[n].coef[c]);
				break;

				case 2: // handle D type

					// skip if printed this p-shell
					if( n>0 && isSameShell(nBasis,n,n-1,gto) &&
					   (gto[n-1].l+gto[n-1].m+gto[n-1].n)==2) continue;
					fprintf(fd, " D%4d 1.00       0.000000000000\n",
					        gto[n].nContract);
					for(c=0; c<gto[n].nContract; c++)
						fprintf(fd, "    %20.10E%20.10E\n",
						        gto[n].exp[c], gto[n].coef[c]);
				break;

				}
			}

		}

		// terminate this entry
		fprintf(fd, " ****\n");
	}

	// output alpha molecular orbital information
	fprintf(fd,"\n" 
	           "     Alpha Molecular Orbital Coefficients\n");
	for(i=0; i<nBasis;){

		fprintf(fd, "                  ");
		for(n=i; n<nBasis && (n-i)<5; n++)
			fprintf(fd, "%10d", n+1);
		fprintf(fd,"\n");

		fprintf(fd, "                     ");
		for(n=i; n<nBasis && (n-i)<5; n++)
			if(n<nEA) fprintf(fd, "%10s", "(A1)--O");
			else      fprintf(fd, "%10s", "(A1)--V");
		fprintf(fd,"\n");

		fprintf(fd, "     EIGENVALUES --  ");
		for(n=i; n<nBasis && (n-i)<5; n++)
			fprintf(fd, "%10.5f",eA[n]);
		fprintf(fd,"\n");

		for(c=0; c<nBasis; c++){
			fprintf(fd, "%4d                 ",c+1);
			for(n=i; n<nBasis && (n-i)<5; n++)
				fprintf(fd, "%10.5f",CA[n*nBasis+c]);
			fprintf(fd,"\n");
		}

		i=n;
	}

	// output beta molecular orbital information
	fprintf(fd,"     Beta Molecular Orbital Coefficients.\n");
	for(i=0; i<nBasis;){

		fprintf(fd, "                  ");
		for(n=i; n<nBasis && (n-i)<5; n++)
			fprintf(fd, "%10d", n+1);
		fprintf(fd,"\n");

		fprintf(fd, "                     ");
		for(n=i; n<nBasis && (n-i)<5; n++)
			if(n<nEB) fprintf(fd, "%10s", "(A1)--O");
			else      fprintf(fd, "%10s", "(A1)--V");
		fprintf(fd,"\n");

		fprintf(fd, "     EIGENVALUES --  ");
		for(n=i; n<nBasis && (n-i)<5; n++)
			fprintf(fd, "%10.5f",eB[n]);
		fprintf(fd,"\n");

		for(c=0; c<nBasis; c++){
			fprintf(fd, "%4d                 ",c+1);
			for(n=i; n<nBasis && (n-i)<5; n++)
				fprintf(fd, "%10.5f",CB[n*nBasis+c]);
			fprintf(fd,"\n");
		}

		i=n;
	}

	// close file
	fclose(fd);

	printf("[Done]  exporting molecular orbital to Gaussian log format\n");
	fflush(stdout);

}

// print out density infomation on a line to 'line.txt'
//
// Jan 4, 2024 - Teepanis Chachiyo
// 		Initial implementation, copied from xsfden
//
void lineden(int nBasis,               // number of basis function
            struct GTOBasis_t * gto,  // pointer to basis set information
            struct Molecule_t * mol,  // pointer to molecular structure
            double *CA,               // calculated molecular orbital coeff.
            double *CB,               // same but for beta spin
            struct option_t *option){ // print according to options

	int i;               // loop index along a line
	double dx, dy, dz;   // displacement vector  [Angstrom]
	double x0, y0, z0;   // centering coordinate [Angstrom]
	double x, y, z;      // current poiont       [Angstrom]
	double ds;           // step size (1 means full displement vector)
	double s;            // current displacement
	double value;        // value to plot
	FILE *fd;            // file pointer
	double *PA;          // density matix
	double *PB;          // density matrix
	int nEA;             // number of occupied alpha electron
	int nEB;             // number of occupied beta electron

	printf("\n[Begin] printing output on a line to 'line.txt'\n");
	fflush(stdout);

	//
	// compute density matrix
	//

	// memory allocation
	PA = calloc(nBasis*nBasis, sizeof(double));
	PB = calloc(nBasis*nBasis, sizeof(double));               
	if(PA==NULL || PB==NULL){                                             
		printf("lineden: Error - Cannot allocate memory\n");
		exit(EXIT_FAILURE);                                  
	}

	// number of occupy orbital
	nEA = get_nEA(mol,option->multiplicity);
	nEB = get_nEB(mol,option->multiplicity);

	// compute densit matrix
	uhf_getDMatrix(nBasis, nEA, CA, PA);
	uhf_getDMatrix(nBasis, nEB, CB, PB);

	// compute displacemnt vector in Angstrom
	dx = option->lineDEF[3] - option->lineDEF[0];
	dy = option->lineDEF[4] - option->lineDEF[1];
	dz = option->lineDEF[5] - option->lineDEF[2];

	// define centering coordinate
	x0 = option->lineDEF[0];
	y0 = option->lineDEF[1];
	z0 = option->lineDEF[2];

	// compute step size
	ds = (option->lineSPAN[1] - option->lineSPAN[0])/(option->lineNP - 1.0);

	// open file for writing
	fd = fopen("line.txt","w");
	if(fd==NULL){
		printf("lineden : Cannot open volume file\n");
		exit(EXIT_FAILURE);
	}

	for(i=0; i < option->lineNP; i++){
		s = option->lineSPAN[0] + i*ds;
		x = x0 + dx * s;
		y = y0 + dy * s;
		z = z0 + dz * s;

		value = uhf_getValue(option->outVolumeType, nBasis,option,gto,mol,PA,PB,x*ANGSTROM2BOHR,
		                                                                        y*ANGSTROM2BOHR,
		                                                                        z*ANGSTROM2BOHR);
		fprintf(fd, "%20.8E %20.8E %20.8E %20.8E\n", x, y, z, value);

	}

	// close file
	fclose(fd);

	// free memory
	free(PA);
	free(PB);

	printf("[Done]  printing output on a line to 'line.txt'\n");
	fflush(stdout);
}


//
// export_gamess: export molecular orbital in $VEC format for GAMESS
//
// Feb 24, 2024 - Teepanis Chaciyo
//   Initial implementation and testing
// 
void export_gamess(int nBasis,      // number of basis functions
            struct GTOBasis_t * gto,  // pointer to basis set information
            struct Molecule_t * mol,  // pointer to molecular structure
            double *CA,               // calculated molecular orbital coeff.
            double *CB,               // save but for beta spin
            double *eA,               // spin up eigen value
            double *eB,               // spin down eigen value
            struct option_t *option){ // print according to options

	int i,n,c;     // generic index
	FILE *fd;      // file pointer
	int nEA, nEB;  // number of spin up and spin down electron

	printf("[Begin] exporting occupied molecular orbital using GAMESS $VEC format\n");
	fflush(stdout);

	// number of occupy orbital
	nEA = get_nEA(mol,option->multiplicity);
	nEB = get_nEB(mol,option->multiplicity);

	// open file for writing
	fd = fopen("gamess.txt","w");
	if(fd==NULL){
		printf("export_gamess : error cannot open gamess output file\n");
		exit(EXIT_FAILURE);
	}

	// output signature of gamess output format
	fprintf(fd," $VEC");

	// output alpha molecular orbital information
	for(i=0; i<nEA; i++){
		for(n=0; n < nBasis; n++){
			if(n%5==0){
				fprintf(fd,"\n");
				fprintf(fd,"%2d %2d", i+1, n/5+1);
			}
			fprintf(fd,"%15.8E",CA[i*nBasis+n]);
		}
	}

	// output beta molecular orbital information
	for(i=0; i<nEB; i++){
		for(n=0; n < nBasis; n++){
			if(n%5==0){
				fprintf(fd,"\n");
				fprintf(fd,"%2d %2d", i+1, n/5+1);
			}
			fprintf(fd,"%15.8E",CB[i*nBasis+n]);
		}
	}

	fprintf(fd, "\n $END\n");

	// close file
	fclose(fd);

	printf("[Done]  exporting occupied molecular orbital using GAMESS $VEC format\n");

	fflush(stdout);

}
