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
#include <stdlib.h>
#include <stdio.h>
#include "lin.h"

//
// LAPACK prototype
//
void dspgvx_(
	long int * ITYPE, char * JOBZ, char * RANGE, char * UPLO, // what to compute
	long int * N,                                             // dimension
	double * AP, double * BP,                                 // input matrices
	double * VL, double * VU, long int * IL, long int * IU,   // eigen range
	double * ABSTOL,                                          // precision
	long int * M,                            // output number of eigen found
	double * W, double * Z, long int * LDZ,  // output eigen value/vectors
	double * WORK,                           // double working space
	long int * IWORK,                        // integer working space
	long int * IFAIL,                        // reported failed eigen
	long int * INFO);                        // exit status

void dspgvd_(
	long int * ITYPE, char * JOBZ, char * UPLO, // what to compute
	long int * N,                               // dimension
	double * AP, double * BP,                   // input matrices
	double * W, double * Z, long int * LDZ,     // output eigen value/vectors
	double * WORK, long int * LWORK,            // double working space
	long int * IWORK, long int * LIWORK,        // integer working space
	long int * INFO);                           // exit status

void dgesv_(long int * N,       // dimension of matrix A
            long int * NRHS,    // number of columns of B
            double   * A,       // input matrix,
            long int * LDA,     // leading dimension of A max(1,N),
            long int * IPIV,    // pivot matrix (output)
            double   * B,       // maxtix dimension LDBxNRHS
            long int * LDB,     // leading dimension of B
            long int * INFO);   // exit status

// gen_sym_eigen : processes eigen vector/value evaluations
// It solves generalize eigen value problem 
// Ac=eSc where A and S is symmetric matrix
//
// Feb 22, 2008 - Teepanis Chachiyo
//     Initial implementation.
//
// Mar  1, 2008 - Teepanis Chachiyo
//     Fix bug the UPLO needs to be 'L'
//
int gen_sym_eigen(
	int nDim,                           // array dimension
	const double *A, const double *S,   // input A and S
	double *e, double *C){              // output value and vector

	int i,j,k;   // looping 

	//
	// LAPACK dspgvd_ variables
	//
	long int ITYPE=1;  // A*x = (lambda)*B*x
	char     JOBZ='V'; // compute eigenvalues and vectors
	char     UPLO='L'; // upper triangle of A and B are stored
	long int N;        // dimension of matrices A and B
	double  *AP;       // upper triangle matrix A
	double  *BP;       // upper triangle matrix B
	double  *W;        // output eigen value
	double  *Z;        // output eigen vectors
	long int LDZ;      // dimension of matrix
	double  *WORK;     // working array size LWORK
	long int LWORK;    // should be >= 1 + 6N + 2N^2
	long int *IWORK;   // array dimension LIWORK
	long int LIWORK;   // shouble be >= 3 + 5N
	long int INFO;     // exit status

	// sanity check
	if(nDim < 1){
		printf("gen_sym_eigen: Error - Invalid dimension range\n");
		exit(-1);
	}

	// allocate memory
	N      = nDim;
	LDZ    = nDim;
	AP     = calloc(N*(N+1)/2, sizeof(double));
	BP     = calloc(N*(N+1)/2, sizeof(double));
	W      = calloc(N,         sizeof(double));
	Z      = calloc(LDZ*N,     sizeof(double));

	// computation of LWORK is from DSPGVD code
	i      = (int)(log(N)/log(2.0));
	if(pow(2.0,i) < (double)N) i++;
	if(pow(2.0,i) < (double)N) i++;
	LWORK  = 1 + 5*N + 2*N*N + 2*N*i; 

	WORK   = calloc(LWORK,     sizeof(double));
	LIWORK = 3 + 5*N;
	IWORK  = calloc(LIWORK,    sizeof(long int));
	INFO   = 0;

	// read lower triangle array A
	// Important!!!!
	// According to LAPACK document, it should be
	// upper. But from emprical experiment, you really
	// need to enter the matrix as lower triangle
	// - Teepanis Chachiyo Dec 23, 2007
	k = 0;
	for(i=0; i < N; i++)
		for(j=i; j < N; j++){
			AP[k] = A[i*N+j];
			k++;
		}

	// read lower triangle array B
	k = 0;
	for(i=0; i < N; i++)
		for(j=i; j < N; j++){
			BP[k] = S[i*N+j];
			k++;
		}

	dspgvd_(&ITYPE, &JOBZ, &UPLO, &N, 
	        AP, BP,
	        W, Z, &LDZ,
	        WORK, &LWORK, IWORK, &LIWORK,
	        &INFO);

	if(INFO!=0){
		printf("gen_sym_eigen: Error - LAPACK subroutine returns error %ld\n"
		       , INFO);
		exit(EXIT_FAILURE);
	}	
	// print eigen values
	//for(i=0; i < N; i++){
	//	fprintf(outFile, "%20.8lE", W[i]);
	//}
	//fprintf(outFile, "\n");
	for(i=0; i < N; i++)
		e[i] = W[i];

	// print eigen vector per row
	//for(i=0; i < N; i++){
	//	for(j=0; j < N; j++){
	//		fprintf(outFile, "%20.8lE", Z[i*N+j]);
	//	}
	//	fprintf(outFile, "\n");
	//}
	for(i=0; i < N; i++)
		for(j=0; j < N; j++)
			C[i*N+j] = Z[i*N+j];

	// clean up memory
	free(AP);
	free(BP);
	free(W);
	free(Z);
	free(WORK);
	free(IWORK);

	// return TRUE on success
	return 1;
}


//
// gen_sym_SomeEigen: compute some of eigen vector and eigen values of
// symmetrix matrix. It solves generalize eigen value problem 
// Ac=eSc where A and S is symmetric matrix
//
//  Dec 17, 2015 - Teepanis Chachiyo
//		Check if number of eigen is zero
//
//  Jan 16, 2014 - Teepanis Chachiyo
//		Bug fix, M needs to be set because sometimes it is not set by 
//		the subroutine
//
//	Jan 15, 2014 - Teepanis Chachiyo
//		Adpated from gen_sym_eigen
//
int gen_sym_SomeEigen(
	int nDim,                           // array dimension
	int nEigen,                         // number of eigen value requested
	const double *A, const double *S,   // input A and S
	double *e, double *C){              // output value and vector

	int i,j,k;   // looping 

	// validate
	if(nEigen==0) return 1;

	//
	// LAPACK dspgvx_ variables
	//
	long int ITYPE=1;     // A*x = (lambda)*B*x
	char     JOBZ='V';    // compute eigenvalues and vectors
	char     RANGE='I';   // compute range of eigen value within index range
	char     UPLO='L';    // upper triangle of A and B are stored
	long int N=nDim;      // dimension of matrices A and B
	double  *AP;          // upper triangle matrix A
	double  *BP;          // upper triangle matrix B
	double VL, VU;        // lower and upper bound of eigenvalue
	long int IL=1;        // lower bound of eigen index
	long int IU=nEigen;   // upper bound of eigen index
	double ABSTOL=1.0E-8; // accuracy
	long int M=nDim;      // number of eigen returned
	double  *W;           // output eigen value
	double  *Z;           // output eigen vectors
	long int LDZ=nDim;    // dimension of matrix
	double  *WORK;        // working array size 8*N
	long int *IWORK;      // array dimension 5*N
	long int *IFAIL;      // array dimension N
	long int INFO=0;      // exit status

	// sanity check
	if(nDim < 1){
		printf("gen_sym_eigen: Error - Invalid dimension range\n");
		exit(-1);
	}

	// allocate memory
	AP     = calloc(N*(N+1)/2, sizeof(double));
	BP     = calloc(N*(N+1)/2, sizeof(double));
	W      = calloc(N,         sizeof(double));
	Z      = calloc(LDZ*N,     sizeof(double));
	WORK   = calloc(8*N,       sizeof(double));
	IWORK  = calloc(5*N,       sizeof(long int));
	IFAIL  = calloc(N,         sizeof(long int));

	// read lower triangle array A
	// Important!!!!
	// According to LAPACK document, it should be
	// upper. But from emprical experiment, you really
	// need to enter the matrix as lower triangle
	// - Teepanis Chachiyo Dec 23, 2007
	k = 0;
	for(i=0; i < N; i++)
		for(j=i; j < N; j++){
			AP[k] = A[i*N+j];
			k++;
		}

	// read lower triangle array B
	k = 0;
	for(i=0; i < N; i++)
		for(j=i; j < N; j++){
			BP[k] = S[i*N+j];
			k++;
		}

	dspgvx_(&ITYPE, &JOBZ, &RANGE, &UPLO, &N, 
	        AP, BP,
	        &VL, &VU, &IL, &IU,
	        &ABSTOL,
	        &M,
	        W, Z, &LDZ,
	        WORK, IWORK, IFAIL,
	        &INFO);

	// 
	// Linux version of lapack seems to have problem with this if I do not set
	// the value of M prior to calling the function. But the WindowsXP version
	// works fine. (Teepanis Jan 16, 2014)
	if(M != nEigen){
		printf("gen_sym_SomeEigen: Error - LAPACK subroutine failed to converge (M=%ld)\n", M);
		exit(EXIT_FAILURE);
	}

	if(INFO!=0){
		printf("gen_sym_SomeEigen: Error - LAPACK subroutine returns error %ld\n"
		       , INFO);
		exit(EXIT_FAILURE);
	}	

	for(i=0; i < N; i++)
		e[i] = W[i];

	for(i=0; i < N; i++)
		for(j=0; j < N; j++)
			C[i*N+j] = Z[i*N+j];

	// clean up memory
	free(AP);
	free(BP);
	free(W);
	free(Z);
	free(WORK);
	free(IWORK);
	free(IFAIL);

	// return TRUE on success
	return 1;
}

//
// linear_solver: solves simple linear equation AX=B , where
// A is a real matrix but need not symmetric. X and B are matrices
// of dimension nRow x nCol 
//
// Oct 4, 2012 - Teepanis Chachiyo
//                Copy code from Siam Solver
// 
// Nov 16, 2009 - Teepanis Chachiyo
//                Copy code from my previous implementation.
//
int linear_solver(int nRow,        // number of row
                  int nCol,        // number column
                  const double *A, // input matrix A
                  double *B){      // input matrix B, output matrix X
	
	int i,j; 

	//
	// LAPACK dgesv_ variables 
	//
	double *LA;    // LAPACK routine matrix
	double *LB;    // LAPACK routine matrix
	long int N;    // number of rows
	long int NRHS; // number of columns
	long int LDA;  // number of rows
	long int *IPIV;// internal usage
	long int LDB;  // number of rows
	long int INFO; // return code

	// validate range
	if(nRow < 1 || nCol < 1){
		printf("linear_solver: Error - invalid dimension range\n");
		exit(-1);
	}
	
	// setup LAPACK variables
	N    = nRow;
	NRHS = nCol;
	LDA  = nRow;
	LDB  = nRow;
	INFO = 0;
	
	// allocate memory
	LA   = calloc(nRow*nRow, sizeof(double));
	LB   = calloc(nCol*nRow, sizeof(double));
	IPIV = calloc(nRow, sizeof(long int));
	if(IPIV==NULL || LA==NULL || LB==NULL){
		printf("linear_solver: Error - cannnot allocate memeory\n");
		exit(-1);
	}

	// read matrix values using transpose due to LAPACK's FORTRAN notation
	for(i=0; i < nRow; i++)
	for(j=0; j < nRow; j++)
		LA[i*nRow+j] = A[j*nRow+i];

	for(i=0; i < nCol; i++)
	for(j=0; j < nRow; j++)
		LB[i*nRow+j] = B[j*nCol+i];
	
	// solve the problem
	dgesv_( &N, &NRHS, LA, &LDA, IPIV, LB, &LDB, &INFO);

	// check for error
	if(INFO != 0){
		printf("linear_solver: Error - LAPACK returns error code %ld\n", INFO);
		exit(-1);
	}

	// output
	for(i=0; i < nRow; i++)
	for(j=0; j < nCol; j++)
 		B[i*nCol+j] = LB[j*nRow+i];
	
	// clean memory
	free(IPIV);
	free(LA);
	free(LB);
	
	// return TRUE on success
	return 1;
}
