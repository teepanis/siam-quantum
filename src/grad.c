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

#include "basis.h"
#include "matrix.h"
#include "lin.h"
#include "int.h"
#include "util.h"
#include "option.h"
#include "uhf.h"
#include "grad.h"
#include "rpc.h"
#include "dft.h"
#include "grid.h"

// GradVnn : compute contribution to the gradient of energy due to
// nuclei-nuclei repulsion energy. 
//
// Dec 28, 2009 - Theerapon Khamla
//   Initial implementation
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//   Added to Siam Quantum source code
//   Chage loop structure
//
void GradVnn(const struct Molecule_t *mol, double *Gx, double *Gy, double *Gz){
	double r;          // distance between two nuclei
	int A,B;           // atom index

	for(A=0; A < mol->nAtom; A++)
	for(B=0; B < A; B++){

		// compute distance
		r = (mol->x[A] - mol->x[B]) * (mol->x[A] - mol->x[B])
		   +(mol->y[A] - mol->y[B]) * (mol->y[A] - mol->y[B])
		   +(mol->z[A] - mol->z[B]) * (mol->z[A] - mol->z[B]);
		r = sqrt(r);

		// gradient in x direction
		Gx[A] -= mol->Z[A]*mol->Z[B]*(mol->x[A] - mol->x[B])/(r*r*r);
		Gx[B] += mol->Z[A]*mol->Z[B]*(mol->x[A] - mol->x[B])/(r*r*r);
		// gradient in y direction
		Gy[A] -= mol->Z[A]*mol->Z[B]*(mol->y[A] - mol->y[B])/(r*r*r);
		Gy[B] += mol->Z[A]*mol->Z[B]*(mol->y[A] - mol->y[B])/(r*r*r);
		// gradient in z directin
		Gz[A] -= mol->Z[A]*mol->Z[B]*(mol->z[A] - mol->z[B])/(r*r*r);
		Gz[B] += mol->Z[A]*mol->Z[B]*(mol->z[A] - mol->z[B])/(r*r*r);
	}
}

//
// GTO_grad_overlap: compute overlap integral
//
//  /  
//  |     d
//  | dV  --  Chi_i * Chi_j
//  |     dXi
//  /
//
// where Chi_i and Chi_j are basis functions; and Xi is the
// center of the basis Chi_i
//
// Jan 4, 2010 - Theerapon Khamla
//    Initial version
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Modified loop structure and added to Siam Quantum source code
// 
#define GRAD_X 0
#define GRAD_Y 1
#define GRAD_Z 2
double GTO_grad_overlap(int i,                          // ith basis
                        int j,                          // jth basis
                        const struct GTOBasis_t *gto,   // basis database
                        int direction){                 // specified direction

	int iCnt, jCnt;   // contracted function 
	double sum;       // integral sum
	int l,m,n;        // direction index
	double s=0.0;

	// determine requsted direction
	switch(direction){
	case GRAD_X: l=1;m=0;n=0; break;
	case GRAD_Y: l=0;m=1;n=0; break;
	case GRAD_Z: l=0;m=0;n=1; break;
	default:
		printf("GTO_grad_overlap: error unrecognized direction flag\n");
		exit(-1);
	}

	// note:
	// The term (l*gto[i].l+m*gto[i].m+n*gto[i].n) will select out only
	// the angular order in specified direction.
	// Teepanis - Oct 20, 2010

#define OVERLAP_SUM(L,M,N) gto[i].coef[iCnt] * gto[j].coef[jCnt] *           \
			               gto[i].norm[iCnt] * gto[j].norm[jCnt] *           \
		        overlap(gto[i].exp[iCnt], gto[i].l L,gto[i].m M, gto[i].n N, \
		                                  gto[i].x0 ,gto[i].y0 , gto[i].z0,  \
		                gto[j].exp[jCnt], gto[j].l  ,gto[j].m  , gto[j].n,   \
			                              gto[j].x0, gto[j].y0, gto[j].z0);

	///////////////////////////
	// handle lowering angular
	///////////////////////////
	sum = 0.0;
	if((l*gto[i].l+m*gto[i].m+n*gto[i].n)>0){
		// looping over contracted functions
		for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
		for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
			sum += OVERLAP_SUM(-l,-m,-n);

		sum = sum*(l*gto[i].l+m*gto[i].m+n*gto[i].n);
	}
	s = -sum;

	///////////////////////////
	// handle raising angular
	///////////////////////////
	sum = 0.0;
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
		sum += gto[i].exp[iCnt] * OVERLAP_SUM(+l,+m,+n);

	s = s + sum + sum;

#undef OVERLAP_SUM

	return s;
}


//
// GTO_grad_moment: compute moment integral
//
//  /  
//  |     d              mx     my     mz
//  | dV  --  Chi_i (x-cx) (y-cy) (z-cz) Chi_j
//  |     dXi
//  /
//
// where Chi_i and Chi_j are basis functions; and Xi is the
// center of the basis Chi_i
//
// March 14, 2016 - Teepanis Chachiyo
//    Initial version
// 
double GTO_grad_moment(int i,                          // ith basis
                       int j,                          // jth basis
                       const struct GTOBasis_t *gto,   // basis database
                       int mx, int my, int mz,         // moment order
                       double cx, double cy, double cz,// center of moment
                       int direction){                 // specified direction

	int iCnt, jCnt;   // contracted function 
	double sum;       // integral sum
	int l,m,n;        // direction index
	double s=0.0;

	// determine requsted direction
	switch(direction){
	case GRAD_X: l=1;m=0;n=0; break;
	case GRAD_Y: l=0;m=1;n=0; break;
	case GRAD_Z: l=0;m=0;n=1; break;
	default:
		printf("GTO_grad_moment: error unrecognized direction flag\n");
		exit(-1);
	}

	// note:
	// The term (l*gto[i].l+m*gto[i].m+n*gto[i].n) will select out only
	// the angular order in specified direction.
	// Teepanis - Oct 20, 2010

#define MOMENT_SUM(L,M,N)  gto[i].coef[iCnt] * gto[j].coef[jCnt] *           \
			               gto[i].norm[iCnt] * gto[j].norm[jCnt] *           \
		         moment(gto[i].exp[iCnt], gto[i].l L,gto[i].m M, gto[i].n N, \
		                                  gto[i].x0 ,gto[i].y0 , gto[i].z0,  \
		                gto[j].exp[jCnt], gto[j].l  ,gto[j].m  , gto[j].n,   \
			                              gto[j].x0, gto[j].y0 , gto[j].z0,  \
		                                  mx,my,mz,cx,cy,cz);

	///////////////////////////
	// handle lowering angular
	///////////////////////////
	sum = 0.0;
	if((l*gto[i].l+m*gto[i].m+n*gto[i].n)>0){
		// looping over contracted functions
		for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
		for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
			sum += MOMENT_SUM(-l,-m,-n);

		sum = sum*(l*gto[i].l+m*gto[i].m+n*gto[i].n);
	}
	s = -sum;

	///////////////////////////
	// handle raising angular
	///////////////////////////
	sum = 0.0;
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
		sum += gto[i].exp[iCnt] * MOMENT_SUM(+l,+m,+n);

	s = s + sum + sum;

#undef MOMENT_SUM

	return s;
}

//
// GTO_grad_kinetic: compute kinetic integral
//
//  /  
//  |     d            1  _2
//  | dV  --  Chi_i (- - \/ )Chi_j
//  |     dXi          2
//  /
//
// where Chi_i and Chi_j are basis functions; and Xi is the
// center of the basis Chi_i
//
// Jan 25, 2010 - Theerapon Khamla
//    Initial version
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Modified loop structure and added to Siam Quantum source code
// 
double GTO_grad_kinetic(int i,                          // ith basis
                        int j,                          // jth basis
                        const struct GTOBasis_t *gto,   // basis database
                        int direction){                 // specified direction

	int iCnt, jCnt;  // contracted function 
	double sum=0.0;  // integral sum
	int l,m,n;       // direction index
	double k;

	// determine requsted direction
	switch(direction){
	case GRAD_X: l=1;m=0;n=0; break;
	case GRAD_Y: l=0;m=1;n=0; break;
	case GRAD_Z: l=0;m=0;n=1; break;
	default:
		printf("GTO_grad_kinetic: error unrecognized direction flag\n");
		exit(-1);
	}

#define KINETIC_SUM(L,M,N)  gto[i].coef[iCnt]  * gto[j].coef[jCnt] *         \
			                gto[i].norm[iCnt]  * gto[j].norm[jCnt] *         \
		        kinetic(gto[i].exp[iCnt], gto[i].l L,gto[i].m M,gto[i].n N,  \
		                                  gto[i].x0, gto[i].y0, gto[i].z0,   \
		                gto[j].exp[jCnt], gto[j].l , gto[j].m,  gto[j].n,    \
			                              gto[j].x0, gto[j].y0, gto[j].z0);

	////////////////////////////
	// handle lowering angular 
	////////////////////////////
	if((l*gto[i].l+m*gto[i].m+n*gto[i].n)>0){
		for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
		for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
			sum += KINETIC_SUM(-l,-m,-n);

		sum = sum*(l*gto[i].l+m*gto[i].m+n*gto[i].n);
	}
	k = -sum;

	///////////////////////////
	// handle raising angular 
	///////////////////////////
	sum = 0.0;
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
		sum += gto[i].exp[iCnt] * KINETIC_SUM(+l,+m,+n);

	k = k + sum + sum;

#undef KINETIC_SUM

	return k;
}


//
// GTO_grad_nai: compute nuclear attraction integral
//
//  /  
//  |     d           ->
//  | dV  --  Chi_i V(r) Chi_j
//  |     dXi        
//  /
//                                              ->
// where Chi_i and Chi_j are basis functions; V(r) is the Coulomb potential
// due to nuclei; and Xi is the center of the basis Chi_i
//
// March 8, 2010 - Theerapon Khamla
//    Initial version
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Modified loop structure and added to Siam Quantum source code
// 
double GTO_grad_nai(int i,                          // ith basis
                    int j,                          // jth basis
                    const struct GTOBasis_t *gto,   // basis database
                    const struct Molecule_t *mol,   // molecule database
                    int direction){                 // specified direction

	int k;           // nucleus index
	int iCnt, jCnt;  // contracted functions
	double sum=0.0;  // integral sum
	int l,m,n;       // direction index
	double dnai;

	// determine requsted direction
	switch(direction){
	case GRAD_X: l=1;m=0;n=0; break;
	case GRAD_Y: l=0;m=1;n=0; break;
	case GRAD_Z: l=0;m=0;n=1; break;
	default:
		printf("GTO_grad_nai: error unrecognized direction flag\n");
		exit(-1);
	}

#define NAI_SUM(L,M,N)       gto[i].coef[iCnt]  * gto[j].coef[jCnt] *  \
				             gto[i].norm[iCnt]  * gto[j].norm[jCnt] *  \
				      (double)mol->Z[k]                             *  \
				        nai(gto[i].x0, gto[i].y0, gto[i].z0,           \
				            1.0,                                       \
				            gto[i].l L,  gto[i].m M,  gto[i].n N,      \
				            gto[i].exp[iCnt],                          \
				            gto[j].x0, gto[j].y0, gto[j].z0,           \
				            1.0,                                       \
				            gto[j].l  ,  gto[j].m,    gto[j].n,        \
				            gto[j].exp[jCnt],                          \
				            mol->x[k], mol->y[k], mol->z[k]);

	//////////////////////
	// handle lowering l
	//////////////////////
	if((l*gto[i].l+m*gto[i].m+n*gto[i].n)>0){
		// looping over contracted functions
		for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
		for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
				
			// looping over nuclei
			for(k=0; k < mol->nAtom; k++)
				sum += NAI_SUM(-l,-m,-n);

		sum = sum*(l*gto[i].l+m*gto[i].m+n*gto[i].n);	
	}
	dnai = -sum;

	//////////////////////
	// handle raising l
	//////////////////////
	sum = 0.0;
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++)
			
		// looping over nuclei
		for(k=0; k < mol->nAtom; k++)
			sum += gto[i].exp[iCnt] * NAI_SUM(+l,+m,+n);

	dnai = dnai + sum + sum;

#undef NAI_SUM

	return dnai;
}


//
// contr_grad_eri: compute two electron integral
//
//  /  
//  |  3   3    d                       1
//  | dr1 dr2  -- Chi_a(r1) Chi_b(r1)  ---   Chi_c(r2) Chi_d(r2) 
//  |          dXa                   |r1-r2|
//  /
//
// where Chi_a, Chi_b, Chi_c, Chi_d are basis functions; and Xa is the center 
// of the basis Chi_a
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Initial version
// 
double contr_grad_eri(
	int lena,double *aexps,double *acoefs,double *anorms, // Chi_a info
	double xa,double ya,double za,int la,int ma,int na,   //
	int lenb,double *bexps,double *bcoefs,double *bnorms, // Chi_b info
	double xb,double yb,double zb,int lb,int mb,int nb,   //
	int lenc,double *cexps,double *ccoefs,double *cnorms, // Chi_c info
	double xc,double yc,double zc,int lc,int mc,int nc,   //
	int lend,double *dexps,double *dcoefs,double *dnorms, // Chi_d info
	double xd,double yd,double zd,int ld,int md,int nd,   //
	int direction){                                       // specified direction

	double val;
	int p,q,r,s;  // contraction counter
	int l,m,n;    // direction index
	double EE;

	// determine requsted direction
	switch(direction){
	case GRAD_X: l=1;m=0;n=0; break;
	case GRAD_Y: l=0;m=1;n=0; break;
	case GRAD_Z: l=0;m=0;n=1; break;
	default:
		printf("contr_grad_eri: error unrecognized direction flag\n");
		exit(-1);
	}

#define ERISUM(L,M,N)  acoefs[p]*bcoefs[q]*ccoefs[r]*dcoefs[s]*       \
			        eri(xa,ya,za,anorms[p],la L,ma M,na N,aexps[p],   \
			            xb,yb,zb,bnorms[q],lb  ,mb  ,nb  ,bexps[q],   \
			            xc,yc,zc,cnorms[r],lc  ,mc  ,nc  ,cexps[r],   \
			            xd,yd,zd,dnorms[s],ld  ,md  ,nd  ,dexps[s]);

	///////////////////////////////
	// handling lowering la term
	///////////////////////////////
	val = 0.0;
	if((l*la+m*ma+n*na)>0){
		// proceed from lowest exponent value
		for (p=lena-1; p>=0; p--)
		for (q=lenb-1; q>=0; q--)
		for (r=lenc-1; r>=0; r--)
		for (s=lend-1; s>=0; s--)
			val += ERISUM(-l,-m,-n);
		val = val*(l*la+m*ma+n*na);
	}
	EE = -val;

	//////////////////////////////
	// handling raising la term
	//////////////////////////////
	val = 0.0;
	// proceed from lowest exponent value
	for (p=lena-1; p>=0; p--)
	for (q=lenb-1; q>=0; q--)
	for (r=lenc-1; r>=0; r--)
	for (s=lend-1; s>=0; s--)
		val += aexps[p] * ERISUM(+l,+m,+n);
	EE = EE + val + val;

#undef ERISUM

	return EE;
}


//
// contr_grad_grad_eri: compute two electron integral of the form
//
//  /  
//  |  3   3    d                       1     d
//  | dr1 dr2  -- Chi_a(r1) Chi_b(r1)  ---   -- Chi_c(r2) Chi_d(r2) 
//  |          dXa                   |r1-r2| dXc
//  /
//
// where Chi_a, Chi_b, Chi_c, Chi_d are basis functions; and Xa is the center 
// of the basis Chi_a, whereas Xc is the center of the basis Chi_c
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Initial version
// 
double contr_grad_grad_eri(
	int lena,double *aexps,double *acoefs,double *anorms, // Chi_a info
	double xa,double ya,double za,int la,int ma,int na,   //
	int lenb,double *bexps,double *bcoefs,double *bnorms, // Chi_b info
	double xb,double yb,double zb,int lb,int mb,int nb,   //
	int lenc,double *cexps,double *ccoefs,double *cnorms, // Chi_c info
	double xc,double yc,double zc,int lc,int mc,int nc,   //
	int lend,double *dexps,double *dcoefs,double *dnorms, // Chi_d info
	double xd,double yd,double zd,int ld,int md,int nd,   //
	int direction){                                       // specified direction

	double val;
	int p,q,r,s;  // contraction counter
	int lA,mA,nA; // direction index for Chi_a
	int lC,mC,nC; // direction index for Chi_c
	double EE=0.0;

	// determine requsted direction
	switch(direction){
	case GRAD_X: lA=1;mA=0;nA=0; lC=1;mC=0;nC=0; break;
	case GRAD_Y: lA=0;mA=1;nA=0; lC=0;mC=1;nC=0; break;
	case GRAD_Z: lA=0;mA=0;nA=1; lC=0;mC=0;nC=1; break;
	default:
		printf("contr_grad_grad_eri: error unrecognized direction flag\n");
		exit(-1);
	}

#define ERISUM(LA,MA,NA,LC,MC,NC)                                        \
                       acoefs[p]*bcoefs[q]*ccoefs[r]*dcoefs[s]*          \
			        eri(xa,ya,za,anorms[p],la LA,ma MA,na NA,aexps[p],   \
			            xb,yb,zb,bnorms[q],lb  ,mb  ,nb  ,bexps[q],      \
			            xc,yc,zc,cnorms[r],lc LC,mc MC,nc NC,cexps[r],   \
			            xd,yd,zd,dnorms[s],ld  ,md  ,nd  ,dexps[s]);

	///////////////////////////////
	// handling lowering la term
	///////////////////////////////
	if((la*lA+ma*mA+na*nA)>0){

		//////////////////////////////
		// handling lowering lc term
		//////////////////////////////
		val = 0.0;
		if((lc*lC+mc*mC+nc*nC)>0)
		// proceed from lowest exponent value
		for (p=lena-1; p>=0; p--)
		for (q=lenb-1; q>=0; q--)
		for (r=lenc-1; r>=0; r--)
		for (s=lend-1; s>=0; s--)
			val += ERISUM(-lA,-mA,-nA,-lC,-mC,-nC);

		val = val*(la*lA+ma*mA+na*nA)*(lc*lC+mc*mC+nc*nC);
		EE  = EE+val;

		/////////////////////////////
		// handling raising lc term
		/////////////////////////////
		val = 0.0;
		for (p=lena-1; p>=0; p--)
		for (q=lenb-1; q>=0; q--)
		for (r=lenc-1; r>=0; r--)
		for (s=lend-1; s>=0; s--)
			val += cexps[r] * ERISUM(-lA,-mA,-nA,+lC,+mC,+nC);

		val = val*(la*lA+ma*mA+na*nA)*2.0;
		EE  = EE-val;

	}

	//////////////////////////////
	// handling raising la term
	//////////////////////////////

	////////////////////////////
	// handle lowering lc term
	////////////////////////////
	val = 0.0;
	if((lc*lC+mc*mC+nc*nC)>0)
	for (p=lena-1; p>=0; p--)
	for (q=lenb-1; q>=0; q--)
	for (r=lenc-1; r>=0; r--)
	for (s=lend-1; s>=0; s--)
		val += aexps[p] * ERISUM(+lA,+mA,+nA,-lC,-mC,-nC);

	val = val*(lc*lC+mc*mC+nc*nC)*2.0;
	EE  = EE-val;

	//////////////////////////////
	// handling raising lc term
	////////////////////////////
	val = 0.0;
	for (p=lena-1; p>=0; p--)
	for (q=lenb-1; q>=0; q--)
	for (r=lenc-1; r>=0; r--)
	for (s=lend-1; s>=0; s--)
		val += aexps[p] * cexps[r] * ERISUM(+lA,+mA,+nA,+lC,+mC,+nC);

	val = val*4.0;
	EE  = EE+val;

#undef ERISUM

	return EE;
}


// create_grad_Schwarz: generates matrix for Schwarz screening
// The matrix element (p,q) in this matrix is defined by
//
//  /  
//  |  3   3    d                       1     d
//  | dr1 dr2  -- Chi_p(r1) Chi_q(r1)  ---   -- Chi_p(r2) Chi_q(r2) 
//  |          dXp                   |r1-r2| dXp
//  /
//
// Oct 18, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Initial implementatoin
//
double * create_grad_Schwarz(
	int nBasis,                     // number of basis functions
	const struct GTOBasis_t *gto,   // basis set info
	int direction){                 // direction flag

	int p,q;            // loop index
	double *ddx_sch;    // Schwarz matrix pointer
	double upBound;     // upper bound from Schwarz inequality

	// allocate memory
	ddx_sch = calloc(nBasis*nBasis, sizeof(double));
	if(ddx_sch==NULL){
		printf("create_grad_Scharz: Cannot allocate Schwarz matrix\n");
		exit(-1);
	}

	//
	// surprisingly I can't seem to assume that ddx_sch(p,q) = ddx_sch(q,p)
	// Teepanis Chachiyo October 15, 2010
	//
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		// evaluate integral
		upBound = contr_grad_grad_eri(
		             gto[p].nContract,
		             gto[p].exp, gto[p].coef, gto[p].norm,
		             gto[p].x0,        gto[p].y0,  gto[p].z0,
		             gto[p].l,         gto[p].m,   gto[p].n,
		             gto[q].nContract,
		             gto[q].exp, gto[q].coef, gto[q].norm,
		             gto[q].x0,        gto[q].y0,  gto[q].z0,
		             gto[q].l,         gto[q].m,   gto[q].n,
		             gto[p].nContract,
		             gto[p].exp, gto[p].coef, gto[p].norm,
		             gto[p].x0,        gto[p].y0,  gto[p].z0,
		             gto[p].l,         gto[p].m,   gto[p].n,
		             gto[q].nContract,
		             gto[q].exp, gto[q].coef, gto[q].norm,	
		             gto[q].x0,        gto[q].y0,  gto[q].z0,
		             gto[q].l,         gto[q].m,   gto[q].n,
		             direction);

		// make sure we have positive value
		ddx_sch[p*nBasis+q] = sqrt(fabs(upBound));
	}

	return ddx_sch;
}


// GradEri: compute contribution of gradition due to electron-electron integral
// In unrestricted hartree-fock formalism, this is equal to
//
//   1                           d          
// + - SUM[pqrs] PT(p,q)*PT(r,s) - (pq|rs)   
//   2                           dXi        
//
//   1                                                 d
// - - SUM[pqrs] ( PT(p,q)*PT(r,s) + PS(p,q)*PS(r,s) ) -  (pr|qs)
//   4                                                 dXi
//
// Note: Gx,Gy,Gz will not be replaced, but "increased" according to the
// contribution calculated.
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Initial implementation
//
void GradEri(
	int nBasis,                     // number of basis function
	const struct GTOBasis_t *gto,   // basis function database
	const int *basis2Atom,          // basis to atom mapping
	const double *PT,               // total density matrix
	const double *PS,               // spin density matrix
	double *Gx,                     // returned gradient in x direction
	double *Gy,                     // returned gradient in y direction
	double *Gz){                    // returned gradient in z direction

	int p,q,i,j;           // basis function index
	double PTPT,PTPTnPSPS; // product of density matrix
	double r;              // generic real number

#define GRAD_ERI_CALL(P,Q,I,J,X) contr_grad_eri(                           \
                                 gto[P].nContract,                         \
                                 gto[P].exp, gto[P].coef, gto[P].norm,     \
                                 gto[P].x0,        gto[P].y0,  gto[P].z0,  \
                                 gto[P].l,         gto[P].m,   gto[P].n,   \
                                 gto[Q].nContract,                         \
                                 gto[Q].exp, gto[Q].coef, gto[Q].norm,     \
                                 gto[Q].x0,        gto[Q].y0,  gto[Q].z0,  \
                                 gto[Q].l,         gto[Q].m,   gto[Q].n,   \
                                 gto[I].nContract,                         \
                                 gto[I].exp, gto[I].coef, gto[I].norm,     \
                                 gto[I].x0,        gto[I].y0,  gto[I].z0,  \
                                 gto[I].l,         gto[I].m,   gto[I].n,   \
                                 gto[J].nContract,                         \
                                 gto[J].exp, gto[J].coef, gto[J].norm,     \
                                 gto[J].x0,        gto[J].y0,  gto[J].z0,  \
                                 gto[J].l,         gto[J].m,   gto[J].n,   \
                                 X);

	// loop thru all permutation of basis
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++) 
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){

		// compute coulomb and exchange multiplier
		PTPT      = 0.50 *   PT[p*nBasis+q]*PT[i*nBasis+j];
		PTPTnPSPS = 0.25 * ( PT[p*nBasis+i]*PT[q*nBasis+j]
		                    +PS[p*nBasis+i]*PS[q*nBasis+j] );
		r         = PTPT - PTPTnPSPS;

		// gradient in x direction
		Gx[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_X);
		Gx[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_X);
		Gx[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_X);
		Gx[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_X);

		// gradient in y direction
		Gy[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_Y);
		Gy[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_Y);
		Gy[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_Y);
		Gy[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_Y);

		// gradient in z direction
		Gz[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_Z);
		Gz[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_Z);
		Gz[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_Z);
		Gz[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_Z);
	}

}


// GradEri_Schwarz: is another variant of GradEri but it uses
// Schwarz screening to elliminate small integral.
//
// Oct 20, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Initial implementation
//
void GradEri_Schwarz(
	int nBasis,                     // number of basis function
	const struct GTOBasis_t *gto,   // basis function database
	const int *basis2Atom,          // basis to atom mapping
	const double *PT,               // total density matrix
	const double *PS,               // spin density matrix
	double *Gx,                     // returned gradient in x direction
	double *Gy,                     // returned gradient in y direction
	double *Gz){                    // returned gradient in z direction

	int p,q,i,j;           // basis function index
	double PTPT,PTPTnPSPS; // product of density matrix
	double r,f;            // generic real number

	double *Schwarz;       // typical Schwarz matrix
	double *Schwarz_dx;    // for derivative in x direction
	double *Schwarz_dy;    // for derivative in y direction
	double *Schwarz_dz;    // for derivative in z direction

	// generate all schwarz matrix
	Schwarz    = create_Schwarz(nBasis, gto);
	Schwarz_dx = create_grad_Schwarz(nBasis, gto, GRAD_X);
	Schwarz_dy = create_grad_Schwarz(nBasis, gto, GRAD_Y);
	Schwarz_dz = create_grad_Schwarz(nBasis, gto, GRAD_Z);

	// loop thru all permutation of basis
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++) 
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){

		// compute coulomb and exchange multiplier
		PTPT      = 0.50 *   PT[p*nBasis+q]*PT[i*nBasis+j];
		PTPTnPSPS = 0.25 * ( PT[p*nBasis+i]*PT[q*nBasis+j]
		                    +PS[p*nBasis+i]*PS[q*nBasis+j] );
		r         = PTPT - PTPTnPSPS;

		// screening 
		f = fabs(r);
		if(f < GRAD_SCHWARZ_CUTOFF) continue;

#define SCREENX(P,Q,I,J) if(f*Schwarz_dx[P*nBasis+Q]*Schwarz[I*nBasis+J] > GRAD_SCHWARZ_CUTOFF)
#define SCREENY(P,Q,I,J) if(f*Schwarz_dy[P*nBasis+Q]*Schwarz[I*nBasis+J] > GRAD_SCHWARZ_CUTOFF)
#define SCREENZ(P,Q,I,J) if(f*Schwarz_dz[P*nBasis+Q]*Schwarz[I*nBasis+J] > GRAD_SCHWARZ_CUTOFF)

		// gradient in x direction
		SCREENX(p,q,i,j) Gx[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_X);
		SCREENX(q,p,i,j) Gx[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_X);
		SCREENX(i,j,p,q) Gx[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_X);
		SCREENX(j,i,p,q) Gx[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_X);

		// gradient in y direction
		SCREENY(p,q,i,j) Gy[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_Y);
		SCREENY(q,p,i,j) Gy[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_Y);
		SCREENY(i,j,p,q) Gy[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_Y);
		SCREENY(j,i,p,q) Gy[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_Y);

		// gradient in z direction
		SCREENZ(p,q,i,j) Gz[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_Z);
		SCREENZ(q,p,i,j) Gz[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_Z);
		SCREENZ(i,j,p,q) Gz[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_Z);
		SCREENZ(j,i,p,q) Gz[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_Z);

	}

	// free memory
	free(Schwarz);
	free(Schwarz_dx);
	free(Schwarz_dy);
	free(Schwarz_dz);
}
#undef SCREENX
#undef SCREENY
#undef SCREENZ


// GradEri_Schwarz_Sym: is also another variant of GradEri but it uses
// Schwarz screening to elliminate small integral and use that fact 
// 2-electron integrals are identical under symmetric exchange of
// r1 <--> r2 and so on.
//
// Oct 21, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Initial implementation
//
void GradEri_Schwarz_Sym(
	int nBasis,                     // number of basis function
	const struct GTOBasis_t *gto,   // basis function database
	const int *basis2Atom,          // basis to atom mapping
	const double *PT,               // total density matrix
	const double *PS,               // spin density matrix
	double *Gx,                     // returned gradient in x direction
	double *Gy,                     // returned gradient in y direction
	double *Gz){                    // returned gradient in z direction

	int p,q,i,j,n;               // basis function index
	double r,f;                  // generic real number

	double *Schwarz;             // typical Schwarz matrix
	double *Schwarz_dx;          // for derivative in x direction
	double *Schwarz_dy;          // for derivative in y direction
	double *Schwarz_dz;          // for derivative in z direction

	// generate all schwarz matrix
	Schwarz    = create_Schwarz(nBasis, gto);
	Schwarz_dx = create_grad_Schwarz(nBasis, gto, GRAD_X);
	Schwarz_dy = create_grad_Schwarz(nBasis, gto, GRAD_Y);
	Schwarz_dz = create_grad_Schwarz(nBasis, gto, GRAD_Z);

	// loop thru all permutation of basis
	for(p=0; p < nBasis; p++)
	for(q=0; q < p+1; q++) 
	for(i=0; i < p+1; i++){ if(i==p) n=q+1; else n=i+1;
	for(j=0; j < n; j++){

#define CONTRIBUTE_GRAD(p,q,i,j) r += 0.50*  PT[p*nBasis+q]*PT[i*nBasis+j];  \
                                 r -= 0.25*( PT[p*nBasis+i]*PT[q*nBasis+j]   \
                                            +PS[p*nBasis+i]*PS[q*nBasis+j] );

		// determine pqij type
		r = 0.0;
		if((p==q)&&(i==j)&&(p==i)){  // all same
			  CONTRIBUTE_GRAD(p,q,i,j);
			//CONTRIBUTE_GRAD(q,p,i,j);
			//CONTRIBUTE_GRAD(p,q,j,i);
			//CONTRIBUTE_GRAD(q,p,j,i);
			//CONTRIBUTE_GRAD(i,j,p,q);
			//CONTRIBUTE_GRAD(j,i,p,q);
			//CONTRIBUTE_GRAD(i,j,q,p);
			//CONTRIBUTE_GRAD(j,i,q,p);
		}else if((p==q)&&(i==j)){    // 2 pairs
			  CONTRIBUTE_GRAD(p,q,i,j);
			//CONTRIBUTE_GRAD(q,p,i,j);
			//CONTRIBUTE_GRAD(p,q,j,i);
			//CONTRIBUTE_GRAD(q,p,j,i);
			  CONTRIBUTE_GRAD(i,j,p,q);
			//CONTRIBUTE_GRAD(j,i,p,q);
			//CONTRIBUTE_GRAD(i,j,q,p);
			//CONTRIBUTE_GRAD(j,i,q,p);
		}else if(p==q){              // pq pair
			  CONTRIBUTE_GRAD(p,q,i,j);
			//CONTRIBUTE_GRAD(q,p,i,j);
			  CONTRIBUTE_GRAD(p,q,j,i);
			//CONTRIBUTE_GRAD(q,p,j,i);
			  CONTRIBUTE_GRAD(i,j,p,q);
			  CONTRIBUTE_GRAD(j,i,p,q);
			//CONTRIBUTE_GRAD(i,j,q,p);
			//CONTRIBUTE_GRAD(j,i,q,p);
		}else if(i==j){              // ij pair
			  CONTRIBUTE_GRAD(p,q,i,j);
			  CONTRIBUTE_GRAD(q,p,i,j);
			//CONTRIBUTE_GRAD(p,q,j,i);
			//CONTRIBUTE_GRAD(q,p,j,i);
			  CONTRIBUTE_GRAD(i,j,p,q);
			//CONTRIBUTE_GRAD(j,i,p,q);
			  CONTRIBUTE_GRAD(i,j,q,p);
			//CONTRIBUTE_GRAD(j,i,q,p);
		}else if((p==i)&&(q==j)){    // pi-qj pair
			  CONTRIBUTE_GRAD(p,q,i,j);
			  CONTRIBUTE_GRAD(q,p,i,j);
			  CONTRIBUTE_GRAD(p,q,j,i);
			  CONTRIBUTE_GRAD(q,p,j,i);
			//CONTRIBUTE_GRAD(i,j,p,q);
			//CONTRIBUTE_GRAD(j,i,p,q);
			//CONTRIBUTE_GRAD(i,j,q,p);
			//CONTRIBUTE_GRAD(j,i,q,p);
		}else{                       // all distinct
			  CONTRIBUTE_GRAD(p,q,i,j);
			  CONTRIBUTE_GRAD(q,p,i,j);
			  CONTRIBUTE_GRAD(p,q,j,i);
			  CONTRIBUTE_GRAD(q,p,j,i);
			  CONTRIBUTE_GRAD(i,j,p,q);
			  CONTRIBUTE_GRAD(j,i,p,q);
			  CONTRIBUTE_GRAD(i,j,q,p);
			  CONTRIBUTE_GRAD(j,i,q,p);
		}

		f = fabs(r);
		if(f < GRAD_SCHWARZ_CUTOFF) continue;

#define SCREENX(P,Q,I,J) if(f*Schwarz_dx[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         > GRAD_SCHWARZ_CUTOFF) goto essemble_grad;
#define SCREENY(P,Q,I,J) if(f*Schwarz_dy[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         > GRAD_SCHWARZ_CUTOFF) goto essemble_grad;
#define SCREENZ(P,Q,I,J) if(f*Schwarz_dz[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         > GRAD_SCHWARZ_CUTOFF) goto essemble_grad;

		SCREENX(p,q,i,j); SCREENX(q,p,i,j); SCREENX(i,j,p,q); SCREENX(j,i,p,q);
		SCREENY(p,q,i,j); SCREENY(q,p,i,j); SCREENY(i,j,p,q); SCREENY(j,i,p,q);
		SCREENZ(p,q,i,j); SCREENZ(q,p,i,j); SCREENZ(i,j,p,q); SCREENZ(j,i,p,q);

		continue;

essemble_grad:

		// gradient in x direction
		Gx[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_X);
		Gx[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_X);
		Gx[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_X);
		Gx[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_X);

		// gradient in y direction
		Gy[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_Y);
		Gy[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_Y);
		Gy[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_Y);
		Gy[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_Y);

		// gradient in z direction
		Gz[basis2Atom[p]] += r * GRAD_ERI_CALL(p,q,i,j,GRAD_Z);
		Gz[basis2Atom[q]] += r * GRAD_ERI_CALL(q,p,i,j,GRAD_Z);
		Gz[basis2Atom[i]] += r * GRAD_ERI_CALL(i,j,p,q,GRAD_Z);
		Gz[basis2Atom[j]] += r * GRAD_ERI_CALL(j,i,p,q,GRAD_Z);

	}
	}

	// free memory
	free(Schwarz);
	free(Schwarz_dx);
	free(Schwarz_dy);
	free(Schwarz_dz);
}
#undef SCREENX
#undef SCREENY
#undef SCREENZ
#undef GRAD_ERI_CALL


//
// shell_prepGrad : prepare shell pointer and variables to be used with
// function GradEri_ShellSet
// It returns the number of shell
//
// Oct 21, 2010 - Teepanis Chachiyo
//    Function prototype document
//
// Oct 17, 2010 - Teepanis Chachiyo
//    Initial implementation
//
int shell_prepGrad(
	int nBasis,                    // number of basis function
	const struct GTOBasis_t *gto,  // basis function data structure
	const double *Schwarz,         // Schwarz matrix for screening
	int **RshellMap,               // returned shellMap
	int **RshellMaxL,              // returned shellMaxL
	double **RSchwarz_Shell){      // returned max Schwarz within shell

	int nShell;                 // the number of shell
	int pS, qS;                 // shell looping indexes
	int p,q;                    // generic basis indexex
	int n;                      // generic index
	double r;                   // generic real
	int *shellMap;
	int *shellMaxL;
	double *Schwarz_Shell;

	// allocate memory
	shellMap  = calloc(nBasis+1, sizeof(int));
	shellMaxL = calloc(nBasis+1, sizeof(int));
	if(shellMap==NULL || shellMaxL==NULL){
		printf("shell_prep: cannot allocate memory\n");
		exit(-1);
	}

	// generate shell mapping
	// first basis contributes first shell always
	nShell       = 1;
	shellMap[0]  = 0;
	shellMaxL[0] = gto[0].l+gto[0].m+gto[0].n;
	for(n=1; n < nBasis; n++){

		// check if we are still in the same shell
		if( isSameShell(nBasis,n,shellMap[nShell-1],gto) ){
			
			// determine maximum angular moment within same shell
			if((gto[n].l+gto[n].m+gto[n].n)>shellMaxL[nShell-1])
				shellMaxL[nShell-1] = gto[n].l+gto[n].m+gto[n].n;

			continue;

		// terminate to form shell
		}else{
			shellMap[nShell] = n;
			nShell++;
		}
	}
	// last item in shellMap must be the nBasis
	shellMap[nShell] = nBasis;

	//
	// allocate and determine schwarz matrix at shell level
	//
	Schwarz_Shell = calloc(nShell*nShell, sizeof(double));
	if(Schwarz_Shell==NULL){
		printf("shell_prep: cannot allocate Schwarz_Shell\n");
		exit(-1);
	}
	// loop for each shell
	for(pS=0; pS < nShell; pS++)
	for(qS=0; qS < nShell; qS++){
		// determine maximum Scharz[p,q] within this shell
		r = 0.0;
		for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
		for(q=shellMap[qS]; q < shellMap[qS+1] && q < (p+1); q++)
			if(Schwarz[p*nBasis+q] > r) r = Schwarz[p*nBasis+q];
		Schwarz_Shell[pS*nShell+qS] = r;
	} 

	// set returned variables
	*RshellMap      = shellMap;
	*RshellMaxL     = shellMaxL;
	*RSchwarz_Shell = Schwarz_Shell;

	return nShell;
}

// GradEri_ShellSet: is also another variant of GradEri but it uses
// Schwarz screening to elliminate small integral and use that fact 
// 2-electron integrals are identical under symmetric exchange of
// r1 <--> r2 and so on; and more importantly, group the integral
// into shell. The code is very similar to GTO_JK_Matrix_ShellSet
//
// Oct 21, 2010 - Teepanis Chachiyo and Theerapon Khamla
//    Initial implementation
//
// Oct 24, 2010 - Teepanis Chachiyo
//    Avoid calculating (pq|idj) using translational invariant
//    Ref: 1) A.Korminicki, K.Ishida, K.Morokuma, R.Ditchfield, and M.Conrad, 
//            Chem.Phys.Lett.45,595(1977).
//         2) Martin-Head Gordon and John A. Pople, "A method for two electron 
//            integral and integral derivative evaluaton us recurrence relations"
//            J.Chem.Phys.89,9(1988)
//
// May 19, 2011 - Teepanis Chachiyo
//    Adjust how to locate index Bx,By,Bz so that it coincide with the new
//    version of genSetBxyzF. Add the "shell_prepGrad" here because the 
//    original one with matrix.c has been changed to the new shell loop
//    structure and is incompatible to the GradEriSet at the moment. 

// support s,p,d,f,g
#define MAXL          4
#define MAXL1         5
#define MAXL1PWR2    25
#define MAXL1PWR3   125
#define MAXL1PWR4   625
#define MAXL1PWR5  3125
#define MAXL1PWR6 15625
#define TWOMAXL1     10
#define FOURMAXL1    20
#define MAXSHELLINT 50625
void GradEri_ShellSet(
	int nBasis,                     // number of basis function
	const struct GTOBasis_t *gto,   // basis function database
	const int *basis2Atom,          // basis to atom mapping
	const double *PT,               // total density matrix
	const double *PS,               // spin density matrix
	double *Gx,                     // returned gradient in x direction
	double *Gy,                     // returned gradient in y direction
	double *Gz){                    // returned gradient in z direction

	int pS,qS,iS,jS;                // looping variables for shells
	int p,q,i,j;                    // looping variables for basis functions
	int pC,qC,iC,jC;                // contracted index
	double EE;                      // coulomb integral
	double *dPQIJx,*dPQIJy,*dPQIJz; // derivative dPQIJ in x,y,z direction
	double *PdQIJx,*PdQIJy,*PdQIJz; // derivative dPQIJ in x,y,z direction
	double *PQdIJx,*PQdIJy,*PQdIJz; // derivative dPQIJ in x,y,z direction
	double *PQIdJx,*PQIdJy,*PQIdJz; // derivative dPQIJ in x,y,z direction
	double *coef;                   // density matrix prefactor
	char *dPQIJxZ,*dPQIJyZ,*dPQIJzZ;// zero flag for dPQIJ in x,y,z direction
	char *PdQIJxZ,*PdQIJyZ,*PdQIJzZ;// zero flag for dPQIJ in x,y,z direction
	char *PQdIJxZ,*PQdIJyZ,*PQdIJzZ;// zero flag for dPQIJ in x,y,z direction
	char *PQIdJxZ,*PQIdJyZ,*PQIdJzZ;// zero flag for dPQIJ in x,y,z direction
	int passScreen;                 // screening test flag
	int nEE;                        // number of electron within shell counter
	double *Bx, *By, *Bz, *F;       // auxillary array THO paper
	int ki,kq,kp;                   // tBx,tBy,tBz index multiplier
	double t;                       // integral scaling
	double *tBx, *tBy, *tBz;        // tempo to B pointer
	int mX, mY, mZ;                 // maximum order in x,y,z direction
	int kX, kY, kZ;                 // summation loop
	int nShell;                     // number of shell
	int *shellMap;                  // mapping between shell and basis function
	int *shellMaxL;                 // maximum angular moment for each shell
	double *Schwarz_Shell;          // schwarz screening at shell level
	int n;                          // generic index
	double r;                       // generic real number

	double *Schwarz;                // typical Schwarz matrix
	double *Schwarz_dx;             // for derivative in x direction
	double *Schwarz_dy;             // for derivative in y direction
	double *Schwarz_dz;             // for derivative in z direction

	double abcd, a1bcd, ab1cd;      // horizontal recurrence variables

	// generate all schwarz matrix
	Schwarz    = create_Schwarz(nBasis, gto);
	Schwarz_dx = create_grad_Schwarz(nBasis, gto, GRAD_X);
	Schwarz_dy = create_grad_Schwarz(nBasis, gto, GRAD_Y);
	Schwarz_dz = create_grad_Schwarz(nBasis, gto, GRAD_Z);

#define ALLOC(array,item,type)                                         \
array=calloc(item,sizeof(type));                                       \
if(array==NULL){                                                       \
	printf("GTO_JK_Matrix_ShellSet - error cannot allocate memory\n"); \
	exit(-1);                                                          \
}

	// memory allocations
	ALLOC(dPQIJx,  MAXSHELLINT,       double);
	ALLOC(dPQIJy,  MAXSHELLINT,       double);
	ALLOC(dPQIJz,  MAXSHELLINT,       double);
	ALLOC(PdQIJx,  MAXSHELLINT,       double);
	ALLOC(PdQIJy,  MAXSHELLINT,       double);
	ALLOC(PdQIJz,  MAXSHELLINT,       double);
	ALLOC(PQdIJx,  MAXSHELLINT,       double);
	ALLOC(PQdIJy,  MAXSHELLINT,       double);
	ALLOC(PQdIJz,  MAXSHELLINT,       double);
	ALLOC(PQIdJx,  MAXSHELLINT,       double);
	ALLOC(PQIdJy,  MAXSHELLINT,       double);
	ALLOC(PQIdJz,  MAXSHELLINT,       double);
	ALLOC(coef,    MAXSHELLINT,       double);
	ALLOC(dPQIJxZ, MAXSHELLINT,         char);
	ALLOC(dPQIJyZ, MAXSHELLINT,         char);
	ALLOC(dPQIJzZ, MAXSHELLINT,         char);
	ALLOC(PdQIJxZ, MAXSHELLINT,         char);
	ALLOC(PdQIJyZ, MAXSHELLINT,         char);
	ALLOC(PdQIJzZ, MAXSHELLINT,         char);
	ALLOC(PQdIJxZ, MAXSHELLINT,         char);
	ALLOC(PQdIJyZ, MAXSHELLINT,         char);
	ALLOC(PQdIJzZ, MAXSHELLINT,         char);
	ALLOC(PQIdJxZ, MAXSHELLINT,         char);
	ALLOC(PQIdJyZ, MAXSHELLINT,         char);
	ALLOC(PQIdJzZ, MAXSHELLINT,         char);
	ALLOC(Bx,      4*MAXL1*MAXL1PWR4, double);
	ALLOC(By,      4*MAXL1*MAXL1PWR4, double);
	ALLOC(Bz,      4*MAXL1*MAXL1PWR4, double);
	ALLOC(F,       4*MAXL1,           double);
#undef ALLOC

	// prepare shell data
	nShell = shell_prepGrad(nBasis, gto, Schwarz, 
	                        &shellMap, &shellMaxL, &Schwarz_Shell);

	// determine maximum Scharz[p,q] within this shell
	for(pS=0; pS < nShell; pS++)
	for(qS=0; qS < nShell; qS++){
		r = 0.0;
		for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
		for(q=shellMap[qS]; q < shellMap[qS+1] && q < (p+1); q++){
			if(Schwarz   [p*nBasis+q] > r) r = Schwarz   [p*nBasis+q];
			if(Schwarz_dx[p*nBasis+q] > r) r = Schwarz_dx[p*nBasis+q];
			if(Schwarz_dy[p*nBasis+q] > r) r = Schwarz_dy[p*nBasis+q];
			if(Schwarz_dz[p*nBasis+q] > r) r = Schwarz_dz[p*nBasis+q];
		}
		Schwarz_Shell[pS*nShell+qS] = r;
	} 

//
// define loop structure for all distinct pqij
//
#define PQIJ_SHELL_BEGIN                                    \
for(p=shellMap[pS]; p < shellMap[pS+1]; p++)                \
for(q=shellMap[qS]; q < shellMap[qS+1] && q < (p+1); q++)   \
for(i=shellMap[iS]; i < shellMap[iS+1] && i < (p+1); i++)   \
for(j=shellMap[jS]; j < shellMap[jS+1]; j++){               \
   if(i==p) n=q+1; else n=i+1;                              \
   if(!(j<n)) break;                                        \


#define PQIJ_SHELL_END                      }               \


	// permute all possible shells
	for(pS=0; pS < nShell; pS++)
	for(qS=0; qS < pS+1; qS++)
	for(iS=0; iS < pS+1; iS++)
	for(jS=0; jS < pS+1; jS++){

		// skip entire shell if possible
		if(Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS]
		   < GRAD_SCHWARZ_CUTOFF) continue;

		// compute density prefactor
		passScreen=0;
		nEE=0;
		PQIJ_SHELL_BEGIN

			// determine pqij type
			r = 0.0;
			if((p==q)&&(i==j)&&(p==i)){  // all same
				  CONTRIBUTE_GRAD(p,q,i,j);
				//CONTRIBUTE_GRAD(q,p,i,j);
				//CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				//CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if((p==q)&&(i==j)){    // 2 pairs
				  CONTRIBUTE_GRAD(p,q,i,j);
				//CONTRIBUTE_GRAD(q,p,i,j);
				//CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if(p==q){              // pq pair
				  CONTRIBUTE_GRAD(p,q,i,j);
				//CONTRIBUTE_GRAD(q,p,i,j);
				  CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				  CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if(i==j){              // ij pair
				  CONTRIBUTE_GRAD(p,q,i,j);
				  CONTRIBUTE_GRAD(q,p,i,j);
				//CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				  CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if((p==i)&&(q==j)){    // pi-qj pair
				  CONTRIBUTE_GRAD(p,q,i,j);
				  CONTRIBUTE_GRAD(q,p,i,j);
				  CONTRIBUTE_GRAD(p,q,j,i);
				  CONTRIBUTE_GRAD(q,p,j,i);
				//CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else{                       // all distinct
				  CONTRIBUTE_GRAD(p,q,i,j);
				  CONTRIBUTE_GRAD(q,p,i,j);
				  CONTRIBUTE_GRAD(p,q,j,i);
				  CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				  CONTRIBUTE_GRAD(j,i,p,q);
				  CONTRIBUTE_GRAD(i,j,q,p);
				  CONTRIBUTE_GRAD(j,i,q,p);
			}

			// store it in array
			coef[nEE] = r;

			// already passed
			if(passScreen) { nEE++; continue; }

			// pre-screening
			r = fabs(r);
			if(r < GRAD_SCHWARZ_CUTOFF){ nEE++; continue; }

#define SCREENX(P,Q,I,J) if(r*Schwarz_dx[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         > GRAD_SCHWARZ_CUTOFF) passScreen=1;
#define SCREENY(P,Q,I,J) if(r*Schwarz_dy[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         > GRAD_SCHWARZ_CUTOFF) passScreen=1;
#define SCREENZ(P,Q,I,J) if(r*Schwarz_dz[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         > GRAD_SCHWARZ_CUTOFF) passScreen=1;

			SCREENX(p,q,i,j); SCREENY(p,q,i,j); SCREENZ(p,q,i,j);
			SCREENX(q,p,i,j); SCREENY(q,p,i,j); SCREENZ(q,p,i,j);
			SCREENX(i,j,p,q); SCREENY(i,j,p,q); SCREENZ(i,j,p,q);
			//SCREENX(j,i,p,q); SCREENY(j,i,p,q); SCREENZ(j,i,p,q);

#undef SCREENX
#undef SCREENY
#undef SCREENZ

			nEE++;
		PQIJ_SHELL_END

		if(!passScreen) continue;

		// check that the number of integral does not exceed maximum 
		if(nEE > MAXSHELLINT){
			printf("GradEri_ShellSet - error too many integral in shell\n");
			exit(-1);
		}

		// preparation
		nEE=0;
		PQIJ_SHELL_BEGIN
			// set zero
			dPQIJx[nEE] = 0.0; dPQIJy[nEE] = 0.0; dPQIJz[nEE] = 0.0;
			PdQIJx[nEE] = 0.0; PdQIJy[nEE] = 0.0; PdQIJz[nEE] = 0.0;
			PQdIJx[nEE] = 0.0; PQdIJy[nEE] = 0.0; PQdIJz[nEE] = 0.0;
			//PQIdJx[nEE] = 0.0; PQIdJy[nEE] = 0.0; PQIdJz[nEE] = 0.0;

			// another round of screening
			
#define SCREENX(P,Q,I,J) if(r*Schwarz_dx[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         < GRAD_SCHWARZ_CUTOFF)  
#define SCREENY(P,Q,I,J) if(r*Schwarz_dy[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         < GRAD_SCHWARZ_CUTOFF)
#define SCREENZ(P,Q,I,J) if(r*Schwarz_dz[P*nBasis+Q]*Schwarz[I*nBasis+J] \
                         < GRAD_SCHWARZ_CUTOFF)
			r = fabs(coef[nEE]);
			SCREENX(p,q,i,j) dPQIJxZ[nEE] = 0; else dPQIJxZ[nEE] = 1;
			SCREENX(q,p,i,j) PdQIJxZ[nEE] = 0; else PdQIJxZ[nEE] = 1;
			SCREENX(i,j,p,q) PQdIJxZ[nEE] = 0; else PQdIJxZ[nEE] = 1;
			//SCREENX(j,i,p,q) PQIdJxZ[nEE] = 0; else PQIdJxZ[nEE] = 1;
			SCREENY(p,q,i,j) dPQIJyZ[nEE] = 0; else dPQIJyZ[nEE] = 1;
			SCREENY(q,p,i,j) PdQIJyZ[nEE] = 0; else PdQIJyZ[nEE] = 1;
			SCREENY(i,j,p,q) PQdIJyZ[nEE] = 0; else PQdIJyZ[nEE] = 1;
			//SCREENY(j,i,p,q) PQIdJyZ[nEE] = 0; else PQIdJyZ[nEE] = 1;
			SCREENZ(p,q,i,j) dPQIJzZ[nEE] = 0; else dPQIJzZ[nEE] = 1;
			SCREENZ(q,p,i,j) PdQIJzZ[nEE] = 0; else PdQIJzZ[nEE] = 1;
			SCREENZ(i,j,p,q) PQdIJzZ[nEE] = 0; else PQdIJzZ[nEE] = 1;
			//SCREENZ(j,i,p,q) PQIdJzZ[nEE] = 0; else PQIdJzZ[nEE] = 1;
#undef SCREENX
#undef SCREENY
#undef SCREENZ
			nEE++;
		PQIJ_SHELL_END

		///////////////////////////////////////////////
		// algorithm:
		//    1) loop through all contraction
		//         1.1) generate Bx,By,Bz,F
		//         1.2) compute (dPQ|IJ) (PdQ|IJ) (PQ|dIJ) (PQ|IdJ)
		//              in all direction x,y,z
		//
		//    2) loop through all basis function in this shell
		//         2.1) add contribution to forces
		///////////////////////////////////////////////

		// set pqij index
		p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		
		/////////////////////////////
		// start contraction loop  //
		/////////////////////////////
		for(pC=0; pC < gto[p].nContract; pC++)
		for(qC=0; qC < gto[q].nContract; qC++)
		for(iC=0; iC < gto[i].nContract; iC++)
		for(jC=0; jC < gto[j].nContract; jC++){

			// generate Bx, By , Bz, and F
			t=genSetBxyzF(
			        gto[p].x0,gto[p].y0,gto[p].z0,shellMaxL[pS]+1,gto[p].exp[pC],
			        gto[q].x0,gto[q].y0,gto[q].z0,shellMaxL[qS]+0,gto[q].exp[qC],
			        gto[i].x0,gto[i].y0,gto[i].z0,shellMaxL[iS]+1,gto[i].exp[iC],
			        gto[j].x0,gto[j].y0,gto[j].z0,shellMaxL[jS]+0,gto[j].exp[jC],
			        Bx,By,Bz,F,MAXL);

			if(t < PRIMITIVE_CUTOFF) continue;

			// reset electron counter
			nEE = 0;
			
			// loop all basis within shell
			PQIJ_SHELL_BEGIN

#define SUM_EE EE=0.0;                                    \
               for(kX=mX;kX>=0;kX--)                      \
               for(kY=mY;kY>=0;kY--)                      \
               for(kZ=mZ;kZ>=0;kZ--)                      \
               EE += tBx[kX]*tBy[kY]*tBz[kZ]*F[kX+kY+kZ];

				// compute prefactor for this primitive set
				r  = t * gto[p].coef[pC]*gto[p].norm[pC]
				       * gto[q].coef[qC]*gto[q].norm[qC]
				       * gto[i].coef[iC]*gto[i].norm[iC]
				       * gto[j].coef[jC]*gto[j].norm[jC];

				// compute index multiplier
				ki =    (shellMaxL[jS]+1);
				kq = ki*(shellMaxL[iS]+2);
				kp = kq*(shellMaxL[qS]+1);
#define IBX(P,Q,I,J) FOURMAXL1*( (gto[p].l P)*kp + (gto[q].l Q)*kq + (gto[i].l I)*ki + (gto[j].l J));
#define IBY(P,Q,I,J) FOURMAXL1*( (gto[p].m P)*kp + (gto[q].m Q)*kq + (gto[i].m I)*ki + (gto[j].m J));
#define IBZ(P,Q,I,J) FOURMAXL1*( (gto[p].n P)*kp + (gto[q].n Q)*kq + (gto[i].n I)*ki + (gto[j].n J));

				// dPQIJx
				if(dPQIJxZ[nEE] || PdQIJxZ[nEE]){
					if(gto[p].l>0){
						tBx = Bx + IBX(-1,  ,  ,  );
						tBy = By + IBY(  ,  ,  ,  );
						tBz = Bz + IBZ(  ,  ,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l-1;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
						SUM_EE;
						EE           = - EE * gto[p].l;
						dPQIJx[nEE] += EE*r;
					}
					tBx = Bx + IBX(+1,  ,  ,  );
					tBy = By + IBY(  ,  ,  ,  );
					tBz = Bz + IBZ(  ,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l+1;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
					SUM_EE;
					a1bcd        = EE;
					EE           = 2.0*EE*gto[p].exp[pC];
					dPQIJx[nEE] += EE*r;
				}

				// PdQIJx
				if(PdQIJxZ[nEE]){
					if(gto[q].l>0){
						tBx = Bx + IBX(  ,-1,  ,  );
						tBy = By + IBY(  ,  ,  ,  );
						tBz = Bz + IBZ(  ,  ,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l-1;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
						SUM_EE;
						EE           = - EE * gto[q].l;
						PdQIJx[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,  ,  );
					tBy = By + IBY(  ,  ,  ,  );
					tBz = Bz + IBZ(  ,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;
					SUM_EE;
					abcd         = EE;
					ab1cd        = a1bcd + (gto[p].x0-gto[q].x0)*abcd;
					EE           = 2.0*ab1cd*gto[q].exp[qC];
					PdQIJx[nEE] += EE*r;
				}

				// PQdIJx
				if(PQdIJxZ[nEE]){
					if(gto[i].l>0){
						tBx = Bx + IBX(  ,  ,-1,  );
						tBy = By + IBY(  ,  ,  ,  );
						tBz = Bz + IBZ(  ,  ,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l-1;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
						SUM_EE;
						EE           = - EE * gto[i].l;
						PQdIJx[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,+1,  );
					tBy = By + IBY(  ,  ,  ,  );
					tBz = Bz + IBZ(  ,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l+1;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
					SUM_EE;
					EE           = 2.0*EE*gto[i].exp[iC];
					PQdIJx[nEE] += EE*r;
				}

				// dPQIJy
				if(dPQIJyZ[nEE] || PdQIJyZ[nEE]){
					if(gto[p].m>0){
						tBx = Bx + IBX(  ,  ,  ,  );
						tBy = By + IBY(-1,  ,  ,  );
						tBz = Bz + IBZ(  ,  ,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m-1;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
						SUM_EE;
						EE           = - EE * gto[p].m;
						dPQIJy[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,  ,  );
					tBy = By + IBY(+1,  ,  ,  );
					tBz = Bz + IBZ(  ,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m+1;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
					SUM_EE;
					a1bcd        = EE;
					EE           = 2.0*EE*gto[p].exp[pC];
					dPQIJy[nEE] += EE*r;
				}

				// PdQIJy
				if(PdQIJyZ[nEE]){
					if(gto[q].m>0){
						tBx = Bx + IBX(  ,  ,  ,  );
						tBy = By + IBY(  ,-1,  ,  );
						tBz = Bz + IBZ(  ,  ,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m-1;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
						SUM_EE;
						EE           = - EE * gto[q].m;
						PdQIJy[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,  ,  );
					tBy = By + IBY(  ,  ,  ,  );
					tBz = Bz + IBZ(  ,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;
					SUM_EE;
					abcd         = EE;
					ab1cd        = a1bcd + (gto[p].y0-gto[q].y0)*abcd;
					EE           = 2.0*ab1cd*gto[q].exp[qC];
					PdQIJy[nEE] += EE*r;
				}

				// PQdIJy
				if(PQdIJyZ[nEE]){
					if(gto[i].m>0){
						tBx = Bx + IBX(  ,  ,  ,  );
						tBy = By + IBY(  ,  ,-1,  );
						tBz = Bz + IBZ(  ,  ,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m-1;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
						SUM_EE;
						EE           = - EE * gto[i].m;
						PQdIJy[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,  ,  );
					tBy = By + IBY(  ,  ,+1,  );
					tBz = Bz + IBZ(  ,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m+1;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n  ;
					SUM_EE;
					EE           = 2.0*EE*gto[i].exp[iC];
					PQdIJy[nEE] += EE*r;
				}

				// dPQIJz
				if(dPQIJzZ[nEE] || PdQIJzZ[nEE]){
					if(gto[p].n>0){
						tBx = Bx + IBX(  ,  ,  ,  );
						tBy = By + IBY(  ,  ,  ,  );
						tBz = Bz + IBZ(-1,  ,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n-1;
						SUM_EE;
						EE           = - EE * gto[p].n;
						dPQIJz[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,  ,  );
					tBy = By + IBY(  ,  ,  ,  );
					tBz = Bz + IBZ(+1,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n+1;
					SUM_EE;
					a1bcd        = EE;
					EE           = 2.0*EE*gto[p].exp[pC];
					dPQIJz[nEE] += EE*r;
				}

				// PdQIJz
				if(PdQIJzZ[nEE]){
					if(gto[q].n>0){
						tBx = Bx + IBX(  ,  ,  ,  );
						tBy = By + IBY(  ,  ,  ,  );
						tBz = Bz + IBZ(  ,-1,  ,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n-1;
						SUM_EE;
						EE           = - EE * gto[q].n;
						PdQIJz[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,  ,  );
					tBy = By + IBY(  ,  ,  ,  );
					tBz = Bz + IBZ(  ,  ,  ,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;
					SUM_EE;
					abcd         = EE;
					ab1cd        = a1bcd + (gto[p].z0-gto[q].z0)*abcd;
					EE           = 2.0*ab1cd*gto[q].exp[qC];
					PdQIJz[nEE] += EE*r;
				}

				// PQdIJz
				if(PQdIJzZ[nEE]){
					if(gto[i].n>0){
						tBx = Bx + IBX(  ,  ,  ,  );
						tBy = By + IBY(  ,  ,  ,  );
						tBz = Bz + IBZ(  ,  ,-1,  );
						mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
						mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
						mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n-1;
						SUM_EE;
						EE           = - EE * gto[i].n;
						PQdIJz[nEE] += EE*r;
					}
					tBx = Bx + IBX(  ,  ,  ,  );
					tBy = By + IBY(  ,  ,  ,  );
					tBz = Bz + IBZ(  ,  ,+1,  );
					mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l  ;
					mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m  ;
					mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n+1;
					SUM_EE;
					EE           = 2.0*EE*gto[i].exp[iC];
					PQdIJz[nEE] += EE*r;
				}

				nEE++;
			PQIJ_SHELL_END

			// reset shell index to proper value in this shell
			p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		}

		// essemble into gradient
		nEE=0;
		PQIJ_SHELL_BEGIN
			// use translational invariant to compute PQIdJ
			PQIdJx[nEE] = -(dPQIJx[nEE]+PdQIJx[nEE]+PQdIJx[nEE]);
			PQIdJy[nEE] = -(dPQIJy[nEE]+PdQIJy[nEE]+PQdIJy[nEE]);
			PQIdJz[nEE] = -(dPQIJz[nEE]+PdQIJz[nEE]+PQdIJz[nEE]);

			// gradient in x direction
			Gx[basis2Atom[p]] += coef[nEE] * dPQIJx[nEE];
			Gx[basis2Atom[q]] += coef[nEE] * PdQIJx[nEE];
			Gx[basis2Atom[i]] += coef[nEE] * PQdIJx[nEE];
			Gx[basis2Atom[j]] += coef[nEE] * PQIdJx[nEE];

			// gradient in y direction
			Gy[basis2Atom[p]] += coef[nEE] * dPQIJy[nEE];
			Gy[basis2Atom[q]] += coef[nEE] * PdQIJy[nEE];
			Gy[basis2Atom[i]] += coef[nEE] * PQdIJy[nEE];
			Gy[basis2Atom[j]] += coef[nEE] * PQIdJy[nEE];
	
			// gradient in z direction
			Gz[basis2Atom[p]] += coef[nEE] * dPQIJz[nEE];
			Gz[basis2Atom[q]] += coef[nEE] * PdQIJz[nEE];
			Gz[basis2Atom[i]] += coef[nEE] * PQdIJz[nEE];
			Gz[basis2Atom[j]] += coef[nEE] * PQIdJz[nEE];

			nEE++;
		PQIJ_SHELL_END
	}

	// clean memory
	free(dPQIJx);
	free(dPQIJy);
	free(dPQIJz);
	free(PdQIJx);
	free(PdQIJy);
	free(PdQIJz);
	free(PQdIJx);
	free(PQdIJy);
	free(PQdIJz);
	free(PQIdJx);
	free(PQIdJy);
	free(PQIdJz);
	free(coef);
	free(dPQIJxZ);
	free(dPQIJyZ);
	free(dPQIJzZ);
	free(PdQIJxZ);
	free(PdQIJyZ);
	free(PdQIJzZ);
	free(PQdIJxZ);
	free(PQdIJyZ);
	free(PQdIJzZ);
	free(PQIdJxZ);
	free(PQIdJyZ);
	free(PQIdJzZ);
	free(shellMap);
	free(shellMaxL);
	free(Bx);
	free(By);
	free(Bz);
	free(F);
	free(Schwarz_Shell);
	free(Schwarz);
	free(Schwarz_dx);
	free(Schwarz_dy);
	free(Schwarz_dz);
}


// GradEri_ShellSet_Parallel: this is the parallel version of the subroutine
//
// May 11, 2016 - Teepanis Chachiyo
//     can now adjust Exchange contribution 
//
// March 8, 2016 - Teepanis Chachiyo
//     Clean up, new loop structure for grad calculation, and better screening
//     Add progress bar report
//
// March 8, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
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
	double *Gz){                    // returned gradient in z direction

	int pqS,pS,qS,iS,jS;            // looping variables for shells
	int p,q,i,j;                    // looping variables for basis functions
	int pC,qC,iC,jC;                // contracted index
	double EE;                      // coulomb integral
	double dPQIJ,PdQIJ,PQdIJ;       // derivative
	double *coef;                   // density matrix prefactor
	int nEE;                        // number of electron within shell counter
	double *Bx, *By, *Bz, *F, *BBF; // auxillary array THO paper
	int ki,kq,kp;                   // tBx,tBy,tBz index multiplier
	double t;                       // integral scaling
	double *tBx, *tBy, *tBz;        // tempo to B pointer
	int mX, mY, mZ;                 // maximum order in x,y,z direction
	int kX, kY, kZ;                 // summation loop
	int nShell;                     // number of shell
	int *shellMap;                  // mapping between shell and basis function
	int *shellMaxL;                 // maximum angular moment for each shell
	double *Schwarz_Shell;          // schwarz screening at shell level
	double *PM;                     // maximum total density within shell
	int n;                          // generic index
	double r;                       // generic real number

	double *Schwarz;                // typical Schwarz matrix
	double *gx,*gy,*gz;             // basis center force

	double abcd, a1bcd, ab1cd;      // horizontal recurrence variables

	// generate schwarz matrix
	Schwarz    = create_Schwarz(nBasis, gto);

#define ALLOC(array,item,type)                                            \
array=calloc(item,sizeof(type));                                          \
if(array==NULL){                                                          \
	printf("GradEri_ShellSet_Parallel - error cannot allocate memory\n"); \
	exit(-1);                                                             \
}

	// memory allocations
	ALLOC(coef,    MAXSHELLINT,       double);
	ALLOC(Bx,      4*MAXL1*MAXL1PWR4, double);
	ALLOC(By,      4*MAXL1*MAXL1PWR4, double);
	ALLOC(Bz,      4*MAXL1*MAXL1PWR4, double);
	ALLOC(F,       4*MAXL1,           double);
	ALLOC(BBF,     4*MAXL1,           double);
	ALLOC(gx,nBasis,                  double);
	ALLOC(gy,nBasis,                  double);
	ALLOC(gz,nBasis,                  double);

	// prepare shell data
	nShell = shell_prepGrad(nBasis, gto, Schwarz, 
	                        &shellMap, &shellMaxL, &Schwarz_Shell);

	ALLOC(PM,nShell*nShell,           double);
#undef ALLOC

	// determine maximum Scharz[p,q] and PT within this shell
	for(pS=0; pS < nShell; pS++)
	for(qS=0; qS < nShell; qS++){
		// schwarz
		r = 0.0;
		for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
		for(q=shellMap[qS]; q < shellMap[qS+1] && q < (p+1); q++){
			if(Schwarz[p*nBasis+q] > r) r = Schwarz[p*nBasis+q];
		}
		Schwarz_Shell[pS*nShell+qS] = r;

		// PT
		r = 0.0;
		for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
		for(q=shellMap[qS]; q < shellMap[qS+1] && q < (p+1); q++){
			if(fabs(PT[p*nBasis+q]) > r) r = fabs(PT[p*nBasis+q]);
		}
		PM[pS*nShell+qS] = r;
	}

//
// define loop structure for all distinct pqij
//
#define PQIJ_SHELL_BEGIN                                    \
for(p=shellMap[pS]; p < shellMap[pS+1]; p++)                \
for(q=shellMap[qS]; q < shellMap[qS+1] && q < (p+1); q++)   \
for(i=shellMap[iS]; i < shellMap[iS+1] && i < (p+1); i++)   \
for(j=shellMap[jS]; j < shellMap[jS+1]; j++){               \
   if(i==p) n=q+1; else n=i+1;                              \
   if(!(j<n)) break;                                        \

#define PQIJ_SHELL_END                      }               \


	// compute maximum number of inside cycle
	int maxPQ=0;
	for(pqS=0,pS=0; pS < nShell; pS++)
	for(qS=0; qS < pS+1; qS++,pqS++) if(pqS%nCPU==childID) maxPQ++;

	// reset progress bar
	stepProgressBar(0,0);

	// permute all possible shells
	for(pqS=0,pS=0; pS < nShell; pS++)
	for(qS=0; qS < pS+1; qS++,pqS++) if(pqS%nCPU==childID){

	// report status
	if(childID==0) stepProgressBar(61,maxPQ);

	for(iS=0; iS < pS+1; iS++)
	for(jS=0; jS < pS+1; jS++){

		// skip entire shell if possible
		t = Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS];
		if(t*PM[pS*nShell+qS]*PM[iS*nShell+jS] < GRAD_SCHWARZ_CUTOFF) continue;

	double Kcoef = Kfactor*0.25;
#undef  CONTRIBUTE_GRAD
#define CONTRIBUTE_GRAD(p,q,i,j) r += 0.50  *  PT[p*nBasis+q]*PT[i*nBasis+j];  \
                                 r -= Kcoef *( PT[p*nBasis+i]*PT[q*nBasis+j]   \
                                              +PS[p*nBasis+i]*PS[q*nBasis+j] );

		// compute density prefactor
		nEE=0;
		PQIJ_SHELL_BEGIN

			// determine pqij type
			r = 0.0;
			if((p==q)&&(i==j)&&(p==i)){  // all same
				  CONTRIBUTE_GRAD(p,q,i,j);
				//CONTRIBUTE_GRAD(q,p,i,j);
				//CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				//CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if((p==q)&&(i==j)){    // 2 pairs
				  CONTRIBUTE_GRAD(p,q,i,j);
				//CONTRIBUTE_GRAD(q,p,i,j);
				//CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if(p==q){              // pq pair
				  CONTRIBUTE_GRAD(p,q,i,j);
				//CONTRIBUTE_GRAD(q,p,i,j);
				  CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				  CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if(i==j){              // ij pair
				  CONTRIBUTE_GRAD(p,q,i,j);
				  CONTRIBUTE_GRAD(q,p,i,j);
				//CONTRIBUTE_GRAD(p,q,j,i);
				//CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				  CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else if((p==i)&&(q==j)){    // pi-qj pair
				  CONTRIBUTE_GRAD(p,q,i,j);
				  CONTRIBUTE_GRAD(q,p,i,j);
				  CONTRIBUTE_GRAD(p,q,j,i);
				  CONTRIBUTE_GRAD(q,p,j,i);
				//CONTRIBUTE_GRAD(i,j,p,q);
				//CONTRIBUTE_GRAD(j,i,p,q);
				//CONTRIBUTE_GRAD(i,j,q,p);
				//CONTRIBUTE_GRAD(j,i,q,p);
			}else{                       // all distinct
				  CONTRIBUTE_GRAD(p,q,i,j);
				  CONTRIBUTE_GRAD(q,p,i,j);
				  CONTRIBUTE_GRAD(p,q,j,i);
				  CONTRIBUTE_GRAD(q,p,j,i);
				  CONTRIBUTE_GRAD(i,j,p,q);
				  CONTRIBUTE_GRAD(j,i,p,q);
				  CONTRIBUTE_GRAD(i,j,q,p);
				  CONTRIBUTE_GRAD(j,i,q,p);
			}

			// store it in array
			coef[nEE] = r;

			nEE++;
		PQIJ_SHELL_END

		// check that the number of integral does not exceed maximum 
		if(nEE > MAXSHELLINT){
			printf("GradEri_ShellSet_Parallel - error too many integral in shell\n");
			exit(-1);
		}

		// compute index multiplier
		ki =    (shellMaxL[jS]+1);
		kq = ki*(shellMaxL[iS]+2);
		kp = kq*(shellMaxL[qS]+1);
		ki *= FOURMAXL1;
		kq *= FOURMAXL1;
		kp *= FOURMAXL1;

		///////////////////////////////////////////////
		// algorithm:
		//    1) loop through all contraction
		//         1.1) generate Bx,By,Bz,F
		//         1.2) compute (dPQ|IJ) (PdQ|IJ) (PQ|dIJ) (PQ|IdJ)
		//              in all direction x,y,z
		//
		//    2) loop through all basis function in this shell
		//         2.1) add contribution to forces
		///////////////////////////////////////////////

		// set pqij index
		p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		
		/////////////////////////////
		// start contraction loop  //
		/////////////////////////////
		for(pC=0; pC < gto[p].nContract; pC++)
		for(qC=0; qC < gto[q].nContract; qC++)
		for(iC=0; iC < gto[i].nContract; iC++)
		for(jC=0; jC < gto[j].nContract; jC++){

			// generate Bx, By , Bz, and F
			t=genSetBxyzF(
				gto[p].x0,gto[p].y0,gto[p].z0,shellMaxL[pS]+1,gto[p].exp[pC],
				gto[q].x0,gto[q].y0,gto[q].z0,shellMaxL[qS]+0,gto[q].exp[qC],
				gto[i].x0,gto[i].y0,gto[i].z0,shellMaxL[iS]+1,gto[i].exp[iC],
				gto[j].x0,gto[j].y0,gto[j].z0,shellMaxL[jS]+0,gto[j].exp[jC],
				Bx,By,Bz,F,MAXL);

			if(t < PRIMITIVE_CUTOFF) continue;

			// reset electron counter
			nEE = 0;

			// loop all basis within shell
			PQIJ_SHELL_BEGIN

				// compute prefactor for this primitive set
				r  = t * gto[p].coef[pC]*gto[p].norm[pC]
				       * gto[q].coef[qC]*gto[q].norm[qC]
				       * gto[i].coef[iC]*gto[i].norm[iC]
				       * gto[j].coef[jC]*gto[j].norm[jC];

				if(fabs(coef[nEE]*r)< GRAD_SCHWARZ_CUTOFF){ nEE++; continue; }

				tBx = Bx + (gto[p].l)*kp + (gto[q].l)*kq + (gto[i].l)*ki + (gto[j].l)*FOURMAXL1;
				tBy = By + (gto[p].m)*kp + (gto[q].m)*kq + (gto[i].m)*ki + (gto[j].m)*FOURMAXL1;
				tBz = Bz + (gto[p].n)*kp + (gto[q].n)*kq + (gto[i].n)*ki + (gto[j].n)*FOURMAXL1;
				mX  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
				mY  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
				mZ  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;

#define GRAD1D(x,X,x0,A,y,Y,z,Z)                                      \
				dPQIJ=PdQIJ=PQdIJ=0.0;                                \
				                                                      \
				for(k##X=m##X+1;k##X>=0;k##X--){                      \
				BBF[k##X]=0.0;                                        \
				for(k##Y=m##Y;k##Y>=0;k##Y--)                         \
				for(k##Z=m##Z;k##Z>=0;k##Z--)                         \
				BBF[k##X]+=tB##y[k##Y]*tB##z[k##Z]*F[k##X+k##Y+k##Z]; \
				}                                                     \
				                                                      \
				if(gto[p].A>0){                                       \
				tB##x-=kp; m##X--;                                                \
				for(EE=0.0,k##X=m##X;k##X>=0;k##X--) EE += tB##x[k##X]*BBF[k##X]; \
				EE = - EE * gto[p].A;                                             \
				dPQIJ += EE*r;                                                    \
				tB##x+=kp; m##X++;                                                \
				}                                                                 \
				tB##x+=kp; m##X++;                                                \
				for(EE=0.0,k##X=m##X;k##X>=0;k##X--) EE += tB##x[k##X]*BBF[k##X]; \
				a1bcd = EE; EE = 2.0*EE*gto[p].exp[pC];                           \
				dPQIJ += EE*r;                                                    \
				tB##x-=kp; m##X--;                                                \
				                                                                  \
				if(gto[q].A>0){                                                   \
				tB##x-=kq; m##X--;                                                \
				for(EE=0.0,k##X=m##X;k##X>=0;k##X--) EE += tB##x[k##X]*BBF[k##X]; \
				                   EE = - EE * gto[q].A;                          \
				                   PdQIJ += EE*r;                                 \
				tB##x+=kq; m##X++;                                                \
				}                                                                 \
				for(EE=0.0,k##X=m##X;k##X>=0;k##X--) EE += tB##x[k##X]*BBF[k##X]; \
				abcd = EE;                                                        \
				ab1cd = a1bcd + (gto[p].x0-gto[q].x0)*abcd;                       \
				EE = 2.0*ab1cd*gto[q].exp[qC];                                    \
				PdQIJ += EE*r;                                                    \
				                                                                  \
				if(gto[i].A>0){                                                   \
				tB##x-=ki; m##X--;                                                \
				for(EE=0.0,k##X=m##X;k##X>=0;k##X--) EE += tB##x[k##X]*BBF[k##X]; \
				EE = - EE * gto[i].A;                                             \
				PQdIJ += EE*r;                                                    \
				tB##x+=ki; m##X++;                                                \
				}                                                                 \
				tB##x+=ki; m##X++;                                                \
				for(EE=0.0,k##X=m##X;k##X>=0;k##X--) EE += tB##x[k##X]*BBF[k##X]; \
				EE = 2.0*EE*gto[i].exp[iC];                                       \
				PQdIJ += EE*r;                                                    \
				tB##x-=ki; m##X--;                                                \
				                                                                  \
				g##x[p] += coef[nEE]*dPQIJ;                                       \
				g##x[q] += coef[nEE]*PdQIJ;                                       \
				g##x[i] += coef[nEE]*PQdIJ;                                       \
				g##x[j] -= coef[nEE]*(dPQIJ+PdQIJ+PQdIJ);

				GRAD1D(x,X,x0,l,y,Y,z,Z);
				GRAD1D(y,Y,y0,m,z,Z,x,X);
				GRAD1D(z,Z,z0,n,x,X,y,Y);

				nEE++;
		PQIJ_SHELL_END

			// reset shell index to proper value in this shell
			p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		}
	}
	}

	// translate to atomic center
	for(p=0; p < nBasis; p++){
		Gx[basis2Atom[p]] += gx[p];
		Gy[basis2Atom[p]] += gy[p];
		Gz[basis2Atom[p]] += gz[p];
	}

	// clean memory
	free(coef);
	free(shellMap);
	free(shellMaxL);
	free(Bx);
	free(By);
	free(Bz);
	free(F);
	free(BBF);
	free(Schwarz_Shell);
	free(PM);
	free(Schwarz);
	free(gx);
	free(gy);
	free(gz);
}
#undef MAXSHELLINT 
#undef MAXL         
#undef MAXL1        
#undef MAXL1PWR2   
#undef MAXL1PWR3   
#undef MAXL1PWR4  
#undef MAXL1PWR5 
#undef MAXL1PWR6 
#undef TWOMAXL1     
#undef FOURMAXL1   


// printForce() : print force data
//
// Apr 6, 2010 - Theerapon Khamla
//   Initial implementation
//
// Oct 20, 2010 - Teepanis Chachiyo
//   Add to Siam Quantum source code
//
void printForce(const struct Molecule_t * mol,       // molecule info
                const double *Fx,                    // force in x direction 
                const double *Fy,                    // force in y direction
                const double *Fz){                   // force in z direction

	int A;      // atom index;

	printf(
"\n"
"-------------------------------------------------------------\n"
" Center  Atomic                Forces (Hartrees/Bohr)        \n"
" Number  Number          Fx             Fy             Fz    \n"
"-------------------------------------------------------------\n");
	for(A=0; A < mol->nAtom; A++){
		printf(
		       "%5d  %5d    %15.9f%15.9f%15.9f\n",A+1,
		                                             mol->Z[A],
		                                             Fx[A],
		                                             Fy[A],
		                                             Fz[A]);
	}
	printf(
"-------------------------------------------------------------\n");
}


// uhf_force : computes net force acting on nuclei. Set pointer Fx,Fy,Fz to
// NULL if there is no need to return value. 
//
// Aug 18, 2019 - Teepanis Chachiyo
//     Add support for LibXC: GGA and Hybrid functionals
//
// May 11, 2016 - Teepanis Chachiyo
//     Add exchange weight Kfactor to support DFT
//
// March 14, 2016 - Teepanis Chachiyo
//     Force due to external electric field
//
// March 8, 2013 - Teepanis Chachiyo
//     Parallel version
//
// Oct 20, 2010 - Teepanis Chachiyo
//     Added to Siam Quantum source code
//
// June 5, 2010 - Theerapon Khamla
//     Iinitial implementation
//
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
	double *Fz){             // returned force in z direction

	double *QT=NULL;         // energy weighted density matrix
	double *PA=NULL;         // alpha electron density matrix
	double *PB=NULL;         // beta electron density matrix
	double *PT=NULL;         // total density matrix
	double *PS=NULL;         // spin density matrix
	
	double *fx=NULL;         // calculated force in x direction
	double *fy=NULL;         // calculated force in y direction
	double *fz=NULL;         // calculated force in z direction

	int *basis2Atom;         // mapping basis 2 atom index
	int p,q;                 // basis looping indexes
	int i;                   // generic looping indexes
	double ex,ey,ez;         // electric field integral

	double Kfactor=1.0;      // for DFT hybrid/LDA/GGA

	// print out title section
	printf(
	"                                                             \n"
	"                                                             \n"
	"-------------------------------------------------------------\n"
	"-----                 FORCE CALCULATIONS                -----\n"
	"-------------------------------------------------------------\n"
	);
	fflush(stdout);

#define ALLOC(array,item,type)                                         \
array=calloc(item,sizeof(type));                                       \
if(array==NULL){                                                       \
	printf("uhf_force - error cannot allocate memory\n");              \
	exit(-1);                                                          \
}

	// memory allocations
	ALLOC(QT, nBasis*nBasis, double);
	ALLOC(PA, nBasis*nBasis, double);
	ALLOC(PB, nBasis*nBasis, double);
	ALLOC(PT, nBasis*nBasis, double);
	ALLOC(PS, nBasis*nBasis, double);
	ALLOC(fx, mol->nAtom,    double);
	ALLOC(fy, mol->nAtom,    double);
	ALLOC(fz, mol->nAtom,    double);
	ALLOC(basis2Atom, nBasis,   int);

#undef ALLOC

	// determine atom center mapping
	for(p=0; p < nBasis; p++){
		for(i=0; i < mol->nAtom; i++)
			if(gto[p].x0==mol->x[i] && 
			   gto[p].y0==mol->y[i] && 
			   gto[p].z0==mol->z[i]){
				basis2Atom[p] = i;
				break;
			}
	}

	// compute various types of density matrix
	uhf_getDMatrix(nBasis, nEA, CA, PA);
	uhf_getDMatrix(nBasis, nEB, CB, PB);
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		PT[p*nBasis+q] = PA[p*nBasis+q] + PB[p*nBasis+q];
		PS[p*nBasis+q] = PA[p*nBasis+q] - PB[p*nBasis+q];
	}

	// compute QT matrix (energy weighted density matrix)
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		QT[p*nBasis+q] = 0.0;
		// add alpha spin contribution
		for(i=0; i < nEA; i++)
			QT[p*nBasis+q] += eA[i]*CA[i*nBasis+p]*CA[i*nBasis+q];
		// add beta spin contribution
		for(i=0; i < nEB; i++)
			QT[p*nBasis+q] += eB[i]*CB[i*nBasis+p]*CB[i*nBasis+q];
	}

	// set force to zero initially
	for(i=0; i < mol->nAtom; i++){ fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0; }

	
	// add nuclei-nuclei contribution
	GradVnn(mol,fx,fy,fz);

	// compute eigen energy weighted density matrix contribution
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		// add contribution
		fx[basis2Atom[p]] -= GTO_grad_overlap(p,q,gto,GRAD_X) * QT[p*nBasis+q];
		fx[basis2Atom[q]] -= GTO_grad_overlap(q,p,gto,GRAD_X) * QT[p*nBasis+q];

		fy[basis2Atom[p]] -= GTO_grad_overlap(p,q,gto,GRAD_Y) * QT[p*nBasis+q];
		fy[basis2Atom[q]] -= GTO_grad_overlap(q,p,gto,GRAD_Y) * QT[p*nBasis+q];

		fz[basis2Atom[p]] -= GTO_grad_overlap(p,q,gto,GRAD_Z) * QT[p*nBasis+q];
		fz[basis2Atom[q]] -= GTO_grad_overlap(q,p,gto,GRAD_Z) * QT[p*nBasis+q];
	}

	// compute derivative of kinetic integral contribution
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		// add contribution
		fx[basis2Atom[p]] += GTO_grad_kinetic(p,q,gto,GRAD_X) * PT[p*nBasis+q];
		fx[basis2Atom[q]] += GTO_grad_kinetic(q,p,gto,GRAD_X) * PT[p*nBasis+q];

		fy[basis2Atom[p]] += GTO_grad_kinetic(p,q,gto,GRAD_Y) * PT[p*nBasis+q];
		fy[basis2Atom[q]] += GTO_grad_kinetic(q,p,gto,GRAD_Y) * PT[p*nBasis+q];

		fz[basis2Atom[p]] += GTO_grad_kinetic(p,q,gto,GRAD_Z) * PT[p*nBasis+q];
		fz[basis2Atom[q]] += GTO_grad_kinetic(q,p,gto,GRAD_Z) * PT[p*nBasis+q];
	}

	// compute derivative of nuclear attraction integral contribution
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		// derivative of the gaussian function
		// p and q contribution is symmetric in the case of
		// summation although each term is not.
		fx[basis2Atom[p]] += 2.0 * GTO_grad_nai(p,q,gto,mol,GRAD_X) * PT[p*nBasis+q];
		//fx[basis2Atom[q]] += GTO_grad_nai(q,p,gto,mol,GRAD_X) * PT[p*nBasis+q];

		fy[basis2Atom[p]] += 2.0 * GTO_grad_nai(p,q,gto,mol,GRAD_Y) * PT[p*nBasis+q];
		//fy[basis2Atom[q]] += GTO_grad_nai(q,p,gto,mol,GRAD_Y) * PT[p*nBasis+q];

		fz[basis2Atom[p]] += 2.0 * GTO_grad_nai(p,q,gto,mol,GRAD_Z) * PT[p*nBasis+q];
		//fz[basis2Atom[q]] += GTO_grad_nai(q,p,gto,mol,GRAD_Z) * PT[p*nBasis+q];
	}
	// derivative of the operator
	for(p=0; p < nBasis; p++)
	for(q=0; q <=p; q++){
		for(i=0; i < mol->nAtom; i++){
			GTO_efi(p,q,gto,mol->x[i],mol->y[i],mol->z[i],&ex,&ey,&ez);
			if(p!=q){ ex = 2.0*ex; ey = 2.0*ey; ez = 2.0*ez; }
			fx[i] += ex * mol->Z[i] * PT[p*nBasis+q];
			fy[i] += ey * mol->Z[i] * PT[p*nBasis+q];
			fz[i] += ez * mol->Z[i] * PT[p*nBasis+q];
		}
	}

	if(opt->Ex != 0.0 || opt->Ey != 0.0 || opt->Ez != 0.0)
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		fx[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,1,0,0,0.0,0.0,0.0,GRAD_X) * PT[p*nBasis+q] * opt->Ex;
		fy[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,1,0,0,0.0,0.0,0.0,GRAD_Y) * PT[p*nBasis+q] * opt->Ex;
		fz[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,1,0,0,0.0,0.0,0.0,GRAD_Z) * PT[p*nBasis+q] * opt->Ex;

		fx[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,0,1,0,0.0,0.0,0.0,GRAD_X) * PT[p*nBasis+q] * opt->Ey;
		fy[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,0,1,0,0.0,0.0,0.0,GRAD_Y) * PT[p*nBasis+q] * opt->Ey;
		fz[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,0,1,0,0.0,0.0,0.0,GRAD_Z) * PT[p*nBasis+q] * opt->Ey;

		fx[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,0,0,1,0.0,0.0,0.0,GRAD_X) * PT[p*nBasis+q] * opt->Ez;
		fy[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,0,0,1,0.0,0.0,0.0,GRAD_Y) * PT[p*nBasis+q] * opt->Ez;
		fz[basis2Atom[p]] += 2.0 * GTO_grad_moment(p,q,gto,0,0,1,0.0,0.0,0.0,GRAD_Z) * PT[p*nBasis+q] * opt->Ez;
	}
	for(i=0; i < mol->nAtom; i++){
		fx[i] -= opt->Ex * mol->Z[i];
		fy[i] -= opt->Ey * mol->Z[i];
		fz[i] -= opt->Ez * mol->Z[i];
	}

	// determine Kfactor based on types DFT method
	switch(opt->method){
	case METHOD_HF:        Kfactor = 1.0; break;
	case METHOD_SVWN5:
	case METHOD_SCHACHIYO:
	case METHOD_CHACHIYO:
	case METHOD_xBECKE88:
		                   Kfactor = 0.0; break;
	case METHOD_HALF:      Kfactor = 0.5; break;
#ifdef LIBXC
	case METHOD_LIBXC:     Kfactor = opt->hyb_hf_coef; break;
#endif
	}

	// Note: There are 4 possible versions. They should give the same results
	// (up to 6 significant figures); but ShellSet is the fastest one.
	// Teepanis Oct 21, 2010
	//

	// add contribution from 2-electron integral
	//GradEri             (nBasis, gto, basis2Atom, PT, PS, fx, fy, fz);
	//GradEri_Schwarz     (nBasis, gto, basis2Atom, PT, PS, fx, fy, fz);
	//GradEri_Schwarz_Sym (nBasis, gto, basis2Atom, PT, PS, fx, fy, fz);
	//GradEri_ShellSet    (nBasis, gto, basis2Atom, PT, PS, fx, fy, fz);

	//
	// parallel version
	//
	double *fxSet, *fySet, *fzSet;   // set of fx,fy,fz for each cpu
	int *status;                     // status for each cpu
	int alldone;                     // all idle flag

	// memory callocation
	status=calloc(opt->nCPU,sizeof(int));
	fxSet =calloc(opt->nCPU * mol->nAtom,sizeof(double));
	fySet =calloc(opt->nCPU * mol->nAtom,sizeof(double));
	fzSet =calloc(opt->nCPU * mol->nAtom,sizeof(double));
	if(status==NULL || fxSet==NULL || fySet==NULL || fzSet==NULL){
		printf("uhf_force - error cannot allocate memory\n");
		exit(-1);
	}

	// reset status to idle
	for(i=(opt->nCPU-1);i>=0;i--) status[i] = RPC_IDLE;

	// loop thru all cpu and compute fx,fy,fz
	do{
		// remote call to all cpu except childID=0
		for(i=(opt->nCPU-1);i>0;i--)
			if(status[i] != RPC_DONE)
			status[i] =
			rpc_GradEri_ShellSet_Parallel(status[i], mol->nAtom, i, opt->nCPU, nBasis,
			                              gto, basis2Atom, PT, PS, Kfactor,
			                              fxSet+i*mol->nAtom,
			                              fySet+i*mol->nAtom,
			                              fzSet+i*mol->nAtom,
			                              opt);

		// local call for childID=0
		if(status[0]==RPC_IDLE){
			GradEri_ShellSet_Parallel(0, opt->nCPU, nBasis,
			                          gto, basis2Atom, PT, PS, Kfactor, fxSet, fySet, fzSet);
			status[0]=RPC_DONE;
		}

		// check if all done
		alldone=1;
		for(i=(opt->nCPU-1);i>=0;i--) if(status[i] != RPC_DONE) alldone=0;

	}while(!alldone);

	// accumulate fx,fy,fz
	for(i=(opt->nCPU-1);i>=0;i--){
		for(p=0; p < mol->nAtom; p++){
			fx[p] += fxSet[mol->nAtom*i + p];
			fy[p] += fySet[mol->nAtom*i + p];
			fz[p] += fzSet[mol->nAtom*i + p];
		}
	}

	// free memory
	free(status);
	free(fxSet);
	free(fySet);
	free(fzSet);

	// force is negative of gradient
	for(i=0; i < mol->nAtom; i++){
		fx[i] = -fx[i];
		fy[i] = -fy[i];
		fz[i] = -fz[i];
	}

	//
	// compute DFT forces
	//
	struct MolGrid_t *grid;
	struct gridBasisCenter_t *gB;
	//struct BGrid_t *bg;
	double *Exc;

	switch(opt->method){
	case METHOD_SVWN5:
	case METHOD_HALF:
	case METHOD_SCHACHIYO:
	case METHOD_CHACHIYO:
	case METHOD_xBECKE88:
#ifdef LIBXC
	case METHOD_LIBXC:
#endif

		// build full molecular grid
		grid = genMolGrid(opt->whichGrid, opt->nradial, opt->nlabedev,
		                  nBasis, gto, mol, BGRID_CUTOFF);
		gB = genGridBasisCenter(nBasis, gto, mol, grid);
		//bg = genBasisGrid(nBasis, gto, grid, BGRID_CUTOFF);

		// allocation
		Exc = calloc(grid->nPoint, sizeof(double));
		if(Exc==NULL){
			printf("uhf: Error - Cannot allocate memory\n");
			exit(EXIT_FAILURE);
		}

		// compute exchange correlation forces
		//getXC(nBasis, PA, PB, grid, bg, opt, gto, mol,
		//      RHO_CUTOFF,0,1,
		//      Exc, NULL, NULL,
		//      fx, fy, fz);
		//getXCDirect(nBasis, nEA, nEB, CA, CB, PA, PB, grid, gB, opt, gto, mol,
		//            RHO_CUTOFF,0,1,
		//            Exc, NULL, NULL,
		//            fx, fy, fz);
		getXCDirect_LibXCcompat(nBasis, nEA, nEB, CA, CB, PA, PB, grid, gB, opt, gto, mol,
		            RHO_CUTOFF,0,1,
		            Exc, NULL, NULL,
		            fx, fy, fz);

		// clean memory
		grid = cleanMolGrid(grid);
		gB = cleanGridBasisCenter(nBasis, gB);
		//bg = cleanBGrid(bg);
		free(Exc);
	break;
	}

	// print out for final results
	printForce(mol,fx,fy,fz);
	fflush(stdout);

	// copy calculated forces if requsted
	for(i=0; i < mol->nAtom; i++){
		if(Fx!=NULL) Fx[i] = fx[i];
		if(Fy!=NULL) Fy[i] = fy[i];
		if(Fz!=NULL) Fz[i] = fz[i];
	}

	// free memory
	free(QT);
	free(PA);
	free(PB);
	free(PT);
	free(PS);
	free(fx);
	free(fy);
	free(fz);
	free(basis2Atom);
}
