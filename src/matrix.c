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
#include <string.h>
#include "int.h"
#include "basis.h"
#include "mol.h"
#include "multipole.h"
#include "matrix.h"

//
// GTO_overlap : computes overlap integral between
// two Gaussian basis functions. It has the following form:
// < GTO_i | GTO_j >
//
// Jan 20, 2008 - Teepanis Chachiyo
//     Iinitial implementation
//
double GTO_overlap(int i,                          // ith basis
                   int j,                          // jth basis
                   const struct GTOBasis_t *gto){  // basis database

	int iCnt, jCnt;   // contracted function 
	double sum=0.0;   // integral sum

	// looping over contracted functions
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++){
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++){
		sum = sum + gto[i].coef[iCnt] * gto[j].coef[jCnt] *
		            gto[i].norm[iCnt] * gto[j].norm[jCnt] *
	        overlap(gto[i].exp[iCnt], gto[i].l,  gto[i].m,  gto[i].n,
	                                  gto[i].x0, gto[i].y0, gto[i].z0,
	                gto[j].exp[jCnt], gto[j].l,  gto[j].m,  gto[j].n,
		                          gto[j].x0, gto[j].y0, gto[j].z0);
	}
	}
	return sum;
}


//
// GTO_moment : computes moment integral between
// two Gaussian basis functions. It has the following form:
//
//                mx     my     mz
// < GTO_i | (x-xc) (y-yc) (z-zc)  | GTO_j >
//
// 2011 - Aniwat Kesorn
//     Adapted from GTO_overlap
//
// Oct 13, 2012 - Teepanis Chachiyo Aniwat Kesorn
//     Add to Siam Quantum
//
double GTO_moment(int i,                            // ith basis
                  int j,                            // jth basis
                  const struct GTOBasis_t *gto,     // basis database
                  int mx, int my, int mz,           // moment order
                  double xc, double yc, double zc){ // moment center

	int iCnt, jCnt;   // contracted function 
	double sum=0.0;   // integral sum

	// looping over contracted functions
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++){
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++){
		sum = sum + gto[i].coef[iCnt] * gto[j].coef[jCnt] *
		            gto[i].norm[iCnt] * gto[j].norm[jCnt] *
	        moment( gto[i].exp[iCnt], gto[i].l,  gto[i].m,  gto[i].n,
	                                  gto[i].x0, gto[i].y0, gto[i].z0,
	                gto[j].exp[jCnt], gto[j].l,  gto[j].m,  gto[j].n,
		                          gto[j].x0, gto[j].y0, gto[j].z0,
		            mx, my, mz, xc, yc, zc);
	}
	}
	return sum;
}


// GTO_kinetic : computes kinetic operator integral between
// two Gaussian basis functions. It has the following form:
// < GTO_i | -1/2 Laplacian | GTO_j >
//
// Feb 22, 2008 - Teepanis Chachiyo
//     Initial implementation
//
double GTO_kinetic(int i,                          // ith basis
                   int j,                          // jth basis
                   const struct GTOBasis_t *gto){  // basis database

	int iCnt, jCnt;  // contracted function 
	double sum=0.0;  // integral sum

	// looping over contracted functions
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++){
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++){
		sum = sum + gto[i].coef[iCnt]  * gto[j].coef[jCnt] *
		            gto[i].norm[iCnt]  * gto[j].norm[jCnt]  *
	        kinetic(gto[i].exp[iCnt], gto[i].l,  gto[i].m,  gto[i].n,
	                                  gto[i].x0, gto[i].y0, gto[i].z0,
	                gto[j].exp[jCnt], gto[j].l,  gto[j].m,  gto[j].n,
		                          gto[j].x0, gto[j].y0, gto[j].z0);
	}
	}
	return sum;
}


// GTO_nuclei : computes nuclei attraction integral between
// two Gaussian basis functions. It has the following form:
// < GTO_i |   SUM Z_k/(r-R_k) |  GTO_j >
// where index k is a sum over all nuclei
//
// Jan 29, 2008 - Teepanis Chachiyo
//     Initial implementation
//
double GTO_nuclei(int i,                          // ith basis
                  int j,                          // jth basis
                  const struct GTOBasis_t *gto,   // basis database
                  const struct Molecule_t *mol){  // molecule database

	int k;           // nucleus index
	int iCnt, jCnt;  // contracted functions
	double sum=0.0;  // integral sum

	// looping over contracted functions
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++){
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++){
			
		// looping over nuclei
		for(k=0; k < mol->nAtom; k++){
			sum = sum +  gto[i].coef[iCnt]  * gto[j].coef[jCnt] *
			             gto[i].norm[iCnt]  * gto[j].norm[jCnt] *
			      (double)mol->Z[k]                             *
			        nai(gto[i].x0, gto[i].y0, gto[i].z0,
			            1.0,
			            gto[i].l,  gto[i].m,  gto[i].n,
			            gto[i].exp[iCnt],
			            gto[j].x0, gto[j].y0, gto[j].z0,
			            1.0,
			            gto[j].l,  gto[j].m,  gto[j].n,
			            gto[j].exp[jCnt],
			            mol->x[k], mol->y[k], mol->z[k]);
		}
	}
	}
	return sum;

}


// GTO_nai: compute nuclear attraction integral between two Gaussian basis
// functions. It has the following form:
//
//                           (i| 1/rC|j)
//
// Oct 26, 2010 - Nanta Sophonrat
//         Initial Implementation     
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//         Added to Siam Quantum source tree
//
double GTO_nai(int i,                            // index of GTO_tho
               int j,                            // index of GTO_tho
               const struct GTOBasis_t *gto,     // basis function info.
               double xc, double yc, double zc){ // position to find potential

	int    iCnt, jCnt;        // index for contracted GTO
	double sum = 0.0;         // potential

	for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++){
		sum += gto[i].coef[iCnt]*gto[i].norm[iCnt]*
		       gto[j].coef[jCnt]*gto[j].norm[jCnt]*
		       nai( gto[i].x0, gto[i].y0, gto[i].z0, 1.0,
		            gto[i].l,  gto[i].m,  gto[i].n,  gto[i].exp[iCnt],
		            gto[j].x0, gto[j].y0, gto[j].z0, 1.0,
		            gto[j].l,  gto[j].m,  gto[j].n,  gto[j].exp[jCnt],
		            xc, yc, zc);
	}
	return sum;
}

// GTO_efi : computes electric field integral between
// two Gaussian basis functions. It has the following form:
//
// (i|ddxC(1/rC)|j), (i|ddyC(1/rC)|j), (i|ddzC(1/rC)|j)
//
// 2010 - Theerapon Khumla
//     Initial implementation
//
void GTO_efi(int i,                               // ith basis
             int j,                               // jth basis
             const struct GTOBasis_t *gto,        // basis database
             double xc, double yc, double zc,     // point to evaluate efi
             double *ex, double *ey, double *ez){ // returned value

	int iCnt, jCnt;   // contracted function 
	double sx,sy,sz;  // integral sum
	double p;         // prefactor

	// looping over contracted functions
	sx = 0.0; sy = 0.0; sz = 0.0;
	for(iCnt=0; iCnt < gto[i].nContract; iCnt++)
	for(jCnt=0; jCnt < gto[j].nContract; jCnt++){

		p = gto[i].coef[iCnt] * gto[j].coef[jCnt] *
		    gto[i].norm[iCnt] * gto[j].norm[jCnt];

		efi(gto[i].x0, gto[i].y0, gto[i].z0,
		    1.0,
		    gto[i].l,  gto[i].m,  gto[i].n,
		    gto[i].exp[iCnt],
		    gto[j].x0, gto[j].y0, gto[j].z0,
		    1.0,
		    gto[j].l,  gto[j].m,  gto[j].n,
		    gto[j].exp[jCnt],
		    xc, yc, zc,
		    ex, ey, ez);

		sx += *ex * p;
		sy += *ey * p;
		sz += *ez * p;
	}

	*ex = sx; *ey = sy; *ez = sz;
}

// contr_coulomb : compute contracted integral using THO.
// The coding structure is from PyQuante.
//
// Feb 19, 2008 - Teepanis Chachiyo
//     Original implementation.
//
double contr_eri(
             int lena,double *aexps,double *acoefs,double *anorms,
             double xa,double ya,double za,int la,int ma,int na,
             int lenb,double *bexps,double *bcoefs,double *bnorms,
             double xb,double yb,double zb,int lb,int mb,int nb,
             int lenc,double *cexps,double *ccoefs,double *cnorms,
             double xc,double yc,double zc,int lc,int mc,int nc,
             int lend,double *dexps,double *dcoefs,double *dnorms,
             double xd,double yd,double zd,int ld,int md,int nd){
			 
	double val = 0.;
	int i,j,k,l;
	double EE;

	// proceed from lowest exponent value
	//for (i=lena-1; i>=0; i--)
	//for (j=lenb-1; j>=0; j--)
	//for (k=lenc-1; k>=0; k--)
	//for (l=lend-1; l>=0; l--){
	//	EE =    acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]*
	//	        eri(xa,ya,za,anorms[i],la,ma,na,aexps[i],
	//	            xb,yb,zb,bnorms[j],lb,mb,nb,bexps[j],
	//	            xc,yc,zc,cnorms[k],lc,mc,nc,cexps[k],
	//	            xd,yd,zd,dnorms[l],ld,md,nd,dexps[l]);
	//	val += EE;
	//}

	// proceed from highest exponent value
	for (i=0; i<lena; i++)
	for (j=0; j<lenb; j++)
	for (k=0; k<lenc; k++)
	for (l=0; l<lend; l++){
		// compute element
		EE  = acoefs[i]*bcoefs[j]*ccoefs[k]*dcoefs[l]
		      * eri(xa,ya,za,anorms[i],la,ma,na,aexps[i],
		            xb,yb,zb,bnorms[j],lb,mb,nb,bexps[j],
		            xc,yc,zc,cnorms[k],lc,mc,nc,cexps[k],
		            xd,yd,zd,dnorms[l],ld,md,nd,dexps[l]);
		val += EE;
	}
	
	return val;
}

// contr_hrr : compute contracted integral using HGP horizontal recursive
// The coding structure is from PyQuante.
//
// July 1, 2010 - Teepanis Chachiyo
//                Original implementation.
//                Initial test suggested that contr_eri is faster (on singlet
//                oxygen)
//
double contr_hrr(
             int lena,double *aexps,double *acoefs,double *anorms,
             double xa,double ya,double za,int la,int ma,int na,
             int lenb,double *bexps,double *bcoefs,double *bnorms,
             double xb,double yb,double zb,int lb,int mb,int nb,
             int lenc,double *cexps,double *ccoefs,double *cnorms,
             double xc,double yc,double zc,int lc,int mc,int nc,
             int lend,double *dexps,double *dcoefs,double *dnorms,
             double xd,double yd,double zd,int ld,int md,int nd){
	if(lb>0){
		return contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la+1, ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb-1, mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc,   mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld,   md, nd)
		       + (xa-xb) *          
		       contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la,   ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb-1, mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc,   mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld,   md, nd);
	}else if(mb>0){
		return contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma+1, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb-1, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc,   nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md,   nd)
		       + (ya-yb) *          
		       contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma,   na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb-1, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc,   nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md,   nd);
	}else if(nb>0){
		return contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma, na+1,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb, nb-1,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md, nd)
		       + (za-zb) *          
		       contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb, nb-1,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md, nd);		
	}else if(ld>0){
		return contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc+1, mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld-1, md, nd)
		       + (xc-xd) *          
		       contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la,   ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb,   mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc,   mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld-1, md, nd);
	}else if(md>0){
		return contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma,   na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb,   nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc+1, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md-1, nd)
		       + (yc-yd) *          
		       contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md-1, nd);
	}else if(nd>0){
		return contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc, nc+1,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md, nd-1)
		       + (zc-zd) *          
		       contr_hrr(lena, aexps, acoefs, anorms, xa, ya, za, la, ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md, nd-1);
	}else
		return contr_eri(lena, aexps, acoefs, anorms, xa, ya, za, la, ma, na,
		                 lenb, bexps, bcoefs, bnorms, xb, yb, zb, lb, mb, nb,
		                 lenc, cexps, ccoefs, cnorms, xc, yc, zc, lc, mc, nc,
		                 lend, dexps, dcoefs, dnorms, xd, yd, zd, ld, md, nd);
}

// create_Schwarz : computes the Schwarz matrix and return pointer to the
// matrix. This technique essentially reduces that size from N**4 problem
// to a managable N**2 problem.
//
// Douglas L. Strout and Gustavo E. Scuseria. "A quantitative study of the
// scaling properties of the Hartree-Fock method." J. Chem. Phys. (1995)
// Vol 102 page 8448
//
// Feb 18, 2008 - Teepanis Chachiyo
//     Initial implementation
//
double * create_Schwarz(
	int nBasis,                     // number of basis functions
	const struct GTOBasis_t *gto){  // basis set info

	int p,q;        // loop index
	double *sch;    // Schwarz matrix pointer
	double upBound; // upper bound from Schwarz inequality

	// allocate memory
	sch = calloc(nBasis*nBasis, sizeof(double));
	if(sch==NULL){
		printf("Cannot allocate Schwarz matrix\n");
		exit(-1);
	}

	// compute the matrix
	// There is a sch[p,q] and sch[q,p] symmetry, so we 
	// need to take advantage of this property to gain
	// speed in the future.
	//
	// Feb 22, 2008 - Teepanis Chachiyo
	//
	for(p=0; p < nBasis; p++){
	for(q=0; q < nBasis; q++){

		// evaluate integral
		upBound =
		      contr_eri(
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
		             gto[q].l,         gto[q].m,   gto[q].n);

		// make sure we have positive value
		sch[p*nBasis+q] = sqrt(fabs(upBound));
	}
	}
	
	return sch;
}

// GTO_2JK_Matrix : compute 2J-K matrix elements taking into account of 
// the symmetry in the Coulomb integrals. This type of matrix is for
// closed-shell calculations.
//
// It computes repulsion integral between molecular orbital
// described by density matrix P and a pair of basis i,j. It has the following
// form: [RETURN]pq = SUM Pij*(2 EE(piqj) - EE(pijq)) 
//
// Note: For closed-shell system: Gpq = [RETURN]pq / 2.0
//
// Important Definition
// --------------------
// 
// (piqj) = Int[  dr1 dr2 phi_p*(r1) phi_i*(r2) (1/r12) phi_q(r1) phi_j(r2) ]
//
// July 10, 2010 - Teepanis Chachiyo
//      Convert name to 2JK
//
// Jan 29, 2008 - Teepanis Chachiyo
//     Original Implementation
//
// Feb 18, 2008 - Teepanis Chachiyo
//     Optimization on G[p,q] multiplications
//
void GTO_2JK_Matrix(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G){                  // return matrix

	int p,q,i,j,maxj;       // looping variables
	double EE;              // coulomb integral
	double upBound;         // Schwarz upper bound

	// algorithm based on my algorithm during graduate years
	// which rooted from Thijssen.
	// This is for looping all unique (piqj)
	// HF.tar.gz  -- Teepanis Chachiyo
	for(p=0; p < nBasis; p++){
		for(q=0; q < p+1; q++){
			for(i=0; i < p+1; i++){
				if(i == p){
					maxj = q+1;
				}else{
					maxj = i+1;
				}
				for(j=0; j < maxj; j++){
#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
					// screen integral using Schwarz inequality
					upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
					if(upBound < SCHWARZ_CUTOFF) continue;

					// compute two-electron integral
					EE =
					contr_eri(
					     gto[p].nContract,
					     gto[p].exp, gto[p].coef, gto[p].norm,
					     gto[p].x0,        gto[p].y0,  gto[p].z0,
					     gto[p].l,         gto[p].m,   gto[p].n,
					     gto[q].nContract,
					     gto[q].exp, gto[q].coef, gto[q].norm,
					     gto[q].x0,        gto[q].y0,  gto[q].z0,
					     gto[q].l,         gto[q].m,   gto[q].n,
					     gto[i].nContract,
					     gto[i].exp, gto[i].coef, gto[i].norm,
					     gto[i].x0,        gto[i].y0,  gto[i].z0,
					     gto[i].l,         gto[i].m,   gto[i].n,
					     gto[j].nContract,
					     gto[j].exp, gto[j].coef, gto[j].norm,
					     gto[j].x0,        gto[j].y0,  gto[j].z0,
					     gto[j].l,         gto[j].m,   gto[j].n);

					// perform matrix symmetry rotations
					if((p==q)&&(i==j)&&(p==i)){  // all same
						G(p,q) = G(p,q) +     EE*P(i,j);
					}else if((p==q)&&(i==j)){    // 2 pairs
						G(p,p) = G(p,p) + 2.0*EE*P(i,i);
						G(i,i) = G(i,i) + 2.0*EE*P(p,p);

						G(p,i) = G(p,i) -     EE*P(i,p);
						G(i,p) = G(i,p) -     EE*P(i,p);
					}else if(p==q){              // pq pair
						G(p,p) = G(p,p) + 4.0*EE*P(i,j);

						G(i,j) = G(i,j) + 2.0*EE*P(p,p);
						G(j,i) = G(j,i) + 2.0*EE*P(p,p);

						G(p,j) = G(p,j) -     EE*P(i,p);
						G(j,p) = G(j,p) -     EE*P(i,p);

						G(i,p) = G(i,p) -     EE*P(p,j);
						G(p,i) = G(p,i) -     EE*P(p,j);
					}else if(i==j){              // ij pair
						G(i,i) = G(i,i) + 4.0*EE*P(p,q);

						G(p,q) = G(p,q) + 2.0*EE*P(i,i);
						G(q,p) = G(q,p) + 2.0*EE*P(i,i);

						G(p,i) = G(p,i) -     EE*P(i,q);
						G(i,p) = G(i,p) -     EE*P(i,q);

						G(q,i) = G(q,i) -     EE*P(i,p);
						G(i,q) = G(i,q) -     EE*P(i,p);
					}else if((p==i)&&(q==j)){    // pi-qj pair
						G(p,q) = G(p,q) + 3.0*EE*P(p,q);
						G(q,p) = G(q,p) + 3.0*EE*P(p,q);

						G(p,p) = G(p,p) -     EE*P(q,q);
						G(q,q) = G(q,q) -     EE*P(p,p);
					}else{                       // all distinct
						G(p,q) = G(p,q) + 4.0*EE*P(i,j);
						G(q,p) = G(q,p) + 4.0*EE*P(i,j);

						G(i,j) = G(i,j) + 4.0*EE*P(p,q);
						G(j,i) = G(j,i) + 4.0*EE*P(p,q);

						G(p,j) = G(p,j) -     EE*P(i,q);
						G(j,p) = G(j,p) -     EE*P(i,q);

						G(p,i) = G(p,i) -     EE*P(j,q);
						G(i,p) = G(i,p) -     EE*P(j,q);

						G(q,i) = G(q,i) -     EE*P(j,p);
						G(i,q) = G(i,q) -     EE*P(j,p);

						G(q,j) = G(q,j) -     EE*P(i,p);
						G(j,q) = G(j,q) -     EE*P(i,p);
					}
				}
			}
		}
	}
#undef G
#undef P
}

// GTO_JK_Matrix_NoSym : compute 2J-K matrix elements without taking 
// advantange of the rotation symmetry.
//
// It computes repulsion integral between molecular orbital
// described by density matrix P and a pair of basis i,j. It has the following
// form: [RETURN]pq = SUM Pij*(2 EE(piqj) - EE(pijq)) 
//
// Note: For closed-shell system: Gpq = [RETURN]pq / 2.0
//
// Important Definition
// --------------------
// 
// (piqj) = Int[  dr1 dr2 phi_p*(r1) phi_i*(r2) (1/r12) phi_q(r1) phi_j(r2) ]
//
// July 10, 2010 - Teepanis Chachiyo
//     Convert name to 2JK
//
// Mar 08, 2010 - Teepanis Chachiyo
//     Original Implementation
//
void GTO_2JK_Matrix_NoSymm(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G){                  // return matrix

	int p,q,i,j;                 // looping variables
	double EE;                   // coulomb integral
	double upBound;              // Schwarz upper bound

	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++) 
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){
#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
		// screen integral using Schwarz inequality
		upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
		if(upBound < SCHWARZ_CUTOFF) continue;

		// compute two-electron integral
		EE =
		contr_eri(
			     gto[p].nContract,
			     gto[p].exp, gto[p].coef, gto[p].norm,
			     gto[p].x0,        gto[p].y0,  gto[p].z0,
			     gto[p].l,         gto[p].m,   gto[p].n,
			     gto[q].nContract,
			     gto[q].exp, gto[q].coef, gto[q].norm,
			     gto[q].x0,        gto[q].y0,  gto[q].z0,
			     gto[q].l,         gto[q].m,   gto[q].n,
			     gto[i].nContract,
			     gto[i].exp, gto[i].coef, gto[i].norm,
			     gto[i].x0,        gto[i].y0,  gto[i].z0,
			     gto[i].l,         gto[i].m,   gto[i].n,
			     gto[j].nContract,
			     gto[j].exp, gto[j].coef, gto[j].norm,
			     gto[j].x0,        gto[j].y0,  gto[j].z0,
			     gto[j].l,         gto[j].m,   gto[j].n);

		// I need to explain how to derive this relation
		// so that other students can do it too.
		// - Teepanis Chachiyo (Mar 08,2010)
		G(p,q) = G(p,q) + 2.0*EE*P(i,j);
		G(p,j) = G(p,j) -     EE*P(i,q);

	}
#undef G
#undef P
}


// GTO_JK_Matrix_Sym : compute 2J-K matrix elements by taking 
// advantange of the rotation symmetry. This is another way to perform
// symmetry checking, and it is suitable for shell permutation.
//
// It computes repulsion integral between molecular orbital
// described by density matrix P and a pair of basis i,j. It has the following
// form: [RETURN]pq = SUM Pij*(2 EE(piqj) - EE(pijq)) 
//
// Note: For closed-shell system: Gpq = [RETURN]pq / 2.0
//
// Important Definition
// --------------------
// 
// (piqj) = Int[  dr1 dr2 phi_p*(r1) phi_i*(r2) (1/r12) phi_q(r1) phi_j(r2) ]
//
// July 10, 2010 - Teepanis Chachiyo
//     Convert name to 2JK
//
// Mar 09, 2010 - Teepanis Chachiyo
//     Original Implementation
//
void GTO_2JK_Matrix_Symm(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G){                  // return matrix

	int p,q,i,j,maxj;            // looping variables
	double EE;                   // coulomb integral
	double upBound;              // Schwarz upper bound

	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){

		// check uniqueness
		if(!(q<p+1)) continue;
		if(!(i<p+1)) continue;
		if(i==p) maxj=q+1; else maxj=i+1;
		if(!(j<maxj)) continue;

#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
		// screen integral using Schwarz inequality
		upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
		if(upBound < SCHWARZ_CUTOFF) continue;

		// compute two-electron integral
		EE =
		contr_eri(
			     gto[p].nContract,
			     gto[p].exp, gto[p].coef, gto[p].norm,
			     gto[p].x0,        gto[p].y0,  gto[p].z0,
			     gto[p].l,         gto[p].m,   gto[p].n,
			     gto[q].nContract,
			     gto[q].exp, gto[q].coef, gto[q].norm,
			     gto[q].x0,        gto[q].y0,  gto[q].z0,
			     gto[q].l,         gto[q].m,   gto[q].n,
			     gto[i].nContract,
			     gto[i].exp, gto[i].coef, gto[i].norm,
			     gto[i].x0,        gto[i].y0,  gto[i].z0,
			     gto[i].l,         gto[i].m,   gto[i].n,
			     gto[j].nContract,
			     gto[j].exp, gto[j].coef, gto[j].norm,
			     gto[j].x0,        gto[j].y0,  gto[j].z0,
			     gto[j].l,         gto[j].m,   gto[j].n);

		// same as funtion GTO_JK_Matrix
		if((p==q)&&(i==j)&&(p==i)){  // all same
			G(p,q) = G(p,q) +     EE*P(i,j);
		}else if((p==q)&&(i==j)){    // 2 pairs
			G(p,p) = G(p,p) + 2.0*EE*P(i,i);
			G(i,i) = G(i,i) + 2.0*EE*P(p,p);

			G(p,i) = G(p,i) -     EE*P(i,p);
			G(i,p) = G(i,p) -     EE*P(i,p);
		}else if(p==q){              // pq pair
			G(p,p) = G(p,p) + 4.0*EE*P(i,j);

			G(i,j) = G(i,j) + 2.0*EE*P(p,p);
			G(j,i) = G(j,i) + 2.0*EE*P(p,p);

			G(p,j) = G(p,j) -     EE*P(i,p);
			G(j,p) = G(j,p) -     EE*P(i,p);

			G(i,p) = G(i,p) -     EE*P(p,j);
			G(p,i) = G(p,i) -     EE*P(p,j);
		}else if(i==j){              // ij pair
			G(i,i) = G(i,i) + 4.0*EE*P(p,q);

			G(p,q) = G(p,q) + 2.0*EE*P(i,i);
			G(q,p) = G(q,p) + 2.0*EE*P(i,i);

			G(p,i) = G(p,i) -     EE*P(i,q);
			G(i,p) = G(i,p) -     EE*P(i,q);

			G(q,i) = G(q,i) -     EE*P(i,p);
			G(i,q) = G(i,q) -     EE*P(i,p);
		}else if((p==i)&&(q==j)){    // pi-qj pair
			G(p,q) = G(p,q) + 3.0*EE*P(p,q);
			G(q,p) = G(q,p) + 3.0*EE*P(p,q);

			G(p,p) = G(p,p) -     EE*P(q,q);
			G(q,q) = G(q,q) -     EE*P(p,p);
		}else{                       // all distinct
			G(p,q) = G(p,q) + 4.0*EE*P(i,j);
			G(q,p) = G(q,p) + 4.0*EE*P(i,j);

			G(i,j) = G(i,j) + 4.0*EE*P(p,q);
			G(j,i) = G(j,i) + 4.0*EE*P(p,q);

			G(p,j) = G(p,j) -     EE*P(i,q);
			G(j,p) = G(j,p) -     EE*P(i,q);

			G(p,i) = G(p,i) -     EE*P(j,q);
			G(i,p) = G(i,p) -     EE*P(j,q);

			G(q,i) = G(q,i) -     EE*P(j,p);
			G(i,q) = G(i,q) -     EE*P(j,p);

			G(q,j) = G(q,j) -     EE*P(i,p);
			G(j,q) = G(j,q) -     EE*P(i,p);
		}


	}
#undef G
#undef P
}

// isSameShell determine that two basis function i and j belong
// to the same shell. A definition of shell is a set of basis 
// function with the same center and same contracted exponent,
// for example, sp shell, d shell.
//
// Mar 08, 2010 - Teepanis Chachiyo
//     Original Implementation
//
int isSameShell(
	int nBasis,                    // number of basis functions
	int i, int j,                  // which basis to test
	const struct GTOBasis_t *gto){ // basis set info

	int k;

	// test invalid range
	if(i >= nBasis || j >= nBasis){
		printf("isSameShell - invalid i,j\n");
		exit(-1);
	}

	// test center
	if(gto[i].x0 != gto[j].x0 || 
	   gto[i].y0 != gto[j].y0 || 
	   gto[i].z0 != gto[j].z0)
		return 0;

	// test contracted exponent
	if(gto[i].nContract != gto[j].nContract) return 0;

	for(k=0; k < gto[i].nContract; k++)
		if(gto[i].exp[k] != gto[j].exp[k]) return 0;

	return 1;
}


// GTO_2JK_Matrix_Sym_Shell : compute 2J-K matrix elements by taking 
// advantange of the rotation symmetry, but split the integral into shells.
//
// It computes repulsion integral between molecular orbital
// described by density matrix P and a pair of basis i,j. It has the following
// form: [RETURN]pq = SUM Pij*(2 EE(piqj) - EE(pijq)) 
//
// Note: For closed-shell system: Gpq = [RETURN]pq / 2.0
//
// Important Definition
// --------------------
// 
// (piqj) = Int[  dr1 dr2 phi_p*(r1) phi_i*(r2) (1/r12) phi_q(r1) phi_j(r2) ]
//
// July 10, 2010 - Teepanis Chachiyo
//      Convert name to 2JK
//
// Mar 08, 2010 - Teepanis Chachiyo
//     Original Implementation
//
#define MAXMEMBER 6
#define MAXINT 1296       // maximum number of integral to store in memory
// NOTE: The number of maximum integral for 1296 is 6*6*6*6 which means
// it only support shell which has maximally 6 basis functions such
// as d orbitals. We need to be careful on how we implement this during
// the generating shell mapping part.
// - Teepanis Chachiyo March 13, 2010.
// 
#define INTTHRESHOLD 2    // threshold number of integral to use eri_list
void GTO_2JK_Matrix_Symm_Shell(
	int nBasis,                  // number of basis functions
	const double *P,             // density matrix
	const struct GTOBasis_t *gto,// basis set info
	const double *Schwarz,       // pointer to schwarz matrix
	double cutoff,               // cutoff to ignore
	double *G){                  // return matrix

	int pS,qS,iS,jS;             // looping variables for shells
	int p,q,i,j,maxj;            // looping variables for basis functions
	double EE;                   // coulomb integral
	double upBound;              // Schwarz upper bound

	int nShell;                  // number of shell
	int *shellMap;               // mapping between shell and basis function

	int *iP, *iQ, *iI, *iJ;      // basis index storage

	int *la,*ma,*na;             // Gaussian angular order
	int *lb,*mb,*nb;             // Gaussian angular order
	int *lc,*mc,*nc;             // Gaussian angular order
	int *ld,*md,*nd;             // Gaussian angular order

	double *EEStore;             // full EE storage
	double *EEPrim;              // primitive EE storage

	int nInt;                    // number of integral
	int n;                       // integral index
	int r;                       // generic counter
	int pC,qC,iC,jC;             // contracted function loop

	int maxSumL;

	// allocating memory
	shellMap = calloc(MAX_BASIS, sizeof(int));
	if(shellMap==NULL){
		printf("GTO_JK_Matrix_Symm_Shell - error allocating memory");
		exit(-1);
	}
#define ALLOC calloc(MAXINT, sizeof(int))
	la = ALLOC; lb = ALLOC; lc = ALLOC; ld = ALLOC;
	ma = ALLOC; mb = ALLOC; mc = ALLOC; md = ALLOC;
	na = ALLOC; nb = ALLOC; nc = ALLOC; nd = ALLOC;
	iP = ALLOC; iQ = ALLOC; iI = ALLOC; iJ = ALLOC;
#undef ALLOC

	EEStore = calloc(MAXINT, sizeof(double));
	EEPrim  = calloc(MAXINT, sizeof(double));

	if(la==NULL || lb==NULL || lc==NULL || ld==NULL ||
	   ma==NULL || mb==NULL || mc==NULL || md==NULL ||
	   na==NULL || nb==NULL || nc==NULL || nd==NULL ||
	   iP==NULL || iQ==NULL || iI==NULL || iJ==NULL ||
	   EEStore==NULL        || EEPrim==NULL){
		printf("GTO_JK_Matrix_Symm_Shell - error allocating memory\n");
		exit(-1);
	}

	// generate shell mapping
	// first basis contributes first shell always
	nShell      = 1;
	shellMap[0] = 0;
	for(i=1; i < nBasis; i++)
		if( isSameShell(nBasis,i,shellMap[nShell-1],gto) &&
		    (i-shellMap[nShell-1]) < MAXMEMBER              ) continue;
		else{
			shellMap[nShell] = i;
			nShell++;
		}
	// last item in shellMap must be the nBasis
	shellMap[nShell] = nBasis;

	// permute all possible shells
	for(pS=0; pS < nShell; pS++)
	for(qS=0; qS < (pS+1); qS++)
	for(iS=0; iS < (pS+1); iS++)
	for(jS=0; jS < (pS+1); jS++){

		// reset counter
		nInt = 0;

		// generate all possible integrals
		for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
		for(q=shellMap[qS]; q < shellMap[qS+1]; q++) if(q<(p+1)) // uniqueness
		for(i=shellMap[iS]; i < shellMap[iS+1]; i++) if(i<(p+1)) // uniqueness
		for(j=shellMap[jS]; j < shellMap[jS+1]; j++){

			// check uniqueness
			if(i==p) maxj=q+1; else maxj=i+1;
			if(!(j<maxj)) continue;

			// screen integral using Schwarz inequality
			upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
			if(upBound < SCHWARZ_CUTOFF) continue;

			// load integral info
			la[nInt] = gto[p].l; ma[nInt] = gto[p].m; na[nInt] = gto[p].n;
			lb[nInt] = gto[q].l; mb[nInt] = gto[q].m; nb[nInt] = gto[q].n;
			lc[nInt] = gto[i].l; mc[nInt] = gto[i].m; nc[nInt] = gto[i].n;
			ld[nInt] = gto[j].l; md[nInt] = gto[j].m; nd[nInt] = gto[j].n;

			// load basis index
			iP[nInt] = p; iQ[nInt] = q; iI[nInt] = i; iJ[nInt] = j;

			// increment the number of integral
			nInt++;
		}

		// if the number of integral is too few call simple contracted eri
		if(nInt < INTTHRESHOLD){

			for(n=0; n < nInt; n++){

				// set basis function index
				p = iP[n]; q = iQ[n]; i = iI[n]; j = iJ[n];

				// compute two-electron integral
				EEStore[n] =

				contr_eri(
					     gto[p].nContract,
					     gto[p].exp, gto[p].coef, gto[p].norm,
					     gto[p].x0,  gto[p].y0,   gto[p].z0,
					     gto[p].l,   gto[p].m,    gto[p].n,
					     gto[q].nContract,
					     gto[q].exp, gto[q].coef, gto[q].norm,
					     gto[q].x0,  gto[q].y0,   gto[q].z0,
					     gto[q].l,   gto[q].m,    gto[q].n,
					     gto[i].nContract,
					     gto[i].exp, gto[i].coef, gto[i].norm,
					     gto[i].x0,  gto[i].y0,   gto[i].z0,
					     gto[i].l,   gto[i].m,    gto[i].n,
					     gto[j].nContract,
					     gto[j].exp, gto[j].coef, gto[j].norm,
					     gto[j].x0,  gto[j].y0,   gto[j].z0,
					     gto[j].l,   gto[j].m,    gto[j].n);			

			}

		}else{

			// prepare before going into contraction loop
			maxSumL = 0;
			for(n=0; n < nInt; n++){

				// determine maxSumL
				r = la[n]+lb[n]+lc[n]+ld[n]+
			        ma[n]+mb[n]+mc[n]+md[n]+
			        na[n]+nb[n]+nc[n]+nd[n];
				if(r > maxSumL) maxSumL = r;

				// reset integral to zero
				EEStore[n] = 0.0;

			}

			// set basis function index
			p = iP[0]; q = iQ[0]; i = iI[0]; j = iJ[0];

			// evaluate integral
			for(pC=0; pC < gto[p].nContract; pC++)
			for(qC=0; qC < gto[q].nContract; qC++)
			for(iC=0; iC < gto[i].nContract; iC++)
			for(jC=0; jC < gto[j].nContract; jC++){
	
				// compute the entire list
				eri_list(
				     gto[p].x0,gto[p].y0,gto[p].z0,la,ma,na,gto[p].exp[pC],
				     gto[q].x0,gto[q].y0,gto[q].z0,lb,mb,nb,gto[q].exp[qC],
				     gto[i].x0,gto[i].y0,gto[i].z0,lc,mc,nc,gto[i].exp[iC],
				     gto[j].x0,gto[j].y0,gto[j].z0,ld,md,nd,gto[j].exp[jC],
				     maxSumL, nInt, EEPrim);

				// adding to the final value
				for(n=0; n < nInt; n++)
					EEStore[n] += EEPrim[n] * gto[iP[n]].coef[pC]
					                        * gto[iQ[n]].coef[qC]
					                        * gto[iI[n]].coef[iC]
					                        * gto[iJ[n]].coef[jC]
					                        * gto[iP[n]].norm[pC]
					                        * gto[iQ[n]].norm[qC]
					                        * gto[iI[n]].norm[iC]
					                        * gto[iJ[n]].norm[jC];

			}
		}

		// convert it to G matrix
		for(n=0; n < nInt; n++){

			// load calculated ee integral
			EE = EEStore[n];

			// set basis function index
			p = iP[n]; q = iQ[n]; i = iI[n]; j = iJ[n];

#define G(x,y) G[nBasis*x+y]
#define P(x,y) P[nBasis*x+y]
			// same as funtion GTO_JK_Matrix
			if((p==q)&&(i==j)&&(p==i)){  // all same
				G(p,q) = G(p,q) +     EE*P(i,j);           // p > q always
			}else if((p==q)&&(i==j)){    // 2 pairs
				G(p,p) = G(p,p) + 2.0*EE*P(i,i);           // diagonal
				G(i,i) = G(i,i) + 2.0*EE*P(p,p);           // diagonal

				G(p,i) = G(p,i) -     EE*P(i,p);           // p > i always
//				G(i,p) = G(i,p) -     EE*P(i,p);           // repetitive
			}else if(p==q){              // pq pair
				G(p,p) = G(p,p) + 4.0*EE*P(i,j);           // diagonal

				G(i,j) = G(i,j) + 2.0*EE*P(p,p);           // i > j always
//				G(j,i) = G(j,i) + 2.0*EE*P(p,p);           // repetitive

				G(p,j) = G(p,j) -     EE*P(i,p);           // p > j always
//				G(j,p) = G(j,p) -     EE*P(i,p);           // repetitive

				G(i,p) = G(i,p) -     EE*P(p,j);           // inconclusive
				G(p,i) = G(p,i) -     EE*P(p,j);           // inconclusive
			}else if(i==j){              // ij pair
				G(i,i) = G(i,i) + 4.0*EE*P(p,q);

				G(p,q) = G(p,q) + 2.0*EE*P(i,i);           // p > q always
//				G(q,p) = G(q,p) + 2.0*EE*P(i,i);           // repetitive

				G(p,i) = G(p,i) -     EE*P(i,q);           // p > i always
//				G(i,p) = G(i,p) -     EE*P(i,q);           // repetitive

				G(q,i) = G(q,i) -     EE*P(i,p);           // inconclusive
				G(i,q) = G(i,q) -     EE*P(i,p);           // inconclusive
			}else if((p==i)&&(q==j)){    // pi-qj pair
				G(p,q) = G(p,q) + 3.0*EE*P(p,q);           // diagonal
//				G(q,p) = G(q,p) + 3.0*EE*P(p,q);           // repetitive

				G(p,p) = G(p,p) -     EE*P(q,q);           // diagonal
				G(q,q) = G(q,q) -     EE*P(p,p);           // diagonal
			}else{                       // all distinct
				G(p,q) = G(p,q) + 4.0*EE*P(i,j);           // p > q always
//				G(q,p) = G(q,p) + 4.0*EE*P(i,j);           // repetitive

				G(i,j) = G(i,j) + 4.0*EE*P(p,q);           // i > j always
//				G(j,i) = G(j,i) + 4.0*EE*P(p,q);           // repetitive

				G(p,j) = G(p,j) -     EE*P(i,q);           // p > j always
//				G(j,p) = G(j,p) -     EE*P(i,q);           // repetitive

				G(p,i) = G(p,i) -     EE*P(j,q);           // inconclusive
				G(i,p) = G(i,p) -     EE*P(j,q);           // inconclusive

				G(q,i) = G(q,i) -     EE*P(j,p);           // inconclusive
				G(i,q) = G(i,q) -     EE*P(j,p);           // inconclusive

				G(q,j) = G(q,j) -     EE*P(i,p);           // inconclusive
				G(j,q) = G(j,q) -     EE*P(i,p);           // inconclusive
			}
		}
	}

	// copy value for symmetry
	for(p=0; p < nBasis; p++)
	for(q=0; q < p; q++)
		G(q,p) = G(p,q);

#undef G
#undef P

	// clean memory
	free(shellMap);
	free(la); free(lb); free(lc); free(ld);
	free(ma); free(mb); free(mc); free(md);
	free(na); free(nb); free(nc); free(nd);
	free(iP); free(iQ); free(iI); free(iJ);
	free(EEStore);
	free(EEPrim);
}

// maximum angular moment of gaussian function in basis set it can handel
// support s,p,d,f
//#define MAXL         3
//#define MAXSHELLINT 10000

// support s,p,d,f,g
#define MAXL          4
#define MAXSHELLINT  50625

// storeJK : store J-K matrix elements by taking advantange of the symmetry,
// for example (pq|ij) = (qp|ij) = ...
//
// for Alpha Spin
// [GB]pq = SUM [PT]ij*(pq|ij) - [PA]ij*(pi|qj) 
//
// for Beta Spin
// [GA]pg = SUM [PT]ij*(pq|ij) - [PB]ij*(pi|qj)
//
// where PT is the toal spin density and is defined by
// PT = PA+PB
//
// Note: Do not forget to manually put G(q,p) = G(p,q) in your code
//
// This is for unrestricted hartree-fock calculation
//
// Feb 2013 - Teepanis Chachiyo
//     No longer omit the "repetitive" case and optimized code has been erased
//
// Jan 7, 2013 - Teepanis Chachiyo
//     Code Optimization
//
// May 18,2011 - Teepanis Chachiyo
//     Take away the grouping for sase PQIJ_ALLSAME so that it can be used
//     with shell index too.
//
// Sep 5, 2010 - Teepanis Chachiyo
//     Sepcified pqij cases before going into the loop
//
// July 11, 2011 - Teepanis Chachiyo
//     Initial implementation and testing
//
#define PQIJ_ALLSAME     5
#define PQIJ_2PAIRS      4
#define PQIJ_PQPAIR      3
#define PQIJ_IJPAIR      2
#define PQIJ_PIQJPAIR    1
#define PQIJ_ALLDISTINCT 0
void storeJK(
       double * __restrict GA, const double * __restrict PT,const double * __restrict PA,
       int nBasis, double EE, int p, int q, int i, int j, char type){
#define GA(x,y)  GA[nBasis*x+y]
#define PT(x,y)  PT[nBasis*x+y]
#define PA(x,y)  PA[nBasis*x+y]

	// go thru all possible cases
	switch(type){
	case PQIJ_ALLSAME: //////// all same /////////
		// there is only one possible case
		// (pq|ij) = ....... = ....... = ....... =
		// ....... = ....... = ....... = .......
		//
			GA(p,q) = GA(p,q) + PT(i,j)*EE;
			GA(p,i) = GA(p,i) - PA(q,j)*EE;
		///	GA(p,q) = GA(p,q) + (PT(i,j)-PA(i,j))*EE; // group (1)
	break;

	case PQIJ_2PAIRS:  //////// 2 pairs ///////////
		// there are 2 identical integrals
		// (pq|ij) = ....... = ....... = ....... =
		// (ij|pq) = ....... = ....... = .......
		//
			GA(p,q) = GA(p,q) + PT(i,j)*EE;
			GA(p,i) = GA(p,i) - PA(q,j)*EE;
	
			GA(i,j) = GA(i,j) + PT(p,q)*EE;
			GA(i,p) = GA(i,p) - PA(j,q)*EE;          // repetitive
	break;

	case PQIJ_PQPAIR:  ///////// pq pair ///////////
		// there are 4 identical integrals
		// (pq|ij) = (pq|ji) = ....... = ....... =
		// (ij|pq) = (ji|pq) = ....... = .......
		//
		///	GA(p,q) = GA(p,q) + PT(i,j)*EE;          // grouped into (1)
			GA(p,i) = GA(p,i) - PA(q,j)*EE;
	
		///	GA(p,q) = GA(p,q) + PT(j,i)*EE;          // grouped into (1)
			GA(p,j) = GA(p,j) - PA(q,i)*EE;
	
			GA(i,j) = GA(i,j) + PT(p,q)*EE;
			GA(i,p) = GA(i,p) - PA(j,q)*EE;
	
			GA(p,q) = GA(p,q) + 2.0*PT(i,j)*EE;      // group(1)
	
		 	GA(j,i) = GA(j,i) + PT(p,q)*EE;          // repetitive
			GA(j,p) = GA(j,p) - PA(i,q)*EE;          // repetitive
	break;

	case PQIJ_IJPAIR: ////////// ij pair ////////////
		// there are 4 identical integrals
		// (pq|ij) = ....... = (qp|ij) = ....... =
		// (ij|pq) = ....... = (ij|qp) = .......
		//
			GA(p,q) = GA(p,q) + PT(i,j)*EE;
			GA(p,i) = GA(p,i) - PA(q,j)*EE;
	
			GA(q,p) = GA(q,p) + PT(i,j)*EE;          // repetitive
			GA(q,i) = GA(q,i) - PA(p,j)*EE;
	
		///	GA(i,j) = GA(i,j) + PT(p,q)*EE;          // grouped into (1)
			GA(i,p) = GA(i,p) - PA(j,q)*EE;          // repetitive
	
		///	GA(i,j) = GA(i,j) + PT(q,p)*EE;          // grouped into (1)
			GA(i,q) = GA(i,q) - PA(j,p)*EE;
	
			GA(i,j) = GA(i,j) + 2.0*PT(p,q)*EE;      // group(1)
	break;

	case PQIJ_PIQJPAIR: //////// pi-qj pair //////////
		// there are 4 identical integrals
		// (pq|ij) = (pq|ji) = (qp|ij) = (qp|ji) =
		//  .....  =  .....  =  .....  =  .....
		//
			GA(p,q) = GA(p,q) + PT(i,j)*EE;
			GA(p,i) = GA(p,i) - PA(q,j)*EE;
	
			GA(p,q) = GA(p,q) + PT(j,i)*EE;
			GA(p,j) = GA(p,j) - PA(q,i)*EE;
	
			GA(q,p) = GA(q,p) + PT(i,j)*EE;          // repetitive
			GA(q,i) = GA(q,i) - PA(p,j)*EE;
	
			GA(q,p) = GA(q,p) + PT(j,i)*EE;          // repetitive
			GA(q,j) = GA(q,j) - PA(p,i)*EE;
	break;

	default:   ///////////////// all distinct  /////////////
		// there are 8 identical integrals
		// (pq|ij) = (pq|ji) = (qp|ij) = (qp|ji) =
		// (ij|pq) = (ji|pq) = (ij|qp) = (ji|qp)
		//
		///	GA(p,q) = GA(p,q) + PT(i,j)*EE;          // grouped into (1)
			GA(p,i) = GA(p,i) - PA(q,j)*EE;
	
		///	GA(p,q) = GA(p,q) + PT(j,i)*EE;          // grouped into (1)
			GA(p,j) = GA(p,j) - PA(q,i)*EE;
	
		///	GA(q,p) = GA(q,p) + PT(i,j)*EE;          // grouped into (2)
			GA(q,i) = GA(q,i) - PA(p,j)*EE;
	
		///	GA(q,p) = GA(q,p) + PT(j,i)*EE;          // grouped into (2)
			GA(q,j) = GA(q,j) - PA(p,i)*EE;
	
		///	GA(i,j) = GA(i,j) + PT(p,q)*EE;          // grouped into (3)
			GA(i,p) = GA(i,p) - PA(j,q)*EE;     
	
		///	GA(j,i) = GA(j,i) + PT(p,q)*EE;          // grouped into (4)        
			GA(j,p) = GA(j,p) - PA(i,q)*EE;
	
		///	GA(i,j) = GA(i,j) + PT(q,p)*EE;          // grouped into (3)
			GA(i,q) = GA(i,q) - PA(j,p)*EE;
	
		///	GA(j,i) = GA(j,i) + PT(q,p)*EE;          // grouped into (4)
		    GA(j,q) = GA(j,q) - PA(i,p)*EE;

		//// effectively multiply by 2
		EE = EE + EE;

			GA(p,q) = GA(p,q) + PT(i,j)*EE;          // group (1)
			GA(q,p) = GA(q,p) + PT(i,j)*EE;          // group (2) repetitive
			GA(i,j) = GA(i,j) + PT(p,q)*EE;          // group (3) 
			GA(j,i) = GA(j,i) + PT(p,q)*EE;          // group (4) repetitive
	break;
	}
}
#undef GA
#undef PT
#undef PA

// storeJ : store only the Coulomb part of the fock matrix.
// Note that is code still omits the "repetitive" part.
// 
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
static void storeJ(
       double *GJ, const double *PT,
       int nBasis, double EE, int p, int q, int i, int j, char type){

#define GJ(p,q) GJ[p*nBasis+q]
#define PT(p,q) PT[p*nBasis+q]
	// go thru all possible cases
	switch(type){
	case PQIJ_ALLSAME:
		GJ(p,q) += PT(i,j)*EE;
	break;

	case PQIJ_2PAIRS:
		GJ(p,q) += PT(i,j)*EE;
		GJ(i,j) += PT(p,q)*EE;
	break;

	case PQIJ_PQPAIR:
		GJ(i,j) += PT(p,q)*EE; 
		EE      += EE;
		GJ(p,q) += PT(i,j)*EE;
	break;

	case PQIJ_IJPAIR:
		GJ(p,q) += PT(i,j)*EE;
		EE      += EE;
		GJ(i,j) += PT(p,q)*EE;
	break;

	case PQIJ_PIQJPAIR:
		EE      += EE;
		GJ(p,q) += PT(i,j)*EE;
	break;

	default:
		EE      += EE;
		GJ(p,q) += PT(i,j)*EE;
		GJ(i,j) += PT(p,q)*EE;
	break;
	}

}
#undef GJ
#undef PT

#define GA(p,q) GA[p*nBasis+q]
#define GB(p,q) GB[p*nBasis+q]
#define PA(p,q) PA[p*nBasis+q]
#define PB(p,q) PB[p*nBasis+q]

// storeK : store only the exchange part of the fock matrix.
// Note that is code still omits the "repetitive" part.
// 
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
static void storeK(
       double *GA, double *GB, const double *PA, const double *PB,
       int nBasis, double EE, int p, int q, int i, int j, char type){

	double r;
	// go thru all possible cases
	switch(type){
	case PQIJ_ALLSAME:
		GA(p,i) -= PA(q,j)*EE;
		GB(p,i) -= PB(q,j)*EE;
	break;

	case PQIJ_2PAIRS:
		GA(p,i) -= PA(q,j)*EE;
		GB(p,i) -= PB(q,j)*EE;
	break;

	case PQIJ_PQPAIR:
		GA(p,i) -= PA(q,j)*EE; GA(p,j) -= PA(q,i)*EE; GA(i,p) -= PA(j,q)*EE;
		GB(p,i) -= PB(q,j)*EE; GB(p,j) -= PB(q,i)*EE; GB(i,p) -= PB(j,q)*EE;
	break;

	case PQIJ_IJPAIR:
		GA(p,i) -= PA(q,j)*EE; GA(q,i) -= PA(p,j)*EE; GA(i,q) -= PA(j,p)*EE;
		GB(p,i) -= PB(q,j)*EE; GB(q,i) -= PB(p,j)*EE; GB(i,q) -= PB(j,p)*EE;
	break;

	case PQIJ_PIQJPAIR:
		GA(p,i) -= PA(q,j)*EE; GA(p,j) -= PA(q,i)*EE; GA(q,i) -= PA(p,j)*EE; GA(q,j) -= PA(p,i)*EE;
		GB(p,i) -= PB(q,j)*EE; GB(p,j) -= PB(q,i)*EE; GB(q,i) -= PB(p,j)*EE; GB(q,j) -= PB(p,i)*EE;
	break;

	default:
		r = PA(q,j)*EE; GA(p,i) -= r; GA(i,p) -= r;
		r = PA(q,i)*EE; GA(p,j) -= r; GA(j,p) -= r;
		r = PA(p,j)*EE; GA(q,i) -= r; GA(i,q) -= r;
		r = PA(p,i)*EE; GA(q,j) -= r; GA(j,q) -= r;

		r = PB(q,j)*EE; GB(p,i) -= r; GB(i,p) -= r;
		r = PB(q,i)*EE; GB(p,j) -= r; GB(j,p) -= r;
		r = PB(p,j)*EE; GB(q,i) -= r; GB(i,q) -= r;
		r = PB(p,i)*EE; GB(q,j) -= r; GB(j,q) -= r;
	break;
	}

}
#undef GA
#undef GB
#undef PA
#undef PB


//
// shell_prep : prepare shell pointer and variables to be used with
// function GTO_JK_Matrix_ShellSet and alike
// It returns the number of shell
//
// May 19, 2011 - Teepanis Chachiyo
//    No longer check for unique basis index q
//
// Oct 21, 2010 - Teepanis Chachiyo
//    Function prototype document
//
// Oct 17, 2010 - Teepanis Chachiyo
//    Initial implementation
//
// Oct 21, 2010 - Teepanis Chachiyo
//    Function prototype document
//
// Oct 17, 2010 - Teepanis Chachiyo
//    Initial implementation
//
int shell_prep(
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
		for(q=shellMap[qS]; q < shellMap[qS+1]; q++)
			if(Schwarz[p*nBasis+q] > r) r = Schwarz[p*nBasis+q];
		Schwarz_Shell[pS*nShell+qS] = r;
	} 

	// set returned variables
	*RshellMap      = shellMap;
	*RshellMaxL     = shellMaxL;
	*RSchwarz_Shell = Schwarz_Shell;

	return nShell;
}


// loop structure to go thru all distinct shell sequence
#define SHELL_LOOP_BEGIN                                      \
	for(pS=nShell-1; pS>=0; pS--)                             \
	for(qS=pS;       qS>=0; qS--)                             \
	for(iS=pS;       iS>=0; iS--){if(iS==pS) n=qS; else n=iS; \
	for(jS=n;        jS>=0; jS--){


#define SHELL_LOOP_END            }                           \
	                              }


// loop structure to go thru all basis function sequence in the shell
#define PQIJ_LOOP_BEGIN                       \
for(p=shellMap[pS]; p < shellMap[pS+1]; p++)  \
for(q=shellMap[qS]; q < shellMap[qS+1]; q++)  \
for(i=shellMap[iS]; i < shellMap[iS+1]; i++)  \
for(j=shellMap[jS]; j < shellMap[jS+1]; j++){

#define PQIJ_LOOP_END                       }


// countErrorBank : returns the number of storage need to keep error profile
// for each shell quartet.
//
// Dec 25, 2012 - Teepanis Chachiyo
//      Initial implementation and testing
//
unsigned long long countErrorBank(
	int nShell,            // number of shell
	int *shellMap,         // shell to basis function mapping
	double *Schwarz_Shell, // schwarz at shell level
	double cutoff){        // cutoff value

	unsigned long long nErrorBank=0;// number of shell quartet
	int pS,qS,iS,jS;                // looping variables for shells
	int n;                          // generic index

	// loop all distinct shell sequence
	SHELL_LOOP_BEGIN

		// screen integral using Schwarz inequality at shell level
		if(Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS] < cutoff)
			continue;

		nErrorBank++;

	SHELL_LOOP_END

	return nErrorBank;
}


// countEEBank : return the number integral needed to calculate when
// the schwarz screen is used with "cutmem" cutoff value
//
// May 22, 2011 - Teepanis Chachiyo
//     Initial implementation and testing
//
// Nov 11, 2012 - Teepanis Chachiyo
//     Use only cutoff value
//
unsigned long long countEEBank(
	int nShell,            // number of shell
	int *shellMap,         // shell to basis function mapping
	double *Schwarz_Shell, // schwarz at shell level
	double cutoff){        // cutoff value

	unsigned long long nEEBank=0;// number of integral
	int pS,qS,iS,jS;             // looping variables for shells
	int p,q,i,j;                 // looping variables for basis functions
	int n;                       // generic index

	// loop all distinct shell sequence
	SHELL_LOOP_BEGIN

		// screen integral using Schwarz inequality at shell level
		if(Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS] < cutoff)
			continue;

		// screen integral for all permutation of basis function in shell
		PQIJ_LOOP_BEGIN
			nEEBank++;
		PQIJ_LOOP_END

	SHELL_LOOP_END

	return nEEBank;
}

//
// GTO_JK_Matrix_ShellSet : calculate G matrix for unrestricted hartree-fock
// calculations. The equations are given in (Szabo and Ostlund, 1989)
// 
// for Alpha Spin
// [GB]pq = SUM [PT]ij*(pq|ij) - [PA]ij*(pi|qj) 
//
// for Beta Spin
// [GA]pg = SUM [PT]ij*(pq|ij) - [PB]ij*(pi|qj)
//
// where PT is the toal spin density and is defined by
// PT = PA+PB
//
// Feb 2013 - Teepanis Chachiyo
//    - This function is now obsolete. Use the Quartet version instead
//
// Jan 23, 2013 - Teepanis Chachiyo
//    - The array Sx,Sy,Sz is done by calling function genSetBxyzSxyzF instead.
//
// Jan 7, 2013  - Teepanis Chachiyo
//    - The density screening is done at a shell level with the variable
//      maxPShell instead of going thru basis at a time
// 
// Nov 11, 2012 - Teepanis Chachiyo
//    - Computing tBz[kZ]*F[kX+kY+kZ] first and put them in tSz[kX+kY]
//    - Bug fix when screening at the entire-shell primitive. We do need
//      to compute when deposit is required
// 
// Nov 7, 2012 - Teepanis Chachiyo
//    - Compute contraction coefficient before entering PQIJ loop and stored
//      in cntCoef array
//    - Optimize the kX,kY,kZ loop
//
// Oct 16, 2012 - Teepanis Chachiyo
//    Accessing EEBank using pointer rather than indexing. 
//
// May 28, 2011 - Teepanis Chachiyo
//    Add another schwarz screening before going in genSetBxyzF
//
// May 24, 2011 - Teepanis Chachiyo
//    Serious! bugfix, no longer trap EE=0 during contraction. The idea was
//    to aviod doing symmetrically zero EE. But this cause error in some
//    basis function like sto3g
//
// May 23, 2011 - Teepanis Chachiyo
//    Adjust corse/fine search when determining the cutoff for memory storage
//
// May 21, 2011 - Teepanis Chachiyo
//    Add memory storage option. Because there is a static array allocated
//    UHF need to call with nBasis=0 to reset every thing
//
// May 19, 2011 - Teepanis Chachiyo
//    The distinct index is now done at shell level
//
// Oct 18, 2010 - Teepanis Chachiyo
//    Now support s,p,d,f,h type orbital
//
// Oct 17, 2010 - Teepanis Chachiyo
//    Put shell handler in a separate function called shell_prep
//
// Oct 16, 2010 - Teepanis Chachiyo
//    Pre-screening using Schwarz inequality at shell level
//    bugfix regarding the case  mZ[nEE]==0 && mX[nEE]==0
//
// Sep 5, 2010 - Teepanis Chachiyo
//    Change from memset(EEStore, 0, ...) to for loop EEStore[nEE] instead.
//    This almost double speed of the calculations in some cases. Remember,
//    block memory operation is very very expensive!
//
// July 11, 2010 - Teepanis Chachiyo
//    Initial implementation and testing
//
void GTO_JK_Matrix_ShellSet(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB,                    // return G for spin down
	struct option_t *opt){         // global option

	double *PT;                    // total spin density matrix
	double *PM;                    // maximum density for screening
	double *GJ;                    // coulomb contribution
	double *shellMaxP;             // for screening with density within shell
	int pS,qS,iS,jS;               // looping variables for shells
	int p,q,i,j;                   // looping variables for basis functions
	int pC,qC,iC,jC;               // contracted index
	double EE;                     // coulomb integral
	int nEE;                       // number of electron within shell counter
	double EEStore[MAXSHELLINT];   // electron integral storage
	int        iBx[MAXSHELLINT];   // index for Bx storage
	int        iBy[MAXSHELLINT];   // index for By storage
	int        iBz[MAXSHELLINT];   // index for Bz storage
	char        mX[MAXSHELLINT];   // maximum in x direction
	char        mY[MAXSHELLINT];   // maximum in y direction
	char        mZ[MAXSHELLINT];   // maximum in z direction
	char        mT[MAXSHELLINT];   // type of mX,mY,mZ for loop optimization
	double  tSum;                  // intermediate variables for kX,kY,kZ loop
	int maxCnt;                    // maximum number of contraction
	double *cntCoef;               // contract function coefficients            
	double Bx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double By[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Bz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sy[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double  F[4*(MAXL+1)];         // Boys functions
	double *tBx,*tBy,*tBz;         // pointer to Bx, By, Bz
	double *tSx,*tSy,*tSz;         // pointer to Sx, Sy, Sz
	int kX, kY, kZ;                // summation loop
	char pqijT;                    // storage type for a set of pqij
	int nShell;                    // number of shell
	int *shellMap;                 // mapping between shell and basis function
	int *shellMaxL;                // maximum angular moment for each shell
	double *Schwarz_Shell;         // schwarz screening at shell level
	static int nShellPrim=0;       // number of shell-primitive
	static double *SchwarzSP=NULL; // schwarz at shell-primitive level
	static int *primMap=NULL;      // shell to primitive mapping
	int n;                         // generic index
	double r;                      // generic real number
	double maxP;                   // maximum value of density matrices
	double maxPShell;              // maximum value of density within shell

	static double *EEBank=NULL;    // global 2e storage
	double *EEBankPtr=NULL;        // pointer to current 2e storage
	static int firstVisit=1;       // first time visit flag
	static double cutmem=0.0;      // schwartz cutoff to store in memory
	unsigned long long nEEBank=0;  // the number of items in EEBank
	int deposit=0;                 // need to deposit to 2E bank flag
	unsigned long long maxEE;      // maximum EE for available memory
#define LNBASEI 2.4663034624
#define ERROR_BASE 1.5
#define ERROR_MIN       -120
#define ERROR_MAX       +120
#define ERROR_UNKNOWN   +121
#define ERROR_NOINFO    +122
	static struct Multipole_t *mpole=NULL; // pointer to multipole moment
	const struct Multipole_t *m1,*m2;      // pointer to both doublet multipole
	double K;                              // inverse of distance

	static int *ErrorBank=NULL;            // global error bank
	int *ErrorBankPtr=NULL;                // pointer to current shell 
	static double errorCutoff=0.0;         // schwarz cutoff for error bank
	signed char errorExp;                  // current shell max error
	unsigned char maxEERatio;              // current shell max EE
	int shellProp;
	unsigned long long nErrorBank=0;       // number of shell bank


	// catch reset call
	if(nBasis==0){
		firstVisit=1;
		cutmem=0.0;
		nShellPrim=0;
		errorCutoff=0.0;
		if(EEBank!=NULL){ free(EEBank); EEBank=NULL; }
		if(ErrorBank!=NULL){ free(ErrorBank); ErrorBank=NULL; }
		if(mpole!=NULL){ free(mpole); mpole=NULL; }
		if(SchwarzSP!=NULL){ free(SchwarzSP); SchwarzSP=NULL;}
		if(primMap!=NULL){ free(primMap); primMap=NULL;}
		return;
	}

#define ALLOC(array,item,type)                                         \
array=calloc(item,sizeof(type));                                       \
if(array==NULL){                                                       \
	printf("GTO_JK_Matrix_ShellSet - error cannot allocate memory\n"); \
	exit(-1);                                                          \
}

	// memory allocations
	ALLOC(GJ,  nBasis*nBasis, double);
	ALLOC(PT,  nBasis*nBasis, double);
	ALLOC(PM,  nBasis*nBasis, double);

	// contraction coefficients
	maxCnt = gto[0].nContract;
	for(p=0; p < nBasis; p++)
		if(gto[p].nContract > maxCnt) maxCnt = gto[p].nContract;
	ALLOC(cntCoef, maxCnt*nBasis, double);
	for(p=0; p < nBasis; p++)
		for(pC=0; pC < gto[p].nContract; pC++)
			cntCoef[pC*nBasis+p] = gto[p].coef[pC]*gto[p].norm[pC];

	// compute total density matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		PT[p*nBasis+q] = PA[p*nBasis+q] + PB[p*nBasis+q];
		                                          PM[p*nBasis+q] = fabs(PT[p*nBasis+q]);
		if(PM[p*nBasis+q] < fabs(PA[p*nBasis+q])) PM[p*nBasis+q] = fabs(PA[p*nBasis+q]);
		if(PM[p*nBasis+q] < fabs(PB[p*nBasis+q])) PM[p*nBasis+q] = fabs(PB[p*nBasis+q]);
	}

	// prepare shell data
	nShell = shell_prep(nBasis,gto,Schwarz,&shellMap,&shellMaxL,&Schwarz_Shell);

	// compute maximum density for screening
	ALLOC(shellMaxP, nShell*nShell, double);
	maxP = 0.0;
	for(pS=0; pS < nShell; pS++)
	for(qS=0; qS < nShell; qS++){
		// determine maximum density within this shell
		r = 0.0;
		for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
		for(q=shellMap[qS]; q < shellMap[qS+1]; q++){
			if(fabs(PT[p*nBasis+q]) > r) r = fabs(PT[p*nBasis+q]);
			if(fabs(PA[p*nBasis+q]) > r) r = fabs(PA[p*nBasis+q]);
			if(fabs(PB[p*nBasis+q]) > r) r = fabs(PB[p*nBasis+q]);
		}
		shellMaxP[pS*nShell+qS] = r;
		if(r > maxP) maxP = r;
	}

	// first time visit
	if(firstVisit){

		// generate multipole moment
		mpole = genMultipole(nBasis, gto);

		// allocate memory for shell-primitive mapping
		ALLOC(primMap, nShell, int);

		// compute primitive-shell mapping
		nShellPrim=0;
		for(pS=0; pS<nShell; pS++){
			primMap[pS] = nShellPrim;
			nShellPrim += gto[shellMap[pS]].nContract;
		}
		
		// allocate schawarz at shell-primitive level
		ALLOC(SchwarzSP, nShellPrim*nShellPrim, double);

		// compute schwarz at shell-primitive level
		for(pS=0; pS<nShell; pS++)
		for(qS=0; qS<nShell; qS++){

			for(pC=0; pC<gto[shellMap[pS]].nContract; pC++)
			for(qC=0; qC<gto[shellMap[qS]].nContract; qC++){

				// search for maximum within this permutation
				r=0.0;
				for(p=shellMap[pS]; p<shellMap[pS+1]; p++)
				for(q=shellMap[qS]; q<shellMap[qS+1]; q++){

					EE = eri(gto[p].x0, gto[p].y0, gto[p].z0, 1.0, 
					         gto[p].l, gto[p].m, gto[p].n, gto[p].exp[pC],
					         gto[q].x0, gto[q].y0, gto[q].z0, 1.0,
					         gto[q].l, gto[q].m, gto[q].n, gto[q].exp[qC],
					         gto[p].x0, gto[p].y0, gto[p].z0, 1.0,
					         gto[p].l, gto[p].m, gto[p].n, gto[p].exp[pC],
					         gto[q].x0, gto[q].y0, gto[q].z0, 1.0,
					         gto[q].l, gto[q].m, gto[q].n, gto[q].exp[qC]);

					EE = EE*gto[p].coef[pC]*gto[p].norm[pC]
					       *gto[q].coef[qC]*gto[q].norm[qC]
					       *gto[p].coef[pC]*gto[p].norm[pC]
					       *gto[q].coef[qC]*gto[q].norm[qC];

					EE = fabs(EE); if(EE>r) r = EE;
				}

				// store value
				SchwarzSP[(primMap[pS]+pC)*nShellPrim+(primMap[qS]+qC)] = sqrt(r);
			}
		}

		// increase cutmem until it is enough for available memory
		cutmem = fixedCutoff;
		maxEE  = (unsigned long long)(opt->maxMem)*1000000/sizeof(double);

		// set corse search flag
		r = 1;
		do{
CutMemSearch:
			nEEBank = 0;
			// loop all distinct shell sequence
			SHELL_LOOP_BEGIN
		
				// screen integral using Schwarz inequality at shell level
				if(Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS] < cutmem)
					continue;
		
				// accumulate the number of shell permutation
				nEE=0;
				PQIJ_LOOP_BEGIN
					nEE++;
				PQIJ_LOOP_END

				nEEBank+=nEE;
				if(nEEBank > maxEE){
					if(r==1) // corse search
						cutmem = cutmem*2.00;
					else     // fine search
						cutmem = cutmem*1.05;
					goto CutMemSearch;
				}

			SHELL_LOOP_END

			// done corse search and begin fine search
			if(r==1 && nEEBank < maxEE){ r = 0; cutmem/=2.0; goto CutMemSearch; }

		}while(nEEBank > maxEE);

		// cutmem should not be smaller than fixedCutoff
		if(cutmem < fixedCutoff) cutmem = fixedCutoff;

		// compute number of integral needed
		nEEBank=countEEBank(nShell, shellMap, Schwarz_Shell, cutmem);

		// allocate memory
		EEBank=calloc(nEEBank,sizeof(double));
		if(EEBank==NULL){
			printf("GTO_JK_Matrix_ShellSet "
			       "- cannot allocate 2e memory, check option -MAXMEM\n");
			exit(-1);
		}

		// allocate memory for error properties
		errorCutoff = opt->SCFCutoff;
		nErrorBank  = countErrorBank(nShell, shellMap, Schwarz_Shell, errorCutoff);

		ErrorBank=calloc(nErrorBank, sizeof(int));
		if(ErrorBank==NULL){
			printf("GTO_JK_Matrix_ShellSet "
			       "- cannot allocate error profile storage\n");
			exit(-1);
		}
		// set value to unkown
		ErrorBankPtr = ErrorBank;
		for(maxEE=0; maxEE < nErrorBank; maxEE++){
			*ErrorBankPtr = ERROR_UNKNOWN;
		     ErrorBankPtr++;
		}

	}

	EEBankPtr = EEBank;
	ErrorBankPtr = ErrorBank;
	// loop all distinct shell sequence
	SHELL_LOOP_BEGIN

		r = Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS];

		if(r < errorCutoff) continue;

		// handle error profile
		shellProp  = *ErrorBankPtr;
		errorExp   = (signed char)(shellProp & 255);
		maxEERatio = (unsigned char)(shellProp >> 8);
		ErrorBankPtr++;

		// screen integral using Schwarz inequality at shell level
		if(r < fixedCutoff) continue;

		// if larger than cutmem we load from memory
		if(r >= cutmem){

			// first time, however, need to evaluate and save them
			if(firstVisit){ deposit=1; maxPShell = maxP; goto computeUniqueShell; }

			// screen with maximum electron density wieghted
			if(r*maxP < fixedCutoff){
				EEBankPtr += (shellMap[pS+1]-shellMap[pS])*
				             (shellMap[qS+1]-shellMap[qS])*
				             (shellMap[iS+1]-shellMap[iS])*
				             (shellMap[jS+1]-shellMap[jS]);
				continue;
			}

			// screen another layer with electron density at shell level weighted
			maxPShell = shellMaxP[pS*nShell+qS] + shellMaxP[iS*nShell+jS] +
			            shellMaxP[pS*nShell+iS] + shellMaxP[pS*nShell+jS] +
			            shellMaxP[qS*nShell+iS] + shellMaxP[qS*nShell+jS];

			if(r*maxPShell < fixedCutoff){
				EEBankPtr += (shellMap[pS+1]-shellMap[pS])*
				             (shellMap[qS+1]-shellMap[qS])*
				             (shellMap[iS+1]-shellMap[iS])*
				             (shellMap[jS+1]-shellMap[jS]);
				continue;
			}

			// screen with individual density
			PQIJ_LOOP_BEGIN
				r = PM[p*nBasis+q] + PM[i*nBasis+j] + 
				    PM[p*nBasis+i] + PM[p*nBasis+j] +
				    PM[q*nBasis+i] + PM[q*nBasis+j];
				if(r*Schwarz[p*nBasis+q]*Schwarz[i*nBasis+j] > fixedCutoff)
					goto loadUniqueEE;
			PQIJ_LOOP_END
			EEBankPtr += (shellMap[pS+1]-shellMap[pS])*
			             (shellMap[qS+1]-shellMap[qS])*
			             (shellMap[iS+1]-shellMap[iS])*
			             (shellMap[jS+1]-shellMap[jS]);
			continue;

loadUniqueEE:
			// determine pS qS iS jS type for storage
			if((pS==qS)&&(iS==jS)&&(pS==iS)){ pqijT = PQIJ_ALLSAME;
			}else if((pS==qS)&&(iS==jS)){     pqijT = PQIJ_2PAIRS;
			}else if(pS==qS){                 pqijT = PQIJ_PQPAIR;
			}else if(iS==jS){                 pqijT = PQIJ_IJPAIR;
			}else if((pS==iS)&&(qS==jS)){     pqijT = PQIJ_PIQJPAIR;
			}else{                            pqijT = PQIJ_ALLDISTINCT;
			}

			// assemble into G matrix
			PQIJ_LOOP_BEGIN
				EE = *EEBankPtr;
				if(fabs(EE*maxPShell)>fixedCutoff){
					storeJ(GJ,PT,nBasis,EE,p,q,i,j,pqijT);
					storeK(GA,GB,PA,PB,nBasis,EE,p,q,i,j,pqijT);
				}
				EEBankPtr++;
			PQIJ_LOOP_END
			continue;
		
		}

		//
		// here we have integral too large to ignore, too small to store
		//

		// reduce the screening to known maximum EE
		if(errorExp!=ERROR_UNKNOWN) r /= maxEERatio;

		// screen with maximum electron density wieghted
		if(r*maxP < fixedCutoff) continue;

		// do not deposit in this case, even if it's the first time
		if(firstVisit) deposit=0;

		// screen another layer with electron density weighted at shell level
		maxPShell = shellMaxP[pS*nShell+qS] + shellMaxP[iS*nShell+jS] +
		            shellMaxP[pS*nShell+iS] + shellMaxP[pS*nShell+jS] +
		            shellMaxP[qS*nShell+iS] + shellMaxP[qS*nShell+jS];

		if(r*maxPShell < fixedCutoff) continue;

		// screen with individual density
		PQIJ_LOOP_BEGIN
			K = PM[p*nBasis+q] + PM[i*nBasis+j] + 
			    PM[p*nBasis+i] + PM[p*nBasis+j] +
			    PM[q*nBasis+i] + PM[q*nBasis+j];
			EE = Schwarz[p*nBasis+q]*Schwarz[i*nBasis+j];
			if(EE > r) EE = r;
			if(K*EE > fixedCutoff)	goto computeUniqueShell;
		PQIJ_LOOP_END
		continue;

computeUniqueShell:
		// perform multipole approximation if possible
		if(errorExp==ERROR_UNKNOWN) goto computeExactUniqueShell;
		EE = pow(2.0,errorExp);
		PQIJ_LOOP_BEGIN
			r = PM[p*nBasis+q] + PM[i*nBasis+j] + 
			    PM[p*nBasis+i] + PM[p*nBasis+j] +
			    PM[q*nBasis+i] + PM[q*nBasis+j];
			if(r*EE > fixedCutoff) goto computeExactUniqueShell;
		PQIJ_LOOP_END

		// determine pS qS iS jS type for storage
		if((pS==qS)&&(iS==jS)&&(pS==iS)){ pqijT = PQIJ_ALLSAME;
		}else if((pS==qS)&&(iS==jS)){     pqijT = PQIJ_2PAIRS;
		}else if(pS==qS){                 pqijT = PQIJ_PQPAIR;
		}else if(iS==jS){                 pqijT = PQIJ_IJPAIR;
		}else if((pS==iS)&&(qS==jS)){     pqijT = PQIJ_PIQJPAIR;
		}else{                            pqijT = PQIJ_ALLDISTINCT;
		}

		// compute and store
		for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
		for(q=shellMap[qS]; q < shellMap[qS+1]; q++){
			m1 = &mpole[p*nBasis+q];
			for(i=shellMap[iS]; i < shellMap[iS+1]; i++){
				m2 = &mpole[i*nBasis+shellMap[jS]];
				for(j=shellMap[jS]; j < shellMap[jS+1]; j++){

					K   = (m1->x - m2->x)*(m1->x - m2->x) + 
					      (m1->y - m2->y)*(m1->y - m2->y) +
					      (m1->z - m2->z)*(m1->z - m2->z);

					if(K < MULTIPOLE_RADII2_CUTOFF)
						EE = 0.0; 
					else
						EE = m1->q * m2->q * sqrt(1/K);

					if(fabs(EE*maxPShell) > fixedCutoff){
						storeJ(GJ,PT,nBasis,EE,p,q,i,j,pqijT);
						storeK(GA,GB,PA,PB,nBasis,EE,p,q,i,j,pqijT);
					}

					m2++;
				}
			}
		}
		continue;

computeExactUniqueShell:
		// index preparation
		nEE = 0;
		PQIJ_LOOP_BEGIN

			// precompute index for efficientcy
			kZ =    (shellMaxL[jS]+1);
			kY = kZ*(shellMaxL[iS]+1);
			kX = kY*(shellMaxL[qS]+1);
			iBx[nEE] = 4*(MAXL+1)*(gto[p].l*kX + gto[q].l*kY + gto[i].l*kZ + gto[j].l);
			iBy[nEE] = 4*(MAXL+1)*(gto[p].m*kX + gto[q].m*kY + gto[i].m*kZ + gto[j].m);
			iBz[nEE] = 4*(MAXL+1)*(gto[p].n*kX + gto[q].n*kY + gto[i].n*kZ + gto[j].n);
			mX[nEE]  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
			mY[nEE]  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
			mZ[nEE]  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;
#define MXZERO 0
#define MYZERO 1
#define MZZERO 2
#define MXMAX  3
#define MYMAX  4
#define MZMAX  5
			if(mX[nEE]==0)                                  mT[nEE] = MXZERO;
			else if(mY[nEE]==0)                             mT[nEE] = MYZERO;
			else if(mZ[nEE]==0)                             mT[nEE] = MZZERO;
			else if(mX[nEE] > mY[nEE] && mX[nEE] > mZ[nEE]) mT[nEE] = MXMAX;
			else if(mY[nEE] > mZ[nEE])                      mT[nEE] = MYMAX;
			else                                            mT[nEE] = MZMAX;

			// reset EEStore to zero
			EEStore[nEE] = 0.0;

			// increase number of EE index
			nEE++;

		PQIJ_LOOP_END

		// check that the number of integral does not exceed maximum 
		if(nEE > MAXSHELLINT){
			printf("GTO_JK_Matrix_ShellSet - error too many integrals\n");
			exit(-1);
		}

		///////////////////////////////////////////////
		// generate (ab|cd)
		//
		// 1) loop all contracted
		// 2) generate Bx,By,Bz for each contracted
		// 3) summation to get (ab|cd)
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

			// schwarz screening for the entire shell-primitive
			if(SchwarzSP[(primMap[pS]+pC)*nShellPrim+(primMap[qS]+qC)] *
			   SchwarzSP[(primMap[iS]+iC)*nShellPrim+(primMap[jS]+jC)] * maxP <
			   fixedCutoff)
				if(!deposit) continue;

			// generate Bx, By , Bz, and F
			//r=genSetBxyzF(
			//        gto[p].x0,gto[p].y0,gto[p].z0,shellMaxL[pS],gto[p].exp[pC],
			//        gto[q].x0,gto[q].y0,gto[q].z0,shellMaxL[qS],gto[q].exp[qC],
			//        gto[i].x0,gto[i].y0,gto[i].z0,shellMaxL[iS],gto[i].exp[iC],
			//        gto[j].x0,gto[j].y0,gto[j].z0,shellMaxL[jS],gto[j].exp[jC],
			//        Bx,By,Bz,F,MAXL);

			r=genSetBxyzSxyzF(
			        gto[p].x0,gto[p].y0,gto[p].z0,shellMaxL[pS],gto[p].exp[pC],
			        gto[q].x0,gto[q].y0,gto[q].z0,shellMaxL[qS],gto[q].exp[qC],
			        gto[i].x0,gto[i].y0,gto[i].z0,shellMaxL[iS],gto[i].exp[iC],
			        gto[j].x0,gto[j].y0,gto[j].z0,shellMaxL[jS],gto[j].exp[jC],
			        Bx,By,Bz,Sx,Sy,Sz,F,MAXL);

			if(r < PRIMITIVE_CUTOFF) continue;
			
			// loop all basis within shell
			nEE = 0;
			PQIJ_LOOP_BEGIN

				// compute two-electron integral using summation
				tBx = Bx + iBx[nEE];
				tBy = By + iBy[nEE];
				tBz = Bz + iBz[nEE];

				tSx = Sx + iBx[nEE];
				tSy = Sy + iBy[nEE];
				tSz = Sz + iBz[nEE];

				EE=0.0;
				switch(mT[nEE]){

				case MXZERO:
					for(kY=mY[nEE];kY>=0;kY--) EE += tBy[kY]*tSz[kY];
				break;

				case MYZERO:
					for(kZ=mZ[nEE];kZ>=0;kZ--) EE += tBz[kZ]*tSx[kZ];
				break;

				case MZZERO:
					for(kX=mX[nEE];kX>=0;kX--) EE += tBx[kX]*tSy[kX];
				break;

				case MXMAX: 
					for(kY=mY[nEE];kY>=0;kY--){
						for(tSum=0.0,kZ=mZ[nEE];kZ>=0;kZ--) tSum+=tBz[kZ]*tSx[kY+kZ];
						EE += tBy[kY]*tSum;
					}
				break;

				case MYMAX:
					for(kZ=mZ[nEE];kZ>=0;kZ--){
						for(tSum=0.0,kX=mX[nEE];kX>=0;kX--) tSum+=tBx[kX]*tSy[kX+kZ];
						EE += tBz[kZ]*tSum;
					}
				break;

				case MZMAX:
					for(kX=mX[nEE];kX>=0;kX--){
						for(tSum=0.0,kY=mY[nEE];kY>=0;kY--) tSum+=tBy[kY]*tSz[kX+kY];
						EE += tBx[kX]*tSum;
					}
				break;
				}

				// multiply contracted coefficient
				//EEStore[nEE] += EE * r
				//                   * gto[p].coef[pC]*gto[p].norm[pC]
				//                   * gto[q].coef[qC]*gto[q].norm[qC]
				//                   * gto[i].coef[iC]*gto[i].norm[iC]
				//                   * gto[j].coef[jC]*gto[j].norm[jC];
				EEStore[nEE] += EE * r
				                   * cntCoef[pC*nBasis+p]
				                   * cntCoef[qC*nBasis+q]
				                   * cntCoef[iC*nBasis+i]
				                   * cntCoef[jC*nBasis+j];

				nEE++;
			
			PQIJ_LOOP_END

			// reset shell index to proper value in this shell
			p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		}

		// determine pS qS iS jS type for storage
		if((pS==qS)&&(iS==jS)&&(pS==iS)){ pqijT = PQIJ_ALLSAME;
		}else if((pS==qS)&&(iS==jS)){     pqijT = PQIJ_2PAIRS;
		}else if(pS==qS){                 pqijT = PQIJ_PQPAIR;
		}else if(iS==jS){                 pqijT = PQIJ_IJPAIR;
		}else if((pS==iS)&&(qS==jS)){     pqijT = PQIJ_PIQJPAIR;
		}else{                            pqijT = PQIJ_ALLDISTINCT;
		}

		// assemble into G matrix
		nEE = 0;
		PQIJ_LOOP_BEGIN
			EE = EEStore[nEE];
			if(fabs(EE*maxPShell) > fixedCutoff){
				storeJ(GJ,PT,nBasis,EE,p,q,i,j,pqijT);
				storeK(GA,GB,PA,PB,nBasis,EE,p,q,i,j,pqijT);
			}
			nEE++;
		PQIJ_LOOP_END

		// compute error profile if currently unknown
		if(errorExp==ERROR_UNKNOWN){

			// compute and compare
			nEE=0;
			r=0.0;
			tSum=0.0;
			for(p=shellMap[pS]; p < shellMap[pS+1]; p++)
			for(q=shellMap[qS]; q < shellMap[qS+1]; q++){
				m1 = &mpole[p*nBasis+q];
				for(i=shellMap[iS]; i < shellMap[iS+1]; i++){
					m2 = &mpole[i*nBasis+shellMap[jS]];
					for(j=shellMap[jS]; j < shellMap[jS+1]; j++){

						K   = (m1->x - m2->x)*(m1->x - m2->x) + 
						      (m1->y - m2->y)*(m1->y - m2->y) +
						      (m1->z - m2->z)*(m1->z - m2->z);

						if(K < MULTIPOLE_RADII2_CUTOFF)
							EE = 0.0; 
						else
							EE = m1->q * m2->q * sqrt(1/K);

						EE -= EEStore[nEE];
						if(fabs(EE) > r) r = fabs(EE);
						if(fabs(EEStore[nEE]) > tSum) tSum = fabs(EEStore[nEE]);
						nEE++;
						m2++;
					}
				}
			}

			// store value
			if(r == 0.0) r = ERROR_MIN;
			else         r = log(r)*LNBASEI;
			if(r < ERROR_MIN) r = ERROR_MIN;
			if(r > ERROR_MAX) r = ERROR_MAX;

			errorExp = (signed char)r;
			r = Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS];
			if(tSum==0.0) r = 100.0; else r = r/tSum;
			if(r < 1.0)   r = 1.0;
			if(r > 100.0) r = 100.0;
			maxEERatio = (unsigned char)r;
			shellProp  = ((int)maxEERatio << 8) |  (int)(errorExp & 255);
			 ErrorBankPtr--;
			*ErrorBankPtr = shellProp;
			 ErrorBankPtr++;

		}

		// handle first visit preparations
		if(firstVisit && deposit){
			nEE=0;
			PQIJ_LOOP_BEGIN
				*EEBankPtr = EEStore[nEE];
				 EEBankPtr++;
				nEE++;
			PQIJ_LOOP_END
		}

	SHELL_LOOP_END

	// symmetrize G matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < p; q++){
		GA[q*nBasis+p] = GA[p*nBasis+q];
		GB[q*nBasis+p] = GB[p*nBasis+q];
		GJ[q*nBasis+p] = GJ[p*nBasis+q];
	}

	// add coulomb contribution and return value 
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		GA[p*nBasis+q] += GJ[p*nBasis+q];
		GB[p*nBasis+q] += GJ[p*nBasis+q];
	}

	// deselect first visit flag
	if(firstVisit) firstVisit=0;

	// clean memory
	free(cntCoef);
	free(PT);
	free(PM);
	free(GJ);
	free(shellMaxP);
	free(shellMap);
	free(shellMaxL);
	free(Schwarz_Shell);
}
#undef ALLOC


// shellCpy : copies data between shell quartet
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void shellCpy(struct GTOShell_t *dest, struct GTOShell_t *src){
	int i,j;
	dest->x         = src->x;
	dest->y         = src->y;
	dest->z         = src->z;
	dest->nBasis    = src->nBasis;
	dest->nContract = src->nContract;
	dest->maxL      = src->maxL;
	dest->min       = src->min;
	dest->max       = src->max;
	dest->type      = src->type;

	for(i=0; i < MAXBASIS; i++){
		dest->l[i]    = src->l[i];
		dest->m[i]    = src->m[i];
		dest->n[i]    = src->n[i];
		dest->exps[i] = src->exps[i];
		for(j=0; j < MAXBASIS; j++)
			dest->coef[i][j] = src->coef[i][j];
	}
}


// sortGTOShell: sort the GTOShell according their shell types
// in order to enhance "instruction-cache" localization
// when calling the computeQuartet and storeQuartet
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void sortGTOShell(int nShell, struct GTOShell_t *shell){
	struct GTOShell_t buffer;
	int P,Q;

	for(P=0; P < nShell; P++)
	for(Q=0; Q < (nShell-1); Q++){
		if(shell[Q].type > shell[Q+1].type){
			shellCpy(&buffer  , shell+Q);
			shellCpy(shell+Q  , shell+Q+1);
			shellCpy(shell+Q+1, &buffer);
		}
	}
}


// printGTOShell : print out information of a particular shell to screen.
// This is used mainly for debuging purpose.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void printGTOShell(int nShell, struct GTOShell_t *shell){
	int i;

	printf("There are %d shells\n", nShell);
	for(i=0; i < nShell; i++){
		printf("%3d Center(%6.1lf,%6.1lf,%6.1lf) Basis(%d-%d)\n", i,
		        shell[i].x, shell[i].y, shell[i].z, shell[i].min, shell[i].max);
		printf("    nContract(%d) maxL(%d)\n",shell[i].nContract, shell[i].maxL);
	}
	printf("\n");
}


// genGTOShell : generate the shell structure from basis set information.
// 
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
// March 8, 2013 - Teepanis Chachiyo
//    Check if the number of contraction does not exceed MAXCONTRACT
//
struct GTOShell_t * genGTOShell(
	int nBasis,                         // number of basis function
	const struct GTOBasis_t *gto,       // basis function data structure
	int *RnShell){                      // returned number of shell created

	int nShell;               // number of shell
	int p,q,i;                // generic counter
	struct GTOShell_t *shell; // pointer to shell data structure


	// check the number of contraction
	for(p=0; p < nBasis; p++)
		if(gto[p].nContract > MAXCONTRACT){
			printf("genGTOShell: too many contractions, increase MAXCONTRACT\n");
			exit(-1);
		}

	// allocate memory
	shell = calloc(nBasis,sizeof(struct GTOShell_t));
	if(shell==NULL){
		printf("genGTOShell: cannot allocate memory\n");
		exit(-1);
	}

	// filling shell data
	nShell=0;
	for(p=0; p < nBasis; p++){

		// copy values
		shell[nShell].min       = p;
		shell[nShell].x         = gto[p].x0;
		shell[nShell].y         = gto[p].y0;
		shell[nShell].z         = gto[p].z0;
		shell[nShell].nContract = gto[p].nContract;
		for(i=0; i < gto[p].nContract; i++)
			shell[nShell].exps[i] = gto[p].exp[i];

		// expand basis to cover this shell
		for(q=p; q < nBasis; q++)

			if( isSameShell(nBasis,p,q,gto) ){

				// determine maxL
				if( shell[nShell].maxL < (gto[q].l+gto[q].m+gto[q].n))
					shell[nShell].maxL = (gto[q].l+gto[q].m+gto[q].n);

				// copy values
				shell[nShell].max = q;
				shell[nShell].l[ shell[nShell].nBasis  ] = gto[q].l;
				shell[nShell].m[ shell[nShell].nBasis  ] = gto[q].m;
				shell[nShell].n[ shell[nShell].nBasis  ] = gto[q].n;
				for(i=0; i < gto[q].nContract; i++)
					shell[nShell].coef[i][ shell[nShell].nBasis ] =
					gto[q].norm[i] * gto[q].coef[i];

				shell[nShell].nBasis++;
			}else				
				break;

		p += (shell[nShell].nBasis-1);
		nShell++;
	}
	shell=realloc(shell,nShell*sizeof(struct GTOShell_t));

	// determine shell type
	for(p=0; p < nShell; p++)
		switch(shell[p].nBasis){
		case  1: shell[p].type = TYPES; break;
		case  3: shell[p].type = TYPEP; break;
		case  4: shell[p].type = TYPEL; break;
		case  6: shell[p].type = TYPED; break;
		case 10: shell[p].type = TYPEF; break;
		default: shell[p].type = TYPEX; break;
		}

	*RnShell = nShell;
	return shell;
}


// getErrExpQuartetMPole : compute the log(x) base 1.5 of maximum error within 
// this shell quartet when using multipole approximation compared to the
// exact value of the integrals.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
signed char getErrExpQuartetMPole(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	const struct Multipole_t *mp,
	int nBasis,
	double *EEStore){

	int p,q,i,j;                            // basis function index
	const struct Multipole_t *m1,*m2;       // pointer to multipole structure
	const struct Multipole_t *m1ptr,*m2ptr; // starting points
	double K;                               // inverse of distance
	double maxErr=0.0;

	m1ptr = mp + P->min * nBasis + Q->min;
	m2ptr = mp + I->min * nBasis + J->min;
	
	// loop thru all basis quartet in the shell
	for(p=0; p < P->nBasis; p++)
	for(q=0; q < Q->nBasis; q++)
	for(i=0; i < I->nBasis; i++)
	for(j=0; j < J->nBasis; j++){

		m1 = m1ptr + p*nBasis + q;
		m2 = m2ptr + i*nBasis + j;

		K   = (m1->x - m2->x)*(m1->x - m2->x) + 
		      (m1->y - m2->y)*(m1->y - m2->y) +
		      (m1->z - m2->z)*(m1->z - m2->z);

		// compute multipole value and its error
		if(K < MULTIPOLE_RADII2_CUTOFF)
			K = *EEStore - 0.0;
		else
			K = *EEStore - m1->q * m2->q * sqrt(1/K);

		// get maximum
		if(fabs(K) > maxErr) maxErr = fabs(K);

		EEStore++;
	}

	// get log(x) base 1.5 of the maximum
	if(maxErr == 0.0) maxErr = ERROR_MIN; else maxErr = log(maxErr)*LNBASEI;
	if(maxErr < ERROR_MIN) maxErr = ERROR_MIN;
	if(maxErr > ERROR_MAX) maxErr = ERROR_MAX;

	return (signed char)maxErr;
}


// computeQuartetMPole : compute the shell quartet integral using 
// multipole expansion
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void computeQuartetMPole(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	const struct Multipole_t *mp,
	int nBasis,
	double *EEStore){

	int p,q,i,j;                            // basis function index
	const struct Multipole_t *m1,*m2;       // pointer to multipole structure
	const struct Multipole_t *m1ptr,*m2ptr; // starting points
	double K;                               // inverse of distance

	m1ptr = mp + P->min * nBasis + Q->min;
	m2ptr = mp + I->min * nBasis + J->min;

	// loop thru all integral in the quartet
	for(p=0; p < P->nBasis; p++)
	for(q=0; q < Q->nBasis; q++)
	for(i=0; i < I->nBasis; i++)
	for(j=0; j < J->nBasis; j++){

		m1 = m1ptr + p*nBasis + q;
		m2 = m2ptr + i*nBasis + j;

		K   = (m1->x - m2->x)*(m1->x - m2->x) + 
		      (m1->y - m2->y)*(m1->y - m2->y) +
		      (m1->z - m2->z)*(m1->z - m2->z);

		// compute only the monopole-monopole term
		if(K < MULTIPOLE_RADII2_CUTOFF)
			*EEStore = 0.0; 
		else
			*EEStore = m1->q * m2->q * sqrt(1/K);

		EEStore++;
	}
}


//
// the global array needed for computeQuartetEE
//
unsigned int iBx[MAXSHELLINT];   // index for Bx storage
unsigned int iBy[MAXSHELLINT];   // index for By storage
unsigned int iBz[MAXSHELLINT];   // index for Bz storage
unsigned int  mX[MAXSHELLINT];   // maximum in x direction
unsigned int  mY[MAXSHELLINT];   // maximum in y direction
unsigned int  mZ[MAXSHELLINT];   // maximum in z direction
unsigned int  mT[MAXSHELLINT];   // type of mX,mY,mZ for loop optimization

double Bx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double By[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Bz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Sx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Sy[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double Sz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
double  F[4*(MAXL+1)];

// computeQuartetEE : compute the integral exactly using THO method.
// For some cases, the "unrolled" version has been created and will
// be much faster.
// 
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void computeQuartetEE(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	double *EEStore){

	double * __restrict tBx, * __restrict tBy, * __restrict tBz;  // pointer to Bx, By, Bz
	double * __restrict tSx, * __restrict tSy, * __restrict tSz;  // pointer to Sx, Sy, Sz
	int kX, kY, kZ;                // summation loop
	double r;                      // generic double precision
	double rpqi;                   // prefactor*p*q*i coefficients
	int pC,qC,iC,jC;               // contraction index
	int p,q,i,j;                   // basis function index
	int nEE;                       // integral index
	register double EE;            // 2e integral
	register double tSum;          // intermediate variables for kX,kY,kZ loop
/* -- Commented out by Teepanis Chachiyo on 16/9/2017. Error detected when computing
   -- Argon atom using ccTZP basis set
#define SS 0
#define LL 1
#define DD 3
	// call pre "unrolled" case if available
	switch(256*P->type + 64*Q->type + 8*I->type + J->type){

	case 256*SS+64*SS+8*SS+SS: computeSSSSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*LL+LL: computeLLLLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*DD+DD: computeDDDDQuartetEE(P,Q,I,J,EEStore); return; break;

	// SL pair
	case 256*SS+64*LL+8*LL+LL: computeSLLLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*LL+LL: computeLSLLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*SS+LL: computeLLSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*LL+SS: computeLLLSQuartetEE(P,Q,I,J,EEStore); return; break;

	case 256*LL+64*SS+8*SS+SS: computeLSSSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*SS+SS: computeSLSSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*SS+8*LL+SS: computeSSLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*SS+8*SS+LL: computeSSSLQuartetEE(P,Q,I,J,EEStore); return; break;

	case 256*SS+64*SS+8*LL+LL: computeSSLLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*SS+LL: computeSLSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*LL+SS: computeSLLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*SS+LL: computeLSSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*LL+SS: computeLSLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*SS+SS: computeLLSSQuartetEE(P,Q,I,J,EEStore); return; break;

	// DS pair
	case 256*DD+64*SS+8*SS+SS: computeDSSSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*DD+8*SS+SS: computeSDSSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*SS+8*DD+SS: computeSSDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*SS+8*SS+DD: computeSSSDQuartetEE(P,Q,I,J,EEStore); return; break;

	case 256*SS+64*DD+8*DD+SS: computeSDDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*DD+8*SS+DD: computeSDSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*SS+8*DD+DD: computeSSDDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*SS+DD: computeDSSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*DD+SS: computeDSDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*SS+SS: computeDDSSQuartetEE(P,Q,I,J,EEStore); return; break;

	case 256*SS+64*DD+8*DD+DD: computeSDDDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*DD+DD: computeDSDDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*SS+DD: computeDDSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*DD+SS: computeDDDSQuartetEE(P,Q,I,J,EEStore); return; break;

	// DL pair
	case 256*DD+64*LL+8*LL+LL: computeDLLLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*DD+8*LL+LL: computeLDLLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*DD+LL: computeLLDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*LL+DD: computeLLLDQuartetEE(P,Q,I,J,EEStore); return; break;

	case 256*LL+64*DD+8*DD+DD: computeLDDDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*DD+DD: computeDLDDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*LL+DD: computeDDLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*DD+LL: computeDDDLQuartetEE(P,Q,I,J,EEStore); return; break;

	case 256*LL+64*DD+8*DD+LL: computeLDDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*DD+8*LL+DD: computeLDLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*DD+DD: computeLLDDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*LL+DD: computeDLLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*DD+LL: computeDLDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*LL+LL: computeDDLLQuartetEE(P,Q,I,J,EEStore); return; break;

	// 2S and DL pair
	case 256*SS+64*SS+8*DD+LL: computeSSDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*SS+8*LL+DD: computeSSLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*DD+8*SS+LL: computeSDSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*SS+DD: computeSLSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*DD+8*LL+SS: computeSDLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*DD+SS: computeSLDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*SS+LL: computeDSSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*SS+DD: computeLSSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*LL+SS: computeDSLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*DD+SS: computeLSDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*SS+SS: computeDLSSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*DD+8*SS+SS: computeLDSSQuartetEE(P,Q,I,J,EEStore); return; break;

	// 2L and SD pair
	case 256*LL+64*LL+8*SS+DD: computeLLSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*LL+8*DD+SS: computeLLDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*LL+DD: computeLSLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*DD+8*LL+SS: computeLDLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*DD+LL: computeLSDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*DD+8*SS+LL: computeLDSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*LL+DD: computeSLLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*LL+SS: computeDLLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*DD+LL: computeSLDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*SS+LL: computeDLSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*DD+8*LL+LL: computeSDLLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*LL+LL: computeDSLLQuartetEE(P,Q,I,J,EEStore); return; break;

	// 2D and SL pair
	case 256*DD+64*DD+8*SS+LL: computeDDSLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*DD+8*LL+SS: computeDDLSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*DD+LL: computeDSDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*DD+SS: computeDLDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*SS+8*LL+DD: computeDSLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*DD+64*LL+8*SS+DD: computeDLSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*DD+8*DD+LL: computeSDDLQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*DD+8*DD+SS: computeLDDSQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*DD+8*LL+DD: computeSDLDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*DD+8*SS+DD: computeLDSDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*SS+64*LL+8*DD+DD: computeSLDDQuartetEE(P,Q,I,J,EEStore); return; break;
	case 256*LL+64*SS+8*DD+DD: computeLSDDQuartetEE(P,Q,I,J,EEStore); return; break;
	}
#undef SS
#undef LL
#undef DD
*/
	// index preparation
	nEE = 0;
	for(p=0; p < P->nBasis; p++)
	for(q=0; q < Q->nBasis; q++)
	for(i=0; i < I->nBasis; i++)
	for(j=0; j < J->nBasis; j++){

		// precompute index for efficientcy
		kZ =    (J->maxL+1);
		kY = kZ*(I->maxL+1);
		kX = kY*(Q->maxL+1);
		iBx[nEE] = 4*(MAXL+1)*(P->l[p]*kX + Q->l[q]*kY + I->l[i]*kZ + J->l[j]);
		iBy[nEE] = 4*(MAXL+1)*(P->m[p]*kX + Q->m[q]*kY + I->m[i]*kZ + J->m[j]);
		iBz[nEE] = 4*(MAXL+1)*(P->n[p]*kX + Q->n[q]*kY + I->n[i]*kZ + J->n[j]);		
		mX[nEE]  = P->l[p] + Q->l[q] + I->l[i] + J->l[j];
		mY[nEE]  = P->m[p] + Q->m[q] + I->m[i] + J->m[j];
		mZ[nEE]  = P->n[p] + Q->n[q] + I->n[i] + J->n[j];
		if(mX[nEE]==0)                                  mT[nEE] = MXZERO;
		else if(mY[nEE]==0)                             mT[nEE] = MYZERO;
		else if(mZ[nEE]==0)                             mT[nEE] = MZZERO;
		else if(mX[nEE] > mY[nEE] && mX[nEE] > mZ[nEE]) mT[nEE] = MXMAX;
		else if(mY[nEE] > mZ[nEE])                      mT[nEE] = MYMAX;
		else                                            mT[nEE] = MZMAX;
		
		// reset EEStore to zero
		EEStore[nEE] = 0.0;

		// increase number of EE index
		nEE++;
	}

	for(pC=0; pC < P->nContract; pC++)
	for(qC=0; qC < Q->nContract; qC++)
	for(iC=0; iC < I->nContract; iC++)
	for(jC=0; jC < J->nContract; jC++){

		r=genSetBxyzSxyzF(P->x,P->y,P->z,P->maxL,P->exps[pC],
			              Q->x,Q->y,Q->z,Q->maxL,Q->exps[qC],
			              I->x,I->y,I->z,I->maxL,I->exps[iC],
			              J->x,J->y,J->z,J->maxL,J->exps[jC],
		                  Bx,By,Bz,Sx,Sy,Sz,F,MAXL);
		if(r < PRIMITIVE_CUTOFF) continue;

		// loop all basis within shell
		nEE = 0;
		for(p=0; p < P->nBasis; p++)
		for(q=0; q < Q->nBasis; q++)
		for(i=0; i < I->nBasis; i++){
		rpqi = r * P->coef[pC][p] * Q->coef[qC][q] * I->coef[iC][i];
		for(j=0; j < J->nBasis; j++){


			// compute two-electron integral using summation
			EE=0.0;

			switch(mT[nEE]){

			case MXZERO:
				tBy = By + iBy[nEE];
				tSz = Sz + iBz[nEE];
				for(kY=mY[nEE];kY>=0;kY--) EE += tBy[kY]*tSz[kY];
			break;

			case MYZERO:
				tBz = Bz + iBz[nEE];
				tSx = Sx + iBx[nEE];
				for(kZ=mZ[nEE];kZ>=0;kZ--) EE += tBz[kZ]*tSx[kZ];
			break;

			case MZZERO:
				tBx = Bx + iBx[nEE];
				tSy = Sy + iBy[nEE];
				for(kX=mX[nEE];kX>=0;kX--) EE += tBx[kX]*tSy[kX];
			break;

			case MXMAX: 
				tBz = Bz + iBz[nEE];
				tSx = Sx + iBx[nEE];
				tBy = By + iBy[nEE];
				for(kY=mY[nEE];kY>=0;kY--){
					for(tSum=0.0,kZ=mZ[nEE];kZ>=0;kZ--) tSum+=tBz[kZ]*tSx[kY+kZ];
					EE += tBy[kY]*tSum;
				}
			break;

			case MYMAX:
				tBx = Bx + iBx[nEE];
				tSy = Sy + iBy[nEE];
				tBz = Bz + iBz[nEE];
				for(kZ=mZ[nEE];kZ>=0;kZ--){
					for(tSum=0.0,kX=mX[nEE];kX>=0;kX--) tSum+=tBx[kX]*tSy[kX+kZ];
					EE += tBz[kZ]*tSum;
				}
			break;

			case MZMAX:
				tBy = By + iBy[nEE];
				tSz = Sz + iBz[nEE];
				tBx = Bx + iBx[nEE];
				for(kX=mX[nEE];kX>=0;kX--){
					for(tSum=0.0,kY=mY[nEE];kY>=0;kY--) tSum+=tBy[kY]*tSz[kX+kY];
					EE += tBx[kX]*tSum;
				}
			break;
			}

			EEStore[nEE] += EE * rpqi * J->coef[jC][j];
			nEE++;
		}
		}
	}
}


// storeQuartetEE : put the integral in the fock matrix by multiplying with
// the electron densities. Pre "unrolled" cases are also available. This
// subroutine is only for the case "PQIJ_ALLDISTINCT".
//
// Feb 2013 - Teepanis Chachiyo
//  - Initial implementation and testing
//
// 24 Feb, 2013 - Teepanis Chachiyo
//  - Put the code in C preprocessor form and now supports
//    all type of shells
//  - Use different code for restricted and unrestricted case
//
// 5 May, 2016 - Teepanis Chachiyo
//   Separate total Coulomb matrix (JT) from the exchange KA,KB
//
void storeQuartetEE(
	const struct GTOShell_t *P,
	const struct GTOShell_t *Q,
	const struct GTOShell_t *I,
	const struct GTOShell_t *J,
	char pqijT,
	int restricted,
	int nBasis,
	double *JT, double *KA, double *KB,
	const double *PT, const double *PA, const double *PB,
	const double *EEStore){

	int p,q,i,j;      // basis index
	int pp,qq,ii,jj;  // multiplier index
	const double *PTptr, *PAptr, *PBptr;
	register double Ga,Gb,EE;
/* -- Commented out by Teepanis Chachiyo on 16/9/2017. Error detected when computing
   -- Argon atom using ccTZP basis set
#define SS 0
#define LL 1
#define DD 3
	if(pqijT==PQIJ_ALLDISTINCT){

	if(restricted)
	switch(256*P->type + 64*Q->type + 8*I->type + J->type){

	case 256*SS+64*SS+8*SS+SS: RstoreSSSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*LL+LL: RstoreLLLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*DD+DD: RstoreDDDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	// SL pair
	case 256*SS+64*LL+8*LL+LL: RstoreSLLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*LL+LL: RstoreLSLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*SS+LL: RstoreLLSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*LL+SS: RstoreLLLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	case 256*LL+64*SS+8*SS+SS: RstoreLSSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*SS+SS: RstoreSLSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*SS+8*LL+SS: RstoreSSLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*SS+8*SS+LL: RstoreSSSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	case 256*SS+64*SS+8*LL+LL: RstoreSSLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*SS+LL: RstoreSLSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*LL+SS: RstoreSLLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*SS+LL: RstoreLSSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*LL+SS: RstoreLSLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*SS+SS: RstoreLLSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	// DS pair
	case 256*DD+64*SS+8*SS+SS: RstoreDSSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*DD+8*SS+SS: RstoreSDSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*SS+8*DD+SS: RstoreSSDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*SS+8*SS+DD: RstoreSSSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	case 256*SS+64*DD+8*DD+SS: RstoreSDDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*DD+8*SS+DD: RstoreSDSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*SS+8*DD+DD: RstoreSSDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*SS+DD: RstoreDSSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*DD+SS: RstoreDSDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*SS+SS: RstoreDDSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	case 256*SS+64*DD+8*DD+DD: RstoreSDDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*DD+DD: RstoreDSDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*SS+DD: RstoreDDSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*DD+SS: RstoreDDDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	// DL pair
	case 256*DD+64*LL+8*LL+LL: RstoreDLLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*DD+8*LL+LL: RstoreLDLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*DD+LL: RstoreLLDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*LL+DD: RstoreLLLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	case 256*LL+64*DD+8*DD+DD: RstoreLDDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*DD+DD: RstoreDLDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*LL+DD: RstoreDDLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*DD+LL: RstoreDDDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	case 256*LL+64*DD+8*DD+LL: RstoreLDDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*DD+8*LL+DD: RstoreLDLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*DD+DD: RstoreLLDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*LL+DD: RstoreDLLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*DD+LL: RstoreDLDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*LL+LL: RstoreDDLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	// 2S and DL pair
	case 256*SS+64*SS+8*DD+LL: RstoreSSDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*SS+8*LL+DD: RstoreSSLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*DD+8*SS+LL: RstoreSDSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*SS+DD: RstoreSLSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*DD+8*LL+SS: RstoreSDLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*DD+SS: RstoreSLDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*SS+LL: RstoreDSSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*SS+DD: RstoreLSSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*LL+SS: RstoreDSLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*DD+SS: RstoreLSDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*SS+SS: RstoreDLSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*DD+8*SS+SS: RstoreLDSSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	// 2L and SD pair
	case 256*LL+64*LL+8*SS+DD: RstoreLLSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*LL+8*DD+SS: RstoreLLDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*LL+DD: RstoreLSLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*DD+8*LL+SS: RstoreLDLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*DD+LL: RstoreLSDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*DD+8*SS+LL: RstoreLDSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*LL+DD: RstoreSLLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*LL+SS: RstoreDLLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*DD+LL: RstoreSLDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*SS+LL: RstoreDLSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*DD+8*LL+LL: RstoreSDLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*LL+LL: RstoreDSLLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;

	// 2D and SL pair
	case 256*DD+64*DD+8*SS+LL: RstoreDDSLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*DD+8*LL+SS: RstoreDDLSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*DD+LL: RstoreDSDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*DD+SS: RstoreDLDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*SS+8*LL+DD: RstoreDSLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*DD+64*LL+8*SS+DD: RstoreDLSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*DD+8*DD+LL: RstoreSDDLQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*DD+8*DD+SS: RstoreLDDSQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*DD+8*LL+DD: RstoreSDLDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*DD+8*SS+DD: RstoreLDSDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*SS+64*LL+8*DD+DD: RstoreSLDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	case 256*LL+64*SS+8*DD+DD: RstoreLSDDQuartetEE(P,Q,I,J,nBasis,JT,KA,PT,PA,EEStore); return; break;
	}

	else
	switch(256*P->type + 64*Q->type + 8*I->type + J->type){

	case 256*SS+64*SS+8*SS+SS: UstoreSSSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*LL+LL: UstoreLLLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*DD+DD: UstoreDDDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	// SL pair
	case 256*SS+64*LL+8*LL+LL: UstoreSLLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*LL+LL: UstoreLSLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*SS+LL: UstoreLLSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*LL+SS: UstoreLLLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	case 256*LL+64*SS+8*SS+SS: UstoreLSSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*SS+SS: UstoreSLSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*SS+8*LL+SS: UstoreSSLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*SS+8*SS+LL: UstoreSSSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	case 256*SS+64*SS+8*LL+LL: UstoreSSLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*SS+LL: UstoreSLSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*LL+SS: UstoreSLLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*SS+LL: UstoreLSSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*LL+SS: UstoreLSLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*SS+SS: UstoreLLSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	// DS pair
	case 256*DD+64*SS+8*SS+SS: UstoreDSSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*DD+8*SS+SS: UstoreSDSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*SS+8*DD+SS: UstoreSSDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*SS+8*SS+DD: UstoreSSSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	case 256*SS+64*DD+8*DD+SS: UstoreSDDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*DD+8*SS+DD: UstoreSDSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*SS+8*DD+DD: UstoreSSDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*SS+DD: UstoreDSSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*DD+SS: UstoreDSDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*SS+SS: UstoreDDSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	case 256*SS+64*DD+8*DD+DD: UstoreSDDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*DD+DD: UstoreDSDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*SS+DD: UstoreDDSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*DD+SS: UstoreDDDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	// DL pair
	case 256*DD+64*LL+8*LL+LL: UstoreDLLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*DD+8*LL+LL: UstoreLDLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*DD+LL: UstoreLLDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*LL+DD: UstoreLLLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	case 256*LL+64*DD+8*DD+DD: UstoreLDDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*DD+DD: UstoreDLDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*LL+DD: UstoreDDLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*DD+LL: UstoreDDDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	case 256*LL+64*DD+8*DD+LL: UstoreLDDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*DD+8*LL+DD: UstoreLDLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*DD+DD: UstoreLLDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*LL+DD: UstoreDLLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*DD+LL: UstoreDLDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*LL+LL: UstoreDDLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	// 2S and DL pair
	case 256*SS+64*SS+8*DD+LL: UstoreSSDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*SS+8*LL+DD: UstoreSSLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*DD+8*SS+LL: UstoreSDSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*SS+DD: UstoreSLSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*DD+8*LL+SS: UstoreSDLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*DD+SS: UstoreSLDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*SS+LL: UstoreDSSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*SS+DD: UstoreLSSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*LL+SS: UstoreDSLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*DD+SS: UstoreLSDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*SS+SS: UstoreDLSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*DD+8*SS+SS: UstoreLDSSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	// 2L and SD pair
	case 256*LL+64*LL+8*SS+DD: UstoreLLSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*LL+8*DD+SS: UstoreLLDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*LL+DD: UstoreLSLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*DD+8*LL+SS: UstoreLDLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*DD+LL: UstoreLSDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*DD+8*SS+LL: UstoreLDSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*LL+DD: UstoreSLLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*LL+SS: UstoreDLLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*DD+LL: UstoreSLDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*SS+LL: UstoreDLSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*DD+8*LL+LL: UstoreSDLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*LL+LL: UstoreDSLLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;

	// 2D and SL pair
	case 256*DD+64*DD+8*SS+LL: UstoreDDSLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*DD+8*LL+SS: UstoreDDLSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*DD+LL: UstoreDSDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*DD+SS: UstoreDLDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*SS+8*LL+DD: UstoreDSLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*DD+64*LL+8*SS+DD: UstoreDLSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*DD+8*DD+LL: UstoreSDDLQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*DD+8*DD+SS: UstoreLDDSQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*DD+8*LL+DD: UstoreSDLDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*DD+8*SS+DD: UstoreLDSDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*SS+64*LL+8*DD+DD: UstoreSLDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	case 256*LL+64*SS+8*DD+DD: UstoreLSDDQuartetEE(P,Q,I,J,nBasis,JT,KA,KB,PT,PA,PB,EEStore); return; break;
	}


	}
#undef SS
#undef LL
#undef DD
*/
	jj = 1;
	ii = J->nBasis;
	qq = I->nBasis * J->nBasis;
	pp = Q->nBasis * I->nBasis * J->nBasis;

	const double *EEptr;

	//
	// the restricted version deals with only KA
	//
#define RSTOREJ_DIDI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)              \
	for(p=0; p < P->nBasis; p++)                               \
	for(q=0; q < Q->nBasis; q++){                              \
		PTptr = PT + I->min * nBasis + J->min;                 \
		EEptr = EEStore + p*pp+q*qq;                           \
		Ga=0.0;                                                \
		for(i=0; i < I->nBasis; i++,EEptr+=ii,PTptr+=nBasis)   \
		for(j=0; j < J->nBasis; j++){                          \
			EE  = EEptr[j*jj]; EE += EE;                       \
			Ga += PTptr[j]*EE;                                 \
		}                                                      \
		JT[(p+P->min)*nBasis+(q+Q->min)] += Ga;                \
		JT[(q+Q->min)*nBasis+(p+P->min)] += Ga;                \
	}

#define RSTOREJ_DISA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)              \
	for(p=0; p < P->nBasis; p++)                               \
	for(q=0; q < Q->nBasis; q++){                              \
		PTptr = PT + I->min * nBasis + J->min;                 \
		EEptr = EEStore + p*pp+q*qq;                           \
		Ga=0.0;                                                \
		for(i=0; i < I->nBasis; i++,EEptr+=ii,PTptr+=nBasis)   \
		for(j=0; j < J->nBasis; j++){                          \
			EE  = EEptr[j*jj];                                 \
			Ga += PTptr[j]*EE;                                 \
		}                                                      \
		JT[(p+P->min)*nBasis+(q+Q->min)] += Ga;                \
		JT[(q+Q->min)*nBasis+(p+P->min)] += Ga;                \
	}

#define RSTOREJ_SADI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)              \
	for(p=0; p < P->nBasis; p++)                               \
	for(q=0; q < Q->nBasis; q++){                              \
		PTptr = PT + I->min * nBasis + J->min;                 \
		EEptr = EEStore + p*pp+q*qq;                           \
		Ga=0.0;                                                \
		for(i=0; i < I->nBasis; i++,EEptr+=ii,PTptr+=nBasis)   \
		for(j=0; j < J->nBasis; j++){                          \
			EE  = EEptr[j*jj]; EE += EE;                       \
			Ga += PTptr[j]*EE;                                 \
		}                                                      \
		JT[(p+P->min)*nBasis+(q+Q->min)] += Ga;                \
	}

#define RSTOREJ_SASA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)              \
	for(p=0; p < P->nBasis; p++)                               \
	for(q=0; q < Q->nBasis; q++){                              \
		PTptr = PT + I->min * nBasis + J->min;                 \
		EEptr = EEStore + p*pp+q*qq;                           \
		Ga=0.0;                                                \
		for(i=0; i < I->nBasis; i++,EEptr+=ii,PTptr+=nBasis)   \
		for(j=0; j < J->nBasis; j++){                          \
			EE  = EEptr[j*jj];                                 \
			Ga += PTptr[j]*EE;                                 \
		}                                                      \
		JT[(p+P->min)*nBasis+(q+Q->min)] += Ga;                \
	}

#define RSTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)              \
	for(p=0; p < P->nBasis; p++)                               \
	for(i=0; i < I->nBasis; i++){                              \
		PAptr = PA + Q->min * nBasis + J->min;                 \
		EEptr = EEStore + p*pp+i*ii;                           \
		Ga=0.0;                                                \
		for(q=0; q < Q->nBasis; q++,EEptr+=qq,PAptr+=nBasis)   \
		for(j=0; j < J->nBasis; j++){                          \
			EE  = EEptr[j*jj];                                 \
			Ga += PAptr[j]*EE;                                 \
		}                                                      \
		KA[(p+P->min)*nBasis+(i+I->min)] -= Ga;                \
		KA[(i+I->min)*nBasis+(p+P->min)] -= Ga;                \
	}

#define RSTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)              \
	for(p=0; p < P->nBasis; p++)                               \
	for(i=0; i < I->nBasis; i++){                              \
		PAptr = PA + Q->min * nBasis + J->min;                 \
		EEptr = EEStore + p*pp+i*ii;                           \
		Ga=0.0;                                                \
		for(q=0; q < Q->nBasis; q++,EEptr+=qq,PAptr+=nBasis)   \
		for(j=0; j < J->nBasis; j++){                          \
			EE  = EEptr[j*jj];                                 \
			Ga += PAptr[j]*EE;                                 \
		}                                                      \
		KA[(p+P->min)*nBasis+(i+I->min)] -= Ga;                \
	}

	//
	// the unrestricted version deals with both KA and KB
	//
#define USTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)                          \
	for(p=0; p < P->nBasis; p++)                                           \
	for(i=0; i < I->nBasis; i++){                                          \
		PAptr = PA + Q->min * nBasis + J->min;                             \
		PBptr = PB + Q->min * nBasis + J->min;                             \
		EEptr = EEStore + p*pp+i*ii;                                       \
		Ga=0.0; Gb=0.0;                                                    \
		for(q=0; q < Q->nBasis; q++,EEptr+=qq,PAptr+=nBasis,PBptr+=nBasis) \
		for(j=0; j < J->nBasis; j++){                                      \
			EE  = EEptr[j*jj];                                             \
			Ga += PAptr[j]*EE;                                             \
			Gb += PBptr[j]*EE;                                             \
		}                                                                  \
		KA[(p+P->min)*nBasis+(i+I->min)] -= Ga;                            \
		KB[(p+P->min)*nBasis+(i+I->min)] -= Gb;                            \
		KA[(i+I->min)*nBasis+(p+P->min)] -= Ga;                            \
		KB[(i+I->min)*nBasis+(p+P->min)] -= Gb;                            \
	}

#define USTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj)                          \
	for(p=0; p < P->nBasis; p++)                                           \
	for(i=0; i < I->nBasis; i++){                                          \
		PAptr = PA + Q->min * nBasis + J->min;                             \
		PBptr = PB + Q->min * nBasis + J->min;                             \
		EEptr = EEStore + p*pp+i*ii;                                       \
		Ga=0.0; Gb=0.0;                                                    \
		for(q=0; q < Q->nBasis; q++,EEptr+=qq,PAptr+=nBasis,PBptr+=nBasis) \
		for(j=0; j < J->nBasis; j++){                                      \
			EE  = EEptr[j*jj];                                             \
			Ga += PAptr[j]*EE;                                             \
			Gb += PBptr[j]*EE;                                             \
		}                                                                  \
		KA[(p+P->min)*nBasis+(i+I->min)] -= Ga;                            \
		KB[(p+P->min)*nBasis+(i+I->min)] -= Gb;                            \
	}

	if(restricted)
	switch(pqijT){
	case PQIJ_ALLSAME:
		RSTOREJ_SASA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
	break;

	case PQIJ_2PAIRS:
		RSTOREJ_SASA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_SASA(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		RSTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
	break;

	case PQIJ_PQPAIR:
		RSTOREJ_SADI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_DISA(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		RSTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREK_DIFF(P,Q,J,I,p,q,j,i,pp,qq,jj,ii);
	break;

	case PQIJ_IJPAIR:
		RSTOREJ_DISA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_SADI(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		RSTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREK_DIFF(Q,P,I,J,q,p,i,j,qq,pp,ii,jj);
	break;

	case PQIJ_PIQJPAIR:
		RSTOREJ_DIDI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREK_SAME(P,Q,J,I,p,q,j,i,pp,qq,jj,ii);
		RSTOREK_SAME(Q,P,I,J,q,p,i,j,qq,pp,ii,jj);
		RSTOREK_SAME(Q,P,J,I,q,p,j,i,qq,pp,jj,ii);
	break;

	case PQIJ_ALLDISTINCT:
		RSTOREJ_DIDI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_DIDI(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		RSTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREK_DIFF(P,Q,J,I,p,q,j,i,pp,qq,jj,ii);
		RSTOREK_DIFF(Q,P,I,J,q,p,i,j,qq,pp,ii,jj);
		RSTOREK_DIFF(Q,P,J,I,q,p,j,i,qq,pp,jj,ii);
	break;
	}

	else
	switch(pqijT){
	case PQIJ_ALLSAME:
		RSTOREJ_SASA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		USTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
	break;

	case PQIJ_2PAIRS:
		RSTOREJ_SASA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_SASA(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		USTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
	break;

	case PQIJ_PQPAIR:
		RSTOREJ_SADI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_DISA(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		USTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		USTOREK_DIFF(P,Q,J,I,p,q,j,i,pp,qq,jj,ii);
	break;

	case PQIJ_IJPAIR:
		RSTOREJ_DISA(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_SADI(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		USTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		USTOREK_DIFF(Q,P,I,J,q,p,i,j,qq,pp,ii,jj);
	break;

	case PQIJ_PIQJPAIR:
		RSTOREJ_DIDI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		USTOREK_SAME(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		USTOREK_SAME(P,Q,J,I,p,q,j,i,pp,qq,jj,ii);
		USTOREK_SAME(Q,P,I,J,q,p,i,j,qq,pp,ii,jj);
		USTOREK_SAME(Q,P,J,I,q,p,j,i,qq,pp,jj,ii);
	break;

	case PQIJ_ALLDISTINCT:
		RSTOREJ_DIDI(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		RSTOREJ_DIDI(I,J,P,Q,i,j,p,q,ii,jj,pp,qq);
		USTOREK_DIFF(P,Q,I,J,p,q,i,j,pp,qq,ii,jj);
		USTOREK_DIFF(P,Q,J,I,p,q,j,i,pp,qq,jj,ii);
		USTOREK_DIFF(Q,P,I,J,q,p,i,j,qq,pp,ii,jj);
		USTOREK_DIFF(Q,P,J,I,q,p,j,i,qq,pp,jj,ii);
	break;
	}
}


// getMemCutoff : get cutoff for memeory of the schwarz matrix so that
// the amount of memory need for storage does not exceed user's defined
// maximum value.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
double getMemCutoff(
	int nShell,
	const struct GTOShell_t *shell,
	double fixedCutoff,
	const double *schwarz_shell,
	const struct option_t *opt){

	int P,Q,I,J,N;                    // shell index
	double cutmem;                    // cutoff value for available memory
	unsigned long maxEE,nEEBank;      // maximum possible number of integrals
	double r;                         // generic double precision

	// increase cutmem until it is enough for available memory
	cutmem = fixedCutoff;
	maxEE  = (unsigned long)(opt->maxMem)*1000000/sizeof(double);

	// set corse search flag
	r = 1;
	do{
CutMemSearch:
		nEEBank = 0;

		// loop all distinct shell sequence
		for(P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
		for(J=0; J <= N; J++){

		
			// screen integral using Schwarz inequality at shell level
			if(schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J] < cutmem)
				continue;
		
			// accumulate the number of shell permutation
			nEEBank += shell[P].nBasis * shell[Q].nBasis * 
			           shell[I].nBasis * shell[J].nBasis;
			if(nEEBank > maxEE){
				if(r==1) // corse search
					cutmem = cutmem*2.00;
				else     // fine search
					cutmem = cutmem*1.05;
				goto CutMemSearch;
			}

		}
		}

		// done corse search and begin fine search
		if(r==1 && nEEBank < maxEE){ r = 0; cutmem/=2.0; goto CutMemSearch; }

	}while(nEEBank > maxEE);

	// cutmem should not be smaller than fixedCutoff
	if(cutmem < fixedCutoff) cutmem = fixedCutoff;

	return cutmem;
}


// GTO_JK_Matrix_Quartet : calculate G matrix for unrestricted hartree-fock
// calculations. The equations are given in (Szabo and Ostlund, 1989)
// 
// for Alpha Spin
// [GB]pq = SUM [PT]ij*(pq|ij) - [PA]ij*(pi|qj) 
//
// for Beta Spin
// [GA]pg = SUM [PT]ij*(pq|ij) - [PB]ij*(pi|qj)
//
// where PT is the toal spin density and is defined by
// PT = PA+PB
//
// As of Feb 2013, this is the faster version of the code the compute
// fock matrix.
//
// Feb 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
// Feb 24, 2013 - Teepanis Chachiyo
//     No longer call storeJK and use storeQuartetEE instead
//
/*
void GTO_JK_Matrix_Quartet(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *schwarz_basis,   // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB,                    // return G for spin down
	struct option_t *opt){         // global option

	static int firstVisit=1;              // first visit flag
	static struct GTOShell_t *shell=NULL; // shell data structure
	static int nShell;                    // the number of shell
	int P,Q,I,J,N;                        // shell loop index
	int p,q;                              // basis loop index
	static double *schwarz_shell=NULL;    // schwarz matrix at shell level
	static double *PM=NULL;               // maximum density at shell level
	static double *PT=NULL;               // total density matrix
	double r;                             // generic double precision
	double maxPShell;                     // maximum density in the shell
	double EEStore[MAXSHELLINT];          // electron integral storage
	char pqijT;                           // type of shell quartet
	int nEE;                              // number of integral within shell

	static struct Multipole_t *mp=NULL;   // multipole data
	static signed char *erBank=NULL;      // multipole error bank at shell level
	signed char *erBankPtr=NULL;          // pointer to current value
	signed char erExp;                    // current value of error

	static double *eeBank=NULL;           // global 2e storage
	double *eeBankPtr=NULL;               // pointer to current 2e storage
	static double memCutoff=0.0;          // schwartz cutoff to store in memory

#define ALLOC(array,item,type)                                        \
array=calloc(item,sizeof(type));                                      \
if(array==NULL){                                                      \
	printf("GTO_JK_Matrix_Quartet - error cannot allocate memory\n"); \
	exit(-1);                                                         \
}

	// reset call
	if(nBasis==0){
		firstVisit = 1;
		if(shell != NULL){ free(shell); shell=NULL; }
		if(schwarz_shell != NULL){ free(schwarz_shell); schwarz_shell=NULL; }
		if(PM != NULL){ free(PM); PM=NULL; }
		if(PT != NULL){ free(PT); PT=NULL; }
		if(mp != NULL){ free(mp); mp=NULL; }
		if(erBank != NULL){ free(erBank); erBank=NULL; }
		if(eeBank != NULL){ free(eeBank); eeBank=NULL; }
		return;
	}

	// first visit preparations
	if(firstVisit){

		// generate shell
		shell = genGTOShell(nBasis, gto, &nShell);
		sortGTOShell(nShell, shell);

		// generate multipole data
		mp = genMultipole(nBasis, gto);

		// allocate memory
		ALLOC(schwarz_shell, nShell*nShell, double);
		ALLOC(PM, nShell*nShell, double);
		ALLOC(PT, nBasis*nBasis, double);

		// compute schwarz at shell level
		for(P=0; P < nShell; P++)
		for(Q=0; Q < nShell; Q++){
			r = 0.0;
			for(p=shell[P].min; p <= shell[P].max; p++)
			for(q=shell[Q].min; q <= shell[Q].max; q++)
				if(r <schwarz_basis[p*nBasis+q]) r = schwarz_basis[p*nBasis+q];
			schwarz_shell[P*nShell+Q] = r;
		}

		// count the number shells needed and allocate memory
		unsigned long nEE=0;
		for(P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
		for(J=0; J <= N; J++){
			r = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
			if(r < opt->SCFCutoff) continue;
			nEE++;
		}
		}
		ALLOC(erBank, nEE, signed char);
		memset(erBank, ERROR_UNKNOWN, nEE);

		// determine memory cutoff, number of integrals, and allocate memory
		memCutoff = getMemCutoff(nShell, shell, fixedCutoff, schwarz_shell, opt);
		nEE=0;
		for(P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
		for(J=0; J <= N; J++){
			if(schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J] < memCutoff)
				continue;

			nEE += shell[P].nBasis * shell[Q].nBasis * 
			       shell[I].nBasis * shell[J].nBasis;
		}
		}
		ALLOC(eeBank, nEE, double);

	}

	// compute total density matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		PT[p*nBasis+q] = PA[p*nBasis+q] + PB[p*nBasis+q];

	// compute schwarz and maximum density at shell level
	for(P=0; P < nShell; P++)
	for(Q=0; Q < nShell; Q++){
		r = 0.0;
		for(p=shell[P].min; p <= shell[P].max; p++)
		for(q=shell[Q].min; q <= shell[Q].max; q++){
			if(r<fabs(PT[p*nBasis+q])) r = fabs(PT[p*nBasis+q]);
			if(r<fabs(PA[p*nBasis+q])) r = fabs(PA[p*nBasis+q]);
			if(r<fabs(PB[p*nBasis+q])) r = fabs(PB[p*nBasis+q]);
		}
		PM[P*nShell+Q] = r;
	}

	// main loop
	erBankPtr = erBank;
	eeBankPtr = eeBank;
	for(P=0; P < nShell; P++)
	for(Q=0; Q <= P; Q++)
	for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
	for(J=0; J <= N; J++){

		// schwarz screening
		r = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
		if(r < opt->SCFCutoff) continue;
		erExp = *erBankPtr; erBankPtr++;

		if(r < fixedCutoff) continue;

		// maximum density in this shell
		maxPShell = PM[P*nShell+Q] + PM[I*nShell+J] +
		            PM[P*nShell+I] + PM[P*nShell+J] +
		            PM[Q*nShell+I] + PM[Q*nShell+J];

		// load from memory if possible
		if(r >= memCutoff){
			nEE = shell[P].nBasis * shell[Q].nBasis *
			      shell[I].nBasis * shell[J].nBasis;
			if(firstVisit){
				computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);
				memcpy(eeBankPtr, EEStore, sizeof(double) * nEE);
				eeBankPtr += nEE;
			}else{
				// screening again with density weighted
				if(r*maxPShell < fixedCutoff){ eeBankPtr += nEE; continue; }

				memcpy(EEStore, eeBankPtr, sizeof(double) * nEE);
				eeBankPtr += nEE;
			}
		}else{
			// screening again with density weighted
			if(r*maxPShell < fixedCutoff) continue;

			// compute the entire shell quartet
			if(erExp!=ERROR_UNKNOWN && pow(ERROR_BASE,erExp)*maxPShell < fixedCutoff)
				computeQuartetMPole(shell+P,shell+Q,shell+I,shell+J,mp,nBasis,EEStore);
			else
				computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);
		}

		// store multipole error
		if(erExp==ERROR_UNKNOWN)
			erBankPtr[-1]=getErrExpQuartetMPole(shell+P,shell+Q,shell+I,shell+J,
			                                    mp,nBasis,EEStore);

		// determine P Q I J type for storage
		if((P==Q)&&(I==J)&&(P==I)){ pqijT = PQIJ_ALLSAME;
		}else if((P==Q)&&(I==J)){   pqijT = PQIJ_2PAIRS;
		}else if(P==Q){             pqijT = PQIJ_PQPAIR;
		}else if(I==J){             pqijT = PQIJ_IJPAIR;
		}else if((P==I)&&(Q==J)){   pqijT = PQIJ_PIQJPAIR;
		}else{                      pqijT = PQIJ_ALLDISTINCT;
		}

		// store in the G matrix
		storeQuartetEE(shell+P,shell+Q,shell+I,shell+J,
		               pqijT,opt->RHF,nBasis,GA,GB,PT,PA,PB,EEStore);

	}
	}

	// handle restricted case
	if(opt->RHF) for(p=0; p < nBasis*nBasis; p++) GB[p] = GA[p];

	// symmetrize G matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < p; q++){
		GA[q*nBasis+p] = GA[p*nBasis+q];
		GB[q*nBasis+p] = GB[p*nBasis+q];
	}

	// reset flag
	if(firstVisit) firstVisit=0;
}
*/

// getMemCutoff_Parallel: this is  the parallel version of the getMemCutoff
//
// Mar 3, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
double getMemCutoff_Parallel(
	int childID,
	int nCPU,
	int nShell,
	const struct GTOShell_t *shell,
	double fixedCutoff,
	const double *schwarz_shell,
	const struct option_t *opt){

	int PQ,P,Q,I,J,N;                 // shell index
	double cutmem;                    // cutoff value for available memory
	unsigned long maxEE,nEEBank;      // maximum possible number of integrals
	double r;                         // generic double precision

	// increase cutmem until it is enough for available memory
	cutmem = fixedCutoff;
	maxEE  = (unsigned long)(opt->maxMem)*1000000/sizeof(double);

	// set corse search flag
	r = 1;
	do{
CutMemSearch:
		nEEBank = 0;

		// loop all distinct shell sequence
		for(PQ=0,P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++,PQ++) if(PQ%nCPU==childID)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
		for(J=0; J <= N; J++){

		
			// screen integral using Schwarz inequality at shell level
			if(schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J] < cutmem)
				continue;
		
			// accumulate the number of shell permutation
			nEEBank += shell[P].nBasis * shell[Q].nBasis * 
			           shell[I].nBasis * shell[J].nBasis;
			if(nEEBank > maxEE){
				if(r==1) // corse search
					cutmem = cutmem*2.00;
				else     // fine search
					cutmem = cutmem*1.05;
				goto CutMemSearch;
			}

		}
		}

		// done corse search and begin fine search
		if(r==1 && nEEBank < maxEE){ r = 0; cutmem/=2.0; goto CutMemSearch; }

	}while(nEEBank > maxEE);

	// cutmem should not be smaller than fixedCutoff
	if(cutmem < fixedCutoff) cutmem = fixedCutoff;

	return cutmem;
}


// GTO_JK_Matrix_Quartet_Parallel: this is a parallel version of the
// GTO_JK_Matrix_Quartet. See the original subroutine for information.
//
// May 5, 2016 - Teepanis Chachiyo
//     Separate total Coulomb matrix (JT) from the exchange KA,KB
//
// Mar 3, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void GTO_JK_Matrix_Quartet_Parallel(
	int childID,                   // child id number
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *schwarz_basis,   // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *JT,                    // returned total Coulomb matrix
	double *KA,                    // returned Exchange for spin up
	double *KB,                    // returned Exchange for spin down
	struct option_t *opt){         // global option

	static int firstVisit=1;              // first visit flag
	static struct GTOShell_t *shell=NULL; // shell data structure
	static int nShell;                    // the number of shell
	int PQ,P,Q,I,J,N;                     // shell loop index
	int p,q;                              // basis loop index
	static double *schwarz_shell=NULL;    // schwarz matrix at shell level
	static double *PM=NULL;               // maximum density at shell level
	static double *PT=NULL;               // total density matrix
	double r;                             // generic double precision
	double maxPShell;                     // maximum density in the shell
	double EEStore[MAXSHELLINT];          // electron integral storage
	char pqijT;                           // type of shell quartet
	int nEE;                              // number of integral within shell

	static struct Multipole_t *mp=NULL;   // multipole data
	static signed char *erBank=NULL;      // multipole error bank at shell level
	signed char *erBankPtr=NULL;          // pointer to current value
	signed char erExp;                    // current value of error

	static double *eeBank=NULL;           // global 2e storage
	double *eeBankPtr=NULL;               // pointer to current 2e storage
	static double memCutoff=0.0;          // schwartz cutoff to store in memory

#define ALLOC(array,item,type)                                                 \
array=calloc(item,sizeof(type));                                               \
if(array==NULL){                                                               \
	printf("GTO_JK_Matrix_Quartet_Parallel - error cannot allocate memory\n"); \
	printf("tried to allocate %.1f MB\n",1.0e-6*item*sizeof(type));            \
	exit(-1);                                                                  \
}

	// reset call
	if(nBasis==0){
		firstVisit = 1;
		if(shell != NULL){ free(shell); shell=NULL; }
		if(schwarz_shell != NULL){ free(schwarz_shell); schwarz_shell=NULL; }
		if(PM != NULL){ free(PM); PM=NULL; }
		if(PT != NULL){ free(PT); PT=NULL; }
		if(mp != NULL){ free(mp); mp=NULL; }
		if(erBank != NULL){ free(erBank); erBank=NULL; }
		if(eeBank != NULL){ free(eeBank); eeBank=NULL; }
		return;
	}

	// first visit preparations
	if(firstVisit){

		// generate shell
		shell = genGTOShell(nBasis, gto, &nShell);
		sortGTOShell(nShell, shell);

		// generate multipole data
		mp = genMultipole(nBasis, gto);

		// allocate memory
		ALLOC(schwarz_shell, nShell*nShell, double);
		ALLOC(PM, nShell*nShell, double);
		ALLOC(PT, nBasis*nBasis, double);

		// compute schwarz at shell level
		for(P=0; P < nShell; P++)
		for(Q=0; Q < nShell; Q++){
			r = 0.0;
			for(p=shell[P].min; p <= shell[P].max; p++)
			for(q=shell[Q].min; q <= shell[Q].max; q++)
				if(r <schwarz_basis[p*nBasis+q]) r = schwarz_basis[p*nBasis+q];
			schwarz_shell[P*nShell+Q] = r;
		}

		// count the number shells needed and allocate memory
		unsigned long nEE=0;
		for(PQ=0,P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++,PQ++) if(PQ%(opt->nCPU)==childID)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
		for(J=0; J <= N; J++){
			r = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
			if(r < opt->SCFCutoff) continue;
			nEE++;
		}
		}
		ALLOC(erBank, nEE, signed char);
		memset(erBank, ERROR_UNKNOWN, nEE);

		// determine memory cutoff, number of integrals, and allocate memory
		memCutoff = getMemCutoff_Parallel(childID,opt->nCPU,nShell, shell,
		                                  fixedCutoff, schwarz_shell, opt);
		nEE=0;
		for(PQ=0,P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++,PQ++) if(PQ%(opt->nCPU)==childID)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
		for(J=0; J <= N; J++){
			if(schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J] < memCutoff)
				continue;

			nEE += shell[P].nBasis * shell[Q].nBasis * 
			       shell[I].nBasis * shell[J].nBasis;
		}
		}
		ALLOC(eeBank, nEE, double);

	}

	// compute total density matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		PT[p*nBasis+q] = PA[p*nBasis+q] + PB[p*nBasis+q];

	// compute schwarz and maximum density at shell level
	for(P=0; P < nShell; P++)
	for(Q=0; Q < nShell; Q++){
		r = 0.0;
		for(p=shell[P].min; p <= shell[P].max; p++)
		for(q=shell[Q].min; q <= shell[Q].max; q++){
			if(r<fabs(PT[p*nBasis+q])) r = fabs(PT[p*nBasis+q]);
			if(r<fabs(PA[p*nBasis+q])) r = fabs(PA[p*nBasis+q]);
			if(r<fabs(PB[p*nBasis+q])) r = fabs(PB[p*nBasis+q]);
		}
		PM[P*nShell+Q] = r;
	}

	// main loop
	erBankPtr = erBank;
	eeBankPtr = eeBank;
	for(PQ=0,P=0; P < nShell; P++)
	for(Q=0; Q <= P; Q++,PQ++) if(PQ%(opt->nCPU)==childID)
	for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
	for(J=0; J <= N; J++){

		// schwarz screening
		r = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
		if(r < opt->SCFCutoff) continue;
		erExp = *erBankPtr; erBankPtr++;

		if(r < fixedCutoff) continue;

		// maximum density in this shell
		maxPShell = PM[P*nShell+Q] + PM[I*nShell+J] +
		            PM[P*nShell+I] + PM[P*nShell+J] +
		            PM[Q*nShell+I] + PM[Q*nShell+J];

		// load from memory if possible
		if(r >= memCutoff){
			nEE = shell[P].nBasis * shell[Q].nBasis *
			      shell[I].nBasis * shell[J].nBasis;
			if(firstVisit){
				computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);
				memcpy(eeBankPtr, EEStore, sizeof(double) * nEE);
				eeBankPtr += nEE;
			}else{
				// screening again with density weighted
				if(r*maxPShell < fixedCutoff){ eeBankPtr += nEE; continue; }

				memcpy(EEStore, eeBankPtr, sizeof(double) * nEE);
				eeBankPtr += nEE;
			}
		}else{
			// screening again with density weighted
			if(r*maxPShell < fixedCutoff) continue;

			// compute the entire shell quartet
			if(erExp!=ERROR_UNKNOWN && pow(ERROR_BASE,erExp)*maxPShell < fixedCutoff)
				computeQuartetMPole(shell+P,shell+Q,shell+I,shell+J,mp,nBasis,EEStore);
			else
				computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);
		}

		// store multipole error
		if(erExp==ERROR_UNKNOWN)
			erBankPtr[-1]=getErrExpQuartetMPole(shell+P,shell+Q,shell+I,shell+J,
			                                    mp,nBasis,EEStore);

		// determine P Q I J type for storage
		if((P==Q)&&(I==J)&&(P==I)){ pqijT = PQIJ_ALLSAME;
		}else if((P==Q)&&(I==J)){   pqijT = PQIJ_2PAIRS;
		}else if(P==Q){             pqijT = PQIJ_PQPAIR;
		}else if(I==J){             pqijT = PQIJ_IJPAIR;
		}else if((P==I)&&(Q==J)){   pqijT = PQIJ_PIQJPAIR;
		}else{                      pqijT = PQIJ_ALLDISTINCT;
		}

		// store in the G matrix
		storeQuartetEE(shell+P,shell+Q,shell+I,shell+J,
		               pqijT,opt->restricted,nBasis,JT,KA,KB,PT,PA,PB,EEStore);

	}
	}

	// handle restricted case
	if(opt->restricted) for(p=0; p < nBasis*nBasis; p++) KB[p] = KA[p];

	// symmetrize G matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < p; q++){
		JT[q*nBasis+p] = JT[p*nBasis+q];
		KA[q*nBasis+p] = KA[p*nBasis+q];
		KB[q*nBasis+p] = KB[p*nBasis+q];
	}

	// reset flag
	if(firstVisit) firstVisit=0;
}
#undef ALLOC


//
// Below is the simple but un-optmized version of the GTO_JK_Matrix_ShellSet.
// I realized that as time progesses, the subroutine will get messy as we
// add memory management, parallelization, and so on. Therefore, the original
// version is given here for students to make sense of it much better.
//
// Warning! It should give the same answer as the optimized version, down
// to the last digit (12-15), so this subroutine serves as a final test 
// of any other implementation of the function.
//
// Do not make the same mistake I did, not testing the ZeroEE[nEE], and it
// is just a time bomb wating to happen with STO3G basis set; where as, the
// sympton does not appear with 321G. See note on May 23, 2011.
//
/*
void GTO_JK_Matrix_ShellSet(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB,                    // return G for spin down
	struct option_t *opt){         // global option

	double *PT;                    // total spin density matrix
	double *PTA, *PAB;             // for screening with density
	int pS,qS,iS,jS;               // looping variables for shells
	int p,q,i,j;                   // looping variables for basis functions
	int pC,qC,iC,jC;               // contracted index
	double EE;                     // coulomb integral
	int nEE;                       // number of electron within shell counter
	double EEStore[MAXSHELLINT];   // electron integral storage
	int        iBx[MAXSHELLINT];   // index for Bx storage
	int        iBy[MAXSHELLINT];   // index for By storage
	int        iBz[MAXSHELLINT];   // index for Bz storage
	char        mX[MAXSHELLINT];   // maximum in x direction
	char        mY[MAXSHELLINT];   // maximum in y direction
	char        mZ[MAXSHELLINT];   // maximum in z direction
	double Bx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double By[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Bz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double  F[4*(MAXL+1)];         // Boys functions
	double *tBx,*tBy,*tBz;         // pointer to Bx, By, Bz
	int kX, kY, kZ;                // summation loop
	char pqijT;                    // storage type for a set of pqij
	int nShell;                    // number of shell
	int *shellMap;                 // mapping between shell and basis function
	int *shellMaxL;                // maximum angular moment for each shell
	double *Schwarz_Shell;         // schwarz screening at shell level
	int n;                         // generic index
	double r;                      // generic real number

	if(nBasis<=0) return;

#define ALLOC(array,item,type)                                         \
array=calloc(item,sizeof(type));                                       \
if(array==NULL){                                                       \
	printf("GTO_JK_Matrix_ShellSet - error cannot allocate memory\n"); \
	exit(-1);                                                          \
}

	// memory allocations
	ALLOC(PT,      nBasis*nBasis,     double);
	ALLOC(PTA,     nBasis*nBasis,     double);
	ALLOC(PAB,     nBasis*nBasis,     double);

#undef ALLOC

	// compute total spin density matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		PT[p*nBasis+q] = PA[p*nBasis+q] + PB[p*nBasis+q];

	// compute absolute value for screening
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++){
		PTA[p*nBasis+q] = fabs(PT[p*nBasis+q]);
		PAB[p*nBasis+q] = fabs(PA[p*nBasis+q]);
		if( PAB[p*nBasis+q]  < fabs(PB[p*nBasis+q]) )
			PAB[p*nBasis+q]  = fabs(PB[p*nBasis+q]);
	}

	// prepare shell data
	nShell = shell_prep(nBasis,gto,Schwarz,&shellMap,&shellMaxL,&Schwarz_Shell);

	// loop all distinct shell sequence
	SHELL_LOOP_BEGIN

		// screen integral using Schwarz inequality at shell level
		r = Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS];
		if(r < fixedCutoff)
			continue;

		// screen another layer with electron density weighted
		PQIJ_LOOP_BEGIN
			r = Schwarz[p*nBasis+q]*Schwarz[i*nBasis+j];
			if(r > fixedCutoff){
				if(r*PTA[p*nBasis+q] > fixedCutoff) goto computeUniqueShell;
				if(r*PTA[i*nBasis+j] > fixedCutoff) goto computeUniqueShell;
				if(r*PAB[p*nBasis+i] > fixedCutoff) goto computeUniqueShell;
				if(r*PAB[p*nBasis+j] > fixedCutoff) goto computeUniqueShell;
				if(r*PAB[q*nBasis+i] > fixedCutoff) goto computeUniqueShell;
				if(r*PAB[q*nBasis+j] > fixedCutoff) goto computeUniqueShell;
			}
		PQIJ_LOOP_END

		continue;

computeUniqueShell:

		// index preparation
		nEE = 0;
		PQIJ_LOOP_BEGIN

			// precompute index for efficientcy
			kZ =    (shellMaxL[jS]+1);
			kY = kZ*(shellMaxL[iS]+1);
			kX = kY*(shellMaxL[qS]+1);
			iBx[nEE] = 4*(MAXL+1)*(gto[p].l*kX + gto[q].l*kY + gto[i].l*kZ + gto[j].l);
			iBy[nEE] = 4*(MAXL+1)*(gto[p].m*kX + gto[q].m*kY + gto[i].m*kZ + gto[j].m);
			iBz[nEE] = 4*(MAXL+1)*(gto[p].n*kX + gto[q].n*kY + gto[i].n*kZ + gto[j].n);
			mX[nEE]  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
			mY[nEE]  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
			mZ[nEE]  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;

			// reset EEStore to zero
			EEStore[nEE] = 0.0;

			// increase number of EE index
			nEE++;

		PQIJ_LOOP_END

		// check that the number of integral does not exceed maximum 
		if(nEE > MAXSHELLINT){
			printf("GTO_JK_Matrix_ShellSet - error too many integrals\n");
			exit(-1);
		}

		///////////////////////////////////////////////
		// generate (ab|cd)
		//
		// 1) loop all contracted
		// 2) generate Bx,By,Bz for each contracted
		// 3) summation to get (ab|cd)
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
			r=genSetBxyzF(
			        gto[p].x0,gto[p].y0,gto[p].z0,shellMaxL[pS],gto[p].exp[pC],
			        gto[q].x0,gto[q].y0,gto[q].z0,shellMaxL[qS],gto[q].exp[qC],
			        gto[i].x0,gto[i].y0,gto[i].z0,shellMaxL[iS],gto[i].exp[iC],
			        gto[j].x0,gto[j].y0,gto[j].z0,shellMaxL[jS],gto[j].exp[jC],
			        Bx,By,Bz,F,MAXL);

			if(r < PRIMITIVE_CUTOFF) continue;
			
			// loop all basis within shell
			nEE = 0;
			PQIJ_LOOP_BEGIN
				
				// compute two-electron integral using summation
				tBx = Bx + iBx[nEE];
				tBy = By + iBy[nEE];
				tBz = Bz + iBz[nEE];

				EE=0.0;
				for(kX=mX[nEE];kX>=0;kX--)
				for(kY=mY[nEE];kY>=0;kY--)
				for(kZ=mZ[nEE];kZ>=0;kZ--)
					EE += tBx[kX]*tBy[kY]*tBz[kZ]*F[kX+kY+kZ];
				
				// multiply contracted coefficient
				EEStore[nEE] += EE * r
				                   * gto[p].coef[pC]*gto[p].norm[pC]
				                   * gto[q].coef[qC]*gto[q].norm[qC]
				                   * gto[i].coef[iC]*gto[i].norm[iC]
				                   * gto[j].coef[jC]*gto[j].norm[jC];

				nEE++;
			
			PQIJ_LOOP_END

			// reset shell index to proper value in this shell
			p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		}

		// determine pS qS iS jS type for storage
		if((pS==qS)&&(iS==jS)&&(pS==iS)){ pqijT = PQIJ_ALLSAME;
		}else if((pS==qS)&&(iS==jS)){     pqijT = PQIJ_2PAIRS;
		}else if(pS==qS){                 pqijT = PQIJ_PQPAIR;
		}else if(iS==jS){                 pqijT = PQIJ_IJPAIR;
		}else if((pS==iS)&&(qS==jS)){     pqijT = PQIJ_PIQJPAIR;
		}else{                            pqijT = PQIJ_ALLDISTINCT;
		}

		// assemble into G matrix
		nEE = 0;
		PQIJ_LOOP_BEGIN
			EE = EEStore[nEE];
			if(fabs(EE) > fixedCutoff){
				storeJK(GA,PT,PA,nBasis,EE,p,q,i,j,pqijT);
				storeJK(GB,PT,PB,nBasis,EE,p,q,i,j,pqijT);
			}
			nEE++;
		PQIJ_LOOP_END

	SHELL_LOOP_END

	// symmetrize G matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < p; q++){
		GA[q*nBasis+p] = GA[p*nBasis+q];
		GB[q*nBasis+p] = GB[p*nBasis+q];
	}

	// clean memory
	free(PT);
	free(PTA);
	free(PAB);
	free(shellMap);
	free(shellMaxL);
	free(Schwarz_Shell);
}
*/

//
// mp2_gen_aqbj_ShellSet : generate 2e integral (aq|bj) for MP2 calculation.
// See the more simpler version "mp2_gen_aqbj_NoSymm to get the basic 
// understanding of what it does.
//
// Oct 20, 2012 - Teepanis Chachiyo
//    Copied from GTO_JK_Matrix_ShellSet given as comments above.
//    Initial implementation and testing
//
double * mp2_gen_aqbj_ShellSet(
	int nBasis,                     // number of basis function
	int nOcc,                       // number of occpupied orbitals
	int nCore,                      // number of core orbitals to ignore
	const double *C,                // molecular orbitals
	const struct GTOBasis_t *gto){  // basis function structure

	// mp2 aqbj generator variables
	int nCorr;                     // number of correlated orbitals to include
	double *aqbjBank;              // pointer to integral storage
	int a,b;                       // orbital index
	double *ptr;                   // pointer to current value
	double *CT;                    // transposed version of molecular orbitals
	double **QJm;                  // QJ index mapping

	// shell set specific variables
	int pS,qS,iS,jS;               // looping variables for shells
	int p,q,i,j;                   // looping variables for basis functions
	int pC,qC,iC,jC;               // contracted index
	double EE;                     // coulomb integral
	int nEE;                       // number of electron within shell counter
	double EEStore[MAXSHELLINT];   // electron integral storage
	int        iBx[MAXSHELLINT];   // index for Bx storage
	int        iBy[MAXSHELLINT];   // index for By storage
	int        iBz[MAXSHELLINT];   // index for Bz storage
	char        mX[MAXSHELLINT];   // maximum in x direction
	char        mY[MAXSHELLINT];   // maximum in y direction
	char        mZ[MAXSHELLINT];   // maximum in z direction
	double Bx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double By[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Bz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double  F[4*(MAXL+1)];         // Boys functions
	double *tBx,*tBy,*tBz;         // pointer to Bx, By, Bz
	int kX, kY, kZ;                // summation loop
	char pqijT;                    // storage type for a set of pqij
	int nShell;                    // number of shell
	int *shellMap;                 // mapping between shell and basis function
	int *shellMaxL;                // maximum angular moment for each shell
	double *Schwarz_Shell;         // schwarz screening at shell level
	double *Schwarz;               // schwarz screening at basis level
	int n;                         // generic index
	double r;                      // generic real number

	// compute number of correlated orbitals
	nCorr = nOcc - nCore;

	// memory allocation
	if((aqbjBank=calloc((nCorr*(nCorr+1))/2*nBasis*nBasis,sizeof(double)))==NULL){
		printf("mp2_gen_aqbj_ShellSet - error cannot allocate memory\n");
		exit(-1); 
	}

	// allocate and compute QJ mapping matrix
	if((QJm=calloc(nBasis*nBasis,sizeof(double*)))==NULL){
		printf("mp2_gen_aqbj_ShellSet - error cannot allocate memory\n");
		exit(-1); 
	}
	for(q=0; q < nBasis; q++)
	for(j=0; j < nBasis; j++)
		QJm[q*nBasis+j] = aqbjBank + (q*nBasis+j)*nCorr*(nCorr+1)/2;

	// allocate and compute transpose matrix
	if((CT=calloc(nBasis*nBasis,sizeof(double)))==NULL){
		printf("mp2_gen_aqbj_ShellSet - error cannot allocate memory\n");
		exit(-1); 
	}
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		CT[p*nBasis+q] = C[q*nBasis+p];

	// prepare shell data
	Schwarz = create_Schwarz(nBasis, gto);
	nShell  = shell_prep(nBasis,gto,Schwarz,&shellMap,&shellMaxL,&Schwarz_Shell);

	// loop all distinct shell sequence
	SHELL_LOOP_BEGIN

		// screen integral using Schwarz inequality at shell level
		r = Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS];
		if(r < MP2_SCHWARZ_CUTOFF)
			continue;

		// index preparation
		nEE = 0;
		PQIJ_LOOP_BEGIN

			// precompute index for efficientcy
			kZ =    (shellMaxL[jS]+1);
			kY = kZ*(shellMaxL[iS]+1);
			kX = kY*(shellMaxL[qS]+1);
			iBx[nEE] = 4*(MAXL+1)*(gto[p].l*kX + gto[q].l*kY + gto[i].l*kZ + gto[j].l);
			iBy[nEE] = 4*(MAXL+1)*(gto[p].m*kX + gto[q].m*kY + gto[i].m*kZ + gto[j].m);
			iBz[nEE] = 4*(MAXL+1)*(gto[p].n*kX + gto[q].n*kY + gto[i].n*kZ + gto[j].n);
			mX[nEE]  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
			mY[nEE]  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
			mZ[nEE]  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;

			// reset EEStore to zero
			EEStore[nEE] = 0.0;

			// increase number of EE index
			nEE++;

		PQIJ_LOOP_END

		// check that the number of integral does not exceed maximum 
		if(nEE > MAXSHELLINT){
			printf("mp2_gen_aqbj_ShellSet - error too many integrals\n");
			exit(-1);
		}

		///////////////////////////////////////////////
		// generate (ab|cd)
		//
		// 1) loop all contracted
		// 2) generate Bx,By,Bz for each contracted
		// 3) summation to get (ab|cd)
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
			r=genSetBxyzF(
			        gto[p].x0,gto[p].y0,gto[p].z0,shellMaxL[pS],gto[p].exp[pC],
			        gto[q].x0,gto[q].y0,gto[q].z0,shellMaxL[qS],gto[q].exp[qC],
			        gto[i].x0,gto[i].y0,gto[i].z0,shellMaxL[iS],gto[i].exp[iC],
			        gto[j].x0,gto[j].y0,gto[j].z0,shellMaxL[jS],gto[j].exp[jC],
			        Bx,By,Bz,F,MAXL);

			if(r < PRIMITIVE_CUTOFF) continue;
			
			// loop all basis within shell
			nEE = 0;
			PQIJ_LOOP_BEGIN
				
				// compute two-electron integral using summation
				tBx = Bx + iBx[nEE];
				tBy = By + iBy[nEE];
				tBz = Bz + iBz[nEE];

				EE=0.0;
				for(kX=mX[nEE];kX>=0;kX--)
				for(kY=mY[nEE];kY>=0;kY--)
				for(kZ=mZ[nEE];kZ>=0;kZ--)
					EE += tBx[kX]*tBy[kY]*tBz[kZ]*F[kX+kY+kZ];
				
				// multiply contracted coefficient
				EEStore[nEE] += EE * r
				                   * gto[p].coef[pC]*gto[p].norm[pC]
				                   * gto[q].coef[qC]*gto[q].norm[qC]
				                   * gto[i].coef[iC]*gto[i].norm[iC]
				                   * gto[j].coef[jC]*gto[j].norm[jC];

				nEE++;
			
			PQIJ_LOOP_END

			// reset shell index to proper value in this shell
			p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		}

		// determine pS qS iS jS type for storage
		if((pS==qS)&&(iS==jS)&&(pS==iS)){ pqijT = PQIJ_ALLSAME;
		}else if((pS==qS)&&(iS==jS)){     pqijT = PQIJ_2PAIRS;
		}else if(pS==qS){                 pqijT = PQIJ_PQPAIR;
		}else if(iS==jS){                 pqijT = PQIJ_IJPAIR;
		}else if((pS==iS)&&(qS==jS)){     pqijT = PQIJ_PIQJPAIR;
		}else{                            pqijT = PQIJ_ALLDISTINCT;
		}

		// assemble into aqbj array
		nEE = 0;
		PQIJ_LOOP_BEGIN
			EE = EEStore[nEE];

#define aqbj(P,Q,I,J)                                            \
ptr = QJm[Q*nBasis+J];                                           \
for(a=0; a < nCorr; a++){                                        \
	r = CT[P*nBasis+nCore+a]*EE;                                 \
	for(b=0; b <= a; b++){                                       \
		*ptr += r*CT[I*nBasis+nCore+b];                          \
		ptr++;                                                   \
	}                                                            \
}

			if(fabs(EE) > MP2_SCHWARZ_CUTOFF)
			switch(pqijT){

			case PQIJ_ALLSAME:
				aqbj(p,q,i,j);
			break;

			case PQIJ_2PAIRS:
				aqbj(p,q,i,j);      aqbj(i,j,p,q);
			break;

			case PQIJ_PQPAIR:
				aqbj(p,q,i,j);      aqbj(i,j,p,q);
				aqbj(p,q,j,i);      aqbj(j,i,p,q);
			break;

			case PQIJ_IJPAIR:
				aqbj(i,j,p,q);      aqbj(p,q,i,j);
				aqbj(i,j,q,p);      aqbj(q,p,i,j);
			break;

			case PQIJ_PIQJPAIR:
				aqbj(p,q,i,j);      //aqbj(i,j,p,q);
				aqbj(p,q,j,i);      //aqbj(i,j,q,p);
				aqbj(q,p,i,j);      //aqbj(j,i,p,q);
				aqbj(q,p,j,i);      //aqbj(j,i,q,p);
			break;

			case PQIJ_ALLDISTINCT:
				aqbj(p,q,i,j);      aqbj(i,j,p,q);
				aqbj(p,q,j,i);      aqbj(i,j,q,p);
				aqbj(q,p,i,j);      aqbj(j,i,p,q);
				aqbj(q,p,j,i);      aqbj(j,i,q,p);
			break;
			}

			nEE++;
		PQIJ_LOOP_END

	SHELL_LOOP_END

	// clean memory
	free(shellMap);
	free(shellMaxL);
	free(Schwarz_Shell);
	free(Schwarz);
	free(CT);
	free(QJm);

	return aqbjBank;
}


//
// mp2_gen_aqij_ShellSet : generate 2e integral (aq|ij) for MP2 calculation.
//
// Nov 5, 2012 - Teepanis Chachiyo
//    Initial implementation and testing. This code is very similar to
//    mp2_gen_aqbj_ShellSet.
//
// Nov 5, 2012 - Teepanis Chachiyo
//    Use the following heirarchy data in the array i->j->a->q. This will
//    slowdown this subroutine a bit, but this subroutine is not the bottleneck
//    so we are not losing much speed.
//
// Nov 7, 2012 - Teepanis Chachiyo
//    Compute contraction coefficient first and stored in cntCoef just like in
//    the subroutine GTO_JK_Matrix_ShellSet
//
// Nov 11, 2012 - Teepanis Chachiyo
//    Use summation tBz[kZ]*F[kX+kY+kZ] --> tSz instead
//
double * mp2_gen_aqij_ShellSet(
	int nBasis,                     // number of basis function
	int nCorr,                      // number of correlated orbitals
	int nCore,                      // number of core orbitals to ignore
	const double *C,                // molecular orbitals
	const struct GTOBasis_t *gto){  // basis function structure

	// mp2 aqij generator variables
	double *aqijBank;              // pointer to integral storage
	int a;                         // orbital index
	double *ptr;                   // pointer to current value
	double *CT;                    // transposed version of molecular orbitals
	double **IJm;                  // ij index mapping

	// shell set specific variables
	int pS,qS,iS,jS;               // looping variables for shells
	int p,q,i,j;                   // looping variables for basis functions
	int pC,qC,iC,jC;               // contracted index
	double EE;                     // coulomb integral
	int nEE;                       // number of electron within shell counter
	double EEStore[MAXSHELLINT];   // electron integral storage
	double  tSum;                  // intermediate variables for kX,kY,kZ loop
	int maxCnt;                    // maximum number of contraction
	double *cntCoef;               // contract function coefficients
	int        iBx[MAXSHELLINT];   // index for Bx storage
	int        iBy[MAXSHELLINT];   // index for By storage
	int        iBz[MAXSHELLINT];   // index for Bz storage
	char        mX[MAXSHELLINT];   // maximum in x direction
	char        mY[MAXSHELLINT];   // maximum in y direction
	char        mZ[MAXSHELLINT];   // maximum in z direction
	char        mT[MAXSHELLINT];   // type of mX,mY,mZ for loop optimization
	double Bx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double By[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Bz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sy[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double  F[4*(MAXL+1)];         // Boys functions
	double *tBx,*tBy,*tBz;         // pointer to Bx, By, Bz
	double *tSx,*tSy,*tSz;         // pointer to Sx, Sy, Sz
	int kX, kY, kZ;                // summation loop
	char pqijT;                    // storage type for a set of pqij
	int nShell;                    // number of shell
	int *shellMap;                 // mapping between shell and basis function
	int *shellMaxL;                // maximum angular moment for each shell
	double *Schwarz_Shell;         // schwarz screening at shell level
	double *Schwarz;               // schwarz screening at basis level
	int n;                         // generic index
	double r;                      // generic real number

	// prepare shell data
	Schwarz = create_Schwarz(nBasis, gto);
	nShell  = shell_prep(nBasis,gto,Schwarz,&shellMap,&shellMaxL,&Schwarz_Shell);

	// compute number of integral needed
	nEE = 0;
	for(i=0; i < nBasis; i++)
	for(j=0; j <= i; j++)
		if(Schwarz[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) continue; else nEE++;

	// contraction coefficients
	maxCnt = gto[0].nContract;
	for(p=0; p < nBasis; p++)
		if(gto[p].nContract > maxCnt) maxCnt = gto[p].nContract;
	if((cntCoef=calloc(maxCnt*nBasis,sizeof(double)))==NULL){
		printf("mp2_gen_aqij_ShellSet - error cannot allocate cntCoef\n");
		exit(-1);
	}
	for(p=0; p < nBasis; p++)
		for(pC=0; pC < gto[p].nContract; pC++)
			cntCoef[pC*nBasis+p] = gto[p].coef[pC]*gto[p].norm[pC];

	// memory allocation
	if((aqijBank=calloc(nEE*nBasis*nCorr,sizeof(double)))==NULL){
		printf("mp2_gen_aqij_ShellSet - error cannot allocate memory\n");
		exit(-1); 
	}

	// allocate and compute transpose matrix
	if(( CT=calloc(nBasis*nBasis,sizeof(double )))==NULL ||
	   (IJm=calloc(nBasis*nBasis,sizeof(double*)))==NULL){
		printf("mp2_gen_aqij_ShellSet - error cannot allocate memory\n");
		exit(-1); 
	}
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		CT[p*nBasis+q] = C[q*nBasis+p];

	// build IJm mapping matrix
	nEE=0;
	for(i=0; i < nBasis; i++)
	for(j=0; j <= i; j++){
		if(Schwarz[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) IJm[i*nBasis+j] = NULL;
		else{
			IJm[i*nBasis+j] = aqijBank + nEE*nBasis*nCorr;
			nEE++;
		}
		IJm[j*nBasis+i] = IJm[i*nBasis+j];
	}

	// loop all distinct shell sequence
	SHELL_LOOP_BEGIN

		// screen integral using Schwarz inequality at shell level
		r = Schwarz_Shell[pS*nShell+qS]*Schwarz_Shell[iS*nShell+jS];
		if(r < MP2_SCHWARZ_CUTOFF)
			continue;

		// index preparation
		nEE = 0;
		PQIJ_LOOP_BEGIN

			// precompute index for efficientcy
			kZ =    (shellMaxL[jS]+1);
			kY = kZ*(shellMaxL[iS]+1);
			kX = kY*(shellMaxL[qS]+1);
			iBx[nEE] = 4*(MAXL+1)*(gto[p].l*kX + gto[q].l*kY + gto[i].l*kZ + gto[j].l);
			iBy[nEE] = 4*(MAXL+1)*(gto[p].m*kX + gto[q].m*kY + gto[i].m*kZ + gto[j].m);
			iBz[nEE] = 4*(MAXL+1)*(gto[p].n*kX + gto[q].n*kY + gto[i].n*kZ + gto[j].n);
			mX[nEE]  = gto[p].l+gto[q].l+gto[i].l+gto[j].l;
			mY[nEE]  = gto[p].m+gto[q].m+gto[i].m+gto[j].m;
			mZ[nEE]  = gto[p].n+gto[q].n+gto[i].n+gto[j].n;

			// determine type of mX,mY,mZ for loop optimization
			if(mX[nEE]==0)                                  mT[nEE] = MXZERO;
			else if(mY[nEE]==0)                             mT[nEE] = MYZERO;
			else if(mZ[nEE]==0)                             mT[nEE] = MZZERO;
			else if(mX[nEE] > mY[nEE] && mX[nEE] > mZ[nEE]) mT[nEE] = MXMAX;
			else if(mY[nEE] > mZ[nEE])                      mT[nEE] = MYMAX;
			else                                            mT[nEE] = MZMAX;

			// reset EEStore to zero
			EEStore[nEE] = 0.0;

			// increase number of EE index
			nEE++;

		PQIJ_LOOP_END

		// check that the number of integral does not exceed maximum 
		if(nEE > MAXSHELLINT){
			printf("mp2_gen_aqij_ShellSet - error too many integrals\n");
			exit(-1);
		}

		///////////////////////////////////////////////
		// generate (ab|cd)
		//
		// 1) loop all contracted
		// 2) generate Bx,By,Bz for each contracted
		// 3) summation to get (ab|cd)
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
			r=genSetBxyzF(
			        gto[p].x0,gto[p].y0,gto[p].z0,shellMaxL[pS],gto[p].exp[pC],
			        gto[q].x0,gto[q].y0,gto[q].z0,shellMaxL[qS],gto[q].exp[qC],
			        gto[i].x0,gto[i].y0,gto[i].z0,shellMaxL[iS],gto[i].exp[iC],
			        gto[j].x0,gto[j].y0,gto[j].z0,shellMaxL[jS],gto[j].exp[jC],
			        Bx,By,Bz,F,MAXL);

			if(r < PRIMITIVE_CUTOFF) continue;

			// compute sum Bx*F, By*F, Bz*F
			tBx = Bx; tSx = Sx;
			tBy = By; tSy = Sy;
			tBz = Bz; tSz = Sz;
			int a,b,c,d;
			int maxL,maxShellL;
			int l,mn;
			maxShellL = shellMaxL[pS]+shellMaxL[qS]+shellMaxL[iS]+shellMaxL[jS];
			for(a=0; a < shellMaxL[pS]+1; a++)
			for(b=0; b < shellMaxL[qS]+1; b++)
			for(c=0; c < shellMaxL[iS]+1; c++)
			for(d=0; d < shellMaxL[jS]+1; d++){
				maxL = a+b+c+d;
				double xSum, ySum, zSum;
				for(mn=0; mn < (maxShellL-maxL+1); mn++){
					xSum = ySum = zSum = 0.0;
					for(l=0;l<(maxL+1);l++){
						xSum+=tBx[l]*F[mn+l];
						ySum+=tBy[l]*F[mn+l];
						zSum+=tBz[l]*F[mn+l];
					}
					tSx[mn] = xSum;
					tSy[mn] = ySum;
					tSz[mn] = zSum;
				}
				tBx += 4*(MAXL+1); tSx += 4*(MAXL+1);
				tBy += 4*(MAXL+1); tSy += 4*(MAXL+1);
				tBz += 4*(MAXL+1); tSz += 4*(MAXL+1);
			}
			
			// loop all basis within shell
			nEE = 0;
			PQIJ_LOOP_BEGIN
				
				// compute two-electron integral using summation
				tBx = Bx + iBx[nEE];
				tBy = By + iBy[nEE];
				tBz = Bz + iBz[nEE];

				tSx = Sx + iBx[nEE];
				tSy = Sy + iBy[nEE];
				tSz = Sz + iBz[nEE];

				EE=0.0;
				//for(kX=mX[nEE];kX>=0;kX--)
				//for(kY=mY[nEE];kY>=0;kY--)
				//for(kZ=mZ[nEE];kZ>=0;kZ--)
				//	EE += tBx[kX]*tBy[kY]*tBz[kZ]*F[kX+kY+kZ];

				switch(mT[nEE]){

				case MXZERO:
					for(kY=mY[nEE];kY>=0;kY--) EE += tBy[kY]*tSz[kY];
				break;

				case MYZERO:
					for(kZ=mZ[nEE];kZ>=0;kZ--) EE += tBz[kZ]*tSx[kZ];
				break;

				case MZZERO:
					for(kX=mX[nEE];kX>=0;kX--) EE += tBx[kX]*tSy[kX];
				break;

				case MXMAX: 
					for(kY=mY[nEE];kY>=0;kY--){
						for(tSum=0.0,kZ=mZ[nEE];kZ>=0;kZ--) tSum+=tBz[kZ]*tSx[kY+kZ];
						EE += tBy[kY]*tSum;
					}
				break;

				case MYMAX:
					for(kZ=mZ[nEE];kZ>=0;kZ--){
						for(tSum=0.0,kX=mX[nEE];kX>=0;kX--) tSum+=tBx[kX]*tSy[kX+kZ];
						EE += tBz[kZ]*tSum;
					}
				break;

				case MZMAX:
					for(kX=mX[nEE];kX>=0;kX--){
						for(tSum=0.0,kY=mY[nEE];kY>=0;kY--) tSum+=tBy[kY]*tSz[kX+kY];
						EE += tBx[kX]*tSum;
					}
				break;
				}
				
				// multiply contracted coefficient
				//EEStore[nEE] += EE * r
				//                   * gto[p].coef[pC]*gto[p].norm[pC]
				//                   * gto[q].coef[qC]*gto[q].norm[qC]
				//                   * gto[i].coef[iC]*gto[i].norm[iC]
				//                   * gto[j].coef[jC]*gto[j].norm[jC];
				EEStore[nEE] += EE * r
				                   * cntCoef[pC*nBasis+p]
				                   * cntCoef[qC*nBasis+q]
				                   * cntCoef[iC*nBasis+i]
				                   * cntCoef[jC*nBasis+j];

				nEE++;
			
			PQIJ_LOOP_END

			// reset shell index to proper value in this shell
			p=shellMap[pS]; q=shellMap[qS]; i=shellMap[iS]; j=shellMap[jS];
		}

		// determine pS qS iS jS type for storage
		if((pS==qS)&&(iS==jS)&&(pS==iS)){ pqijT = PQIJ_ALLSAME;
		}else if((pS==qS)&&(iS==jS)){     pqijT = PQIJ_2PAIRS;
		}else if(pS==qS){                 pqijT = PQIJ_PQPAIR;
		}else if(iS==jS){                 pqijT = PQIJ_IJPAIR;
		}else if((pS==iS)&&(qS==jS)){     pqijT = PQIJ_PIQJPAIR;
		}else{                            pqijT = PQIJ_ALLDISTINCT;
		}

		// assemble into aqbj array
		nEE = 0;
		PQIJ_LOOP_BEGIN
			EE = EEStore[nEE];

#define aqij(P,Q,I,J)                                        \
if(J<=I && (ptr=IJm[I*nBasis+J])){                           \
	for(a=0; a < nCorr; a++){                                \
		ptr[a*nBasis+Q] += CT[P*nBasis+nCore+a]*EE;          \
	}                                                        \
}

			if(fabs(EE) > MP2_SCHWARZ_CUTOFF)
			switch(pqijT){

			case PQIJ_ALLSAME:
				aqij(p,q,i,j);
			break;

			case PQIJ_2PAIRS:
				aqij(p,q,i,j);      aqij(i,j,p,q);
			break;

			case PQIJ_PQPAIR:
				aqij(p,q,i,j);      aqij(i,j,p,q);
				aqij(p,q,j,i);      aqij(j,i,p,q);
			break;

			case PQIJ_IJPAIR:
				aqij(i,j,p,q);      aqij(p,q,i,j);
				aqij(i,j,q,p);      aqij(q,p,i,j);
			break;

			case PQIJ_PIQJPAIR:
				aqij(p,q,i,j);      //aqij(i,j,p,q);
				aqij(p,q,j,i);      //aqij(i,j,q,p);
				aqij(q,p,i,j);      //aqij(j,i,p,q);
				aqij(q,p,j,i);      //aqij(j,i,q,p);
			break;

			case PQIJ_ALLDISTINCT:
				aqij(p,q,i,j);      aqij(i,j,p,q);
				aqij(p,q,j,i);      aqij(i,j,q,p);
				aqij(q,p,i,j);      aqij(j,i,p,q);
				aqij(q,p,j,i);      aqij(j,i,q,p);
			break;
			}

			nEE++;
		PQIJ_LOOP_END

	SHELL_LOOP_END

	// clean memory
	free(cntCoef);
	free(shellMap);
	free(shellMaxL);
	free(Schwarz_Shell);
	free(Schwarz);
	free(CT);
	free(IJm);

	return aqijBank;
}


//
// mp2_gen_aqij_Quartet : generate 2e integral (aq|ij) for MP2 calculation.
//
// Feb 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
double * mp2_gen_aqij_Quartet(
	int nBasis,                     // number of basis function
	int nCorr,                      // number of correlated orbitals
	int nCore,                      // number of core orbitals to ignore
	const double *C,                // molecular orbitals
	const struct GTOBasis_t *gto,   // basis function structure
	const double *schwarz_basis){   // schwarz matrix at basis level

	// mp2 aqij generator variables
	double *aqijBank;              // pointer to integral storage
	int a;                         // orbital index
	double *ptr;                   // pointer to current value
	double *CT;                    // transposed version of molecular orbitals
	double **IJm;                  // ij index mapping

	// quartet set specific variables
	static int firstVisit=1;              // first visit flag
	static struct GTOShell_t *shell=NULL; // shell data structure
	static int nShell;                    // the number of shell
	int P,Q,I,J,N;                        // shell loop index
	int p,q,i,j;                          // basis loop index
	static double *schwarz_shell=NULL;    // schwarz matrix at shell level
	double EE;                            // generic double precision
	int nEE;                              // generic int
	double EEStore[MAXSHELLINT];          // electron integral storage
	double *EEStorePtr;                   // pointer to the calculated values
	char pqijT;                           // type of shell quartet

	static struct Multipole_t *mp=NULL;   // multipole data
	static signed char *erBank=NULL;      // multipole error bank at shell level
	signed char *erBankPtr=NULL;          // pointer to current value
	signed char erExp;                    // current value of error

	// reset call
	if(nBasis==0){
		firstVisit = 1;
		if(shell != NULL){ free(shell); shell=NULL; }
		if(schwarz_shell != NULL){ free(schwarz_shell); schwarz_shell=NULL; }
		if(mp != NULL){ free(mp); mp=NULL; }
		if(erBank != NULL){ free(erBank); erBank=NULL; }
		return NULL;
	}

	// first visit preparations
	if(firstVisit){

		// generate shell
		shell = genGTOShell(nBasis, gto, &nShell);
		sortGTOShell(nShell, shell);

		// generate multipole data
		mp = genMultipole(nBasis, gto);

		// allocate memory
		if((schwarz_shell=calloc(nShell*nShell, sizeof(double)))==NULL){
			printf("mp2_gen_aqij_Quartet - error cannot allocate memory\n");
			exit(-1); 
		}

		// compute schwarz at shell level
		for(P=0; P < nShell; P++)
		for(Q=0; Q < nShell; Q++){
			EE = 0.0;
			for(p=shell[P].min; p <= shell[P].max; p++)
			for(q=shell[Q].min; q <= shell[Q].max; q++)
				if(EE <schwarz_basis[p*nBasis+q]) EE = schwarz_basis[p*nBasis+q];
			schwarz_shell[P*nShell+Q] = EE;
		}

		// count the number shells needed and allocate memory
		unsigned long nEE=0;
		for(P=0; P < nShell; P++)
		for(Q=0; Q <= P; Q++)
		for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
		for(J=0; J <= N; J++){
			EE = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
			if(EE < MP2_SCHWARZ_CUTOFF) continue;
			nEE++;
		}
		}
		erBank=calloc(nEE,sizeof(signed char));
		if(erBank==NULL){
			printf("mp2_gen_aqij_Quartet - error cannot allocate memory\n");
			exit(-1); 
		}
		memset(erBank, ERROR_UNKNOWN, nEE);
	}

	// compute number of integral needed
	nEE = 0;
	for(i=0; i < nBasis; i++)
	for(j=0; j <= i; j++)
		if(schwarz_basis[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) continue; else nEE++;

	// memory allocation
	if((aqijBank=calloc(nEE*nBasis*nCorr,sizeof(double)))==NULL){
		printf("mp2_gen_aqij_Quartet - error cannot allocate memory\n");
		exit(-1); 
	}

	// allocate and compute transpose matrix
	if(( CT=calloc(nBasis*nBasis,sizeof(double )))==NULL ||
	   (IJm=calloc(nBasis*nBasis,sizeof(double*)))==NULL){
		printf("mp2_gen_aqij_Quartet - error cannot allocate memory\n");
		exit(-1); 
	}
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		CT[p*nBasis+q] = C[q*nBasis+p];

	// build IJm mapping matrix
	nEE=0;
	for(i=0; i < nBasis; i++)
	for(j=0; j <= i; j++){
		if(schwarz_basis[i*nBasis+j] < MP2_SCHWARZ_CUTOFF) IJm[i*nBasis+j] = NULL;
		else{
			IJm[i*nBasis+j] = aqijBank + nEE*nBasis*nCorr;
			nEE++;
		}
		IJm[j*nBasis+i] = IJm[i*nBasis+j];
	}

	erBankPtr = erBank;
	for(P=0; P < nShell; P++)
	for(Q=0; Q <= P; Q++)
	for(I=0; I <= P; I++){ if(I==P) N=Q; else N=I; 
	for(J=0; J <= N; J++){

		// schwarz screening
		EE = schwarz_shell[P*nShell+Q]*schwarz_shell[I*nShell+J];
		if(EE < MP2_SCHWARZ_CUTOFF) continue;
		erExp = *erBankPtr; erBankPtr++;

		// compute the entire shell quartet
		if(erExp!=ERROR_UNKNOWN && pow(ERROR_BASE,erExp) < MP2_SCHWARZ_CUTOFF)
			computeQuartetMPole(shell+P,shell+Q,shell+I,shell+J,mp,nBasis,EEStore);
		else
			computeQuartetEE(shell+P,shell+Q,shell+I,shell+J,EEStore);

		// store multipole error
		if(erExp==ERROR_UNKNOWN)
			erBankPtr[-1]=getErrExpQuartetMPole(shell+P,shell+Q,shell+I,shell+J,
			                                    mp,nBasis,EEStore);

		// determine P Q I J type for storage
		if((P==Q)&&(I==J)&&(P==I)){ pqijT = PQIJ_ALLSAME;
		}else if((P==Q)&&(I==J)){   pqijT = PQIJ_2PAIRS;
		}else if(P==Q){             pqijT = PQIJ_PQPAIR;
		}else if(I==J){             pqijT = PQIJ_IJPAIR;
		}else if((P==I)&&(Q==J)){   pqijT = PQIJ_PIQJPAIR;
		}else{                      pqijT = PQIJ_ALLDISTINCT;
		}

		// store in the G matrix
		EEStorePtr = EEStore;
		for(p=shell[P].min; p <= shell[P].max; p++)
		for(q=shell[Q].min; q <= shell[Q].max; q++)
		for(i=shell[I].min; i <= shell[I].max; i++)
		for(j=shell[J].min; j <= shell[J].max; j++){
			EE = *EEStorePtr; EEStorePtr++;

			if(fabs(EE) > MP2_SCHWARZ_CUTOFF)
			switch(pqijT){

			case PQIJ_ALLSAME:
				aqij(p,q,i,j);
			break;

			case PQIJ_2PAIRS:
				aqij(p,q,i,j);      aqij(i,j,p,q);
			break;

			case PQIJ_PQPAIR:
				aqij(p,q,i,j);      aqij(i,j,p,q);
				aqij(p,q,j,i);      aqij(j,i,p,q);
			break;

			case PQIJ_IJPAIR:
				aqij(i,j,p,q);      aqij(p,q,i,j);
				aqij(i,j,q,p);      aqij(q,p,i,j);
			break;

			case PQIJ_PIQJPAIR:
				aqij(p,q,i,j);      //aqij(i,j,p,q);
				aqij(p,q,j,i);      //aqij(i,j,q,p);
				aqij(q,p,i,j);      //aqij(j,i,p,q);
				aqij(q,p,j,i);      //aqij(j,i,q,p);
			break;

			case PQIJ_ALLDISTINCT:
				aqij(p,q,i,j);      aqij(i,j,p,q);
				aqij(p,q,j,i);      aqij(i,j,q,p);
				aqij(q,p,i,j);      aqij(j,i,p,q);
				aqij(q,p,j,i);      aqij(j,i,q,p);
			break;
			}
		}

	}
	}

	// clean memory
	free(CT);
	free(IJm);

	// reset flag
	if(firstVisit) firstVisit=0;

	return aqijBank;
}
#undef MAXSHELLINT 
#undef MAXL         

//
// GTO_JK_Matrix_NoSymm : calculate G matrix for unrestricted hartree-fock
// calculations. The equations are given in (Szabo and Ostlund, 1989)
// This version is slow but because it is very simple it is suitable
// for testing.
//
// for Alpha Spin
// [GA]pq = SUM [PT]ij*(pq|ij) - [PA]ij*(pi|qj) 
//
// for Beta Spin
// [GB]pg = SUM [PT]ij*(pq|ij) - [PB]ij*(pi|qj)
//
// where PT is the toal spin density and is defined by
// PT = PA+PB
//
// Sep 1, 2010 - Teepanis Chachiyo
//    Initial implementation and testing
//
// May 13, 2011 - Teepanis Chachiyo
//    Put the density matrix multipler back in
//
void GTO_JK_Matrix_NoSymm(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double cutoff,                 // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB){                   // return G for spin down

	double *PT;                    // total spin density matrix
	int p,q,i,j;                   // looping variables for basis functions
	double EE;                     // coulomb integral
	double upBound;                // Schwarz inequaliity upper bound

	// allocating memory
	PT        = calloc(nBasis*nBasis, sizeof(double));
	if(PT==NULL){
		printf("GTO_JK_Matrix_NoSymm - error allocating memory");
		exit(-1);
	}

	// compute total spin density matrix
	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++)
		PT[p*nBasis+q] = PA[p*nBasis+q] + PB[p*nBasis+q];

	for(p=0; p < nBasis; p++)
	for(q=0; q < nBasis; q++) 
	for(i=0; i < nBasis; i++)
	for(j=0; j < nBasis; j++){

		// screen integral using Schwarz inequality
		upBound = Schwarz[i*nBasis+j]*Schwarz[p*nBasis+q];
		if(upBound < cutoff) continue;

		// I have tested and re-tested. Multiplying density matrix
		// to the Schwarz inequality prescreening does not reduce
		// scaling properties. With Scharz alone, it is already N^2.
		// Multiplying density matrix will only make it tad bit faster
		// not worth the trouble.
		// Teepanis Chachiyo Sep 2, 2010
		//
		// Now that I use density difference, the density matrix
		// multiplier is important
		// Teepanis Chachiyo May 13, 2011
		//
		if(fabs(upBound*PT[i*nBasis+j]) > cutoff) goto computeEENoSym;
		if(fabs(upBound*PA[q*nBasis+j]) > cutoff) goto computeEENoSym;
		if(fabs(upBound*PB[q*nBasis+j]) > cutoff) goto computeEENoSym;

		continue;

computeEENoSym:

		// compute two-electron integral
		EE =
		contr_eri(
			     gto[p].nContract,
			     gto[p].exp, gto[p].coef, gto[p].norm,
			     gto[p].x0,        gto[p].y0,  gto[p].z0,
			     gto[p].l,         gto[p].m,   gto[p].n,
			     gto[q].nContract,
			     gto[q].exp, gto[q].coef, gto[q].norm,
			     gto[q].x0,        gto[q].y0,  gto[q].z0,
			     gto[q].l,         gto[q].m,   gto[q].n,
			     gto[i].nContract,
			     gto[i].exp, gto[i].coef, gto[i].norm,
			     gto[i].x0,        gto[i].y0,  gto[i].z0,
			     gto[i].l,         gto[i].m,   gto[i].n,
			     gto[j].nContract,
			     gto[j].exp, gto[j].coef, gto[j].norm,
			     gto[j].x0,        gto[j].y0,  gto[j].z0,
			     gto[j].l,         gto[j].m,   gto[j].n);

		GA[p*nBasis+q] = GA[p*nBasis+q] + PT[i*nBasis+j]*EE;
		GA[p*nBasis+i] = GA[p*nBasis+i] - PA[q*nBasis+j]*EE;
	
		GB[p*nBasis+q] = GB[p*nBasis+q] + PT[i*nBasis+j]*EE;
		GB[p*nBasis+i] = GB[p*nBasis+i] - PB[q*nBasis+j]*EE;

	}

	// free memory
	free(PT);
}


struct GTOPrime_t{
	double x,y,z;        // primitive center
	double exp;          // exponent
	int l,m,n;           // angular index
};

// GTO_JK_Matrix_PrimeSpace: generate JK matrix using primitive space
//
// Jan 19, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void GTO_JK_Matrix_PrimeSpace(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double cutoff,                 // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB){                   // return G for spin down

	int s,t,u,v,n;                 // looping variables for primitive functions
	double EE;                     // coulomb integral

	int nPrime;                    // number of primitive
	struct GTOPrime_t *prime;      // primitive data structure
	int *p2b;                      // primitive --> basis mapping
	double *Cfp;                   // coefficient for each primitive
	double *GAp,*GBp;              // G in primitive space
	double *PTp,*PAp,*PBp;         // density matrix in primitive space
	double *Fcp;                   // ERI prefactor in primitive space

	/*
 	 * manage primitive space
	 */

	// compute number of primitive function
	for(nPrime=0, u=0; u < nBasis; u++) nPrime += gto[u].nContract;

	// allocate memory
	GAp   = calloc(nPrime*nPrime, sizeof(double));
	GBp   = calloc(nPrime*nPrime, sizeof(double));
	PTp   = calloc(nPrime*nPrime, sizeof(double));
	PAp   = calloc(nPrime*nPrime, sizeof(double));
	PBp   = calloc(nPrime*nPrime, sizeof(double));
	Fcp   = calloc(nPrime*nPrime, sizeof(double));
	Cfp   = calloc(nPrime,        sizeof(double));
	p2b   = calloc(nPrime,        sizeof(int));
	prime = calloc(nPrime,        sizeof(struct GTOPrime_t)); 
	if(GAp==NULL || GBp==NULL || PTp==NULL || PAp==NULL || PBp==NULL ||
	   Fcp==NULL || p2b==NULL || Cfp==NULL || prime==NULL){
		printf("GTO_JK_Matrix_PrimeSpace - error allocating memory");
		exit(-1);
	}

	// generate primitive to basis mapping and coefficients values
	for(s=0, u=0; u < nBasis; u++)
		for(v=0; v < gto[u].nContract; v++){

			// generate primitive to basis mapping
			p2b[s] = u; 

			// generate coefficient value
			Cfp[s] = gto[u].norm[v] * gto[u].coef[v];

			// fill in primitive data 
			prime[s].x   = gto[u].x0;
			prime[s].y   = gto[u].y0;
			prime[s].z   = gto[u].z0;
			prime[s].l   = gto[u].l;
			prime[s].m   = gto[u].m;
			prime[s].n   = gto[u].n;
			prime[s].exp = gto[u].exp[v];

			s++; 
		}

	for(s=0; s < nPrime; s++)
	for(t=0; t <= s; t++){

		// generate density matrix in primitive space
		PAp[s*nPrime+t] = PA[ p2b[s]*nBasis + p2b[t] ] * Cfp[s] * Cfp[t];
		PBp[s*nPrime+t] = PB[ p2b[s]*nBasis + p2b[t] ] * Cfp[s] * Cfp[t];
		PTp[s*nPrime+t] = PAp[s*nPrime+t]   + PBp[s*nPrime+t];

		// generate ERI prefactor
#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+\
                                  (y1-y2)*(y1-y2)+\
                                  (z1-z2)*(z1-z2))
#define SQRT_TWOPI_5HALF 5.9149671728
		Fcp[s*nPrime+t] = SQRT_TWOPI_5HALF *
		                  exp( -  prime[s].exp * prime[t].exp 
		                       / (prime[s].exp + prime[t].exp)
		                       *  DIST2(prime[s].x,prime[s].y,prime[s].z,
		                                prime[t].x,prime[t].y,prime[t].z) )
		                  /(prime[s].exp + prime[t].exp);
#undef DIST2

		PTp[t*nPrime+s] = PTp[s*nPrime+t];
		PAp[t*nPrime+s] = PAp[s*nPrime+t];
		PBp[t*nPrime+s] = PBp[s*nPrime+t];
		Fcp[t*nPrime+s] = Fcp[s*nPrime+t];
	}

	/*
	 * main loop: G matrix alignment
	 */

	/*
	for(s=0; s < nPrime; s++)
	for(t=0; t <= s; t++)
	for(u=0; u < nPrime; u++)
	for(v=0; v < nPrime; v++){

		// coulomb part
		Pf = PTp[u*nPrime+v]*Fcp[s*nPrime+t]*Fcp[u*nPrime+v];
		if(fabs(Pf) > cutoff){
			nA++;
			EE = primeERI(
				 prime[s].x,prime[s].y,prime[s].z,prime[s].l,prime[s].m,prime[s].n,prime[s].exp,
			     prime[t].x,prime[t].y,prime[t].z,prime[t].l,prime[t].m,prime[t].n,prime[t].exp,
			     prime[u].x,prime[u].y,prime[u].z,prime[u].l,prime[u].m,prime[u].n,prime[u].exp,
			     prime[v].x,prime[v].y,prime[v].z,prime[v].l,prime[v].m,prime[v].n,prime[v].exp);
			GAp[s*nPrime+t] += Pf*EE;
			GBp[s*nPrime+t] += Pf*EE;
		}

		// exchange part
		PAf = PAp[u*nPrime+v]*Fcp[s*nPrime+u]*Fcp[t*nPrime+v];
		PBf = PBp[u*nPrime+v]*Fcp[s*nPrime+u]*Fcp[t*nPrime+v];
		if(fabs(PAf) > cutoff || fabs(PBf) > cutoff){
			nB++;
			EE = primeERI(
			     prime[s].x,prime[s].y,prime[s].z,prime[s].l,prime[s].m,prime[s].n,prime[s].exp,
			     prime[u].x,prime[u].y,prime[u].z,prime[u].l,prime[u].m,prime[u].n,prime[u].exp,
			     prime[t].x,prime[t].y,prime[t].z,prime[t].l,prime[t].m,prime[t].n,prime[t].exp,
			     prime[v].x,prime[v].y,prime[v].z,prime[v].l,prime[v].m,prime[v].n,prime[v].exp);
			GAp[s*nPrime+t] -= PAf*EE;
			GBp[s*nPrime+t] -= PBf*EE;
		}
	}

	// convert back to basis space
	for(s=0; s < nPrime; s++)
	for(t=0; t <= s; t++){
		GA[ p2b[s]*nBasis+p2b[t] ] += GAp[s*nPrime+t]*Cfp[s]*Cfp[t];
		GB[ p2b[s]*nBasis+p2b[t] ] += GBp[s*nPrime+t]*Cfp[s]*Cfp[t];
		if(s!=t){
			GA[ p2b[t]*nBasis+p2b[s] ] += GAp[s*nPrime+t]*Cfp[s]*Cfp[t];
			GB[ p2b[t]*nBasis+p2b[s] ] += GBp[s*nPrime+t]*Cfp[s]*Cfp[t];
		}
	}
	*/

	/*
	 * main loop: integrals alignment
	 */

	double Pfuv,Pfst;
	double PAftv,PBftv,PAfsv,PBfsv;
	double PAftu,PBftu,PAfsu,PBfsu;
	int stuvT;

	for(s=0; s < nPrime; s++)
	for(t=0; t <= s; t++)
	for(u=0; u <= s; u++){ if(u==s) n=t; else n=u;
	for(v=0; v <= n; v++){

		// compute prefactor and screen at the same time
		EE=Fcp[s*nPrime+t]*Fcp[u*nPrime+v];
		if(fabs(EE) < cutoff) continue;

		// screen with density
		Pfuv  = PTp[u*nPrime+v]*EE; Pfst  = PTp[s*nPrime+t]*EE;
		PAftv = PAp[t*nPrime+v]*EE; PAfsv = PAp[s*nPrime+v]*EE;
		PBftv = PBp[t*nPrime+v]*EE; PBfsv = PBp[s*nPrime+v]*EE; 
		PAftu = PAp[t*nPrime+u]*EE; PAfsu = PAp[s*nPrime+u]*EE;
		PBftu = PBp[t*nPrime+u]*EE; PBfsu = PBp[s*nPrime+u]*EE; 

		if(fabs(Pfuv)  < cutoff && fabs(Pfst)  < cutoff &&
		   fabs(PAftv) < cutoff && fabs(PBftv) < cutoff &&
		   fabs(PAfsv) < cutoff && fabs(PBfsv) < cutoff &&
		   fabs(PAftu) < cutoff && fabs(PBftu) < cutoff &&
		   fabs(PAfsu) < cutoff && fabs(PBfsu) < cutoff ) continue;

		// compute integral
		EE = primeERI(
		     prime[s].x,prime[s].y,prime[s].z,prime[s].l,prime[s].m,prime[s].n,prime[s].exp,
		     prime[t].x,prime[t].y,prime[t].z,prime[t].l,prime[t].m,prime[t].n,prime[t].exp,
		     prime[u].x,prime[u].y,prime[u].z,prime[u].l,prime[u].m,prime[u].n,prime[u].exp,
		     prime[v].x,prime[v].y,prime[v].z,prime[v].l,prime[v].m,prime[v].n,prime[v].exp);

#define STUV_ALLSAME     5
#define STUV_2PAIRS      4
#define STUV_STPAIR      3
#define STUV_UVPAIR      2
#define STUV_SUTVPAIR    1
#define STUV_ALLDISTINCT 0

		// determine pS qS iS jS type for storage
		if((s==t)&&(u==v)&&(s==u)){ stuvT = STUV_ALLSAME;
		}else if((s==t)&&(u==v)){   stuvT = STUV_2PAIRS;
		}else if(s==t){             stuvT = STUV_STPAIR;
		}else if(u==v){             stuvT = STUV_UVPAIR;
		}else if((s==u)&&(t==v)){   stuvT = STUV_SUTVPAIR;
		}else{                      stuvT = STUV_ALLDISTINCT;
		}

		switch(stuvT){

		case STUV_ALLSAME:
			GAp[s*nPrime+u] -= PAftv*EE; GBp[s*nPrime+u] -= PBftv*EE;
		//	GAp[t*nPrime+u] -= PAfsv*EE; GBp[t*nPrime+u] -= PBfsv*EE;	
		//	GAp[s*nPrime+v] -= PAftu*EE; GBp[s*nPrime+v] -= PBftu*EE;
		//	GAp[t*nPrime+v] -= PAfsu*EE; GBp[t*nPrime+v] -= PBfsu*EE;

		//	GAp[u*nPrime+s] -= PAftv*EE; GBp[u*nPrime+s] -= PBftv*EE;
		//	GAp[v*nPrime+s] -= PAftu*EE; GBp[v*nPrime+s] -= PBftu*EE;	
		//	GAp[u*nPrime+t] -= PAfsv*EE; GBp[u*nPrime+t] -= PBfsv*EE;
		//	GAp[v*nPrime+t] -= PAfsu*EE; GBp[v*nPrime+t] -= PBfsu*EE;

		//	EE+=EE;
			GAp[s*nPrime+t] += Pfuv*EE; GBp[s*nPrime+t] += Pfuv*EE;
		//	GAp[u*nPrime+v] += Pfst*EE; GBp[u*nPrime+v] += Pfst*EE;
		break;

		case STUV_2PAIRS:
			GAp[s*nPrime+u] -= PAftv*EE; GBp[s*nPrime+u] -= PBftv*EE;
		//	GAp[t*nPrime+u] -= PAfsv*EE; GBp[t*nPrime+u] -= PBfsv*EE;	
		//	GAp[s*nPrime+v] -= PAftu*EE; GBp[s*nPrime+v] -= PBftu*EE;
		//	GAp[t*nPrime+v] -= PAfsu*EE; GBp[t*nPrime+v] -= PBfsu*EE;

			GAp[u*nPrime+s] -= PAftv*EE; GBp[u*nPrime+s] -= PBftv*EE;
		//	GAp[v*nPrime+s] -= PAftu*EE; GBp[v*nPrime+s] -= PBftu*EE;	
		//	GAp[u*nPrime+t] -= PAfsv*EE; GBp[u*nPrime+t] -= PBfsv*EE;
		//	GAp[v*nPrime+t] -= PAfsu*EE; GBp[v*nPrime+t] -= PBfsu*EE;

		//	EE+=EE;
			GAp[s*nPrime+t] += Pfuv*EE; GBp[s*nPrime+t] += Pfuv*EE;
			GAp[u*nPrime+v] += Pfst*EE; GBp[u*nPrime+v] += Pfst*EE;
		break;

		case STUV_STPAIR:
			GAp[s*nPrime+u] -= PAftv*EE; GBp[s*nPrime+u] -= PBftv*EE;
		//	GAp[t*nPrime+u] -= PAfsv*EE; GBp[t*nPrime+u] -= PBfsv*EE;	
			GAp[s*nPrime+v] -= PAftu*EE; GBp[s*nPrime+v] -= PBftu*EE;
		//	GAp[t*nPrime+v] -= PAfsu*EE; GBp[t*nPrime+v] -= PBfsu*EE;

			GAp[u*nPrime+s] -= PAftv*EE; GBp[u*nPrime+s] -= PBftv*EE;
			GAp[v*nPrime+s] -= PAftu*EE; GBp[v*nPrime+s] -= PBftu*EE;	
		//	GAp[u*nPrime+t] -= PAfsv*EE; GBp[u*nPrime+t] -= PBfsv*EE;
		//	GAp[v*nPrime+t] -= PAfsu*EE; GBp[v*nPrime+t] -= PBfsu*EE;

			GAp[u*nPrime+v] += Pfst*EE; GBp[u*nPrime+v] += Pfst*EE;
			EE+=EE;
			GAp[s*nPrime+t] += Pfuv*EE; GBp[s*nPrime+t] += Pfuv*EE;
		break;

		case STUV_UVPAIR:
			GAp[s*nPrime+u] -= PAftv*EE; GBp[s*nPrime+u] -= PBftv*EE;
			GAp[t*nPrime+u] -= PAfsv*EE; GBp[t*nPrime+u] -= PBfsv*EE;	
		//	GAp[s*nPrime+v] -= PAftu*EE; GBp[s*nPrime+v] -= PBftu*EE;
		//	GAp[t*nPrime+v] -= PAfsu*EE; GBp[t*nPrime+v] -= PBfsu*EE;

			GAp[u*nPrime+s] -= PAftv*EE; GBp[u*nPrime+s] -= PBftv*EE;
		//	GAp[v*nPrime+s] -= PAftu*EE; GBp[v*nPrime+s] -= PBftu*EE;	
			GAp[u*nPrime+t] -= PAfsv*EE; GBp[u*nPrime+t] -= PBfsv*EE;
		//	GAp[v*nPrime+t] -= PAfsu*EE; GBp[v*nPrime+t] -= PBfsu*EE;

			GAp[s*nPrime+t] += Pfuv*EE; GBp[s*nPrime+t] += Pfuv*EE;
			EE+=EE;
			GAp[u*nPrime+v] += Pfst*EE; GBp[u*nPrime+v] += Pfst*EE;
		break;

		case STUV_SUTVPAIR:
			GAp[s*nPrime+u] -= PAftv*EE; GBp[s*nPrime+u] -= PBftv*EE;
			GAp[t*nPrime+u] -= PAfsv*EE; GBp[t*nPrime+u] -= PBfsv*EE;	
			GAp[s*nPrime+v] -= PAftu*EE; GBp[s*nPrime+v] -= PBftu*EE;
			GAp[t*nPrime+v] -= PAfsu*EE; GBp[t*nPrime+v] -= PBfsu*EE;

		//	GAp[u*nPrime+s] -= PAftv*EE; GBp[u*nPrime+s] -= PBftv*EE;
		//	GAp[v*nPrime+s] -= PAftu*EE; GBp[v*nPrime+s] -= PBftu*EE;	
		//	GAp[u*nPrime+t] -= PAfsv*EE; GBp[u*nPrime+t] -= PBfsv*EE;
		//	GAp[v*nPrime+t] -= PAfsu*EE; GBp[v*nPrime+t] -= PBfsu*EE;

		//	EE+=EE;
			GAp[s*nPrime+t] += Pfuv*EE; GBp[s*nPrime+t] += Pfuv*EE;
			GAp[u*nPrime+v] += Pfst*EE; GBp[u*nPrime+v] += Pfst*EE;
		break;

		case STUV_ALLDISTINCT:
			GAp[s*nPrime+u] -= PAftv*EE; GBp[s*nPrime+u] -= PBftv*EE;
			GAp[t*nPrime+u] -= PAfsv*EE; GBp[t*nPrime+u] -= PBfsv*EE;	
			GAp[s*nPrime+v] -= PAftu*EE; GBp[s*nPrime+v] -= PBftu*EE;
			GAp[t*nPrime+v] -= PAfsu*EE; GBp[t*nPrime+v] -= PBfsu*EE;

			GAp[u*nPrime+s] -= PAftv*EE; GBp[u*nPrime+s] -= PBftv*EE;
			GAp[v*nPrime+s] -= PAftu*EE; GBp[v*nPrime+s] -= PBftu*EE;	
			GAp[u*nPrime+t] -= PAfsv*EE; GBp[u*nPrime+t] -= PBfsv*EE;
			GAp[v*nPrime+t] -= PAfsu*EE; GBp[v*nPrime+t] -= PBfsu*EE;

			EE+=EE;
			GAp[s*nPrime+t] += Pfuv*EE; GBp[s*nPrime+t] += Pfuv*EE;
			GAp[u*nPrime+v] += Pfst*EE; GBp[u*nPrime+v] += Pfst*EE;
		break;
		}

	}
	}

	// convert back to basis space
	for(s=0; s < nPrime; s++)
	for(t=0; t <= s; t++){
		GA[ p2b[s]*nBasis+p2b[t] ] += GAp[s*nPrime+t]*Cfp[s]*Cfp[t];
		GB[ p2b[s]*nBasis+p2b[t] ] += GBp[s*nPrime+t]*Cfp[s]*Cfp[t];
		if(s!=t){
			GA[ p2b[t]*nBasis+p2b[s] ] += GAp[s*nPrime+t]*Cfp[s]*Cfp[t];
			GB[ p2b[t]*nBasis+p2b[s] ] += GBp[s*nPrime+t]*Cfp[s]*Cfp[t];
		}
	}


	// free memory
	free(PTp);
	free(PAp);
	free(PBp);
	free(GAp);
	free(GBp);
	free(Fcp);
	free(p2b);
	free(Cfp);
	free(prime);
}

#define MAXL 4
#define DIST2(x1,y1,z1,x2,y2,z2) ((x1-x2)*(x1-x2)+\
                                  (y1-y2)*(y1-y2)+\
                                  (z1-z2)*(z1-z2))
#define MAXPRIMESHELL 16
struct PrimeShell_t{
	double x,y,z;         // center
	double exp;           // exponent
	int nPrime;           // number of included primitives
	int min,max;          // minimum and maximum primitive index
	int maxL;             // maximum angular index
	int l[MAXPRIMESHELL]; // x-direction angular index
	int m[MAXPRIMESHELL]; // y-direction angular index
	int n[MAXPRIMESHELL]; // z-direction angular index
};

// GTO_JK_Matrix_PrimeSpaceShellSet: generate JK matrix using primitive space
//
// Jan 19, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void GTO_JK_Matrix_PrimeSpaceShellSet(
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double cutoff,                 // cutoff to ignore
	double *GA,                    // return G for spin up
	double *GB){                   // return G for spin down

	int s,t,u,v;                   // looping variables for primitive functions
	int S,T,U,V,N;                 // looping variables for shell
	double EE;                     // coulomb integral

	double Bx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double By[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Bz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sx[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sy[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double Sz[4*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)*(MAXL+1)];
	double  F[4*(MAXL+1)];         // Boys functions
	double *tBx,*tBy,*tBz;         // pointer to Bx, By, Bz
	double *tSx,*tSy,*tSz;         // pointer to Sx, Sy, Sz
	int kX, kY, kZ;                // summation loop
	int mX, mY, mZ;                // maximum angular parts
	double  tSum;                  // intermediate variables for kX,kY,kZ loop

	static int nPrime;                  // number of primitive
	struct GTOPrime_t *prime;           // primitive data structure
	static int firstVisit=1;            // first visit flag
	static int nShell;                  // number of shell
	static struct PrimeShell_t *shell;  // primitive shell data structure
	static int *p2b;                    // primitive --> basis mapping
	static double *FcS;                 // ERI prefactor in shell space
	static double *Cfp;                 // coefficient for each primitive
	static double *GAp,*GBp;            // G in primitive space
	static double *PTp,*PAp,*PBp;       // density matrix in primitive space
	static double *PMS;                 // maximum density in shell space
	double maxP=0;                      // maximum density
	double maxFcS=0;                    // maximum prefactor

	// reset call
	if(nBasis==0){
		firstVisit=1;
		free(PTp);
		free(PMS);
		free(PAp);
		free(PBp);
		free(GAp);
		free(GBp);
		free(FcS);
		free(p2b);
		free(Cfp);
		free(shell);
	}

	/*
	 * first visit preparations
	 */
	if(firstVisit){
		firstVisit=0;

		// compute number of primitive function
		for(nPrime=0, u=0; u < nBasis; u++) nPrime += gto[u].nContract;

		// allocate memory
		GAp   = calloc(nPrime*nPrime, sizeof(double));
		GBp   = calloc(nPrime*nPrime, sizeof(double));
		PTp   = calloc(nPrime*nPrime, sizeof(double));
		PAp   = calloc(nPrime*nPrime, sizeof(double));
		PBp   = calloc(nPrime*nPrime, sizeof(double));
		Cfp   = calloc(nPrime,        sizeof(double));
		p2b   = calloc(nPrime,        sizeof(int));
		prime = calloc(nPrime,        sizeof(struct GTOPrime_t)); 
		if(GAp==NULL || GBp==NULL || PTp==NULL || PAp==NULL || PBp==NULL ||
		   p2b==NULL || Cfp==NULL || prime==NULL){
			printf("GTO_JK_Matrix_PrimeSpaceShellSet - error allocating memory");
			exit(-1);
		}

		// generate primitive to basis mapping and coefficients values
		for(s=0, u=0; u < nBasis; u++)
			for(v=0; v < gto[u].nContract; v++){

				// generate primitive to basis mapping
				p2b[s] = u; 

				// generate coefficient value
				Cfp[s] = gto[u].norm[v] * gto[u].coef[v];

				// fill in primitive data 
				prime[s].x   = gto[u].x0;
				prime[s].y   = gto[u].y0;
				prime[s].z   = gto[u].z0;
				prime[s].l   = gto[u].l;
				prime[s].m   = gto[u].m;
				prime[s].n   = gto[u].n;
				prime[s].exp = gto[u].exp[v];

				s++; 
			}

		// re-ordering of primitive in order to group into shell
		for(s=0; s < nPrime; s++){
			for(t=s+1; t < nPrime; t++){
				// different center stop searching
				if(prime[t].x != prime[s].x ||
				   prime[t].y != prime[s].y ||
				   prime[t].z != prime[s].z) break;
				// swap if matched
				if(prime[t].exp == prime[s].exp){
					EE = prime[t].exp; prime[t].exp = prime[s+1].exp; prime[s+1].exp = EE;
					u  = prime[t].l;   prime[t].l   = prime[s+1].l;   prime[s+1].l   = u;
					u  = prime[t].m;   prime[t].m   = prime[s+1].m;   prime[s+1].m   = u;
					u  = prime[t].n;   prime[t].n   = prime[s+1].n;   prime[s+1].n   = u;
					u  =   p2b[t];       p2b[t]     =   p2b[s+1];       p2b[s+1]     = u;
					EE =   Cfp[t];       Cfp[t]     =   Cfp[s+1];       Cfp[s+1]     = EE;
				}
			}
		}

		// allocate and copy data into shell
		shell = calloc(nPrime, sizeof(struct PrimeShell_t));
		if(shell==NULL){
			printf("GTO_JK_Matrix_PrimeSpaceShellSet - error allocating memory");
			exit(-1);
		}
		for(S=0, s=0; s < nPrime; s++){
			// save values
			shell[S].x      = prime[s].x;
			shell[S].y      = prime[s].y;
			shell[S].z      = prime[s].z;
			shell[S].exp    = prime[s].exp;
			shell[S].min    = s;
			shell[S].nPrime = 0;
			shell[S].maxL   = 0;
	
			// until shell terminates
			for(t=s; t < nPrime; t++)
				if(prime[t].x   == prime[s].x &&
				   prime[t].y   == prime[s].y &&
				   prime[t].z   == prime[s].z &&
				   prime[t].exp == prime[s].exp){
					shell[S].l[shell[S].nPrime] = prime[t].l;
					shell[S].m[shell[S].nPrime] = prime[t].m;
					shell[S].n[shell[S].nPrime] = prime[t].n;
					if( shell[S].maxL < prime[t].l+prime[t].m+prime[t].n)
						shell[S].maxL = prime[t].l+prime[t].m+prime[t].n;
					shell[S].nPrime++;
					if(shell[S].nPrime > MAXPRIMESHELL){
						printf("GTO_JK_Matrix_PrimeSpaceShellSet - MAXPRIMESHELL too small");
						exit(-1);
					}
				}
				else
					break;
	
			shell[S].max = shell[S].min + (shell[S].nPrime-1);
			s           += (shell[S].nPrime-1);
			S++; 
		}
		free(prime);
		nShell = S;
		shell  = realloc(shell, sizeof(struct PrimeShell_t)*nShell);
		FcS    = calloc(nShell*nShell, sizeof(double));
		PMS    = calloc(nShell*nShell, sizeof(double));
		if(FcS==NULL || PMS==NULL){
			printf("GTO_JK_Matrix_PrimeSpaceShellSet - error allocating memory");
			exit(-1);
		}
		for(S=0; S < nShell; S++)
		for(T=0; T <= S; T++){
#define SQRT_TWOPI_5HALF 5.9149671728
			// generate ERI prefactor
			FcS[S*nShell+T] = SQRT_TWOPI_5HALF *
			                  exp( -  shell[S].exp * shell[T].exp 
			                       / (shell[S].exp + shell[T].exp)
			                       *  DIST2(shell[S].x,shell[S].y,shell[S].z,
			                                shell[T].x,shell[T].y,shell[T].z) )
			                  /(shell[S].exp + shell[T].exp);
	
			FcS[T*nShell+S] = FcS[S*nShell+T];
		}

	}

	/*
	 * handle density matrix in primitive space
	 */

	for(s=0; s < nPrime; s++)
	for(t=0; t <= s; t++){

		// reset GAp GBp
		GAp[s*nPrime+t] = 0.0;
		GBp[s*nPrime+t] = 0.0;
		
		// generate density matrix in primitive space
		PAp[s*nPrime+t] = PA[ p2b[s]*nBasis + p2b[t] ] * Cfp[s] * Cfp[t];
		PBp[s*nPrime+t] = PB[ p2b[s]*nBasis + p2b[t] ] * Cfp[s] * Cfp[t];
		PTp[s*nPrime+t] = PAp[s*nPrime+t]   + PBp[s*nPrime+t];

		PTp[t*nPrime+s] = PTp[s*nPrime+t];
		PAp[t*nPrime+s] = PAp[s*nPrime+t];
		PBp[t*nPrime+s] = PBp[s*nPrime+t];
		GAp[t*nPrime+s] = GAp[s*nPrime+t];
		GBp[t*nPrime+s] = GAp[s*nPrime+t];
	}

	// compute maximum density and prefactors
	for(S=0; S < nShell; S++)
		for(T=0; T <= S; T++){
		PMS[S*nShell+T] = 0;
		for(s=shell[S].min; s <= shell[S].max; s++)
		for(t=shell[T].min; t <= shell[T].max; t++){
			if(PMS[S*nShell+T] < fabs(PTp[s*nPrime+t])) PMS[S*nShell+T] = fabs(PTp[s*nPrime+t]);
			if(PMS[S*nShell+T] < fabs(PAp[s*nPrime+t])) PMS[S*nShell+T] = fabs(PAp[s*nPrime+t]);
			if(PMS[S*nShell+T] < fabs(PBp[s*nPrime+t])) PMS[S*nShell+T] = fabs(PBp[s*nPrime+t]);
		}
		if(maxP < PMS[S*nShell+T]) maxP = PMS[S*nShell+T];
		if(maxFcS < FcS[S*nShell+T]) maxFcS = FcS[S*nShell+T];
		PMS[T*nShell+S] = PMS[S*nShell+T];
	}

	/*
	 * main loop: integrals alignment
	 */

	double K;
	double Pfuv,Pfst;
	double PAftv,PBftv,PAfsv,PBfsv;
	double PAftu,PBftu,PAfsu,PBfsu;
	int stuvT,sI,tI,uI,vI;
	struct PrimeShell_t *sS,*tS,*uS,*vS;

	for(S=0; S < nShell; S++)
	for(T=0; T <= S; T++){
	if(FcS[S*nShell+T]*maxP*maxFcS < cutoff) continue;
	sS = &shell[S]; tS = &shell[T]; 
	for(U=0; U <= S; U++){ if(U==S) N=T; else N=U;
	for(V=0; V <= N; V++){

		// compute prefactor and screen at the same time
		K=FcS[S*nShell+T]*FcS[U*nShell+V];
		if(K*maxP < cutoff) continue;

		uS = &shell[U]; vS = &shell[V];
		K /= sqrt(sS->exp + tS->exp + uS->exp + vS->exp);

		if(K*PMS[S*nShell+T] < cutoff && K*PMS[U*nShell+V] < cutoff &&
		   K*PMS[S*nShell+U] < cutoff && K*PMS[S*nShell+V] < cutoff &&
		   K*PMS[T*nShell+U] < cutoff && K*PMS[T*nShell+V] < cutoff) continue;

		genPrimeSetBxyzSxyzF(sS->x,sS->y,sS->z,sS->maxL,sS->exp,
		                     tS->x,tS->y,tS->z,tS->maxL,tS->exp,
		                     uS->x,uS->y,uS->z,uS->maxL,uS->exp,
		                     vS->x,vS->y,vS->z,vS->maxL,vS->exp,
		                     Bx,By,Bz,Sx,Sy,Sz,F);
		
		// determine S T U V type for storage
		if((S==T)&&(U==V)&&(S==U)){ stuvT = STUV_ALLSAME;
		}else if((S==T)&&(U==V)){   stuvT = STUV_2PAIRS;
		}else if(S==T){             stuvT = STUV_STPAIR;
		}else if(U==V){             stuvT = STUV_UVPAIR;
		}else if((S==U)&&(T==V)){   stuvT = STUV_SUTVPAIR;
		}else{                      stuvT = STUV_ALLDISTINCT;
		}

		for(s=sS->max,sI=sS->nPrime-1; sI >=0 ; s--,sI--)
		for(t=tS->max,tI=tS->nPrime-1; tI >=0 ; t--,tI--)
		for(u=uS->max,uI=uS->nPrime-1; uI >=0 ; u--,uI--)
		for(v=vS->max,vI=vS->nPrime-1; vI >=0 ; v--,vI--){

			// screen with density
			Pfuv  = PTp[u*nPrime+v]*K; Pfst  = PTp[s*nPrime+t]*K;
			PAftv = PAp[t*nPrime+v]*K; PAftu = PAp[t*nPrime+u]*K;
			PAfsv = PAp[s*nPrime+v]*K; PAfsu = PAp[s*nPrime+u]*K;
			PBftv = PBp[t*nPrime+v]*K; PBftu = PBp[t*nPrime+u]*K;
			PBfsv = PBp[s*nPrime+v]*K; PBfsu = PBp[s*nPrime+u]*K;

			if(fabs(Pfuv)  < cutoff && fabs(Pfst)  < cutoff &&
			   fabs(PAftv) < cutoff && fabs(PBftv) < cutoff &&
			   fabs(PAfsv) < cutoff && fabs(PBfsv) < cutoff &&
			   fabs(PAftu) < cutoff && fabs(PBftu) < cutoff &&
			   fabs(PAfsu) < cutoff && fabs(PBfsu) < cutoff ) continue;

			kZ  =    (vS->maxL+1);
			kY  = kZ*(uS->maxL+1);
			kX  = kY*(tS->maxL+1);
			mX  = 4*(MAXL+1)*(sS->l[sI]*kX + tS->l[tI]*kY + uS->l[uI]*kZ + vS->l[vI]);
			mY  = 4*(MAXL+1)*(sS->m[sI]*kX + tS->m[tI]*kY + uS->m[uI]*kZ + vS->m[vI]);
			mZ  = 4*(MAXL+1)*(sS->n[sI]*kX + tS->n[tI]*kY + uS->n[uI]*kZ + vS->n[vI]);
			tBx = Bx + mX;
			tBy = By + mY;
			tBz = Bz + mZ; 
			tSx = Sx + mX;
			tSy = Sy + mY;
			tSz = Sz + mZ;

			mX = sS->l[sI] + tS->l[tI] + uS->l[uI] + vS->l[vI];
			mY = sS->m[sI] + tS->m[tI] + uS->m[uI] + vS->m[vI];
			mZ = sS->n[sI] + tS->n[tI] + uS->n[uI] + vS->n[vI];

			// determine type of mX,mY,mZ for loop optimization
			EE=0.0;
			if(mX==0)       			for(kY=mY;kY>=0;kY--) EE += tBy[kY]*tSz[kY];

			else if(mY==0)  			for(kZ=mZ;kZ>=0;kZ--) EE += tBz[kZ]*tSx[kZ];

			else if(mZ==0)  			for(kX=mX;kX>=0;kX--) EE += tBx[kX]*tSy[kX];

			else if(mX > mY && mX > mZ)	for(kY=mY;kY>=0;kY--){
											for(tSum=0.0,kZ=mZ;kZ>=0;kZ--)
												tSum+=tBz[kZ]*tSx[kY+kZ];
											EE += tBy[kY]*tSum;
										}

			else if(mY > mZ)			for(kZ=mZ;kZ>=0;kZ--){
											for(tSum=0.0,kX=mX;kX>=0;kX--)
												tSum+=tBx[kX]*tSy[kX+kZ];
											EE += tBz[kZ]*tSum;
										}

			else						for(kX=mX;kX>=0;kX--){
											for(tSum=0.0,kY=mY;kY>=0;kY--)
												tSum+=tBy[kY]*tSz[kX+kY];
											EE += tBx[kX]*tSum;
										}

			//EE = 0.0; 
			//for(kX=mX;kX>=0;kX--)
			//for(kY=mY;kY>=0;kY--)
			//for(kZ=mZ;kZ>=0;kZ--)
			//	EE += tBx[kX]*tBy[kY]*tBz[kZ]*F[kX+kY+kZ];

			// storage
			switch(stuvT){
	
			case STUV_ALLSAME:
				GAp[s*nPrime+u] -= PAftv*EE; GBp[s*nPrime+u] -= PBftv*EE;
				GAp[s*nPrime+t] += Pfuv*EE;  GBp[s*nPrime+t] += Pfuv*EE;
			break;
	
			case STUV_2PAIRS:
				GAp[s*nPrime+u] -= PAftv*EE; GAp[u*nPrime+s] -= PAftv*EE;
				GBp[s*nPrime+u] -= PBftv*EE; GBp[u*nPrime+s] -= PBftv*EE;

				GAp[s*nPrime+t] += Pfuv*EE;  GBp[s*nPrime+t] += Pfuv*EE;
				GAp[u*nPrime+v] += Pfst*EE;  GBp[u*nPrime+v] += Pfst*EE;
			break;
	
			case STUV_STPAIR:
				GAp[s*nPrime+u] -= PAftv*EE; GAp[s*nPrime+v] -= PAftu*EE;
				GAp[u*nPrime+s] -= PAftv*EE; GAp[v*nPrime+s] -= PAftu*EE;
				GBp[s*nPrime+u] -= PBftv*EE; GBp[s*nPrime+v] -= PBftu*EE;
				GBp[u*nPrime+s] -= PBftv*EE; GBp[v*nPrime+s] -= PBftu*EE;	

				GAp[u*nPrime+v] += Pfst*EE;  GBp[u*nPrime+v] += Pfst*EE; EE+=EE;
				GAp[s*nPrime+t] += Pfuv*EE;  GBp[s*nPrime+t] += Pfuv*EE;
			break;

			case STUV_UVPAIR:
				GAp[s*nPrime+u] -= PAftv*EE; GAp[t*nPrime+u] -= PAfsv*EE;
				GAp[u*nPrime+s] -= PAftv*EE; GAp[u*nPrime+t] -= PAfsv*EE;
				GBp[s*nPrime+u] -= PBftv*EE; GBp[t*nPrime+u] -= PBfsv*EE;	
				GBp[u*nPrime+s] -= PBftv*EE; GBp[u*nPrime+t] -= PBfsv*EE;

				GAp[s*nPrime+t] += Pfuv*EE;  GBp[s*nPrime+t] += Pfuv*EE; EE+=EE;
				GAp[u*nPrime+v] += Pfst*EE;  GBp[u*nPrime+v] += Pfst*EE;
			break;

			case STUV_SUTVPAIR:
				GAp[s*nPrime+u] -= PAftv*EE; GAp[s*nPrime+v] -= PAftu*EE;
				GAp[t*nPrime+u] -= PAfsv*EE; GAp[t*nPrime+v] -= PAfsu*EE;
				GBp[s*nPrime+u] -= PBftv*EE; GBp[s*nPrime+v] -= PBftu*EE;
				GBp[t*nPrime+u] -= PBfsv*EE; GBp[t*nPrime+v] -= PBfsu*EE;

				GAp[s*nPrime+t] += Pfuv*EE;  GAp[u*nPrime+v] += Pfst*EE;
				GBp[s*nPrime+t] += Pfuv*EE;  GBp[u*nPrime+v] += Pfst*EE;
			break;
	
			case STUV_ALLDISTINCT:
				GAp[s*nPrime+u] -= PAftv*EE; GAp[s*nPrime+v] -= PAftu*EE;
				GAp[t*nPrime+u] -= PAfsv*EE; GAp[t*nPrime+v] -= PAfsu*EE;
				GAp[u*nPrime+s] -= PAftv*EE; GAp[u*nPrime+t] -= PAfsv*EE;
				GAp[v*nPrime+s] -= PAftu*EE; GAp[v*nPrime+t] -= PAfsu*EE;
				GBp[s*nPrime+u] -= PBftv*EE; GBp[s*nPrime+v] -= PBftu*EE;
				GBp[t*nPrime+u] -= PBfsv*EE; GBp[t*nPrime+v] -= PBfsu*EE;
				GBp[u*nPrime+s] -= PBftv*EE; GBp[u*nPrime+t] -= PBfsv*EE;
				GBp[v*nPrime+s] -= PBftu*EE; GBp[v*nPrime+t] -= PBfsu*EE;

				EE+=EE;
				GAp[s*nPrime+t] += Pfuv*EE;  GAp[u*nPrime+v] += Pfst*EE;
				GBp[s*nPrime+t] += Pfuv*EE;  GBp[u*nPrime+v] += Pfst*EE;
			break;
			}
		}
	}
	}
	}

	// convert back to basis space
	for(s=0; s < nPrime; s++)
	for(t=0; t <= s; t++){
		GA[ p2b[s]*nBasis+p2b[t] ] += GAp[s*nPrime+t]*Cfp[s]*Cfp[t];
		GB[ p2b[s]*nBasis+p2b[t] ] += GBp[s*nPrime+t]*Cfp[s]*Cfp[t];
		if(s!=t){
			GA[ p2b[t]*nBasis+p2b[s] ] += GAp[s*nPrime+t]*Cfp[s]*Cfp[t];
			GB[ p2b[t]*nBasis+p2b[s] ] += GBp[s*nPrime+t]*Cfp[s]*Cfp[t];
		}
	}

}
#undef MAXL
