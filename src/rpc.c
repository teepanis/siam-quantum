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
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/file.h>

#ifdef WIN32
#include <windows.h>
#endif

#include "matrix.h"
#include "rpc.h"
#include "mp2.h"
#include "grad.h"
//#include "excite.h"

// pack_Molecule: packs the Molecule_t data structure into a long string
// so that it can be transfered across network/file system
//
// Mar 6, 2013 - Teepanis Chachiyo
//      Initial implementation and testing
//
void * pack_Molecule(const struct Molecule_t *mol, int *nBytes){
	int len=0;    // size of buffer
	void *buffer; // buffer to store the data
	char *ptr;    // generic pointer

	// adding nAtom, Q
	len += 2*sizeof(int);

	// adding Z,x,y,z
	len += mol->nAtom*sizeof(int) + mol->nAtom*3*sizeof(double);

	// set length of buffer
	*nBytes = len;

	// allocate buffer
	if((buffer=calloc(1,len))==NULL){
		printf("pack_Molecule - error cannot allocate memory\n");
		exit(-1);
	}

	//
	// pack data
	//
	ptr = (char *)buffer;

	// copy nAton, Q
	len = 2*sizeof(int);             memcpy(ptr,mol,len);    ptr += len;

	// copy Z
	len = mol->nAtom*sizeof(int);    memcpy(ptr,mol->Z,len); ptr += len;

	// copy x
	len = mol->nAtom*sizeof(double); memcpy(ptr,mol->x,len); ptr += len;

	// copy y
	len = mol->nAtom*sizeof(double); memcpy(ptr,mol->y,len); ptr += len;

	// copy z
	len = mol->nAtom*sizeof(double); memcpy(ptr,mol->z,len); ptr += len;

	return buffer;
}


// unpack_GTOBasis: unpacks the buffer and turned it into the original
// Molecule_t data structure. This is meant to be used with the pack_Molecule
// subroutine.
//
// Mar 6, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
struct Molecule_t * unpack_Molecule(void *buffer, int *nBytes){

	int len;                  // generic counter
	struct Molecule_t *mol;   // molecule data structure 
	char *ptr;                // generic pointer

	// allocate memory
	if((mol=calloc(1,sizeof(struct Molecule_t)))==NULL){
		printf("unpack_Molecule - error cannot allocate memory\n");
		exit(-1);
	}

	// unpack data
	ptr = (char *)buffer;

	// copy nAtom, Q
	len = 2*sizeof(int); memcpy(mol,ptr,len); ptr += len;

	// allocate Z,x,y,z
	mol->Z = calloc(mol->nAtom,sizeof(int));
	mol->x = calloc(mol->nAtom,sizeof(double));
	mol->y = calloc(mol->nAtom,sizeof(double));
	mol->z = calloc(mol->nAtom,sizeof(double));

	// copy Z
	len = mol->nAtom*sizeof(int);    memcpy(mol->Z,ptr,len); ptr += len;

	// copy x
	len = mol->nAtom*sizeof(double); memcpy(mol->x,ptr,len); ptr += len;

	// copy y
	len = mol->nAtom*sizeof(double); memcpy(mol->y,ptr,len); ptr += len;

	// copy z
	len = mol->nAtom*sizeof(double); memcpy(mol->z,ptr,len); ptr += len;

	// set the number of bytes
	*nBytes = (ptr-(char *)buffer);

	return mol;
}


// pack_GTOBasis: packs the GTOBasis_t data structure into a long string
// so that it can be transfered across network/file system.
//
// Mar 1, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void * pack_GTOBasis(int nBasis, const struct GTOBasis_t *gto, int *nBytes){
	int len=0;    // size of buffer
	int i;        // generic counter
	void *buffer; // buffer to store the data	
	char *ptr;    // generic pointer

	// count size of the buffer
	for(i=0; i < nBasis; i++){

		// adding nContract,l,m,n,x0,y0,z0
		len += 4*sizeof(int) + 3*sizeof(double);

		// adding coef,exp,norm
		len += 3*gto[i].nContract*sizeof(double);
	}
	*nBytes = len;

	// allocate buffer
	if((buffer=calloc(1,len))==NULL){
		printf("pack_GTOBasis - error cannot allocate memory\n");
		exit(-1);
	}

	// pack data
	ptr = (char *)buffer;
	for(i=0; i < nBasis; i++){

		// copy nContract,l,m,n,x0,y0,z0
		len = 4*sizeof(int) + 3*sizeof(double);
		memcpy(ptr,&gto[i],len);
		ptr += len;

		// copy coef
		len = gto[i].nContract*sizeof(double);
		memcpy(ptr,gto[i].coef,len);
		ptr += len;

		// copy exp
		len = gto[i].nContract*sizeof(double);
		memcpy(ptr,gto[i].exp,len);
		ptr += len;

		// copy norm
		len = gto[i].nContract*sizeof(double);
		memcpy(ptr,gto[i].norm,len);
		ptr += len;
	}

	return buffer;
}


// unpack_GTOBasis: unpacks the buffer and turned it into the original
// GTOBasis_t data structure. This is meant to be used with the pack_GTOBasis
// subroutine.
//
// Mar 1, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
struct GTOBasis_t * unpack_GTOBasis(int nBasis, void *buffer, int *nBytes){

	int i,len;                // generic counter
	struct GTOBasis_t *gto;   // basis data structure 
	char *ptr;                // generic pointer

	// allocate memory
	if((gto=calloc(nBasis,sizeof(struct GTOBasis_t)))==NULL){
		printf("unpack_GTOBasis - error cannot allocate memory\n");
		exit(-1);
	}

	// unpack data
	ptr = (char *)buffer;
	for(i=0; i < nBasis; i++){
		
		// copy nContract,l,m,n,x0,y0,z0
		len = 4*sizeof(int) + 3*sizeof(double);
		memcpy(&gto[i],ptr,len);
		ptr += len;

		// allocate coef, exp, norm
		gto[i].coef = calloc(gto[i].nContract,sizeof(double));
		gto[i].exp  = calloc(gto[i].nContract,sizeof(double));
		gto[i].norm = calloc(gto[i].nContract,sizeof(double));
		if(gto[i].coef==NULL || gto[i].exp==NULL || gto[i].norm==NULL){
			printf("unpack_GTOBasis - error cannot allocate memory\n");
			exit(-1);
		}

		// copy coef
		len = gto[i].nContract*sizeof(double);
		memcpy(gto[i].coef,ptr,len);
		ptr += len;

		// copy exp
		len = gto[i].nContract*sizeof(double);
		memcpy(gto[i].exp,ptr,len);
		ptr += len;

		// copy norm
		len = gto[i].nContract*sizeof(double);
		memcpy(gto[i].norm,ptr,len);
		ptr += len;
	}

	// set the number of bytes
	*nBytes = (ptr-(char *)buffer);

	return gto;
}


// pack_GTO_JK_Matrix_Quartet_Parallel: packs the arguments for this subroutine
// into a long string so that it can be transferred over network/file system.
//
// May 6, 2016 - Teepanis Chachiyo
//     Use JT,KA,KB instead of GA,GB
//
// Mar 1, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void * pack_GTO_JK_Matrix_Quartet_Parallel(
	int childID,                   // child id number
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *schwarz_basis,   // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *JT,                    // return J total
	double *KA,                    // return G for spin up
	double *KB,                    // return G for spin down
	struct option_t *opt,          // global option
	int *nBytes){                  // size of the buffer
	
	int len=0;       // size of the buffer
	void *gtoBuffer; // buffer for gto
	int gtoLen;      // size of the gto buffer
	void *buffer;    // buffer
	char *ptr;       // generic pointer

	// add childID,nBasis,PA,PB
	len += 2*sizeof(int) + 2*nBasis*nBasis*sizeof(double);

	// add gto
	gtoBuffer = pack_GTOBasis(nBasis,gto,&gtoLen);
	len += gtoLen;

	// add schwarz_basis,fixedCutoff,JT,KA,KB,opt
	len += nBasis*nBasis*sizeof(double) 
	       + sizeof(double) 
	       + 3*nBasis*nBasis*sizeof(double)
	       + sizeof(struct option_t);

	// set the number of bytes
	*nBytes = len;

	// allocate memory
	if((buffer=calloc(1,len))==NULL){
		printf("pack_GTO_JK_Matrix_Quartet - error cannot allocate memory\n");
		exit(-1);
	}

	// initialize pointer
	ptr = (char *)buffer;

	// copy childID
	len = sizeof(int); memcpy(ptr,&childID,len); ptr += len;

	// copy nBasis
	len = sizeof(int); memcpy(ptr,&nBasis,len); ptr += len;

	// copy PA and PB
	len = nBasis*nBasis*sizeof(double);
	memcpy(ptr,PA,len); ptr += len;
	memcpy(ptr,PB,len); ptr += len;

	// copy gto
	len = gtoLen; memcpy(ptr,gtoBuffer,len); ptr += len;

	// copy schwarz_basis
	len = nBasis*nBasis*sizeof(double);
	memcpy(ptr,schwarz_basis,len); ptr += len;

	// copy fixedCutoff
	len = sizeof(double); memcpy(ptr,&fixedCutoff,len); ptr += len;

	// copy GA and GB
	len = nBasis*nBasis*sizeof(double);
	memcpy(ptr,JT,len); ptr += len;
	memcpy(ptr,KA,len); ptr += len;
	memcpy(ptr,KB,len); ptr += len;

	// copy opt
	len = sizeof(struct option_t); memcpy(ptr,opt,len); ptr += len;

	// clean memory
	free(gtoBuffer);

	return buffer;
}


// lpc_GTO_JK_Matrix_Quartet_Parallel: calls the function GTO_JK_Matrix_Quartet locally.
// It then packs the returned GA and GB and the reply.
//
// May 6, 2016 - Teepanis Chachiyo
//     Use JT,KA,KB instead of GA,GB
//
// Mar 2, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void * lpc_GTO_JK_Matrix_Quartet_Parallel(
	void *buffer,                  // buffer containing arguments
	int *nBytes){                  // size of the buffer

	// arguments
	int childID;             // child id number
	int nBasis;              // number of basis functions
	double *PA;              // density matrix for spin up
	double *PB;              // density matrix for spin down 
	struct GTOBasis_t *gto;  // basis set info
	double *schwarz_basis;   // pointer to schwarz matrix
	double fixedCutoff;      // cutoff to ignore
	double *JT;              // return J total
	double *KA;              // return G for spin up
	double *KB;              // return G for spin down
	struct option_t *opt;    // global option

	char *ptr;               // generic pointer
	int len;                 // generic int
	void *reBuffer;          // returned buffer

	// initialize pointer
	ptr = buffer;

	// unpack childID
	len = sizeof(int); memcpy(&childID,ptr,len); ptr += len;

	// unpack nBasis
	len = sizeof(int); memcpy(&nBasis,ptr,len); ptr += len;

	// unpack PA,PB
	len = nBasis*nBasis*sizeof(double);
	PA = (double *)ptr; ptr += len;
	PB = (double *)ptr; ptr += len;

	// unpack gto
	gto = unpack_GTOBasis(nBasis,ptr,&len); ptr += len;

	// unpack schwarz_basis
	len = nBasis*nBasis*sizeof(double); schwarz_basis = (double *)ptr; ptr += len;

	// unpack fixedCutoff
	len = sizeof(double); memcpy(&fixedCutoff,ptr,len); ptr += len;

	// unpack JT,KA,KB,opt
	len = nBasis*nBasis*sizeof(double); JT = (double *)ptr;          ptr += len;
	                                    KA = (double *)ptr;          ptr += len;
	                                    KB = (double *)ptr;          ptr += len;
	len = sizeof(struct option_t);      opt= (struct option_t *)ptr; ptr += len;

	// call the actual function
	GTO_JK_Matrix_Quartet_Parallel(childID,nBasis,PA,PB,gto,schwarz_basis,fixedCutoff,JT,KA,KB,opt);

	// pack the returned JT, KA and KB value
	if((reBuffer=calloc(3*nBasis*nBasis,sizeof(double)))==NULL){
		printf("lpc_GTO_JK_Matrix_Quartet - error cannot allocate memory\n");
		exit(-1);
	}
	ptr = (char *)reBuffer;
	len = nBasis*nBasis*sizeof(double); memcpy(ptr,JT,len); ptr += len;
	                                    memcpy(ptr,KA,len); ptr += len;
	                                    memcpy(ptr,KB,len); ptr += len;

	// set the size of reBuffer
	*nBytes = 3*nBasis*nBasis*sizeof(double);

	// clean memory
	cleanGTOBasis(gto,nBasis);

	return reBuffer;
}


// delay: delay execution for stability
//
// Mar 2, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void delay(){
#ifdef WIN32
	Sleep(RPC_DELAY);
#else
	usleep(RPC_DELAY*1000);
#endif
}


// rpcWriteMsg: writes the message into the file
//
// Mar 3, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void rpcWriteMsg(char *fname, struct rpcHdr_t *hdr, void *buffer){
	FILE *fd;      // file pointer
	int len;       // generic int

	// open file for writing
	if((fd=fopen(fname,"wb"))==NULL){
		printf("rpcWriteMsg - error cannot open file\n");
		exit(-1);
	}

	// lock file
	//#ifndef WIN32
	//flock(fileno(fd), LOCK_EX);
	//#endif

	// write header to file
	len = sizeof(struct rpcHdr_t);
	if(fwrite(hdr,1,len,fd) != len){
		printf("rpcWriteMsg - error writing header\n");
		exit(-1);
	}

	// write output to file
	len = hdr->len;
	if(fwrite(buffer,1,len,fd) != len){
		printf("rpcWriteMsg - error writing buffer\n");
		exit(-1);
	}

	// flush
	if(fflush(fd)!=0){
		printf("rpcWriteMsg - error cannot flush file\n");
		exit(-1);
	}

	// release lock
	//#ifndef WIN32
	//flock(fileno(fd), LOCK_UN);
	//#endif

	// close file
	if(fclose(fd) != 0){
		printf("rpcWriteMsg - error cannot close file\n");
		exit(-1);
	}
}


// rpcReadMsg: read the message and return the pointer
// to the buffer. Also, the hdr value is set.
//
// Mar 3, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void * rpcReadMsg(char *fname, struct rpcHdr_t *hdr){
	FILE *fd;     // file pointer
	int len;      // generic int
	void *buffer; // buffer

	// open file for reading
	if((fd=fopen(fname,"rb"))==NULL){
		printf("rpcReadMsg - error cannot open file\n");
		exit(-1);
	}

	// lock file
	//#ifndef WIN32
	//flock(fileno(fd), LOCK_EX);
	//#endif

	// read header
	len = sizeof(struct rpcHdr_t);
	if(fread(hdr,len,1,fd) != 1){
		printf("rpcReadMsg - error reading header message\n");
		exit(-1);
	}

	// allocate memory for arguments
	if((buffer=calloc(1,hdr->len))==NULL){
		printf("rpcReadMsg - error cannot allocate memory\n");
		exit(-1);
	}

	// read from file
	if(fread(buffer,1,hdr->len,fd) != hdr->len){
		printf("rpcReadMsg - error reading reply\n");
		exit(-1);
	}

	// release lock
	//#ifndef WIN32
	//flock(fileno(fd), LOCK_UN);
	//#endif

	// close file
	if(fclose(fd) != 0){
		printf("rpcReadMsg - error cannot close file\n");
		exit(-1);
	}

	return buffer;
}


// rpc_GTO_JK_Matrix_Quartet: calls the function GTO_JK_Matrix_Quartet remotely.
//
// If status is RPC_IDLE it sends message to child to perform calculation
//                       then return the message id for future message matching.
// If status is INT > 0  it waits for the reply message from the child
//                       and return RPC_DONE if successful
//
// May 6, 2016 - Teepanis Chachiyo
//     Use JT,KA,KB instead of GA,GB
// 
// Mar 2, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
int rpcMsgID=0;
int rpc_GTO_JK_Matrix_Quartet_Parallel(
	int status,                    // current rpc status
	int childID,                   // childID to call
	int nBasis,                    // number of basis functions
	const double *PA,              // density matrix for spin up
	const double *PB,              // density matrix for spin down 
	const struct GTOBasis_t *gto,  // basis set info
	const double *schwarz_basis,   // pointer to schwarz matrix
	double fixedCutoff,            // cutoff to ignore
	double *JT,                    // return J total
	double *KA,                    // return K for spin up
	double *KB,                    // return K for spin down
	struct option_t *opt){         // global option

	void *buffer;   // argument buffer
	void *reBuffer; // returned buffer
	int len;        // size of buffer
	char *ptr;      // generic pointer

	char msgParent2Child[256];    // message from parent to child file name
	char msgChild2Parent[256];    // message from child to parent file name 
	struct rpcHdr_t hdr;          // message header

#define BUILD_MSG_FILENAME(prefix,id)                    \
	sprintf(msgParent2Child,"%s_P2C%.4d.msg",prefix,id); \
	sprintf(msgChild2Parent,"%s_C2P%.4d.msg",prefix,id);

	// construct file name
	BUILD_MSG_FILENAME(opt->prefixStr,childID);

	// send message to call remotely
	if(status==RPC_IDLE){

		// construct arguments buffer
		buffer = pack_GTO_JK_Matrix_Quartet_Parallel(childID,nBasis,PA,PB,
		                                             gto,schwarz_basis,
		                                             fixedCutoff,JT,KA,KB,opt,&len);

		// compose header
		rpcMsgID++;
		hdr.id   = rpcMsgID;
		hdr.type = RPC_GTOJKMATRIXQUARTET;
		hdr.len  = len;

		// status becomes the message id
		status = hdr.id;

		// send message to child
		rpcWriteMsg(msgParent2Child,&hdr,buffer);

		// clean memory
		free(buffer);
	}

	// wait for reply until timeout
	if(status != RPC_IDLE && status != RPC_DONE){

		// wait until timeout
		delay();

		// check incoming message
		if(access(msgChild2Parent,F_OK)==0){
		
			// delay to avoid parent/child conflict
			delay();

			// read message from child
			reBuffer = rpcReadMsg(msgChild2Parent,&hdr);

			// delete the file
			delay();
			if(unlink(msgChild2Parent) != 0){
				printf("rpc_GTO_JK_Matrix_Quartet - error cannot delete %s\n",msgChild2Parent);
				exit(-1);
			}
			delay();

			// check if we get the reply
			if(hdr.id != status){
				free(reBuffer);
				return status;
			}else
				status=RPC_DONE;

			// copy the returned values
			ptr = (char *)reBuffer;
			len = nBasis*nBasis*sizeof(double); memcpy(JT,ptr,len); ptr += len;
			                                    memcpy(KA,ptr,len); ptr += len;
			                                    memcpy(KB,ptr,len); ptr += len;

			// clean memory
			free(reBuffer);
		}
	}

	return status;
}


// pack_mp2_rhf_aqij_Parallel_Contribute: packs the arguments for this subroutine
// into a long string so that it can be transferred over network/file system.
//
// Mar 6, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void * pack_mp2_rhf_aqij_Parallel_Contribute(
	int a,                         // starting orbital index
	int maxCorr,                   // number of correlated orbitals
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const double *Schwarz,         // schwarz inequality
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt,    // global options
	int *nBytes){                  // size of the buffer

	int len=0;       // size of the buffer
	void *gtoBuffer; // buffer for gto
	int gtoLen;      // size of the gto buffer
	void *molBuffer; // buffer to mol
	int molLen;      // size of the mol buffer
	void *buffer;    // buffer
	char *ptr;       // generic pointer

	// add a,maxCorr,nBasis,nOcc,  e,C,  Schwarz
	len += 4*sizeof(int)
	       + nBasis*sizeof(double) + nBasis*nBasis*sizeof(double)
	       + nBasis*nBasis*sizeof(double);

	// add gto
	gtoBuffer = pack_GTOBasis(nBasis,gto,&gtoLen);
	len += gtoLen;

	// add mol
	molBuffer = pack_Molecule(mol,&molLen);
	len += molLen;

	// add opt
	len += sizeof(struct option_t);

	// set the number of bytes
	*nBytes = len;

	// allocate memory
	if((buffer=calloc(1,len))==NULL){
		printf("pack_mp2_rhf_aqij_Parallel_Contribute - error cannot allocate memory\n");
		exit(-1);
	}

	// initialize pointer
	ptr = (char *)buffer;

	// copy data
	len = sizeof(int);                  memcpy(ptr,&a,len);       ptr += len; // a
	len = sizeof(int);                  memcpy(ptr,&maxCorr,len); ptr += len; // maxCorr
	len = sizeof(int);                  memcpy(ptr,&nBasis,len);  ptr += len; // nBasis
	len = sizeof(int);                  memcpy(ptr,&nOcc,len);    ptr += len; // nOcc
	len = nBasis*sizeof(double);        memcpy(ptr,e,len);        ptr += len; // e
	len = nBasis*nBasis*sizeof(double); memcpy(ptr,C,len);        ptr += len; // C
	len = nBasis*nBasis*sizeof(double); memcpy(ptr,Schwarz,len);  ptr += len; // Schwarz
	len = gtoLen;                       memcpy(ptr,gtoBuffer,len);ptr += len; // gto
	len = molLen;                       memcpy(ptr,molBuffer,len);ptr += len; // mol
	len = sizeof(struct option_t);      memcpy(ptr,opt,len);      ptr += len; // opt

	// clean memory
	free(gtoBuffer);
	free(molBuffer);

	return buffer;
}


// lpc_mp2_rhf_aqij_Parallel_Contribute: calls the function  
// mp2_rhf_aqij_Parallel_Contribute locally. It then packs the returned mp2
// contribution as the reply
//
// Mar 6, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void * lpc_mp2_rhf_aqij_Parallel_Contribute(
	void *buffer,                  // buffer containing arguments
	int *nBytes){                  // size of the buffer

	// arguments
	int a;                   // starting orbital index
	int maxCorr;             // number of correlated orbitals
	int nBasis;              // number of basis function
	int nOcc;                // number of occupied orbitals
	double *e;               // eigen values
	double *C;               // molecular orbitals
	double *Schwarz;         // schwarz inequality
	struct GTOBasis_t *gto;  // basis function structure
	struct Molecule_t *mol;  // molecule information
	struct option_t *opt;    // global options

	char *ptr;               // generic pointer
	int len;                 // generic int
	void *reBuffer;          // returned buffer
	double mp2;              // mp2 contribution

	// initialize pointer
	ptr = buffer;

	// unpack a
	len = sizeof(int); memcpy(&a,ptr,len); ptr += len;

	// unpack maxCorr
	len = sizeof(int); memcpy(&maxCorr,ptr,len); ptr += len;

	// unpack nBasis
	len = sizeof(int); memcpy(&nBasis,ptr,len); ptr += len;

	// unpack nOcc
	len = sizeof(int); memcpy(&nOcc,ptr,len); ptr += len;

	// unpack e
	len = nBasis*sizeof(double); e = (double *)ptr; ptr += len;

	// unpack C
	len = nBasis*nBasis*sizeof(double); C = (double *)ptr; ptr += len;

	// unpack Schwarz
	len = nBasis*nBasis*sizeof(double); Schwarz = (double *)ptr; ptr += len;

	// unpack gto
	gto = unpack_GTOBasis(nBasis,ptr,&len); ptr += len;

	// unpack mol
	mol = unpack_Molecule(ptr,&len); ptr += len;

	// unpack opt
	len = sizeof(struct option_t); opt= (struct option_t *)ptr; ptr += len;

	// call the actual function
	mp2 = mp2_rhf_aqij_Parallel_Contribute(a,maxCorr,nBasis,nOcc,e,C,Schwarz,gto,mol,opt);

	// allocate returned buffer
	reBuffer = calloc(1,sizeof(double));
	if(reBuffer==NULL){
		printf("lpc_mp2_rhf_aqij_Parallel_Contribute - cannot allocate memory\n");
		exit(-1);
	}
	*nBytes = sizeof(double);

	// pack mp2 value
	ptr = (char *)reBuffer;
	len = sizeof(double); memcpy(ptr,&mp2,len); ptr += len;

	// clean memory
	cleanGTOBasis(gto,nBasis);
	cleanMolecule(mol);

	return reBuffer;
}


// rpc_mp2_rhf_aqij_Parallel_Contribute: calls the mp2_rhf_aqij_Parallel_Contribute
// remotely. The rpc calling mechanism is similar to the
// rpc_GTO_JK_Matrix_Quartet_Parallel subroutine
//
// Mar 6, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
int rpc_mp2_rhf_aqij_Parallel_Contribute(
	int status,                    // current rpc status
	int childID,                   // childID to call
	double *mp2,                   // returned contribution
	int a,                         // starting orbital index
	int maxCorr,                   // number of correlated orbitals
	int nBasis,                    // number of basis function
	int nOcc,                      // number of occupied orbitals
	const double *e,               // eigen values
	const double *C,               // molecular orbitals
	const double *Schwarz,         // schwarz inequality
	const struct GTOBasis_t *gto,  // basis function structure
	const struct Molecule_t *mol,  // molecule information
	const struct option_t *opt){   // global options

	void *buffer;
	void *reBuffer;
	int len;

	char msgParent2Child[256];    // message from parent to child file name
	char msgChild2Parent[256];    // message from child to parent file name 
	struct rpcHdr_t hdr;          // message header

#define BUILD_MSG_FILENAME(prefix,id)                    \
	sprintf(msgParent2Child,"%s_P2C%.4d.msg",prefix,id); \
	sprintf(msgChild2Parent,"%s_C2P%.4d.msg",prefix,id);

	// construct file name
	BUILD_MSG_FILENAME(opt->prefixStr,childID);

	// send message to call remotely
	if(status==RPC_IDLE){

		// construct arguments buffer
		buffer = pack_mp2_rhf_aqij_Parallel_Contribute(a,maxCorr,nBasis,nOcc,
		                                               e,C,Schwarz,gto,mol,opt, &len);

		// compose header
		rpcMsgID++;
		hdr.id   = rpcMsgID;
		hdr.type = RPC_MP2RHFAQIJCONTRIBUTE;
		hdr.len  = len;

		// status becomes the message id
		status = hdr.id;

		// send message to child
		rpcWriteMsg(msgParent2Child,&hdr,buffer);

		// clean memory
		free(buffer);
	}

	// wait for reply until timeout
	if(status != RPC_IDLE && status != RPC_DONE){

		// wait until timeout
		delay();

		// check incoming message
		if(access(msgChild2Parent,F_OK)==0){
		
			// delay to avoid parent/child conflict
			delay();

			// read message from child
			reBuffer = rpcReadMsg(msgChild2Parent,&hdr);

			// delete the file
			delay();
			if(unlink(msgChild2Parent) != 0){
				printf("rpc_mp2_rhf_aqij_Parallel_Contribute - error cannot delete %s\n",msgChild2Parent);
				exit(-1);
			}
			delay();

			// check if we get the reply
			if(hdr.id != status){
				free(reBuffer);
				return status;
			}else
				status=RPC_DONE;

			// copy the returned values
			*mp2 = ((double *)reBuffer)[0];

			// clean memory
			free(reBuffer);
		}
	}

	return status;
}

 
// pack_GradEri_ShellSet_Parallel: packs the arguments for this subroutine
// into a long string so that it can be transferred over network/file system.
//
// May 11, 2016 - Teepanis Chachiyo
//     Add Kfactor to adjust exchange weight
//
// Mar 8, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void * pack_GradEri_ShellSet_Parallel(
	int nAtom,                      // number of atoms
	int childID,                    // this child id
	int nCPU,                       // number of CPUs
	int nBasis,                     // number of basis function
	const struct GTOBasis_t *gto,   // basis function database
	const int *basis2Atom,          // basis to atom mapping
	const double *PT,               // total density matrix
	const double *PS,               // spin density matrix
	const double Kfactor,          // exchange weight [0,1]
	double *Gx,                     // returned gradient in x direction
	double *Gy,                     // returned gradient in y direction
	double *Gz,                     // returned gradient in z direction
	int *nBytes){                   // size of the buffer

	int len=0;       // size of the buffer
	void *gtoBuffer; // buffer for gto
	int gtoLen;      // size of the gto buffer
	void *buffer;    // buffer
	char *ptr;       // generic pointer

	// add nAtom,childID,nCPU,nBasis
	len += 4*sizeof(int);

	// add gto
	gtoBuffer = pack_GTOBasis(nBasis,gto,&gtoLen);
	len += gtoLen;

	// add basis2Atom,  PT,PS, Kfactor, Gx,Gy,Gz
	len += nBasis*sizeof(int)
	       + 2*nBasis*nBasis*sizeof(double)
	       + sizeof(double)
	       + 3*nAtom*sizeof(double);

	// set the number of bytes
	*nBytes = len;

	// allocate memory
	if((buffer=calloc(1,len))==NULL){
		printf("pack_GradEri_ShellSet_Parallel - error cannot allocate memory\n");
		exit(-1);
	}

	// initialize pointer
	ptr = (char *)buffer;

	// copy data
	len = sizeof(int);                  memcpy(ptr,&nAtom,len);    ptr += len; // nAtom
	len = sizeof(int);                  memcpy(ptr,&childID,len);  ptr += len; // childID
	len = sizeof(int);                  memcpy(ptr,&nCPU,len);     ptr += len; // nCPU
	len = sizeof(int);                  memcpy(ptr,&nBasis,len);   ptr += len; // nBasis
	len = gtoLen;                       memcpy(ptr,gtoBuffer,len); ptr += len; // gto
	len = nBasis*sizeof(int);           memcpy(ptr,basis2Atom,len);ptr += len; // basis2Atom
	len = nBasis*nBasis*sizeof(double); memcpy(ptr,PT,len);        ptr += len; // PT
	len = nBasis*nBasis*sizeof(double); memcpy(ptr,PS,len);        ptr += len; // PS
	len = sizeof(double);               memcpy(ptr,&Kfactor,len);  ptr += len; // Kfactor
	len = nAtom*sizeof(double);         memcpy(ptr,Gx,len);        ptr += len; // Gx
	len = nAtom*sizeof(double);         memcpy(ptr,Gy,len);        ptr += len; // Gy
	len = nAtom*sizeof(double);         memcpy(ptr,Gz,len);        ptr += len; // Gz

	// clean memory
	free(gtoBuffer);

	return buffer;
}


// lpc_GradEri_ShellSet_Parallel: calls the function  
// GradEri_ShellSet_Parallel locally. It then packs the returned Gx,Gy,Gz
// as the reply
//
// May 11, 2016 - Teepanis Chachiyo
//     Add Kfactor to adjust exchange weight
//
// Mar 8, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void * lpc_GradEri_ShellSet_Parallel(
	void *buffer,                  // buffer containing arguments
	int *nBytes){                  // size of the buffer

	// arguments
	int nAtom;                // number of atoms
	int childID;              // this child id
	int nCPU;                 // number of CPUs
	int nBasis;               // number of basis function
	struct GTOBasis_t *gto;   // basis function database
	int *basis2Atom;          // basis to atom mapping
	double *PT;               // total density matrix
	double *PS;               // spin density matrix
	double Kfactor;           // exchange weight [0,1]
	double *Gx;               // returned gradient in x direction
	double *Gy;               // returned gradient in y direction
	double *Gz;               // returned gradient in z direction

	char *ptr;               // generic pointer
	int len;                 // generic int
	void *reBuffer;          // returned buffer

	// initialize pointer
	ptr = buffer;

	// unpack nAtom
	len = sizeof(int); memcpy(&nAtom,ptr,len); ptr += len;

	// unpack childID
	len = sizeof(int); memcpy(&childID,ptr,len); ptr += len;

	// unpack nCPU
	len = sizeof(int); memcpy(&nCPU,ptr,len); ptr += len;

	// unpack nBasis
	len = sizeof(int); memcpy(&nBasis,ptr,len); ptr += len;

	// unpack gto
	gto = unpack_GTOBasis(nBasis,ptr,&len); ptr += len;

	// unpack basis2Atom
	len = nBasis*sizeof(int); basis2Atom = (int *)ptr; ptr += len;

	// unpack PT
	len = nBasis*nBasis*sizeof(double); PT = (double *)ptr; ptr += len;

	// unpack PS
	len = nBasis*nBasis*sizeof(double); PS = (double *)ptr; ptr += len;

	// unpack Kfactor
	len = sizeof(double); memcpy(&Kfactor,ptr,len); ptr += len;

	// unpack Gx
	len = nAtom*sizeof(double); Gx = (double *)ptr; ptr += len;

	// unpack Gy
	len = nAtom*sizeof(double); Gy = (double *)ptr; ptr += len;

	// unpack Gz
	len = nAtom*sizeof(double); Gz = (double *)ptr; ptr += len;

	// call the actual function
	GradEri_ShellSet_Parallel(childID,nCPU,nBasis,gto,basis2Atom,PT,PS,Kfactor,Gx,Gy,Gz);

	// allocate returned buffer
	reBuffer = calloc(3*nAtom,sizeof(double));
	if(reBuffer==NULL){
		printf("lpc_GradEri_ShellSet_Parallel - cannot allocate memory\n");
		exit(-1);
	}
	*nBytes = 3*nAtom*sizeof(double);

	// set pointer preparing to pack values
	ptr = (char *)reBuffer;

	// pack Gx,Gy,Gz value
	len = nAtom*sizeof(double); memcpy(ptr,Gx,len); ptr += len;
	len = nAtom*sizeof(double); memcpy(ptr,Gy,len); ptr += len;
	len = nAtom*sizeof(double); memcpy(ptr,Gz,len); ptr += len;

	// clean memory
	cleanGTOBasis(gto,nBasis);

	return reBuffer;
}


// rpc_GradEri_ShellSet_Parallel: calls the subroutine GradEri_ShellSet_Parallel
// remotely.
//
// March 8, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
int rpc_GradEri_ShellSet_Parallel(
	int status,                     // current rpc status
	int nAtom,                      // number of atoms
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
	double *Gz,                     // returned gradient in z direction
	const struct option_t *opt){    // global options

	void *buffer;   // argument buffer
	void *reBuffer; // returned buffer
	int len;        // size of buffer
	char *ptr;      // generic pointer

	char msgParent2Child[256];    // message from parent to child file name
	char msgChild2Parent[256];    // message from child to parent file name 
	struct rpcHdr_t hdr;          // message header

#define BUILD_MSG_FILENAME(prefix,id)                    \
	sprintf(msgParent2Child,"%s_P2C%.4d.msg",prefix,id); \
	sprintf(msgChild2Parent,"%s_C2P%.4d.msg",prefix,id);

	// construct file name
	BUILD_MSG_FILENAME(opt->prefixStr,childID);

	// send message to call remotely
	if(status==RPC_IDLE){

		// construct arguments buffer
		buffer = pack_GradEri_ShellSet_Parallel(nAtom,childID,nCPU,nBasis,gto,
		                                        basis2Atom,PT,PS,Kfactor,Gx,Gy,Gz,&len);

		// compose header
		rpcMsgID++;
		hdr.id   = rpcMsgID;
		hdr.type = RPC_GRADERISHELLSET;
		hdr.len  = len;

		// status becomes the message id
		status = hdr.id;

		// send message to child
		rpcWriteMsg(msgParent2Child,&hdr,buffer);

		// clean memory
		free(buffer);
	}

	// wait for reply until timeout
	if(status != RPC_IDLE && status != RPC_DONE){

		// wait until timeout
		delay();

		// check incoming message
		if(access(msgChild2Parent,F_OK)==0){
		
			// delay to avoid parent/child conflict
			delay();

			// read message from child
			reBuffer = rpcReadMsg(msgChild2Parent,&hdr);

			// delete the file
			delay();
			if(unlink(msgChild2Parent) != 0){
				printf("rpc_GradEri_ShellSet_Parallel - error cannot delete %s\n",msgChild2Parent);
				exit(-1);
			}
			delay();

			// check if we get the reply
			if(hdr.id != status){
				free(reBuffer);
				return status;
			}else
				status=RPC_DONE;

			// copy the returned values
			ptr = (char *)reBuffer;
			len = nAtom*sizeof(double); memcpy(Gx,ptr,len); ptr += len;
			len = nAtom*sizeof(double); memcpy(Gy,ptr,len); ptr += len;
			len = nAtom*sizeof(double); memcpy(Gz,ptr,len); ptr += len;

			// clean memory
			free(reBuffer);
		}
	}

	return status;
}


void * lpc_exciteSubMatrix_Parallel(
	void *buffer,                  // buffer containing arguments
	int *nBytes){                  // size of the buffer

	// arguments
	int childID;             // child id number
	int isSinglet;           // singlet flag
	int Ni, Nf;              // starting and ending excitation
	int nBasis;              // number of basis function
	int *occ;                // occ index of size nBasis**2
	int *vir;                // vir index of size nBasis**2
	double *C;               // molecular orbital
	struct GTOBasis_t *gto;  // basis set info
	double *Schwarz;         // pointer to schwarz matrix
	double *subH;            // returned subH of size (Nf-Ni+1)**2
	struct option_t *opt;    // global option

	char *ptr;               // generic pointer
	int len;                 // generic int
	void *reBuffer;          // returned buffer

	// initialize pointer
	ptr = buffer;

	// unpack childID
	len = sizeof(int); memcpy(&childID,ptr,len); ptr += len;

	// unpack isSinglet
	len = sizeof(int); memcpy(&isSinglet,ptr,len); ptr += len;

	// unpack Ni
	len = sizeof(int); memcpy(&Ni,ptr,len); ptr += len;

	// unpack Nf
	len = sizeof(int); memcpy(&Nf,ptr,len); ptr += len;

	// unpack nBasis
	len = sizeof(int); memcpy(&nBasis,ptr,len); ptr += len;

	// unpack occ, vir
	len = nBasis*nBasis*sizeof(int);
	occ = (int *)ptr; ptr += len;
	vir = (int *)ptr; ptr += len;

	// unpack C
	len = nBasis*nBasis*sizeof(double);
	C = (double *)ptr; ptr += len;

	// unpack gto
	gto = unpack_GTOBasis(nBasis,ptr,&len); ptr += len;

	// unpack Schwarz
	len = nBasis*nBasis*sizeof(double); Schwarz = (double *)ptr; ptr += len;

	// unpack subH,opt
	len = (Nf-Ni+1)*(Nf-Ni+1)*sizeof(double); subH = (double *)ptr;          ptr += len;
	len = sizeof(struct option_t);            opt  = (struct option_t *)ptr; ptr += len;

	// call the actual function
	//exciteSubMatrix_Parallel(childID, isSinglet, Ni, Nf, 
	//                         nBasis, occ, vir, C, gto, Schwarz, subH, opt);

	// pack the returned subH value
	if((reBuffer=calloc((Nf-Ni+1)*(Nf-Ni+1),sizeof(double)))==NULL){
		printf("lpc_exciteSubMatrix_Parallel - error cannot allocate memory\n");
		exit(-1);
	}
	ptr = (char *)reBuffer;
	len = (Nf-Ni+1)*(Nf-Ni+1)*sizeof(double); memcpy(ptr,subH,len); ptr += len;

	// set the size of reBuffer
	*nBytes = (Nf-Ni+1)*(Nf-Ni+1)*sizeof(double);

	// clean memory
	cleanGTOBasis(gto,nBasis);

	return reBuffer;
}


void * pack_exciteSubMatrix_Parallel(
	int childID,                   // child id number
	int isSinglet,                 // singlet flag
	int Ni, int Nf,                // starting and ending excitation
	int nBasis,                    // number of basis function
	const int *occ,                // occ index of size nBasis**2
	const int *vir,                // vir index of size nBasis**2
	const double *C,               // molecular orbital
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double *subH,                  // returned subH of size (Nf-Ni+1)**2
	const struct option_t *opt,    // global option
	int *nBytes){                  // size of the buffer
	
	int len=0;       // size of the buffer
	void *gtoBuffer; // buffer for gto
	int gtoLen;      // size of the gto buffer
	void *buffer;    // buffer
	char *ptr;       // generic pointer

	// add childID, isSinglet, Ni, Nf, nBasis, occ, vir, C
	len += 5*sizeof(int) 
	       + 2*nBasis*nBasis*sizeof(int) 
	       + nBasis*nBasis*sizeof(double);

	// add gto
	gtoBuffer = pack_GTOBasis(nBasis,gto,&gtoLen);
	len += gtoLen;

	// add Schwarz, subH, opt
	len += nBasis*nBasis*sizeof(double) 
	       + (Nf-Ni+1)*(Nf-Ni+1)*sizeof(double)
	       + sizeof(struct option_t);

	// set the number of bytes
	*nBytes = len;

	// allocate memory
	if((buffer=calloc(1,len))==NULL){
		printf("exciteSubMatrix_Parallel - error cannot allocate memory\n");
		exit(-1);
	}

	// initialize pointer
	ptr = (char *)buffer;

	// copy childID
	len = sizeof(int); memcpy(ptr,&childID,len); ptr += len;

	// copy isSinglet
	len = sizeof(int); memcpy(ptr,&isSinglet,len); ptr += len;

	// copy Ni
	len = sizeof(int); memcpy(ptr,&Ni,len); ptr += len;

	// copy Nf
	len = sizeof(int); memcpy(ptr,&Nf,len); ptr += len;

	// copy nBasis
	len = sizeof(int); memcpy(ptr,&nBasis,len); ptr += len;

	// copy occ vir
	len = nBasis*nBasis*sizeof(int);
	memcpy(ptr,occ,len); ptr += len;
	memcpy(ptr,vir,len); ptr += len;

	// copy C
	len = nBasis*nBasis*sizeof(double);
	memcpy(ptr,C,len); ptr += len;

	// copy gto
	len = gtoLen; memcpy(ptr,gtoBuffer,len); ptr += len;

	// copy Schwarz
	len = nBasis*nBasis*sizeof(double);
	memcpy(ptr,Schwarz,len); ptr += len;

	// copy subH
	len = (Nf-Ni+1)*(Nf-Ni+1)*sizeof(double);
	memcpy(ptr,subH,len); ptr += len;

	// copy opt
	len = sizeof(struct option_t); memcpy(ptr,opt,len); ptr += len;

	// clean memory
	free(gtoBuffer);

	return buffer;
}


int rpc_exciteSubMatrix_Parallel(
	int status,                    // current rpc status
	int childID,                   // childID to call
	int isSinglet,                 // singlet flag
	int Ni, int Nf,                // starting and ending excitation
	int nBasis,                    // number of basis function
	const int *occ,                // occ index of size nBasis**2
	const int *vir,                // vir index of size nBasis**2
	const double *C,               // molecular orbital
	const struct GTOBasis_t *gto,  // basis set info
	const double *Schwarz,         // pointer to schwarz matrix
	double *subH,                  // returned subH of size (Nf-Ni+1)**2
	const struct option_t *opt){   // global option

	void *buffer;   // argument buffer
	void *reBuffer; // returned buffer
	int len;        // size of buffer
	char *ptr;      // generic pointer

	char msgParent2Child[256];    // message from parent to child file name
	char msgChild2Parent[256];    // message from child to parent file name 
	struct rpcHdr_t hdr;          // message header

#define BUILD_MSG_FILENAME(prefix,id)                    \
	sprintf(msgParent2Child,"%s_P2C%.4d.msg",prefix,id); \
	sprintf(msgChild2Parent,"%s_C2P%.4d.msg",prefix,id);

	// construct file name
	BUILD_MSG_FILENAME(opt->prefixStr,childID);

	// send message to call remotely
	if(status==RPC_IDLE){

		// construct arguments buffer
		buffer = pack_exciteSubMatrix_Parallel(childID,isSinglet, Ni, Nf,
		                                       nBasis, occ, vir, C, gto,
		                                       Schwarz,subH,opt,&len);

		// compose header
		rpcMsgID++;
		hdr.id   = rpcMsgID;
		hdr.type = RPC_EXCITESUBMATRIX;
		hdr.len  = len;

		// status becomes the message id
		status = hdr.id;

		// send message to child
		rpcWriteMsg(msgParent2Child,&hdr,buffer);

		// clean memory
		free(buffer);
	}

	// wait for reply until timeout
	if(status != RPC_IDLE && status != RPC_DONE){

		// wait until timeout
		delay();

		// check incoming message
		if(access(msgChild2Parent,F_OK)==0){
		
			// delay to avoid parent/child conflict
			delay();

			// read message from child
			reBuffer = rpcReadMsg(msgChild2Parent,&hdr);

			// delete the file
			delay();
			if(unlink(msgChild2Parent) != 0){
				printf("rpc_exciteSubMatrix_Parallel - error cannot delete %s\n",msgChild2Parent);
				exit(-1);
			}
			delay();

			// check if we get the reply
			if(hdr.id != status){
				free(reBuffer);
				return status;
			}else
				status=RPC_DONE;

			// copy the returned values
			ptr = (char *)reBuffer;
			len = (Nf-Ni+1)*(Nf-Ni+1)*sizeof(double);
			memcpy(subH,ptr,len); ptr += len;

			// clean memory
			free(reBuffer);
		}
	}

	return status;
}


// trapChild: trap the arguments passed to Siam Quantum for the child 
// process. If it is so, the processs will be stuck here indefinitely.
//
// Mar 2, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void trapChild(int argc, char *argv[]){

	int childID;                  // child identification number
	char prefixStr[256];          // job prefix string
	char msgParent2Child[256];    // message from parent to child file name
	char msgChild2Parent[256];    // message from child to parent file name 
	struct rpcHdr_t hdr;          // message header
	int len;                      // generic int
	void *buffer, *reBuffer;      // buffer
	int lastID=-1;                // previous message id

	// two arguments are required here: -ID=INT and -PREFIX=STR
	if(argc==3){

		// load child id
		if(strncmp(argv[1],"-ID=",4)==0){
			childID = atoi(argv[1]+4);
		}else return;

		// load prefix string
		if(strncmp(argv[2],"-PREFIX=",8)==0){
			strcpy(prefixStr, argv[2]+8);
		}else return;

		// validation
		if(childID <= 0){
			printf("trapChild - error invalid child ID\n");
			exit(-1);
		}

		////////////////////////////
		// the trap is successful
		////////////////////////////

		// construct file name
		BUILD_MSG_FILENAME(prefixStr,childID);

		while(1){
			// check incoming message
			if(access(msgParent2Child,F_OK)==0){

				// delay to avoid parent/child conflict
				delay();

				// read message from parent
				buffer = rpcReadMsg(msgParent2Child, &hdr);

				// delete the file
				delay();
				if(unlink(msgParent2Child) != 0){
					printf("trapChild - error cannot delete %s\n",msgParent2Child);
					exit(-1);
				}
				delay();

				// check if the same message
				if(hdr.id == lastID) { free(buffer); continue; } else lastID = hdr.id;

				// handle message
				switch(hdr.type){

				case RPC_CHILDEXIT: exit(0); break;

				case RPC_GTOJKMATRIXQUARTET:

					// call the subroutine locally
					reBuffer = lpc_GTO_JK_Matrix_Quartet_Parallel(buffer,&len);

					// compose header
					hdr.id   = lastID;
					hdr.type = RPC_GTOJKMATRIXQUARTET;
					hdr.len  = len;

					// send message to parent
					rpcWriteMsg(msgChild2Parent,&hdr,reBuffer);

					// clean memory
					free(buffer);
					free(reBuffer);

				break;

				case RPC_MP2RHFAQIJCONTRIBUTE:

					// call the subroutine locally
					reBuffer = lpc_mp2_rhf_aqij_Parallel_Contribute(buffer, &len);

					// compose header
					hdr.id   = lastID;
					hdr.type = RPC_MP2RHFAQIJCONTRIBUTE;
					hdr.len  = len;

					// send message to parent
					rpcWriteMsg(msgChild2Parent,&hdr,reBuffer);

					// clean memory
					free(buffer);
					free(reBuffer);

				break;

				case RPC_GRADERISHELLSET:

					// call the subroutine locally
					reBuffer = lpc_GradEri_ShellSet_Parallel(buffer,&len);

					// compose header
					hdr.id   = lastID;
					hdr.type = RPC_GRADERISHELLSET;
					hdr.len  = len;

					// send message to parent
					rpcWriteMsg(msgChild2Parent,&hdr,reBuffer);

					// clean memory
					free(buffer);
					free(reBuffer);

				break;


				case RPC_EXCITESUBMATRIX:

					// call the subroutine locally
					reBuffer = lpc_exciteSubMatrix_Parallel(buffer,&len);

					// compose header
					hdr.id   = lastID;
					hdr.type = RPC_EXCITESUBMATRIX;
					hdr.len  = len;

					// send message to parent
					rpcWriteMsg(msgChild2Parent,&hdr,reBuffer);

					// clean memory
					free(buffer);
					free(reBuffer);

				break;


				default:
					printf("trapChild - error recieved unknown message\n");
					exit(-1);
				}
			}
			delay();
		}
	}
}


// rpcSpawnChildren: spawn children in preparation for parallel run.
//
// Mar 3, 2013 - Teepanis Chachiyo
//    Initial implementation and testing
//
void rpcSpawnChildren(int argc, char *argv[], struct option_t *opt){
	
	int n;                        // generic int
	char cmd[256];                // command string
	char msgParent2Child[256];    // message from parent to child file name
	char msgChild2Parent[256];    // message from child to parent file name 

	// trap single cpu run
	if(opt->nCPU < 2) return;

	// do not spawn childID=0
	for(n=1; n < opt->nCPU; n++){

		// delete left-over messaging files
		BUILD_MSG_FILENAME(opt->prefixStr,n);
		delay();
		unlink(msgParent2Child);
		unlink(msgChild2Parent);
		delay();

#ifdef WIN32
		sprintf(cmd,"START %s -ID=%d -PREFIX=%s",argv[0],n,opt->prefixStr);
#else
		sprintf(cmd,"%s -ID=%d -PREFIX=%s &",argv[0],n,opt->prefixStr);
#endif
		system(cmd);
	}
}


// rpcExitChildren: send exit signal to all children
//
// Mar 3, 2013 - Teepanis Chachiyo
//     Initial implementation and testing
//
void rpcExitChildren(struct option_t *opt){

	int n;                        // generic int
	struct rpcHdr_t hdr;          // message header
	char msgParent2Child[256];    // message from parent to child file name
	char msgChild2Parent[256];    // message from child to parent file name 

	// trap single cpu run
	if(opt->nCPU < 2) return;

	// send exit signal
	for(n=1; n < opt->nCPU; n++){

		// construct file name
		BUILD_MSG_FILENAME(opt->prefixStr,n);

		// compose header
		rpcMsgID++;
		hdr.id   = rpcMsgID;
		hdr.type = RPC_CHILDEXIT;
		hdr.len  = 0;

		// send message to child
		rpcWriteMsg(msgParent2Child,&hdr,NULL);
		delay();
	}
}


