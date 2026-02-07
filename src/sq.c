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
#include "util.h"
#include "mol.h"
#include "basis.h"
#include "rhf.h"
#include "uhf.h"
#include "option.h"
#include "xsf.h"
#include "grad.h"
#include "hessian.h"
#include "pop.h"
#include "optimize.h"
#include "check.h"
#include "mp2.h"
#include "rpc.h"
#include "mecp.h"
//#include "excite.h"
#include "qmd.h"

#define VERSION_NUMBER "2.0.0"

// greeting : Printout initial message at the beginning of the
// program. The message should be related to copyright issue
// and the credit to various contributors.
//
// Feb 21, 2008 - Teepanis Chachiyo
//     Initial implementation.
static void greeting(FILE *outFile){

	char time[256];

	time_str(256, time);

	fprintf(outFile,
	"                                                            \n"
	"+----------------------------------------------------------+\n"
	"|| Siam Quantum 2                                          |\n"
	"|| Teepanis Chachiyo <teepanisc@nu.ac.th>                  |\n"
	"|| Hathaithip Chachiyo <hathaithip.chachiyo@gmail.com>     |\n"
	"||                                                         |\n"
	"|| Please cite this work as:                               |\n"
	"||                                                         |\n"
	"|| Teepanis Chachiyo, Hathaithip Chachiyo;                 |\n" 
	"|| Siam Quantum 2: An open-source C toolbox for quantum    |\n"
	"|| modeling and electronic structure development.          |\n"
	"|| J. Chem. Phys. 7 February 2026; 164 (5): 052501.        |\n"
	"|| https://doi.org/10.1063/5.0310183                       |\n"
	"||                                                         |\n"
	"|| Version: %7s                                        |\n"
	"||                                                         |\n"
	"|| Supported by                                            |\n"
	"||      - Thailand Center of Excellence in Physics,        |\n"
	"||        Ministry of Higher Education, Science, Research  |\n"
	"||        and Innovation, Thailand                         |\n"
	"||      - Department of Physics, Faculty of Science,       |\n"
	"||        Naresuan University, Phitsanulok, Thailand       |\n"
	"||      - Department of Physics, Faculty of Science,       |\n"
	"||        Khon Kaen University, Khon Kaen, Thailand        |\n"
	"||                                                         |\n"
	"|| Author                                                  |\n"
	"|| [] Teepanis Chachiyo                                    |\n"
	"||      First prototype RHF                                |\n"
	"||      Unrestricted Hartree-Fock     (UHF)                |\n"
	"||      Density Functional Theory     (DFT)                |\n"
	"||      Moller-Plesset Perturbation   (MP2)                |\n"
	"||      Minimum Energy Crossing Point (MECP)               |\n"
	"||      Convergence Acceleration                           |\n"
	"||      Excited States                                     |\n"
	"||      Parallelization                                    |\n"
	"||      Geometry Optimization                              |\n"
	"||      Integral Memory Storage Scheme                     |\n"
	"||      Hessian and Molecular Vibrations                   |\n"
	"||      Fast Energy and Gradients Evaluation               |\n"
	"||      Chachiyo correlation (uniform electron gas)        |\n"
	"|| [] Hathaithip Chachiyo                                  |\n"
	"||      Chachiyo correlation (general electron gas)        |\n"
	"||      Chachiyo exchange    (general electron gas)        |\n"
	"||                                                         |\n"
	"|| Contributors                                            |\n"
	"|| [] Theerapon Khamla                                     |\n"
	"||      First prototype gradient of energy                 |\n"
	"||      First prototype geometry optimization              |\n"
	"||      Exporting densities for visualization              |\n"
	"|| [] Keerati Maneesai  - Boys gamma function              |\n"
	"|| [] Nanta Sophonrat   - Electric field                   |\n"
	"|| [] Chutchawan Jaisuk - Molecular dynamics               |\n"
	"|| [] Aniwat Kesorn     - Dipole moment                    |\n"
	"||                                                         |\n"
	"|| Acknowledgement                                         |\n"
	"|| [] Narong Pootaddouwng                                  |\n"
	"|| [] Nawee  Jaroonchokanan                                |\n"
	"|| [] Labedev grid: from a C/C++ code by John Burkardt     |\n"
	"|| [] Becke grid: adapted from DFTLIBS by Radovan Bast     |\n"
	"|| [] DFT functional formulae: from LibXC                  |\n"
	"|| [] DFT functional formulae: monologue by Timothy Giesea |\n"
	"|| [] Integrating subroutine (first prototype only): from  |\n"
	"||    PyQuante written by Richard Muller                   |\n"
	"||                                                         |\n"
	"+----------------------------------------------------------+\n"
	"                                                            \n"
	" Copyright (C) 2008-2026 Teepanis Chachiyo and Hathaithip Chachiyo  \n"
	"                                                                    \n"
	"This program is free software; you can redistribute it and/or modify\n"
	"it under the terms of the GNU General Public License as published by\n"
	"the Free Software Foundation; either version 2, or (at your option) \n"
	"any later version.                                                  \n"
	"                                                                    \n"
	"                                                                    \n"
	"Begin: %s\n",
	VERSION_NUMBER, time);
}



// main : SQ main subroutine
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//  Add CUBE and XSF options
//
// Oct 20, 2010 - Theerapon Khumla
//  Add force section
//
// June 4, 2011 - Teepanis Chachiyo and Nanta Sophonrat
//  Add electrostatic properties
//
// Oct 13, 2012 - Teepanis Chachiyo and Aniwat Kesorn
//  Add electric multipole moment
//
// Feb 25, 2016 - Teepanis Chachiyo
//  switch to UHF for excited state
//
int main(int argc, char *argv[]){

	struct Molecule_t *mol=NULL;       // pointer to molecular structure info
	struct GTOBasis_t *gto=NULL;       // pointer to basis function
	struct option_t    opt;            // store user specifed options
	struct GTOBasisSet_t *dbSet=NULL;  // pointer to basis set information
	double *CA, *eA, *CB, *eB;         // molecular orbitals and eigen values
	int nBasis;                        // number of basis function
	int dbSetItem;                     // number of basis set item
	FILE *inFile;                      // xyz file pointer
	FILE *basisFile;                   // basis function file pointer
	double Etot;                       // total energy
	char time[256];

	// trap child process
	trapChild(argc, argv);

	greeting(stdout);

	if(argc < 3){
		printf("\n"
		       "usage: %s <xyz file>  <basis set file> [OPTIONS]\n"
		       "<xyz file>       -  Molecular structure in XYZ file format\n"
		       "<basis set file> -  Basis set database in GAMESS-US format\n"
		       "                    Many types of basis set are available \n"
		       "                    at http://bse.pnl.gov                 \n", 
		       argv[0]);
		option_help();
		exit(EXIT_FAILURE);
	}

	parse_option( &opt, argc, argv);

	// need to kill all left-over processes
	rpcExitChildren(&opt);	

	// spawn children for parallel run
	rpcSpawnChildren(argc, argv, &opt);

	printf(
	"                                                             \n"
	"                                                             \n"
	"*************************************************************\n"
	"*****                                                   *****\n"
	"*****         SECTION: STRUCTURE SPECIFICATION          *****\n"
	"*****                                                   *****\n"
	"*************************************************************\n"
	);

	// read molecule from checkpoint if requested
	if(opt.loadCheck){

		mol = load_checkpoint_molecule(&opt);

		printf("Reading molecule from checkpoint\n");
		printf("Detected %d atoms\n", mol->nAtom);
		printMolecule_XYZ(mol, stdout);

	// otherwise read the info from xyz file
	}else{

		inFile = fopen(argv[1],"r");
		if(inFile == NULL){
			printf("Cannot open file %s\n", argv[1]);
			exit(EXIT_FAILURE);
		}
		mol    = readMolecule_XYZ(inFile);
		mol->Q = opt.molCharge;
	}

	////////////////////////////////////
	// handle default for multiplicity
	////////////////////////////////////
	// even number of electron set multiplicity to one by default
	if(opt.multiplicity==0 && get_nElectron(mol)%2==0)
		opt.multiplicity = 1;
	// odd number of electron set multiplicity to two by default
	if(opt.multiplicity==0 && get_nElectron(mol)%2==1)
		opt.multiplicity = 2;

	// read multiplicity from checkpoint if requested
	if(opt.loadCheck)
		opt.multiplicity = load_checkpoint_multiplicity(&opt);

	// validate mecp spin multiplicity
	if(opt.MECP){
		get_nEA(mol,opt.mecpMA);
		get_nEB(mol,opt.mecpMA);
		get_nEA(mol,opt.mecpMB);
		get_nEB(mol,opt.mecpMB);
	}

	/////////////////////////////////////
	// handle default for RHF and UHF
	/////////////////////////////////////
	if(opt.multiplicity>1) opt.restricted = 0;
	if(opt.MECP)           opt.restricted = 0;

	// handle excited state
	if(opt.excite != 0){

		//validate
		if(opt.excite < 0){
			printf("main - error the number of excited states less than zero\n");
			exit(-1);
		}

		// restricted only 
		if(! opt.restricted){
			printf("main - error excitation only supports RHF\n");
			exit(-1);
		}

		if(opt.force){
			printf("main - error excitation does not support force calculation\n");
			exit(-1);
		}

		if(opt.opt){
			printf("main - error geometry optimization does not support excited states\n");
			exit(-1);
		}

		if(opt.qmd){
			printf("main - error QMD does not support excited states\n");
			exit(-1);
		}

		if(opt.MECP){
			printf("main - error excitation does not support MECP\n");
			exit(-1);
		}

		if(opt.method != METHOD_HF){
			printf("main - error excitation only support Hartree-Fock at the moment\n");
			exit(-1);
		}
	}

	// validate restricted mp2
	//if(opt.MP2 && opt.RHF==0){
	//	printf("main - error MP2 only supports RHF\n");
	//	exit(-1);
	//}

	printf(
	"                                                             \n"
	"                                                             \n"
	"*************************************************************\n"
	"*****                                                   *****\n"
	"*****         SECTION: BASIS SET REPRESENTATION         *****\n"
	"*****                                                   *****\n"
	"*************************************************************\n"
	);

	basisFile = fopen(argv[2],"r");
	if(basisFile == NULL){
		printf("Cannot open file %s\n", argv[2]);
		exit(EXIT_FAILURE);
	}
	dbSet = read_GAMESS_BasisSet(basisFile, argv[2], &dbSetItem);

	/////////////////////////////////////////////
	// handle geometry optimization if requested
	/////////////////////////////////////////////
	if(opt.opt){
		printf(
		"                                                             \n"
		"                                                             \n"
		"*************************************************************\n"
		"*****                                                   *****\n"
		"*****         SECTION: GEOMETRY OPTIMIZATION            *****\n"
		"*****                                                   *****\n"
		"*************************************************************\n"
		);
		fflush(stdout);
		optimize(dbSetItem, dbSet, mol, &opt);
	}

	/////////////////////////////////////////////
	// handle quantum molecular dynamics
	/////////////////////////////////////////////
	if(opt.qmd){
		printf(
		"                                                             \n"
		"                                                             \n"
		"*************************************************************\n"
		"*****                                                   *****\n"
		"*****         SECTION: QUANTUM MOLECULAR DYNAMICS       *****\n"
		"*****                                                   *****\n"
		"*************************************************************\n"
		);
		fflush(stdout);
		qmd(dbSetItem, dbSet, mol, &opt);

		// do not perform post SCF
		goto Clean2Exit;
	}


	////////////////////////////////////////////////////////////////////////
	// handle Minimum Energy Crossing Point (MECP) if requested
	// See. Chachiyo, Teepanis, and Jorge H. Rodriguez. 
	// "A direct method for locating minimum-energy crossing points (MECPs) 
	// in spin-forbidden transitions and nonadiabatic reactions." 
	// The Journal of chemical physics 123 (2005): 094711.
	////////////////////////////////////////////////////////////////////////
	if(opt.MECP){
		printf(
		"                                                             \n"
		"                                                             \n"
		"*************************************************************\n"
		"*****                                                   *****\n"
		"*****         SECTION: MECP CALCULATIONS                *****\n"
		"*****                                                   *****\n"
		"*************************************************************\n"
		);
		fflush(stdout);
		mecp(dbSetItem, dbSet, mol, &opt);

		// do not perform SCF nor Post-SCF because there are two states here
		goto Clean2Exit;
	}


	printf(
	"                                                             \n"
	"                                                             \n"
	"*************************************************************\n"
	"*****                                                   *****\n"
	"*****         SECTION: SCF CALCULATIONS                 *****\n"
	"*****                                                   *****\n"
	"*************************************************************\n"
	);
	fflush(stdout);

	// load basis from checkpoint if requested
	if(opt.loadCheck){

		gto = load_checkpoint_basis(&opt, &nBasis);
		printf("There are %d Cartesian-Gaussian contracted basis functions\n", nBasis);
		printf("There are %d Primitive GTO functions\n", get_nPrim(nBasis, gto));

	// otherwise generate basis function
	}else
		gto = genBasis(mol, &nBasis, dbSetItem, dbSet);

	// allocate meory for molecular orbitals and their eigen values
	CA = calloc(nBasis*nBasis,sizeof(double));
	eA = calloc(nBasis,sizeof(double));
	CB = calloc(nBasis*nBasis,sizeof(double));
	eB = calloc(nBasis,sizeof(double));
	if(CA==NULL || eA==NULL || CB==NULL || eB==NULL){
		printf("main: error - cannot allocate memory\n");
		exit(-1);
	}

	// calling restricted calculation is depreciated, use uhf with
	// option_t->RHF instead
	// rhf(nBasis, gto, mol, get_nElectron(mol), CA, eA, &opt);
	if(opt.loadCheck)
		Etot = load_checkpoint_orbital(nBasis, CA, CB, eA, eB, &opt);
	else
	    Etot = uhf(nBasis, gto, mol, 
	               get_nEA(mol,opt.multiplicity), get_nEB(mol,opt.multiplicity), 
	               CA, CB, eA, eB, &opt);

	// skip post-scf if convergence is not reached
	if(Etot==0.0) goto Clean2Exit;

	printf(
	"                                                             \n"
	"                                                             \n"
	"*************************************************************\n"
	"*****                                                   *****\n"
	"*****         SECTION: POST-SCF CALCULATIONS            *****\n"
	"*****                                                   *****\n"
	"*************************************************************\n"
	);

	if(opt.method==METHOD_MP2){
		if(opt.restricted)
		mp2_rhf_aqij_Parallel(nBasis, get_nEA(mol,opt.multiplicity), eA, CA, gto, mol, &opt);
		//mp2_rhf_aqij(nBasis, get_nEA(mol,opt.multiplicity), eA, CA, gto, mol, &opt);
		//mp2_rhf_aqbj(nBasis, get_nEA(mol,opt.multiplicity), eA, CA, gto, mol, &opt);
		//mp2_rhf_direct(nBasis, get_nEA(mol,opt.multiplicity), eA, CA, gto, mol, &opt);
		//mp2_rhf_semi_direct_aqbj(nBasis, get_nEA(mol,opt.multiplicity), eA, CA, gto, mol, &opt);
		//mp2_rhf_semi_direct_aqij(nBasis, get_nEA(mol,opt.multiplicity), eA, CA, gto, mol, &opt);
		else
		//mp2_uhf_direct(nBasis,
		//               get_nEA(mol,opt.multiplicity),
		//               get_nEB(mol,opt.multiplicity),
		//               eA, eB, CA, CB, gto, mol, &opt);
		//mp2_uhf_semi_direct_aqij(nBasis,
		//               get_nEA(mol,opt.multiplicity),
		//               get_nEB(mol,opt.multiplicity),
		//               eA, eB, CA, CB, gto, mol, &opt);
		mp2_uhf_aqij(nBasis,
		             get_nEA(mol,opt.multiplicity),
		             get_nEB(mol,opt.multiplicity),
		             eA, eB, CA, CB, gto, mol, &opt);


	}

	mulliken(nBasis, gto, mol, 
	         get_nEA(mol,opt.multiplicity), get_nEB(mol,opt.multiplicity),
	         CA, CB, eA, eB, &opt);

	electrostatic(nBasis, gto, mol, 
	         get_nEA(mol,opt.multiplicity), get_nEB(mol,opt.multiplicity),
	         CA, CB, eA, eB, &opt);

	electric_multipole(nBasis, gto, mol, 
	         get_nEA(mol,opt.multiplicity), get_nEB(mol,opt.multiplicity),
	         CA, CB, eA, eB, &opt);

	// extra info

	if(opt.outVolumeType != VOLUME_NONE){
		if(opt.outVolume)
		switch(opt.outVolumeFormat){
		case VOLUME_XSF:	xsfden(nBasis, gto, mol, CA, CB, 
							       0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
							       0, 0, 0,
							       &opt);
		break;

		case VOLUME_CUBE:	cubeden(nBasis, gto, mol, CA, CB, 
							        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
							        0, 0, 0,
							        &opt);
		break;
		}

		if(opt.outLine) lineden(nBasis, gto, mol, CA, CB, &opt);
	}

	if(opt.outGAUSSIAN)
		export_gaussian(nBasis, gto, mol, CA, CB, eA, eB, &opt);

	if(opt.outGAMESS)
		export_gamess(nBasis, gto, mol, CA, CB, eA, eB, &opt);

	// force calculation
	if(opt.force){
		uhf_force(nBasis, gto, mol, 
		      get_nEA(mol,opt.multiplicity), get_nEB(mol,opt.multiplicity), 
		      CA, CB, eA, eB, &opt, NULL, NULL, NULL);
	}

	// hessian calculation
	if(opt.hessian){
		hessian_numerical(
		      dbSetItem, dbSet,
		      nBasis, mol, 
		      get_nEA(mol,opt.multiplicity), get_nEB(mol,opt.multiplicity), 
		      CA, CB, eA, eB, &opt);
	}


	if(opt.excite != 0){
		printf("Excitation has been disabled pending publication\n");
		exit(-1);
		//MP1Excite_RHF(nBasis, gto, mol,
		//              get_nEA(mol,opt.multiplicity), CA, eA, &opt);
	}

Clean2Exit:

	// cleaning memory and allocations
	if(mol!=NULL) cleanMolecule(mol);
	if(gto!=NULL) cleanGTOBasis(gto, nBasis);

	// exit all children for parallel run
	rpcExitChildren(&opt);	

	time_str(256, time);
	printf("\nEnd: %s\n", time);

	// successful caculations
	return EXIT_SUCCESS;
}
