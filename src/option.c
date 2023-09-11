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

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "option.h"
#include "util.h"
#include "grid.h"
#include "dft.h"

#ifdef LIBXC
#include <xc.h>
#endif

// option_help : print out all available options
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// Oct 20, 2010 - Theerapon Khumla
//  Add force option
//
// Oct 22, 2010 - Teepanis Chachiyo
//  Add geometry optimization and SCFConv option
//
// May 21, 2011 - Teepanis Chachiyo
//  Add DIIS option and set it to the default option
//
// May 21, 2011 - Teepanis Chachiyo
//  Add maxmem option
//
// June 4, 2011 - Nanta Sophonrat
//  Add potential option
//
// Oct 4, 2012 - Teepanis Chachiyo
//  Add and DIIS checkpoint options
//
// Oct 20, 2012 - Teepanis Chachiyo
//  Add MP2 options
//
// Jan 26, 2013 - Teepanis Chachiyo
//  Add DIIS3 options and set default to DIIS4
//
// March, 17, 2013 - Teepanis Chachiyo
//  Change default of SCFDrag to 0.25
//
// Jan 2, 2016 - Teepanis Chachiyo
//  Correct typo for GabEdit
//
// Feb 25, 2016 - Teepanis Chachiyo
//  Add excitation options
//
// Mar 11, 2016 - Teepanis Chachiyo
//  Add uniform electric field and QMD options
//
// August 15, 2016 - Teepanis Chachiyo
//  Add Chachiyo's correlation as alternative to VWN5
//  Take Excitation out of source tree pending publication
//
// August 18, 2019 - Teepanis Chachiyo
//  Add options for Chachiyo GGA and LibXC supports
//  Add options for DFT integration grid
//
// March, 22, 2020 - Teepanis Chachiyo
//  Add Hessian options
//
void option_help(){
	printf(
	"[OPTIONS] :                                                    \n"
	"                                                               \n"
	"[x] Ab Initio Method:                                          \n"
	"-HF             Hartree-Fock (default)                         \n"
/*
	"-DFT=X:C        Density Functional Theory                      \n"
	"                                                               \n"
	"                    Available exchange keywords X are:         \n"
	"                    S        - (LDA) Dirac exchange            \n"
	"                    PBE      - (GGA) Perdew-B                  \n"
	"                    PY86     - (GGA) Perdew-Yue  1986 exchange \n"
	"                    PW91     - (GGA) Perdew-Wang 1991 exchange \n"
	"                    BECKE    - (GGA) Becke 1988 exchange       \n"
	"                    CHACHIYO - (GGA) The Chachiyo exchange     \n"
	"                    B3       - (HYB) Dirac+B88+HF exchange     \n"
	"                    HALF     - (HYB) 0.5*Dirac + 0.5*HF        \n"
	"                                                               \n"
	"                    Available correlation keywords C are:      \n"
	"                    VWN5     - (LDA) Volko-Wilk-Nusair         \n"
	"                    PW92     - (LDA) Perdew-Wang 1992          \n"
	"                    CHACHIYO - (LDA) Chachiyo 2016             \n"
	"                    PBE      - (GGA) Perdew-B ...              \n"
	"                    PW91     - (GGA) Perdew-Wang 1991          \n"
	"                                                               \n"
*/
	"                Density Functional Theory\n"
	"-LDA=S:VWN5       ex: Slater          corr: VWN5            \n"
	"                  use -VWN5_PAR1_A=REAL  -VWN5_PAR1_x=REAL  \n"
	"                      -VWN5_PAR1_b=REAL  -VWN5_PAR1_c=REAL  \n"
	"                  to override VWN5 parameters               \n"
	"                  (change PAR1 to PAR2 or PAR3 if needed)   \n"

	"-LDA=S:CHACHIYO   ex: Slater          corr: Chachiyo        \n"
	"-DFT=HALF         ex: 0.5*(HF+Slater) corr: VWN5            \n"
#ifdef LIBXC
	"-DFT=PBE          ex: PBE             corr: PBE             \n"
	"-DFT=BLYP         ex: Becke88         corr: LYP             \n"
#endif
	"-DFT=CHACHIYO     ex: Chachiyo        corr: Chachiyo        \n"
#ifdef LIBXC
	"-DFT=B3LYP        The famous Becke's three parameters hybrid\n"
	"-LIBXC=STR        Use LibXC with specific ex. and corr.\n"
	"                  See its website for the avialable functionals   \n"
	"                  Example, -LIBXC=GGA_X_CHACHIYO+GGA_C_CHACHIYO   \n"
	"                           -LIBXC=HYB_GGA_XC_B3LYP                \n"
#endif

	//"-DFT=BECKE88:          exchange: Becke-88 correlation: none    \n"

	"-Q=INT          Set total molecular charge (default=0)         \n"
	"-M=INT          Set molecular spin multiplicity (M=2S+1)       \n"
	"-R              Use restricted   orbitals (default if singlet) \n"
	"-U              Use unrestricted orbitals (default if  M > 1 ) \n"
	"                                                               \n"
	"                                                               \n"
	"[x] Compute DFT exchange from Hartree-Fock densities:          \n"
	"-xSlater        Slater LDA exchange\n"
	"-xPerdewYue     J.P.Perdew and W.Yue 1986 exchange\n"
	"-xPBE           Perdew–Burke–Ernzerhof exchange\n"
	"-xBecke88       A.D.Becke 1988 exchange\n"
	"-xMVS           meta-GGA Made Very Simple exchange\n"
	"-xChachiyo      T.Chachiyo and H.Chachiyo exchange\n"
	"                                                               \n"
	"                                                               \n"
	"[x] Post SCF:                                                  \n"
	"-FORCE          Compute forces acting on nuclei\n"
	"-HESSIAN        Compute hessian matrix and molecular vibrations\n"
	"-OPT            Request geometry optimization\n"
	"-MP2            Request MP2 energy calculations\n"
	"-MECP=INT,INT   Request MECP between two spin multiplicities\n"
	"-EXCITE         Request excited state calculations\n"
	"-QMD            Request quantum molecular dynamics simulation\n"
	"                                                                      \n"
	"                                                                      \n"
	"[x] SCF Cycle:                                                        \n"
	"-GUESS=DIAG     Use identity density matrix as initial guess (default)\n"
	"-GUESS=CORE     Use density from core hamiltonian as initial guess\n"
	"-GUESS=CHECK    Use density from checkpoint as initial guess\n"
	"-SCFDIIS        Use 4-point DIIS method for convergence (default)\n"
	"-SCFDIIS3       Use 3-point DIIS method for convergence\n"
	"-SCFDIIS2       Use 2-point DIIS method for convergence\n"
	"-SCFDAMP        Use simple weighting method for convergence\n"
	"-SCFDRAG=REAL   Set SCF drag coefficient between 0 to 1 (default=0.25)\n"
	"-SCFCONV=REAL   Set SCF convergence threshold (default=1.0E-6)\n"
	"-SCFMAX=INT     Set maximum number of scf cycle (default=80)\n"
	"-SCFACC=3STEP   Use increasing integral accuracy in 3 steps (default)\n"
	"-SCFACC=1STEP   Use fixed integral accuracy\n"
	"-MAXMEM=INT     Set maximum memory per CPU in Mbyte (default=250)\n"
	"                                                               \n"
	"                                                               \n"
	"[x] Grid:\n"
	"-GRIDSIZE=S     Good enough for a few milli-hartree accuracy\n"
	"-GRIDSIZE=M     Below milli-hartree (default for energy run)\n"
	"-GRIDSIZE=L     A few micro hartree (default for optimization)\n"
	"-GRIDSIZE=XL    Good for testing (96/590 radial/anuglar point)\n"
	"-GRID=BALLS     Use Becke partition + uniform radius (default)\n"
	"-GRID=BECKE     Use Becke partition + bragg radius\n"
	"-GRIDINIT=S     Use size S for early SCF cycles (default)\n"
	"-GRIDINIT=M     Use size M for early SCF cycles\n"

	"                                                               \n"
	"                                                               \n"
	"[x] Checkpoint File:\n"
	"-LCHECK         Do not perform SCF but load info from checkpoint\n"
	"-SCHECK         Save checkpoint file at the end (default=no)\n"
	"-SCHECK=ALL     Save checkpoint file every scf cycle (default=no)\n"
	"-FCHECK=STR     Set file name for checkpoint (default=checkpoint.txt)\n"
	"-LDMATRIX       Load density matrix at the beginning (default=no)\n"
	"-SDMATRIX       Save density matrix at the end (default=no)\n"
	"-FDMATRIX=STR   Set file name for density matrix (default=dmatrix.txt)\n"
	"                                                               \n"
	"                                                               \n"
	"[x] Output:                                           \n"
	"-DENSITY        Print electron density   volume information\n"
	"-POTENTIAL      Print electric potential volume information\n"
	"-GRADOVER43     Gradient over rho^4/3    volume information\n"
	"-GRADRS         Gradient electron radius volume information\n"
	"-MOUP=INT       Print spin up mo. volume info (index starts at 1)\n"
	"-MODN=INT       Print spin dn mo. volume info (index starts at 1)\n"
	"-VOLCUT=REAL    Set accuracy for computing volume info (default=1.0E-4)\n"
	"-VOLGRID=INT    Set the number of grid points per angstrom (default=10)\n"
	"-XSF            Volume info. will be in XSF  format to 'volume.xsf'\n"
	"-CUBE           Volume info. will be in CUBE format to 'volume.cube'\n"
	"-GAUSSIAN       Emulate Gaussian output to 'gaussian.log' (for GabEdit)\n"
	"                                                             \n"
	"                                                             \n"
	"[x] Parallel Run:                                            \n"
	"-NCPU=INT       Set the number of CPUs (default=1)           \n"
	"-PREFIX=STR     Set prefix string for the job (default=SQ)   \n"
	"                                                             \n"
	"                                                             \n"
	"[x] Geometry Optimization:                                   \n"
	"-OPTMAX=INT     Maximum number of iterations (default=30)    \n"
	"                                                             \n"
	"                                                             \n"
	"[x] Minimum Energy Crossing Point (MECP):                    \n"
	"-MECPMAX=INT    Maximum number of iterations (default=30) \n"
	"-FCHECKA=STR    State A checkpoint file name (default=checkpointA.txt)\n"
	"-FCHECKB=STR    State B checkpoint file name (default=checkpointB.txt)\n"
	"-FDMATRIXA=STR  State A density matrix file name (default=dmatrixA.txt)\n"
	"-FDMATRIXB=STR  State B density matrix file name (default=dmatrixB.txt)\n"
	"-GAUSSINA=STR   State A Gaussian input file name (excluding .com)\n"
	"-GAUSSINB=STR   State B Gaussian input file name (excluding .com)\n"
	"-GAUSSEXE=STR   Gaussian program execution string\n"
	//"                                                      \n"
	//"Excited State:                                        \n"
	//"-EXCITE=INT     Set the number of excitation energies (default=5)\n"
	"                                                                     \n"
	"                                                                     \n"
	"[x] External Field:                                                  \n"
	"-EF=EX,EY,EZ    Uniform electric field in AU (default=0.0,0.0,0.0)   \n"
	"                Electric field 1 AU = 51.4220652 Volt/Angstrom       \n"
	"                                                                     \n"
	"                                                                     \n"
	"[x] Quantum Molecular Dynamics (QMD):                                \n"
	"-INITVEL=STR    Initial velocity from the file in xyz format  (nm/ps)\n"
	"-INITTEMP=REAL  Initial velocity at random using the temp.   (kelvin)\n"
	"-KEEPTEMP=REAL  Rescale velocity to maintain the temp. (default=none)\n"
	"-TRAJ=STR       Set output trajectory file         (default=traj.xyz)\n"
	"-DT=REAL        Set time step in pico-sec             (default=0.001)\n"
	"-QMDMAX=INT     Set maximum number of steps           (default=25)\n"
	"-EFREQ=REAL     Set electric field frequency in THz   (default=0.0)\n"
	);
}


// parse_option : parse argv program argument into options
//
// 2008 - Teepanis Chachiyo
// 	Initial implementation
//
// Dec 31, 2009 - Teepanis Chachiyo
//  Add printing molecular orbital option
//
// July 11, 2010 - Teepanis Chachiyo
//  Unrestricted calculation options
//
// Oct 20, 2010 - Theerapon Khumla
//  Add force option
//
// Oct 21, 2010 - Teepanis Chachiyo
//  Add geoemtry optimization and optin SCFConv
//
// May 21, 2011 - Teepanis Chachiyo
//  Add DIIS and make it a default option
//  Bugfix, compare strcmp the whole string if the option needs no argument
//
// May 21, 2011 - Teepanis Chachiyo
//  Add MAXMEM option
//
// June 4, 2011 - Nanta Sophonrat
//  Add potential option
//
// Sep 30, 2012 - Teepanis Chachiyo
//  Add DIIS2 option
//
// Oct 4, 2012 - Teepanis Chachiyo
//  Add checkpoint options
//
// Nov 26, 2012 - Teepanis Chachiyo
//  Add checkpoint at every scf and multiple integral accuracy cutoff
//
// Jan 26, 2013 - Teepanis Chachiyo
//  Add DIIS-4 convergence method
//
// Feb 24, 2013 - Teepanis Chachiyo
//  SCFMax is now 80 by default
//
// March 17, 2013 - Teepanis Chachiyo
//  Add MECP
//
// Jan 16, 2014 - Teepanis Chachiyo
//  Handle whichMO index (which the users think that the index start at 1)
//
// Feb 25, 2016 - Teepanis Chachiyo
//  Add excitation options
//
// Mar 11, 2016 - Teepanis Chachiyo
//  Add uniform electric field and QMD options
//
// April 11, 2016 - Teepanis Chachiyo
//  Print current directory and full execution command for logging
//
// May 12, 2016 - Teepanis Chachiyo
//  Add DFT
//
// May 20, 2017 - Teepanis Chachiyo
//  Add Chachiyo exchange
//
// Nov 18, 2017 - Teepanis Chachiyo
//  Add outVolumeGrid to select the number of points per Angstrom
//
// Aug 18, 2019 - Teepanis Chachiyo
//  Add default for DFT integratio grid
//  Add DFT functional and LibXC parsing
//
// March 22, 2020 - Teepanis Chachiyo
//  Add Hessian option
//
void parse_option(struct option_t *opt, int argc, char *argv[]){

	int i;

	// set default
	opt->molCharge       = 0;
	opt->multiplicity    = 0;
	opt->restricted      = 1;
	opt->method          = METHOD_HF;
	opt->force           = 0;
	opt->hessian         = 0;
	opt->opt             = 0;
	opt->outVolumeGrid   = 10;
	opt->outVolumeType   = VOLUME_NONE;
	opt->outWhichMO      = 0;
	opt->outFormat       = VOLUME_XSF;
	opt->outVolumeCut    = 1.0E-4;
	opt->outGAUSSIAN     = 0;
	opt->SCFGuess        = SCFGUESS_DIAG;
	opt->SCFConv         = 1.0E-6;
	opt->SCFCutoff       = 1.0E-15;
	opt->SCFDrag         = 0.25;
	opt->SCFMax          = 80;
	opt->convMethod      = CONVMETHOD_DIIS4;
	opt->SCFAccuracy     = SCFACCURACY_3STEP;
	opt->maxMem          = 250;
	opt->loadDMatrix     = 0;
	opt->saveDMatrix     = 0;
	strcpy(opt->DMatrixFile,"dmatrix.txt");
	opt->saveCheck       = 0;
	opt->saveCheckAll    = 0;
	opt->loadCheck       = 0;
	strcpy(opt->CheckFile,"checkpoint.txt");
	opt->nCPU            = 1;
	strcpy(opt->prefixStr,"SQ");
	opt->MECP            = 0;
	opt->mecpMax         = 30;
	opt->mecpMA          = 0;
	opt->mecpMB          = 0;
	opt->excite          = 0;
	strcpy(opt->DMatrixFileA,"dmatrixA.txt");
	strcpy(opt->DMatrixFileB,"dmatrixB.txt");
	strcpy(opt->CheckFileA,"checkpointA.txt");
	strcpy(opt->CheckFileB,"checkpointB.txt");
	strcpy(opt->gaussEXE,"\0");
	strcpy(opt->gaussINA,"\0");
	strcpy(opt->gaussINB,"\0");
	opt->optMax          = 30;

	opt->Ex              = 0.0;
	opt->Ey              = 0.0;
	opt->Ez              = 0.0;

	opt->qmd             = 0;
	opt->initTemp        = 0.0;
	opt->keepTemp        = 0.0;
	opt->qmdFreq         = 0.0;
	opt->tStep           = 0.001;
	opt->maxQMDStep      = 25;
	strcpy(opt->initVel,"\0");
	strcpy(opt->traj,"traj.xyz");

	opt->whichHF2xDFT = HF2xDFT_xSLATER;

#ifdef LIBXC
	opt->func1Id = 0;
	opt->func2Id = 0;
	opt->hyb_hf_coef = 0.0;
#endif

	opt->whichGrid    = GRID_BALLS;
	opt->nradialInit  = GRIDSIZE_S_NRADIAL;
	opt->nlabedevInit = GRIDSIZE_S_NLABEDEV;
	opt->nradial      = GRIDSIZE_M_NRADIAL;
	opt->nlabedev     = GRIDSIZE_M_NLABEDEV;

	// print execution command
	char str[1024];
	if (getcwd(str, sizeof(str)) != NULL){
		printf("\nCurrent Directory: %s\n", str);
	}else{
		printf("parse_option: error cannot get current directory\n");
		exit(-1);
	}
	printf("Execution Command: ");
	for(i=0; i < argc; i++){
		printf("%s ",argv[i]);
	}
	printf("\n"); fflush(stdout);

	// loop throu all options
	for(i=3; i < argc; i++){
		if(strncmp(argv[i],"-Q=",3)==0){
			opt->molCharge = atoi(argv[i]+3);
			continue;
		}
		if(strncmp(argv[i],"-M=",3)==0){
			opt->multiplicity = atoi(argv[i]+3);
			continue;
		}
		if(strcmp(argv[i],"-HF")==0){
			opt->method = METHOD_HF;
			continue;
		}
		if(strcmp(argv[i],"-LDA=S:VWN5")==0){
			opt->method = METHOD_SVWN5;
			continue;
		}

		// override VWN5 parameters
		if(strncmp(argv[i],"-VWN5_PAR1_A=",13)==0){
			VWN5_A[0] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR1_x=",13)==0){
			VWN5_x[0] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR1_b=",13)==0){
			VWN5_b[0] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR1_c=",13)==0){
			VWN5_c[0] = atof(argv[i]+13);
			continue;
		}

		if(strncmp(argv[i],"-VWN5_PAR2_A=",13)==0){
			VWN5_A[1] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR2_x=",13)==0){
			VWN5_x[1] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR2_b=",13)==0){
			VWN5_b[1] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR2_c=",13)==0){
			VWN5_c[1] = atof(argv[i]+13);
			continue;
		}

		if(strncmp(argv[i],"-VWN5_PAR3_A=",13)==0){
			VWN5_A[2] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR3_x=",13)==0){
			VWN5_x[2] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR3_b=",13)==0){
			VWN5_b[2] = atof(argv[i]+13);
			continue;
		}
		if(strncmp(argv[i],"-VWN5_PAR3_c=",13)==0){
			VWN5_c[2] = atof(argv[i]+13);
			continue;
		}

		if(strcmp(argv[i],"-HALF")==0){
			opt->method = METHOD_HALF;
			continue;
		}
		if(strcmp(argv[i],"-LDA=S:CHACHIYO")==0){
			opt->method = METHOD_SCHACHIYO;
			continue;
		}
		if(strcmp(argv[i],"-DFT=CHACHIYO")==0){
			opt->method = METHOD_CHACHIYO;
			continue;
		}
		if(strcmp(argv[i],"-DFT=BECKE88:")==0){
			opt->method = METHOD_xBECKE88;
			continue;
		}
		if(strcmp(argv[i],"-xSlater")==0){
			opt->method = METHOD_HF2xDFT;
			opt->whichHF2xDFT = HF2xDFT_xSLATER;
			continue;
		}
		if(strcmp(argv[i],"-xPerdewYue")==0){
			opt->method = METHOD_HF2xDFT;
			opt->whichHF2xDFT = HF2xDFT_xPY86;
			continue;
		}
		if(strcmp(argv[i],"-xPBE")==0){
			opt->method = METHOD_HF2xDFT;
			opt->whichHF2xDFT = HF2xDFT_xPBE;
			continue;
		}
		if(strcmp(argv[i],"-xBecke88")==0){
			opt->method = METHOD_HF2xDFT;
			opt->whichHF2xDFT = HF2xDFT_xBECKE88;
			continue;
		}
		if(strcmp(argv[i],"-xMVS")==0){
			opt->method = METHOD_HF2xDFT;
			opt->whichHF2xDFT = HF2xDFT_xMVS;
			continue;
		}
		if(strcmp(argv[i],"-xChachiyo")==0){
			opt->method = METHOD_HF2xDFT;
			opt->whichHF2xDFT = HF2xDFT_xCHACHIYO;
			continue;
		}
		if(strcmp(argv[i],"-MP2")==0){
			opt->method = METHOD_MP2;
			continue;
		}
		if(strcmp(argv[i],"-R")==0){
			opt->restricted = 1;
			continue;
		}
		if(strcmp(argv[i],"-U")==0){
			opt->restricted = 0;
			continue;
		}
		// deprecated, for backward compatibility only 
		if(strcmp(argv[i],"-RHF")==0){
			opt->restricted = 1;
			opt->method     = METHOD_HF;
			continue;
		}
		// deprecated, for backward compatibility only 
		if(strcmp(argv[i],"-UHF")==0){
			opt->restricted = 0;
			opt->method     = METHOD_HF;
			continue;
		}
		if(strcmp(argv[i],"-FORCE")==0){
			opt->force = 1;
			opt->nradial  = GRIDSIZE_L_NRADIAL;
			opt->nlabedev = GRIDSIZE_L_NLABEDEV;
			continue;
		}
		if(strcmp(argv[i],"-HESSIAN")==0){
			opt->hessian = 1;
			opt->nradial  = GRIDSIZE_L_NRADIAL;
			opt->nlabedev = GRIDSIZE_L_NLABEDEV;
			continue;
		}
		if(strcmp(argv[i],"-OPT")==0){
			opt->opt = 1;
			opt->nradial  = GRIDSIZE_L_NRADIAL;
			opt->nlabedev = GRIDSIZE_L_NLABEDEV;
			continue;
		}
		if(strncmp(argv[i],"-OPTMAX=",8)==0){
			opt->optMax = atoi(argv[i]+8);
			continue;
		}
		if(strcmp(argv[i],"-DENSITY")==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_DENSITY_TOTAL;
			continue;
		}
		if(strcmp(argv[i],"-POTENTIAL")==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_POTENTIAL;
			continue;
		}
		if(strcmp(argv[i],"-GRADOVER43")==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_GRADRHO_OVER_RHO43;
			continue;
		}
		if(strcmp(argv[i],"-GRADRS")==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_GRAD_RS;
			continue;
		}
		if(strncmp(argv[i],"-VOLGRID=",9)==0){
			opt->outVolumeGrid = atoi(argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-MOUP=",6)==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_MO_ALPHA;
			opt->outWhichMO = atoi(argv[i]+6);

			// users think that the index starts at 1
			opt->outWhichMO = opt->outWhichMO - 1;
			continue;
		}
		if(strncmp(argv[i],"-MODN=",6)==0){
			if(opt->outVolumeType != VOLUME_NONE){
				printf("parse_option - error multiple volume types requested\n");
				exit(-1);
			}
			opt->outVolumeType = VOLUME_MO_BETA;
			opt->outWhichMO = atoi(argv[i]+6);

			// users think that the index starts at 1
			opt->outWhichMO = opt->outWhichMO - 1;
			continue;
		}
		if(strncmp(argv[i],"-VOLCUT=",8)==0){
			opt->outVolumeCut = atof(argv[i]+8);
			continue;
		}
		if(strcmp(argv[i],"-GAUSSIAN")==0){
			opt->outGAUSSIAN  = 1;
			continue;
		}
		if(strcmp(argv[i],"-XSF")==0){
			opt->outFormat  = VOLUME_XSF;
			continue;
		}
		if(strcmp(argv[i],"-CUBE")==0){
			opt->outFormat = VOLUME_CUBE;
			continue;
		}
		if(strcmp(argv[i],"-GUESS=DIAG")==0){
			opt->SCFGuess = SCFGUESS_DIAG;
			continue;
		}
		if(strcmp(argv[i],"-GUESS=CORE")==0){
			opt->SCFGuess = SCFGUESS_CORE;
			continue;
		}
		if(strcmp(argv[i],"-GUESS=CHECK")==0){
			opt->SCFGuess = SCFGUESS_CHECK;
			continue;
		}
		if(strcmp(argv[i],"-SCFACC=3STEP")==0){
			opt->SCFAccuracy = SCFACCURACY_3STEP;
			continue;
		}
		if(strcmp(argv[i],"-SCFACC=1STEP")==0){
			opt->SCFAccuracy = SCFACCURACY_1STEP;
			continue;
		}
		if(strcmp(argv[i],"-SCFDIIS")==0){
			opt->convMethod = CONVMETHOD_DIIS4;
			continue;
		}
		if(strcmp(argv[i],"-SCFDIIS3")==0){
			opt->convMethod = CONVMETHOD_DIIS3;
			continue;
		}
		if(strcmp(argv[i],"-SCFDIIS2")==0){
			opt->convMethod = CONVMETHOD_DIIS2;
			continue;
		}
		if(strcmp(argv[i],"-SCFDAMP")==0){
			opt->convMethod = CONVMETHOD_DAMPING;
			continue;
		}
		if(strncmp(argv[i],"-MAXMEM=",8)==0){
			opt->maxMem = atoi(argv[i]+8);
			continue;
		}
		if(strncmp(argv[i],"-SCFCONV=",9)==0){
			opt->SCFConv = atof(argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-SCFDRAG=",9)==0){
			opt->SCFDrag = atof(argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-SCFMAX=",8)==0){
			opt->SCFMax = atoi(argv[i]+8);
			continue;
		}
		if(strcmp(argv[i],"-LDMATRIX")==0){
			opt->loadDMatrix = 1;
			continue;
		}
		if(strcmp(argv[i],"-SDMATRIX")==0){
			opt->saveDMatrix = 1;
			continue;
		}
		if(strncmp(argv[i],"-FDMATRIX=",10)==0){
			strcpy(opt->DMatrixFile, argv[i]+10);
			continue;
		}
		if(strcmp(argv[i],"-SCHECK")==0){
			opt->saveCheck = 1;
			continue;
		}
		if(strcmp(argv[i],"-SCHECK=ALL")==0){
			opt->saveCheckAll = 1;
			continue;
		}
		if(strcmp(argv[i],"-LCHECK")==0){
			opt->loadCheck = 1;
			continue;
		}
		if(strncmp(argv[i],"-FCHECK=",8)==0){
			strcpy(opt->CheckFile, argv[i]+8);
			continue;
		}
		if(strncmp(argv[i],"-NCPU=",6)==0){
			opt->nCPU = atoi(argv[i]+6);
			continue;
		}
		if(strncmp(argv[i],"-PREFIX=",8)==0){
			strcpy(opt->prefixStr, argv[i]+8);
			continue;
		}
		if(strncmp(argv[i],"-MECP=",6)==0){
			if(sscanf(argv[i]+6,"%d,%d",&opt->mecpMA,&opt->mecpMB)!=2){
				printf("parse_option - error cannot recognize option %s\n",argv[i]);
				exit(-1);
			}
			opt->MECP = 1;
			opt->nradial  = GRIDSIZE_L_NRADIAL;
			opt->nlabedev = GRIDSIZE_L_NLABEDEV;
			continue;
		}
		if(strncmp(argv[i],"-MECPMAX=",9)==0){
			opt->mecpMax = atoi(argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-FDMATRIXA=",11)==0){
			strcpy(opt->DMatrixFileA, argv[i]+11);
			continue;
		}
		if(strncmp(argv[i],"-FDMATRIXB=",11)==0){
			strcpy(opt->DMatrixFileB, argv[i]+11);
			continue;
		}
		if(strncmp(argv[i],"-FCHECKA=",9)==0){
			strcpy(opt->CheckFileA, argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-FCHECKB=",9)==0){
			strcpy(opt->CheckFileB, argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-GAUSSEXE=",10)==0){
			strcpy(opt->gaussEXE, argv[i]+10);
			continue;
		}
		if(strncmp(argv[i],"-GAUSSINA=",10)==0){
			strcpy(opt->gaussINA, argv[i]+10);
			continue;
		}
		if(strncmp(argv[i],"-GAUSSINB=",10)==0){
			strcpy(opt->gaussINB, argv[i]+10);
			continue;
		}
		if(strcmp(argv[i],"-EXCITE")==0){
			opt->excite = EXCITED_STATES;
			continue;
		}
		if(strncmp(argv[i],"-EXCITE=",8)==0){
			if(sscanf(argv[i]+8,"%d",&opt->excite)!=1){
				printf("parse_option - error cannot recognize option %s\n",argv[i]);
				exit(-1);
			}
			continue;
		}

		if(strcmp(argv[i],"-QMD")==0){
			opt->qmd = 1;
			opt->nradial  = GRIDSIZE_L_NRADIAL;
			opt->nlabedev = GRIDSIZE_L_NLABEDEV;
			continue;
		}
		if(strncmp(argv[i],"-INITTEMP=",10)==0){
			opt->initTemp = atof(argv[i]+10);
			continue;
		}
		if(strncmp(argv[i],"-KEEPTEMP=",10)==0){
			opt->keepTemp = atof(argv[i]+10);
			continue;
		}
		if(strncmp(argv[i],"-INITVEL=",9)==0){
			strcpy(opt->initVel, argv[i]+9);
			continue;
		}
		if(strncmp(argv[i],"-FREQ=",6)==0){
			opt->qmdFreq = atof(argv[i]+6);
			continue;
		}
		if(strncmp(argv[i],"-EF=",4)==0){
			if(sscanf(argv[i]+4,"%lf,%lf,%lf",&opt->Ex,&opt->Ey,&opt->Ez)!=3){
				printf("parse_option - error cannot recognize option %s\n",argv[i]);
				exit(-1);
			}
			continue;
		}
		if(strncmp(argv[i],"-TRAJ=",6)==0){
			strcpy(opt->traj, argv[i]+6);
			continue;
		}
		if(strncmp(argv[i],"-DT=",4)==0){
			opt->tStep = atof(argv[i]+4);
			continue;
		}
		if(strncmp(argv[i],"-QMDMAX=",8)==0){
			if(sscanf(argv[i]+8,"%d",&opt->maxQMDStep)!=1){
				printf("parse_option - error cannot recognize option %s\n",argv[i]);
				exit(-1);
			}
			continue;
		}
#ifdef LIBXC
		if(strcmp(argv[i],"-DFT=PBE")==0){
			opt->method  = METHOD_LIBXC;
			opt->func1Id = XC_GGA_X_PBE;
			opt->func2Id = XC_GGA_C_PBE;
			opt->hyb_hf_coef = 0.0;
			continue;
		}
		if(strcmp(argv[i],"-DFT=BLYP")==0){
			opt->method  = METHOD_LIBXC;
			opt->func1Id = XC_GGA_X_B88;
			opt->func2Id = XC_GGA_C_LYP;
			opt->hyb_hf_coef = 0.0;
			continue;
		}
		if(strcmp(argv[i],"-DFT=B3LYP")==0){
			opt->method  = METHOD_LIBXC;
			opt->func1Id = XC_HYB_GGA_XC_B3LYP;
			opt->func2Id = 0;

			// get hartree-fock coefficient
			xc_func_type func1;
			xc_func_init(&func1, opt->func1Id, XC_POLARIZED);
			opt->hyb_hf_coef = xc_hyb_exx_coef(&func1);
			xc_func_end(&func1);
			continue;
		}
		if(strncmp(argv[i],"-LIBXC=",7)==0){

			opt->method = METHOD_LIBXC;

			// input: -LIBXC=str1   OR  -LIBXC=str1+str2
			// extract functional names to str1 and str2
			char *str1, *str2;
			char buffer[128];
			int nStr;
			strncpy(buffer, argv[i]+7,128);
			str1 = buffer;
			str2 = strstr(buffer,"+");
			if(str2==NULL){ nStr = 1; }
			else{
				nStr = 2;
				*str2 = '\0';
				 str2++;
				if(strlen(str2)==0){
					printf("parse_option - error cannot parse option %s\n",argv[i]);
					exit(-1);
				}
			}

			// load functional id from the names
			opt->func1Id = xc_functional_get_number(str1);
			if(opt->func1Id <= 0){
				printf("parse_option - error, no functional %s in LibXC\n",str1);
				exit(-1);
			}
			if(nStr==2){
				opt->func2Id = xc_functional_get_number(str2);
				if(opt->func2Id <= 0){
					printf("parse_option - error, no functional %s in LibXC\n",str2);
					exit(-1);
				}
			}

			// check if func1Id type is supported
			xc_func_type func;
			xc_func_init(&func, opt->func1Id, XC_POLARIZED);
			switch(func.info->family){
				case XC_FAMILY_LDA:
				case XC_FAMILY_GGA:
				case XC_FAMILY_HYB_GGA:
					// read hf coefficient
					opt->hyb_hf_coef = xc_hyb_exx_coef(&func);
				break;
				case XC_FAMILY_MGGA:
				case XC_FAMILY_HYB_MGGA:
					printf("parse_option - error, Meta-GGA is not supported yet\n");
					exit(-1);
				break;
				default:
					printf("parse_option - error, the functional is not supported\n");
					exit(-1);
			}
			xc_func_end(&func);

			// check if fundId2 is supported
			if(nStr==2){
			xc_func_init(&func, opt->func2Id, XC_POLARIZED);
			switch(func.info->family){
				case XC_FAMILY_LDA:
				case XC_FAMILY_GGA:
				case XC_FAMILY_HYB_GGA:
					// read hf coefficient and add to the first one
					opt->hyb_hf_coef += xc_hyb_exx_coef(&func);
				break;
				case XC_FAMILY_MGGA:
				case XC_FAMILY_HYB_MGGA:
					printf("parse_option - error, Meta-GGA is not supported yet\n");
					exit(-1);
				break;
				default:
					printf("parse_option - error, the functional is not supported\n");
					exit(-1);
			}
			xc_func_end(&func);
			}

			continue;
		}
#endif

		if(strcmp(argv[i],"-GRIDSIZE=S")==0){
			opt->nradial  = GRIDSIZE_S_NRADIAL;
			opt->nlabedev = GRIDSIZE_S_NLABEDEV;
			continue;
		}
		if(strcmp(argv[i],"-GRIDSIZE=M")==0){
			opt->nradial  = GRIDSIZE_M_NRADIAL;
			opt->nlabedev = GRIDSIZE_M_NLABEDEV;
			continue;
		}
		if(strcmp(argv[i],"-GRIDSIZE=L")==0){
			opt->nradial  = GRIDSIZE_L_NRADIAL;
			opt->nlabedev = GRIDSIZE_L_NLABEDEV;
			continue;
		}
		if(strcmp(argv[i],"-GRIDSIZE=XL")==0){
			opt->nradial  = GRIDSIZE_XL_NRADIAL;
			opt->nlabedev = GRIDSIZE_XL_NLABEDEV;
			continue;
		}
		if(strcmp(argv[i],"-GRID=BECKE")==0){
			opt->whichGrid  = GRID_BECKE;
			continue;
		}
		if(strcmp(argv[i],"-GRID=BALLS")==0){
			opt->whichGrid  = GRID_BALLS;
			continue;
		}
		if(strcmp(argv[i],"-GRIDINIT=M")==0){
			opt->nradialInit  = GRIDSIZE_M_NRADIAL;
			opt->nlabedevInit = GRIDSIZE_M_NLABEDEV;
			continue;
		}
		if(strcmp(argv[i],"-GRIDINIT=S")==0){
			opt->nradialInit  = GRIDSIZE_S_NRADIAL;
			opt->nlabedevInit = GRIDSIZE_S_NLABEDEV;
			continue;
		}

		// cannot recognize parameter
		printf("parse_option - error cannot recognize option %s\n", argv[i]);
		exit(-1);
	}

	// validate 
	if(opt->outVolumeType==VOLUME_MO_ALPHA || opt->outVolumeType==VOLUME_MO_BETA)
	if(opt->outWhichMO <= 0){
		printf("parse_option - error molecular orbital index should be greater than zero\n");
		exit(-1);
	}

	// validate SCFConv
	if(opt->SCFConv <= 0.0){
		printf("parse_option - error invalid SCFConv range\n");
		exit(-1);
	}

	// validate volumeCut
	if(opt->outVolumeCut <= 0.0){
		printf("parse_option - error invalid volumeCut range\n");
		exit(-1);
	}

	// validate SCFDrag
	if(opt->SCFDrag <= 0.0 || opt->SCFDrag > 1){
		printf("parse_option - error invalid SCFDrag range\n");
		exit(-1);
	}

	// validate SCFMax
	if(opt->SCFMax < 0){
		printf("parse_option - error invalid SCFMax range\n");
		exit(-1);
	}

	// validate maxMem
	if(opt->maxMem < 0){
		printf("parse_option - error invalid MAXMEM range\n");
		exit(-1);
	}

	// validate checkpoint related options
	if(opt->loadCheck && opt->opt){
		printf("parse_option - error OPT cannot be used with LCHECK\n");
		exit(-1);
	}
	if(opt->loadCheck && (opt->saveCheck || opt->saveCheckAll)){
		printf("parse_option - error LCHECK cannot be used with saving checkpoints\n");
		exit(-1);
	}
	if(opt->saveCheck && opt->saveCheckAll){
		printf("parse_option - error SCHECK cannot be used with SCHECK=ALL");
		exit(-1);
	}

	// optimization does not support mp2 yet
	if(opt->opt && opt->method==METHOD_MP2){
		printf("parse_option - error OPT does not support MP2 at the moment\n");
		exit(-1);
	}

	// validate nCPU
	if(opt->nCPU <= 0){
		printf("parse_option - error invalid number of cpus\n");
		exit(-1);
	}

	// validate MECP
	if(opt->MECP){
		if(opt->opt){
			printf("parse_option - error MECP cannot be used with OPT\n");
			exit(-1);
		}
		if(opt->method==METHOD_MP2){
			printf("parse_option - error MECP cannot be used with MP2\n");
			exit(-1);
		}
	}

	// validate QMD
	if(opt->qmd){

		if(opt->opt){
			printf("parse_option - error choose either OPT or QMD\n");
			exit(-1);
		}

		if(opt->initTemp < 0.0){
			printf("parse_option - error INITTEMP less than zero\n");
			exit(-1);
		}
		if(opt->keepTemp < 0.0){
			printf("parse_option - error KEEPTEMP less than zero\n");
			exit(-1);
		}
		if(opt->maxQMDStep < 0){
			printf("parse_option - error QMDMAX less than zero\n");
			exit(-1);
		}
		if(opt->tStep < 0.0){
			printf("parse_option - error DT less than zero\n");
			exit(-1);
		}
		if(opt->initTemp>0.0 && opt->initVel[0]!='\0'){
			printf("parse_option - choose either INITTEMP or INITVEL\n");
			exit(-1);
		}
	}

}

