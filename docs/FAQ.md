# Frequently Asked Questions

- [Where can I get some help running Siam Quantum](#where-can-i-get-some-help-running-siam-quantum)
- [How do I run a simple calculation](#how-do-i-run-a-simple-calculation)
- [How do I specify molecular charge and spin](#how-do-i-specify-molecular-charge-and-spin)
- [Are there other options I can use](#are-there-other-options-i-can-use)
- [How do I compute the equilibrium structure of a molecule](#how-do-i-compute-the-equilibrium-structure-of-a-molecule)
- [How can I speed up the calculation](#how-can-i-speed-up-the-calculation)
- [Does Siam Quantum support parallel runs](#does-siam-quantum-support-parallel-runs)
- [How do the results compare to that of other programs](#how-do-the-results-compare-to-that-of-other-programs)
- [How can I compute MP2 energy](#how-can-i-compute-mp2-energy)
- [How do I visualize molecular orbital or electron density](#how-do-i-visualize-molecular-orbital-or-electron-density)
- [How to make use of the checkpoint file](#how-to-make-use-of-the-checkpoint-file)
- [What if the SCF does not converge](#what-if-the-scf-does-not-converge)
- [Can I apply uniform electric field](#can-i-apply-uniform-electric-field)
- [How to compute the MECP structure](#how-to-compute-the-mecp-structure)
- [How do I find the MECP for other quantum model](#how-do-i-find-the-mecp-for-other-quantum-model)



## Where can I get some help running Siam Quantum

We are eager to hear from you and use your feedback to develop Siam Quantum. Please feel free to contact us via email  **teepanisc@nu.ac.th** or **hathaithip.chachiyo@gmail.com**

## How do I run a simple calculation

Siam Quantum needs two files to start a calculation: 1) molecular coordinate, and 2) basis set information. The first file contains coordinates of atoms in a molecule which must be in XYZ format. There are a few examples given in the directory ```examples```. The second file contains basis set information which is used to construct the wavefunction of the molecule. Please see the directory ```basis``` for available basis sets such as STO-3G, 3-21G, 6-31G, and 6-31G*.
For example, if we want to compute Hartree-Fock energy of a benzene molecule using 3-21G basis set, we can execute the following command:

```bash
# ./sq  benzene.xyz  321g.txt
```

Here, the file ```benzene.xyz``` contains the atomic coordinates of the molecule, and ```321g.txt``` contains the basis set information. The preceding ```./``` symbol tells the shell prompt that the executable is in the current directory. If the executable or the two files are not located at the current directory, <u>be sure to provide the appropriate path for the files</u> so that Siam Quantum can find and load them. An alternative would be to copy the files to the current directory.


## How do I specify molecular charge and spin

By default, Siam Quantum sets the molecular charge to zero. The spin is set to zero (S=0) if the molecule contains even number of electrons, or to one-half (S=1/2) otherwise. Users can set the charge and the spin by using the options ```-Q=INTEGER``` and ```-M=INTEGER```, for example,

```bash
# ./sq  benzene.xyz  321g.txt  -Q=2  -M=3
```

Here, the molecular charge is set to +2; and the spin multiplicity is 3. Note that the multiplicity M is related to the spin via M = 2S+1. Therefore, in this example, the benzene molecule has lost 2 electrons and is in the spin triplet state S=1.

## Are there other options I can use

Siam Quantum will print out a list of available options if we simply execute it without any argument. For example, using the following command

```bash
# ./sq
```

will produce the output:

```
[OPTIONS] :  
[x] Ab Initio Method: 
-Q=INT       Set total molecular charge (default=0) 
-M=INT       Set molecular spin multiplicity (M=2S+1) 
-RHF         Restricted Hartree-Fock (default for singlet state) 
-UHF         Unrestricted Hartree-Fock (default if M > 1) 
-FORCE       Calculate force acting on nuclei 
-OPT         Request geometry optimization 
-MP2         Request MP2 energy calculations  
... continue ...  
```

## How do I compute the equilibrium structure of a molecule

Aside from calculating the total energy of the molecule using Hartree-Fock theory, Siam Quantum can compute the three dimensional structure of the molecule such that the energy is minimum. This process is typically called **Geometry Optimization**, which can be done using option ```-OPT```. For example,

```bash
# ./sq  benzene.xyz  321g.txt  -OPT
```

will use the molecular structure contained in the file ```benzene.xyz``` as a starting point, and iteratively update the structure until the energy comes to a minimum.
Every time SQ updates geometry, it has to recalculate the wave function; and you can ask SQ to use the wave function from the previous step as a starting point so that it finds the solution easier. Please look at how to use checkpoint file for details.

## How can I speed up the calculation

In many cases, we can speed up the calculation significantly by storing the integrals in the memory as supposed to repeatedly evaluate them. The amount of memory available to Siam Quantum is specified by the option ```-MAXMEM=INTEGER```. For example, if we want Siam Quantum to use up to 500 megabytes of memory, simply execute:

```bash
# ./sq  benzene.xyz  321g.txt  -MAXMEM=500
```

By default, Siam Quantum goes through a cycle of calculations called **SCF Cycle** until the total energy converges to within 10-6 Hartrees. However, if the accuracy is too high for your purpose, you can loosen up the convergence criteria by using option ```-SCFCONV=REAL```, and as a result, making it faster. For example, the following command will direct Siam Quantum to set the convergence criteria to 10-4 Hartrees instead. For example,

```bash
# ./sq  benzene.xyz  321g.txt  -SCFCONV=1.0E-4
```

The command above will direct Siam Quantum to set the convergence criteria to 10<sup>-4</sup> Hartrees instead.

## Does Siam Quantum support parallel runs

Yes, with a few precautions. First, if your calculation is taking only a few minutes to complete, then adding more CPUs will not help. But, if your tasks take longer than 10 - 15 minutes, then it's a good idea to run in parallel. For example, in the following command, we will compute the Hartree-Fock energy of the C60 molecule using 4 CPUs, each of which is allowed to allocate up to 1024 MBytes of memory.

```bash
# ./sq   c60.xyz   321g.txt  -MAXMEM=1024  -NCPU=4
```

Second, if you decide to run multiple jobs under the same directory. You need to specify the prefix name for each job to avoid confusing Siam Quantum when it tries to communicate with each CPU during parallel runs. For instance, suppose you need to perform two calculations at the same time under the same directory, you could execute the following commands:

```bash
# cd   c60_calculations  
# ./sq   c60.xyz   321g.txt  -NCPU=4  -PREFIX=JOB1_NAME & 
# ./sq   c60.xyz   631g.txt  -NCPU=4  -PREFIX=JOB2_NAME &
```

## How do the results compare to that of other programs

They should be identical to that of the well-established software such as GAMESS or Gaussian, but with a few precautions.

First, Siam Quantum does not support symmetry. It uses the input geometry as is; whereas Gaussian or GAMESS tend to rotate the input geometry before going into SCF cycles. To turn off symmetry in Gaussian, use the keyword **NoSymm** in the route section. For GAMESS, in the $DATA section, use the keyword **C1** right after the comment line.

Second, Siam Quantum uses full Cartesian basis functions; whereas Gaussian uses 5D/7F or sometimes 6D/10F depending on types of basis function. To be sure, if you wish to compare the results between Siam Quantum and Gaussian, specify the keyword 6D/10F in the Gaussian's route section.

## How can I compute MP2 energy

When it comes to MP2 or Moller-Plesset Perturbation Theory, Siam Quantum is one of the fastest and requires very small storage. Simply use the ```-MP2``` option. For example, if you need to compute MP2 energy of the C60 molecule, simply execute:

```bash
# ./sq  c60.xyz  321g.txt  -MAXMEM=3000  -NCPU=2  -MP2
```

The above calculation requires the total of 6GB of memory (3G each CPU). On an i7 @3.4GHz machine, this took only 47 minutes to complete. Siam Quantum does not use scratch files. It stores all integrals in memory.


## How do I visualize molecular orbital or electron density

Use option ```-GAUSSIAN``` to direct Siam Quantum to generate a file called **gaussian.log** which contains the information about the molecular orbital and basis function. The content has the format similar to the one from Gaussian program when using GFInput IOP(6/7=3) option.

Once you have the file **gaussian.log**, use a visualization software called **Gabedit** to load the molecular orbital. This method might work for other software like Molden or Jmol, but so far, only with Gabedit has been tested to work.

## How to make use of the checkpoint file

Siam Quantum can save the calculated wave function and other related information so that the user can restart the job or can use the previous wave function as a starting point for further calculations. Consider the following sequence of calculations.

```bash
# ./sq  benzene.xyz  321g.txt  -SCHECK 
# ./sq  benzene.xyz  631g_star.txt  -GUESS=CHECK

```

In the first step we use a small basis set 3-21G to compute the wave function of a benzene molecule. The option ```-SCHECK``` instructs Siam Quantum to save the converged wave function to a checkpoint file. By default, the checkpoint file name is **checkpoint.txt**. However, you can use the option ```-FCHECK=STR``` to change to file name to fit your need, for example, ```-FCHECK=benzene_hf_321g.chk```.

In the second step, with the option ```-GUESS=CHECK```, Siam Quantum looks for the checkpoint file, reads it, and uses the molecular orbitals to construct an initial guess for the larger basis set 6-31G*. Having an already good initial guess, it becomes faster to find the converged solution.

Another advantage of the checkpoint file is to avoid recalculating the SCF cycle and skip ahead to the post-SCF steps. Consider the following sequence of runs.

```bash
# ./sq  benzene.xyz  321g.txt  -SCHECK 
# ./sq  benzene.xyz  321g.txt  -LCHECK -DENSITY 
# ./sq  benzene.xyz  321g.txt  -LCHECK -POTENTIAL

```

In the first step, we perform SCF cycle to compute the wave function of the molecule and then save it to a checkpoint file. In the second and the third step, we simply invoke the option ```-LCHECK``` which instructs Siam Quantum to simply load the calculated wave function from the checkpoint file and perform post-SCF tasks, which in this example, compute the 3D volume of the electron density and the electrostatic potential.

Typically with the option ```-SHCECK```, Siam Quantum only save the data at the end of the calculation even if the SCF cycle is not converged. However, you could ask the program to save the checkpoint file every SCF cycle with the option ```-SCHECK -SCHECK=ALL```.
Furthermore, we can take advantage of the checkpoint file during geometry optimization as well. For example,

```bash
# ./sq   benzene.xyz  631g.txt   -SCHECK 
# ./sq   benzene.xyz  631g.txt   -SCHECK  -GUESS=CHECK -OPT
```

In the first step, we perform a single point calculation to simply get the converged wave function of the initial molecular structure and save it to the checkpoint file. The second step is the geometry optimization using the option ```-OPT```. In this step, with additional option ```-SCHECK -GUESS=CHECK```, Siam Quantum will read the initial guess from a checkpoint file from the previous iterative cycle. Once it is done computing the SCF energy for the cycle, it saves the checkpoint file for the next cycle.


## What if the SCF does not converge

Although you could increase the maximum number of SCF iteration simply by the option ```-SCFMAX=INTEGER```, this is usually not the root cause of the problem. As a rule of thumb, you will find the wave function difficult to converge due to several reasons.

1) The molecular charge or spin makes no chemical sense. Or, the molecular geometry is too far off from where it’s supposed to be.
2) The energy levels of the molecule are very close together that it becomes too difficult for the program to settle down and pick which one is actually the lower energy configuration. This usually happens when there are metal atoms in the molecule such as Fe, Co, or Cu.
3) The basis set is very large; and hence, there are so many degree of freedom for the wave function that it becomes difficult to find the right one.
4) The SCF cycles starts off with a poor initial guess. As of version 1.2.8, SQ is using a very poor initial guess because it basically assumes the electron density is evenly distributed among the basis functions. We are working on this problem and hope to use a better scheme in the next version.

As a good practice when running quantum chemistry program, you should start with a small basis set and increasingly step up to bigger ones. For example, suppose your intended basis is 6-31G*, here is what I recommend you always do.

```bash
# ./sq    molecule.xyz   sto3g.txt       -SCHECK 
# ./sq    molecule.xyz   321g.txt        -SCHECK    -GUESS=CHECK 
# ./sq    molecule.xyz   631g.txt        -SCHECK    -GUESS=CHECK 
# ./sq    molecule.xyz   631g_star.txt   -SCHECK    -GUESS=CHECK 
```

In the example above, we start with a small basis set and work out way up to the target basis set. The option ```-SCHECK``` asks Siam Quantum to save the checkpoint file after the calculation is done; and ```-GUESS=CHECK``` ask SQ to read the current checkpoint file in the current directory and use the stored molecular orbital as an initial guess. You might even find that running in steps as shown above actually takes shorter time than starting with the 6-31G* right away.


## Can I apply uniform electric field

Yes, with the option ```-EF=EX,EY,EZ``` For example, the command

```bash
# ./sq benzene.xyz 321g.txt  -EF=0.01,0,0
```

adds an external electric field 0.01 AU in the x-direction. Be careful, the field is in atomic unit 1 AU = 5.14x10<sup>11</sup> V/m. A typical value should be on the order of 0.01 AU for a laser beam.

The option ```-EF-EX,EY,EZ``` is compatible with geometry optimization ```-OPT``` and quantum molecular dynamics ```-QMD```.



## How to compute the MECP structure

Consider a phenyl cation molecule which can have two possible spin states: single S=0 and triplet S=1. For a typical three-dimensional structure of this molecule, the non-relativistic electronic energies of both spin states are different. But for a special set of 3D molecular structures, the energies happen to be the same. This is called the "Crossing" of the energy surfaces. Out of these set in which they cross, they is a minimum called "Minimum Energy Crossing Point" or MECP.

Siam Quantum can compute the MECP with the option ```-MECP=INT,INT```. For example, with the provided phenyl structure in the "examples" directory, we can execute the following command

```bash
# ./sq phenyl.xyz 321g.txt -Q=1 -MECP=1,3 -GUESS=CORE -SCFMAX=120
``` 

The option ```-MECP=1,3``` means that the crossing is between the state with multiplicity M=1 (singlet S=0) and the multiplicity M=3 (triplet S=1). Unfortunately, the phenyl is a very tricky molecule in that the triplet state wave function is difficult to converge and is very sensitive an initial guess. Here, we add option ```-SCFMAX=120``` to increase the number of maximum SCF cycle. The option  ```-GUESS=CORE``` is to use the wave function from core-hamiltonian (without electron repulsion) as an initial guess. Generally this is not a good initial guess, but happens to give the right wave function for this tricky case.

The algorithm is based on Teepanis Chachiyo, and Jorge H. Rodriguez. <em>A direct method for locating minimum-energy crossing points (MECPs) in spin-forbidden transitions and nonadiabatic reactions</em>. J. Chem. Phys. 123 (2005): 094711. Note that there is a typo in equation (11) in the paper: the last term should have been a minus sign. A more detailed tutorial can be found [here.](https://www.researchgate.net/publication/278244924_Minimum_Energy_Crossing_Point_-_Tutorial_and_Theory)


## How do I find the MECP for other quantum model

Energy gradient (e.g. force) is required to compute the MECP. Presently, Siam Quantum can compute Hartree-Fock and DFT energy gradient, and so the ```-MECP``` option is fully compatible with the Hartree-Fock and DFT. 

However, for other methods such as MP2 or CI you would need other software to compute the energy gradient and use Siam Quantum to advance the search into the next step. The calculation of MECP consists of two parts: 1) computing the forces acting on nuclei and 2) updating the geometry of the molecule so that it becomes closer to the MECP, based on the algorithm by <em>T. Chachiyo, and J. H. Rodriguez. J. Chem. Phys. 123 (2005): 094711.</em>. 

In the first step, Siam Quantum can call an external program such as Gaussian, or GAMESS and retrieve the calculated forces to proceed in the updating steps. To find the MECP in this way, first the user has to prepare two input files, one for each electronic state of interest. 

For example, we want to consider phenyl cation molecule and find the MECP between the spin singlet and spin triplet state. We then prepare two following Gaussian input files, named phenyl_m1.com and phenyl_m3.com accordingly.

```
#P HF/3-21G NOSYMM FORCE SCF=Tight  MECP  1 1 SQ_GEOMETRY

#P HF/3-21G NOSYMM FORCE SCF=Tight  MECP  1 3 SQ_GEOMETRY  
```

Both files look almost exactly like an ordinary Gaussian input file with two exceptions: 1) there must be a capitalized option **NOSYMM FORCE**, and 2) replace the molecular structure with the keyword **SQ_GEOMETRY**.

The option **NOSYMM** is needed to prevent Gaussian from rotating the molecule into a standard orientation and confuses Siam Quantum in the process. The option **FORCE** asks Gaussian to compute the forces on nuclei so that Siam Quantum and retrieve them.

As Siam Quantum proceeds in an iterative fashion to compute the MECP, it continuously updates the geometry and calls Gaussian program to compute the forces. When calling the program, Siam Quantum will replace the **SQ_GEOMETRY** with the current molecular geometry.

Even though in the examples, the file <em>phenyl_m1.com</em> and <em>phenyl_m3.com</em> use the option **HF/3-21G** which instructs Gaussian program to use Hartree-Fock method, the users can choose other methods as well such as B3LYP/3-21G, BLYP/6-31G*, or whichever is supported by the Gaussian program.

Once the preparation is done, simply execute the following command in a <u>single line</u>:

```bash
# ./sq phenyl.xyz 321g.txt -Q=1 -MECP=1,3 -GAUSSEXE=g09 -GAUSSINA=phenyl_m1 -GAUSSINB=phenyl_m3
```

The option ```-GAUSSEXE=STR``` tells Siam Quantum how to execute the Gaussian program. In this example, we use Gaussian version 09. In addition, when specifying the two prepared input files, do not include the .com extension.







