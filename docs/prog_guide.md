# Programming Guides

## 1. The Big Picture - Standard Execution

Below is the map of the source code as SQ undergoes a typical calculation, for example, computing the Hartree-Fock energy and simple population analysis. Please see the details of the theoretical and computational implementation in the Markdown document for each source code.


```
[sq.c] - main subroutine, greetings
 |
 +- 1) call [option.c] - parse user inputs
 +- 2) call [mol.c]    - builds 3D molecular structure in memory
 +- 3) call [basis.c]  - prepares basis functions for the molecules
 +- 4) call [uhf.c]    - solves for the wave function using Hartree-Fock
 |            |
 |            +- 4.1) call [matrix.c] - compute various matrix elements
 |            |                         such as overlap, kinetic, nuclear,
 |            |                         and 2-electron integral
 |            +- 4.2) call [lin.c]    - front-end for LAPACK interface
 |            |                         for solving eigen equation or
 |            |                         an inverse of a matrix
 |            +- 4.3) call [conv.c]   - accelerate convergence using
 |            |                         methods like DIIS and damping
 |            +- 4.4) call [check.c]  - save the wave function and execution
 |                                      status into a checkpoint file
 |
 +- 5) call [pop.c]    - use the wave function (in the form of molecular
                         orbitals) to compute electric properties such as
                         electric field, potential, and Mulliken population
```


## 2. Optional Executions

Many optional calculations (e.g. geometry optimization, MP2, printing density for visualization) are available in Siam Quantum. The set of molecular orbitals $\{\phi_i(\vec{r})\}$ from the previously completed single-point calculations are needed to perform further executions. The map of the source code is as follow.

```
[sq.c] - main subroutine
 |
 +- perform single-point calculation (e.g. HF, DFT)
 |
 +- call [mp2.c]      - compute 2nd order Moller-Plesset perturbation energy
 +- call [optimize.c] - perform geometry optimization
 +- call [mecp.c]     - compute minimum energy crossing point
 +- call [xsf.c]      - save 3D volume information for visualization
```

### 3. Handling of the Integrals

Evaluating integral between molecular orbitals or between basis functions is the essential part of solving the electronic structure problem. These are typically called from **[uhf.c]** during the construction of the Fock matrix.

```
[matrix.c] - compute on-demand matrix elements
 |
 +- calls [int.c]       - for simple integrals such as overlap, kinetic,
 |         |              nuclei, and 2-electron integrals
 |         |
 |         +- calls [fgamma.c] - compute Boys’s F-gamma function
 |
 +- calls [grad.c]      - compute gradients of 2-electron integral for geometry
 |                        optimization or MECP calculation
 +- calls [quartet.c]   - the faster version of 2-electron integral which
 |                        unrolls all the loops
 +- calls [multipole.c] - the fastest but approximated version of 2-electron
                          integrals, used only the error is well below cut-off
                          value (1E-12)
```


### 4. Parallel Execution

Siam Quantum implements its own version of parallel execution without relying on standard platform such as MPI because we want the program to be very simple to install and because the calculations require very little data exchange so that parallelization is very easy to implement in this case.

Basically, when the parent process of Siam Quantum is executed, it prepares a text file for each child, and spawn the children with a specific option so that each child knows which text file it should read. The file contains a task description along with necessary initial data. Then, the parent goes into an indefinite idle loop waiting for all the children to complete the tasks. When each child completes its task, it writes output to a file with a specific name, the same filename that the parent is expecting to see. 

**[rpc.c]** is the source code that handles parallelization. Currently, only three tasks are available. They are 1) computing Coulomb-Exchange matrix, 2) evaluating the MP2 energy, and 3) calculating the gradient of 2-electron integral for geometry optimization.

