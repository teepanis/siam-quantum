# hessian.h and hessian.c

## Introduction

In addition to the geometry optimization during the Newton method, the Hessian is needed for frequency analysis of a molecule. Unlike the geometry optimization where the iterative Hessian from the BFGS scheme is sufficient to carry out the optimization, the frequency analysis requires a more accurate version of the Hessian.

Consider a molecule with its nuclei positions grouped into a parameter vector space:

$$
\vec{R} = \begin{bmatrix} x_1 \\ y_1 \\ z_1 \\ x_2 \\ y_2 \\ z_2 \\ \vdots \\ x_N \\ y_N \\ z_N \end{bmatrix}
$$

Here, $N$ is the number of atoms in the molecule. Using this notation, the Hessian matrix is the second derivative of the total energy with respect to the change of the nuclei positions.

$$
H_{ij} \equiv \frac{\partial^2 E}{\partial R_i \partial R_j}
$$

From this, there is an obvious symmeric property that $H_{ij} = H_{ji}$, so it can be used to reduce the computation by a factor of two. The total energy $E$, also includes the nuclei Coulomb repulsion energy, can be evaluated using electronic structure theory such as Hartree–Fock or DFT. In principle, the first and the second derivative are also analytically possible; but the second derivative can be quite time consuming, more so in large molecules.

There is also complexity regarding the analytical expression for the Hessian. Hence, presently Siam Quantum uses numerical evaluation of the second derivative. The values are tested against the analytical ones from GAMESS to a satisfactory agreement.

## subroutine **hessian_numerical(...)**

```c
void hessian_numerical(
        int dbSize,                         // number of record in basis set db
        const struct GTOBasisSet_t *basisDB,// pointer to basis set database
        int nBasis,              // number of basis functions
        struct Molecule_t * mol, // pointer to molecule structure
        int nEA,                 // total number of spin up electrons
        int nEB,                 // total number of spin down electrons
        const double *CA,        // molecular alpha spin orbital
        const double *CB,        // molecular beta spin orbital 
        const double *eA,        // eigen values
        const double *eB,        // eigen values
        struct option_t *opt){   // options
```

This subroutine computes the Hessian numerically using center finite-difference method, giving $O(h^2)$ accuracy, where $h$ is the step size. For example,

$$
H_{ij} \equiv \frac{\partial^2 E}{\partial R_i \partial R_j} \approx \frac{1}{2 h}\bigg( \frac{\partial E}{\partial R_i}\Big|_{R_j + h} -  \frac{\partial E}{\partial R_i}\Big|_{R_j - h} \bigg)
$$

The gradient $g_i = \frac{\partial E}{\partial R_i}$ is evaluated _analytically_ and readily available during routine geometry optimization. The step size $h=$`HESSIAN_STEPSIZE`, however, requires calibration.

If $h$ is too large, it obviously gives poor estimation. But if $h$ is too small, it produces numerical instability even though the gradient $\vec{g}$ are computed analytically. Various factors such as integration cut-off and approximation can contribute to this instability. During the test, we found the following behavior:

```c
// numerical step in Bohr
#define HESSIAN_STEPSIZE        0.008
// 0.005 causes error for the lowest frequency of CH3 HF/6-31G*
//              the analytical freq is 308, but numerical one is 314
//              setting it 0.01 brings the numerical freq to 308
//
// 0.01  causes non-degeneraacy of some modes in CH4/6-31G* at
//              frequency around 1703 >> 1703 and 1702
//              setting it to 0.008, restores the degeneracy
//
// 0.008 seems to be the sweet spot
//
// Teepanis Chachiyo 27 March, 2020
```

As a result, the optimal value for $h = 0.008$ Bohr is used presently.

## subroutine **hessian_output(...)**

```c
void hessian_output(
        int nDim,
        const double *H,
        const struct Molecule_t *mol){
```

This subroutine prints out formatted Hessian so a simple script can extracted the information from the output. In addition, frequency analysis of the molecule is performed and the Zero-Point vibration energy is reported.

To perform normal mode analysis, first the _mass-weighted_ Hessian is computed as follows:

$$
\mathcal{H}^{}_{ij} \equiv \frac{1}{\sqrt{m_i m_j}} H_{ij} \quad \times \quad \frac{1}{\sqrt{(\text{AMU2KG})^2}}\frac{\text{HARTREE2J}}{(\text{BOHR2ANGSTROM}\times 10^{-10})^2}
$$

Here, $m_i, m_j$ are the mass of the most abundant isotope of the nuclei corresponding to the coordinate parameter $i, j$ respectively. Care must be taken to handle the unit properly from the typical atomic unit in the code to the SI unit.

Once the mass-weighted Hessian $\widetilde{\mathcal{H}}$ matrix is calculated, we can use the package `LAPACK` to solve for this eigenvalue problem.

$$
\widetilde{\mathcal{H}} \vec{c}_n = \lambda_n \vec{c}_n
$$

In general, there are as many eigenvalues and eigenvectors as the dimension of the parameter space, which is `nDim`$=3N$. 

The oscillation frequency in $\text{cm}^{-1}$ can be calculated as follows:

$$
f_n = \frac{1}{2\pi} \sqrt{\lambda_n} \quad \times \quad 10^{-6}\cdot\text{MHZ2CM}
$$

Summing all the frequencies would give the _Zero-Point Energy_ in Hartree.

$$
E_{ZPE} =  \sum_n{ \frac{1}{2} f_n} \quad \times \quad \text{CM2HARTREE}
$$

Here, all conversion factors are defined in `hessian.h` as follows:

```c
 // CRC Handbook of Chemistry and Physics 85th Edition
#define HARTREE2J 4.35974417E-18  // conversion hartree to joule
#define AMU2KG    1.66053886E-27  // conversion u to kilogram
#define AMU2ME    1.82288848E+3   // conversion u to atomic unit
#define ME_KG     9.10938262E-31  // electron mass in kg
#define MHZ2CM    3.33564E-5      // conversion MHz to cm^-1
#define CM2HARTREE 4.556335E-6    // conversion cm^-1 to hartree
#define CM2KCALMOL 2.85914E-3     // conversion cm^-1 to kcal/mol
```

In addition to the electronic energy, the Zero-Point Energy is an important factor when analyzing chemical reaction because the bond energy can change appreciably due to the vibration energy.