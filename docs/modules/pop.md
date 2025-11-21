# pop.h and pop.c

## Introduction

This module computes and prints simple molecular properties such as Mulliken population, dipole moment, and other electrostatic properties. It is mostly based on the electron density and its relation to the density matrix as follows:

$$
\rho(\vec{r}) = \sum_{ij}P_{ij}\chi_i{(\vec{r})}\chi_j(\vec{r})
$$

## subroutine **mulliken(...)**

```c
void mulliken(
        int nBasis,               // number of basis function
        struct GTOBasis_t * gto,  // pointer to basis set information
        struct Molecule_t * mol,  // pointer to molecular structure
        int nEA,                  // number of spin up electron
        int nEB,                  // number of spin down electron
        double *CA,               // calculated molecular orbital coeff.
        double *CB,               // save but for beta spin
        double *eA,               // spin up eigen value
        double *eB,               // spin down eigen value
        struct option_t *option){ // print according to options
```

This subroutine takes the molecular orbital coefficients $C_{ij}^{(\alpha)}$ and $C_{ij}^{(\beta)}$ can compute the density matrix for each spin as follows:

$$
\begin{aligned}
P^{(\alpha)}_{ij} & = \sum_{k \in \text{occ}} C_{ki}^{(\alpha)}C_{kj}^{(\alpha)} \\
P^{(\beta)}_{ij} & = \sum_{k \in \text{occ}} C_{ki}^{(\beta)}C_{kj}^{(\beta)}
\end{aligned}
$$

Then lead to two types of density matrices, the total and the spin.


$$
\begin{aligned}
P^{(T)}_{ij} & = P^{(\alpha)}_{ij} + P^{(\beta)}_{ij} \\
P^{(S)}_{ij} & = P^{(\alpha)}_{ij} - P^{(\beta)}_{ij} 
\end{aligned}
$$

These two types of density matrices are needed to calculate Mulliken _charge_ population ans Mulliken _spin_ population.

The Mulliken charge population is defined as follows:

$$
M^{(\text{charge})}_A = \sum_{i = A} \sum_{j} P^{(T)}_{ij} S_{ij} + Z_A
$$

On the left-hand-side, $M^{(\text{charge})}_A$ is the Mulliken charge population for the atom at position $A$. On the right-hand-side, we see summations over basis function indices $i$ and $j$.

Note that each basis function is specific to its own center, typically assigned to its corresponding atom it is representing. The summation condition $i \in A$ means that only the basis functions whose centers are located at $A$ are to participate in the summation. The last term on the right-and-side is the nuclei contribution of the charge.

The Mulliken spin population is quite similar, except that it has no nuclei contribution.

$$
M^{(\text{spin})}_A = \sum_{i = A} \sum_{j} P^{(S)}_{ij} S_{ij}
$$


## subroutine **electric_multipole(...)**

```c
void electric_multipole(
        int nBasis,               // number of basis function
        struct GTOBasis_t * gto,  // pointer to basis set information
        struct Molecule_t * mol,  // pointer to molecular structure
        int nEA,                  // number of spin up electron
        int nEB,                  // number of spin down electron
        double *CA,               // calculated molecular orbital coeff.
        double *CB,               // save but for beta spin
        double *eA,               // spin up eigen value
        double *eB,               // spin down eigen value
        struct option_t *option){ // print according to options
```

Electric multipole moment (e.g., dipole, quadrupole) reveals important chemical properties of molecules. For example, the fact that water is polar can be seen clearly from its electric dipole moment; whereas benzene is non-polar, hence its electric dipole moment is zero.

At a higher level of computation, electric dipole moment is quite easy to compute, leaving the tedious mathematical matrix element calculations in the `matrix.c` module.

$$
\begin{aligned}
p_x & = -\sum_{ij} ( P^{(\alpha)}_{ij}+P^{(\beta)}_{ij})(\chi_i|\hat{x}|\chi_j) + \sum_{A} Z_A x_A \\
p_y & = -\sum_{ij} (P^{(\alpha)}_{ij}+P^{(\beta)}_{ij})(\chi_i|\hat{y}|\chi_j) + \sum_{A} Z_A y_A  \\
p_z & = -\sum_{ij} (P^{(\alpha)}_{ij}+P^{(\beta)}_{ij})(\chi_i|\hat{z}|\chi_j) + \sum_{A} Z_A z_A 
\end{aligned}
$$

Here, $\vec{p} \equiv (p_x, p_y, p_z)$ is the dipole moment vector. The matrix element, for example, $(\chi_i|\hat{x}|\chi_j)$ is the expectation value of the position $\hat{x}$.

Notice how the electron contribution (e.g. the first summation) is negative and the nuclei contribution is positive. The summation index $A$ runs through all atoms in the molecule.

