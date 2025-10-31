# matrix.h and matrix.c

## Definitions

Here, we lists most often used variables to represent matrix elements. The naming convention for accessing elements in the matrix is `[i*nBasis+j]`, meaning the i<sup>th</sup> row, and j<sup>th</sup> columns.

Also, the chemist's notation for electron integral is particularly useful for writing an integral in a compact form.

$$(a|\hat{O}|b) \equiv \int d^3 r \, a(\vec{r}) b(\vec{r}) \hat{O}  b(\vec{r})$$

Here, $a(\vec{r})$ and $b(\vec{r})$ are any functions, and $\hat{O}$ is an operator acting on $b(\vec{r})$. In some cases, like for electron repulsion energy, a double integral is needed. The chemist's notation reads,

$$(ab|cd) \equiv  \iint d^3 r_1 d^3 r_2 \, a(\vec{r}_1)b(\vec{r}_1)\frac{1}{r_{12}}c(\vec{r}_2)d(\vec{r}_2)$$


### Density Matrix 

```PA[i*nBasis+j] and PB[i*nBasis+j]```  is the density matrix for alpha and beta electron, denoting $P^{(\alpha)}_{ij}, P^{(\beta)}_{ij}$.

A density matrix is the sum over occupied orbitals of the product between the molecular coefficients.

$$P_{ij} = \sum_{k \in \text{occ}} C_{ki} C_{kj}$$

Its direct usage is the computation of the electron density,

$$\rho(\vec{r}) = \sum_{k \in \text{occ}}|\phi_i(\vec{r})|^2 = \sum_{ij}P_{ij}\chi_i(\vec{r})\chi_j(\vec{r})$$

Because the density matrix represent the density of electron, it is also used predominantly to compute electronic energies.

### Core Hamiltonian Matrix

```H[i*nBasis+j]``` the the matrix element of the basis function embracing over the core hamiltonian operator, denoting $H_{ij}$.

$$H_{ij} = (\chi_i|-\frac{1}{2}\nabla^2+\hat{v}_\text{ext}|\chi_j)$$

The **external** potential means potential coming for interaction **other than** the electron themselves. Generally, they are the nuclei attraction potential, but external electric field is also possible.

$$\hat{v}_\text{ext}(\vec{r}) = \sum_{A \in \text{nuclei}} - \frac{Z_A}{|\vec{r}-\vec{R}_A|}$$

Also, the naming convention for denoting the nuclei index is the upper roman letter (i.e., A, B).

### Coulomb Potential Matrix

```JT[i*nBasis+j]```  is the matrix element over Coulomb potential $\hat{v}_J$ due to electron repulsion, denoteing $J^{(T)}_{ij}$

The ```T``` signifies that **total** electron density because the interaction does not distinguish between the alpah and beta elecron, using only the total is sufficient.

$$\hat{v}_J(\vec{r_1}) = \int d^3 r_2 \frac{\rho(\vec{r}_2)}{|\vec{r}_1 - \vec{r}_2|} $$

Therefore, the matrix element is,

$$J^{(T)}_{ij} = (\chi_i|\hat{v}_J|\chi_j)$$

If we expand the density $\rho(\vec{r}_2)$ in terms of the density matrix, we have:

$$J^{(T)}_{ij} = \sum_{pq}P^{(T)}_{pq}(\chi_i\chi_j | \chi_p \chi_q)$$

Note that we are using the chemist's notation to write the integrals.

### Hartree-Fock Exchange Matrix

```KA[i*nBasis+j] and KB[i*nBasis]``` are the matrix element for the HF exchange, denoting $K^{(\alpha)}_{ij},K^{(\beta)}_{ij}$.  

The elements are different between the alpha and beta electrons, so we need the two variants for computations.

$$\begin{align*}
K^{(\alpha)}_{ij} & = -\sum_{pq}P^{(\alpha)}_{pq}(ip|jq) \\
K^{(\beta)}_{ij} & = - \sum_{pq}P^{(\beta)}_{pq}(ip|jq)
\end{align*}
$$


