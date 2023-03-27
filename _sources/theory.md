# Theory of unfolding


Notation: $\vec{k}$ is the kpoint of the supercell and $\vec{k}$ is that of the primitive cell. 
Reference: [https://link.aps.org/doi/10.1103/PhysRevB.85.085201](https://link.aps.org/doi/10.1103/PhysRevB.85.085201)  

Each $k$ in the primitive cell can be mapped to that in the supercell where:

$$
\vec{K} = \vec{k} - \vec{G}_0.
$$

with $\vec{G}_0$ being an reciprocal lattice vector.
Naturally, each $\vec{K}$ in the supercell can be unfolded into a set of $\vec{k}_i$:

$$
\vec{k_i} = \vec{K} + \vec{G}_i, i=1,...,N_\vec{K}.
$$

The key implication is that for a given $\vec{k}$, there is a unique $\vec{K}$ that it folds to.

The goal of the band folding procedure is to obtain the $E(\vec{k})$ from the complicated $E(\vec{K})$, where E is the energy of the Kohn-Sham states. 
This can be achieved by projecting $\ket{\vec{K}m}$ on all of the primitive cell Block states $\ket{\vec{k}_i}$ of a set of $\vec{k_i}$ and compute the spectral weight:

$$
P_{\vec{K}m}(\vec{k}_i) = \sum_n |\braket{\vec{K}m |\vec{k}_i n}|^2.
$$

where $P$ represents the probability of finding a set of primitive cell stats $\ket{\vec{k}_in}$ contributing to the supercell state $\ket{\vec{K}m}$, or the amount of Bloch character $\vec{k}_i$ preserved in $\ket{\vec{K}m}$ at the same energy. 
Based on this, one can further dervice the spectral function of $E$:

$$
A(\vec{k}_i, E) = \sum_m P_{\vec{K}m}(\vec{k}_i)\delta(E_m - E).
$$

In practice, the $\delta$ function is replaced with a Gaussian or Lorentzian function to smear the contribution with discretised $E$. 

Hence, the central quantity to be calculated is the $P_{\vec{K}m}(\vec{k}_i)$.
For plane-wave basis, it can shown that (equation 15 in Popescu et al.):

$$
P_{\vec{K}m}(\vec{k}_j) = \sum_{\vec{g}} |C_{\vec{Km}}(\vec{g} + \vec{G_j})|^2,
$$

where $C_{\vec{Km}}(\vec{g} + \vec{G_j})$ the plane wave coefficients of the supercell. 


## Symmetry considerations

In practice, the set of primitive cell kpoints $\vec{k}_i$ is taken from a given path in the Brouline zone going through multiple high symmetry points. 
This selection of the path depends on the space group of the primitive cell.
Only a limited set of paths are needed because of the presence of point group symmetry, as symmetrically equivalent kpoints contains the same eigenstates. 

However, the supercell to be unfolded often contains reduced point group symmetry compared to that of the primitive cell due to imperfections (presence of defects, strain, distorations e.t.c).
This broken symmetry means that previously equivalent kpoints are no longer equivalent. 

We account of this by including images of the primitive kpoints that are no longer equivalent under the point group of that of the supercell. 
Weighted contributions of these extra kpoints are added to the spectral function. 
