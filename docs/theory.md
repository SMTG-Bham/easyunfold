# Theory of unfolding

Theoretical background of the methodology used in this package is available in the literature[^1][^2].
Here, only a brief summary is given for brevity.

Notation: $\vec{k}$ is the kpoint of the supercell and $\vec{k}$ is that of the primitive cell. 

Reference: [https://link.aps.org/doi/10.1103/PhysRevB.85.085201](https://link.aps.org/doi/10.1103/PhysRevB.85.085201)  

Each $k$ in the primitive cell can be mapped to that in the supercell where:

$$
\vec{K} = \vec{k} - \vec{G}_0.
$$

with $\vec{G}_0$ being an reciprocal lattice vector.
Naturally, each $\vec{K}$ in the supercell can be unfolded into a set of $\vec{k}_i$:

$$
\vec{k_i} = \vec{K} + \vec{G}_i  i=1,...,N_{\vec{K}}
$$

The key implication is that for a given $\vec{k}$, there is a unique $\vec{K}$ that it folds to.

The goal of the band folding procedure is to obtain the $E(\vec{k})$ from the complicated $E(\vec{K})$, where E is the energy of the Kohn-Sham states. 
This can be achieved by projecting $\ket{\vec{K}m}$ on all of the primitive cell Block states $\ket{\vec{k}_i}$ of a set of $\vec{k_i}$ and compute the spectral weight:

$$
P_{\vec{K}m}(\vec{k}_i) = \sum_n |\braket{\vec{K}m |\vec{k}_i n}|^2.
$$

where $P$ represents the probability of finding a set of primitive cell stats $\ket{\vec{k}_in}$ contributing to the supercell state $\ket{\vec{K}m}$, or the amount of Bloch character $\vec{k}_i$ preserved in $\ket{\vec{K}m}$ at the same energy. 

Based on this, one can further derive the spectral function of $E$:

$$
A(\vec{k}_i, E) = \sum_m P_{\vec{K}m}(\vec{k}_i)\delta(E_m - E).
$$

In practice, the $\delta$ function is replaced with a Gaussian or Lorentzian function to smear the contribution with discretised $E$. 

Hence, the central quantity to be calculated is the $P_{\vec{K}m}(\vec{k}_i)$.
For plane-wave basis, it can shown that (equation 15 in Popescu et al.):

$$
P_{\vec{K}m}(\vec{k}_j) = \sum_{\vec{g}} |C_{\vec{Km}}(\vec{g} + \vec{G_j})|^2,
$$

where $C_{\vec{Km}}(\vec{g} + \vec{G_j})$ are the plane wave coefficients of the supercell.
Because the supercell is commensurate with the primitive cell, the vectors $\vec{g}$ are included in the basis set of the supercell calculations. 
Since $\vec{G_{j}}$ is the reciprocal lattice vector that wraps $\vec{k}$ to $\vec{K}$,
all plane wave coefficients that are needed to obtain $P_{\vec{K}_m}(\vec{k}_j)$ are present in the wave function output of the supercell calculation. 


## Symmetry considerations

The set of primitive cell kpoints $\vec{k}_i$ are located on a given path in the Brillouin zone, going through multiple high symmetry points. 
The selection of the path is not unique, but often depends on the space group of the primitive cell.
Only a limited set of paths are needed because of the presence of point group symmetry, as symmetrically equivalent $\vec{k}$ contains the same eigenstates. 

However, the supercell to be unfolded does not necessarily contain the same point group symmetry due to imperfections (presence of defects, strain, distortions).
The broken symmetry means that previously equivalent $\vec{k}$ are no longer equivalent. 

We address this by first expanding each $\vec{k}$ based on the symmetry operations of the primitive cell, followed by a reduction using the symmetry of the supercell.
The spectral weight at each $\vec{k}$ is then a weighted combination of that from a set of $\vec{k_s^\prime}$ points that are inequivalent under the symmetry of the supercell.

[^1]: Popescu, V.; Zunger, A. Effective Band Structure of Random Alloys. Phys. Rev. Lett. 2010, 104 (23), 236403. https://doi.org/10.1103/PhysRevLett.104.236403.
[^2]: Popescu, V.; Zunger, A. Extracting $E$ versus $\vec{k}$ Effective Band Structure from Supercell Calculations on Alloys and Impurities. Phys. Rev. B 2012, 85 (8), 085201. https://doi.org/10.1103/PhysRevB.85.085201.