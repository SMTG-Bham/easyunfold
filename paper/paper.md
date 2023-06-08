---
title: 'easyunfold: A Python package for unfolding electronic band structures'
tags:
  - Python
  - materials science
  - first-principles
  - density functional theory
authors:
  - name: Bonan Zhu
    orcid: 0000-0001-5601-6130
    equal-contrib: false
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: Seán R. Kavanagh
    orcid: 0000-0003-4577-9647
    equal-contrib: false
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
  - name: David Scanlon
    orcid: 0000-0001-9174-8601
    equal-contrib: false
    affiliation: "1, 2, 3" # (Multiple affiliations must be quoted)
affiliations:
 - name: Department of Chemistry, University College London, London, United Kingdom
   index: 1
 - name: Thomas Young Centre, University College London, London, United Kingdom
   index: 2
 - name: The Faraday Institution, Didcot, United Kingdom
   index: 3
date: 6 June 2023
bibliography: paper.bib
---

# Summary

The electronic band structure is an important quantity to calculate in order
to analyse, understand and design solid crystalline materials in many fields
such as photovoltaic, catalysis, thermoelectric materials and transparent conducting
devices. Obtaining the band structure for a given crystal through first-principles
density functional theory calculation is a well-established routine operation.
However, the materials of interest are often complex and contains multiple primitive
cell of the ideal simple archetypical structure, for example, when site-occupation disordering,
cation/anion mutated structures and defects. The use of supercells in the real
space causes the band structure to be folded in the reciprocal space making it difficult
to analyse directly. Band structure unfolding facilitates mapping the electronic structure
back to that of the primitive cell, helping researchers understand
structure-property relationships and compare the effect of various crystal
imperfections on an equal footing.

# Statement of need

There are existing package that provide similar band structure unfolding capability, such as
[@bandup], [@vaspbandunfolding] using outputs of density functional theory calculation.
`easyunfold` is written in Python with a focus of easy-to-use, data provenance and reproducibility.
It also includes plotting functionalities and can generate publication-ready figures.
A key feature of `easyunfold` is to provide data serialization complying the FAIR principle.
Both the input settings and calculated spectral weights are stored as a single JSON file.
This enables the unfolded band structure to be re-plotted and further analysed without reprocessing the wave function data.
In addition, multi-file wave function support allows DFT calculations to be performed in a flexible way,
which is essential for improving parallel efficiency using resource-heavy hybrid functionals.

We chose Python as the programming language due to its low entry-barrier, flexibility and popularity in the materials modelling field.
An object-oriented approach is used when designing the package to allow abstractions when reading and process wave function data.
The code current supports two DFT codes, VASP and CASTEP, and others can be added without much effort.
`easyunfold` depends on common scientific computing packages such as `numpy` [@numpy] and `matplotlib` [@matplotlib].
The Atomic Simulation Environment [@ase] is used for reading input crystal structures.

`easyunfold` is designed for researchers with or without prior knowledge of Python.
A command-line interface is provided as the primary way of use, and power users can use the Python
API directly for advanced analysis or customise the plotting further.
The package has been already been used in several scientific publications [@nicolson:2023, @wang:2022, @huang:2022], and as well as graduate-student research projects.
The combination of easy of use, flexibility, and efficiency will improve the accessibility of
band structure unfolding technique for materials modelling and help training new researchers.

# Theory

The math of band structure unfolding is discussed in detail in @Popescu:2009,
and here we only give a brief summary.

Each $k$ point in the primitive cell's Brounline zone can be mapped to a $K$ in the supercell's Brouline where:

$$
\vec{K} = \vec{k} - \vec{G}_0.
$$

with $\vec{G}_0$ is a reciprocal lattice vector of the supercell.
Conversely, each $\vec{K}$ in the supercell can be unfolded into a set of $\vec{k}_i$:

$$
\vec{k_i} = \vec{K} + \vec{G}_i, i=1,...,N
$$

The key here is that for a given $\vec{k}$, there is a unique $\vec{K}$ that it folds to,
but the reverse it not true, e.g. a single $K$ may map to multiple $k$s, as the Brouline zone of the primitive cell covers that of the supercell.

The goal of the band folding is therefore to obtain the energy of Kohn-Sham states as a function of $\vec{k}$ from on the more complicated $E(\vec{K})$.

To achieve this, we project  $\langle \vec{K}m|$ on all of the primitive cell Bloch states $\langle \vec{k}_i|$ of a set of $\vec{k_i}$ and compute the spectral weight:

$$
P_{\vec{K}m}(\vec{k}_i) = \sum_n |\langle \vec{K}m |\vec{k}_i n \rangle |^2.
$$

where $P$ represents the probability of finding a set of primitive cell states $\langle \vec{k}_in$ contributing to the supercell state $\langle \vec{K}m |$,
or the amount of Bloch character $\vec{k}_i$ preserved in $\langle \vec{K}m \rangle$ at the same energy.

The spectral function of $E$ is defined as:

$$
A(\vec{k}_i, E) = \sum_m P_{\vec{K}m}(\vec{k}_i)\delta(E_m - E).
$$

In practice, the $\delta$ function is replaced with a Gaussian or Lorentzian function which smear the contribution with discretised $E$.

For plane wave basis, the $P_{\vec{K}m(\vec{k}_i)}$ can be shown to be (equation 15 @Popescu:2009):

$$
P_{\vec{K}m}(\vec{k}_j) = \sum_{\vec{g}} |C_{\vec{Km}}(\vec{g} + \vec{G_j})|^2,
$$

where $C_{\vec{Km}}(\vec{g} + \vec{G_j})$ is the plane wave coefficients of the supercell wave function.

A further complication in the process is that the imperfect supercell and the ideal primitive cells do not always match in full after geometry optimisation.
In addition, the reduction in symmetry means the $k$ points along the original path may map to multiple inequivalent points in the supercell.
`easyunfold` handles such symmetry breaking effect by including additional supercell $k$ points.
This involves first expanding the original kpoint with symmetry of the primitive cell, and  then have them reduced with the symmetry of the supercell.

# Acknowledgements

We acknowledge contributions from Adair Nicolson and help from these people testing the code and providing feedback: Joe Willis, Sabrine Hachmioune.

# References
