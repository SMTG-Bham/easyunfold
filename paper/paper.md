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
to analyse directly. Band structure unfolding enables mapping the electronic structure
of such structure back to that of the primitive cell, helping research to gain
insights of the structure-property relationships and comparing effects of various crystal
imperfections with an equal ground.

# Statement of need

There are existing package that provide similar band structure unfolding capability, such as
[@bandup], [@vaspbandunfolding] using outputs of density functional theory calculation.
`easyunfold` is written in python with a focus on ease of use, data provenance and reproducibility.
It also includes plotting functionality for generating publication-ready figures.
A key feature `easyunfold` is to provide data serialization complying the FAIR principle, as input settings
and calculated spectral weights are stored as a single JSON file.
This enables the unfolded band structure to be re-plotted and further analysis without reprocessing the wave function data.
The multiple-file for wave function input allows the underlying DFT calculations to be performed in a flexible way, which is essential for improving parallel efficiency using resource-heavy hybrid functionals.

We chose Python as the programming language due to its ease of use, flexibility and popularity in the materials modelling field.
An object-oriented approach is used when designing the package to allows reading and handling wave function data to be abstracted away.
The code current supports two DFT codes, VASP and CASTEP, but others can be added without much effort.
`easyunfold` depends on common scientific computing packages such as `numpy` [@numpy] and `matplotlib` [@matplotlib].
The Atomic Simulation Environment [@ase] is used for reading input crystal structures.

`easyunfold` is designed to be used by users with or without prior knowledge of Python.
A command-line interface is provided as the primary way of use, and power users can use the Python
API directly for advanced analysis or customise the plotting further.
The package has been already been used in several scientific publications [@nicolson:2023, @wang:2022, @huang:2022], and as well as graduate-student research projects.
The combination of easy of use, flexibility, and efficiency will improve the accessibility of
band structure unfolding technique for materials modelling and help training new researchers.

# Theory

The math of band structure unfolding can be found in `@Popescu:2009`.

Each $k$ in the primitive cell can be mapped to that in the supercell where:

$$
\vec{K} = \vec{k} - \vec{G}_0.
$$

with $\vec{G}_0$ being an reciprocal lattice vector.
Naturally, each $\vec{K}$ in the supercell can be unfolded into a set of $\vec{k}_i$:

$$
\vec{k_i} = \vec{K} + \vec{G}_i, i=1,...,N
$$

The key implication is that for a given $\vec{k}$, there is a unique $\vec{K}$ that it folds to.

The goal of the band folding procedure is to obtain the $E(\vec{k})$ from the complicated $E(\vec{K})$, where E is the energy of the Kohn-Sham states.

This can be achieved by projecting $\langle \vec{K}m$ on all of the primitive cell Block states $\langle \vec{k}_i$ of a set of $\vec{k_i}$ and compute the spectral weight:

$$
P_{\vec{K}m}(\vec{k}_i) = \sum_n |\langle \vec{K}m |\vec{k}_i n \rangle |^2.
$$

where $P$ represents the probability of finding a set of primitive cell stats $\langle \vec{k}_in$ contributing to the supercell state $\langle \vec{K}m$, or the amount of Bloch character $\vec{k}_i$ preserved in $\langle \vec{K}m \rangle$ at the same energy.

The spectral function of $E$ is defined as:

$$
A(\vec{k}_i, E) = \sum_m P_{\vec{K}m}(\vec{k}_i)\delta(E_m - E).
$$

In practice, the $\delta$ function is replaced with a Guassian or Lorentzian function which smear the contribution with discretised $E$.

For plane wave basis, the $P_{\vec{K}m(\vec{k}_i)}$ can be shown to equal (equation 15 `@Popescu:2009`):

$$
P_{\vec{K}m}(\vec{k}_j) = \sum_{\vec{g}} |C_{\vec{Km}}(\vec{g} + \vec{G_j})|^2,
$$

where $C_{\vec{Km}}(\vec{g} + \vec{G_j})$ the plane wave coefficients of the supercell supercell calculation.

A further complication is that the supercell and the ideal primitive cells do not always match in full after geometry optimisation.
In addition, the reduction in symmetry means the $k$ points along the original path should map to multiple inequivalent points in the supercell.
`easyunfold` handles such symmetry breaking effect by including additional supercell $k$ points.
This involves first expanding the original kpoint with symmetry of the primitive cell, and  then have them reduced with the symmetry of the supercell.

# Acknowledgements

We acknowledge contributions from Adair Nicolson and help from these people testing the code and providing feedback: Joe Willis, Sabrine Hachmioune.

# References
