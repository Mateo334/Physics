# Exact quantum Hamiltonian geometry versus variational flow

Verified 2026-11-24. Source-focused companion to `../theorems/observational-entropy-qfi-pesin.md`. Two open PDFs downloaded without overwriting files. PDF text read with installed `pypdf`; Python 3.11 has no `fitz`/PyMuPDF, and nothing was installed. Page numbers below are PDF pages, also the printed page numbers in these copies. This is a targeted reading, not a full review of either paper.

## Verified sources and locations

### Ashtekar–Schilling

Abhay Ashtekar and Troy A. Schilling, *Geometrical Formulation of Quantum Mechanics*, arXiv:gr-qc/9706069v1 (1997), https://arxiv.org/abs/gr-qc/9706069. Later a chapter in *On Einstein's Path*, DOI https://doi.org/10.1007/978-1-4612-1422-9_3.

Local PDF: `../../papers/Ashtekar_Schilling_1997_Geometrical_Formulation_Quantum_Mechanics.pdf`; download https://arxiv.org/pdf/gr-qc/9706069v1.

- **Pp. 6–8, Eqs. (2.1)–(2.5):** real Hilbert metric and symplectic form are 2ħ times the real and imaginary inner product. Multiplication by i is J. Convention is contraction of the Hamiltonian vector into the symplectic form equals dh. Schrödinger evolution is precisely the Hamiltonian flow of the expectation value.
- **Pp. 11–12, Definition II.1 and Eq. (2.14):** the construction descends to projective Hilbert space; the Poisson bracket reproduces the operator commutator. The globally smooth observable discussion explicitly uses bounded self-adjoint operators.
- **Pp. 13–14, Eqs. (2.16)–(2.21), Corollary 1:** projective metric; quantum-observable Hamiltonian vector fields are Killing. Eq. (2.18) relates their covariant derivative to the Hessian. Corollary 1 characterizes observable functions by the Killing condition. Thus “Hamiltonian” alone is weaker than “ordinary quantum evolution.”
- **Pp. 33–34, Eqs. (4.7)–(4.10), §IV B:** expectation coordinates following classical oscillator trajectories do not imply invariance of a chosen coherent-state section. The energy includes uncertainty terms. For the oscillator the ground-state-generated standard coherent section is invariant; general displaced fiducial-state sections need not be.
- **P. 7, immediately after (2.3):** unbounded generators have only densely defined vector fields although their unitary flow exists on all Hilbert space. Do not use the smooth bounded-observable statements indiscriminately for an inverted oscillator.

### Hackl et al.

Lucas Hackl, Tommaso Guaita, Tao Shi, Jutho Haegeman, Eugene Demler and J. Ignacio Cirac, *Geometry of variational methods: dynamics of closed quantum systems*, **SciPost Physics 9**, 048 (2020), DOI https://doi.org/10.21468/SciPostPhys.9.4.048; arXiv:2004.01015.

Local published PDF: `../../papers/Hackl_et_al_2020_Geometry_Variational_Methods_Closed_Quantum_Systems.pdf`; download https://scipost.org/SciPostPhys.9.4.048/pdf.

- **Pp. 26–28, Eqs. (81)–(96), Propositions 4–6:** normalized projective variational action, tangent projection, Hamiltonian parameter flow, energy conservation, and equivalence of Lagrangian/McLachlan/Dirac–Frenkel prescriptions for Kähler families. For non-Kähler families they generally differ; degeneracy of the pulled-back symplectic form needs separate treatment.
- **P. 27, Eq. (88):** parameter velocity is minus the inverse symplectic matrix times the energy gradient. The authors explicitly warn about the symplectic-sign convention. Their units set ħ=1.
- **P. 35, Eqs. (122)–(123) and footnote 15:** linearization at a stationary point is minus inverse symplectic matrix times energy Hessian. Connection dependence drops out *at the fixed point*. The discarded derivative-of-symplectic-matrix term vanishes because the energy gradient vanishes there.
- **P. 36, Eqs. (124)–(125) and following paragraphs:** the full unitary flow preserves all Kähler structures, whereas the variational projected flow generally preserves only the symplectic structure, even on a Kähler variational manifold. This directly verifies the distinction needed here.
- **Pp. 36–37, Eqs. (126)–(128):** differentiating the projected field includes the derivative of the tangent projector. Linearizing TDVP is not merely projecting the Hamiltonian onto tangent space.

**Version pitfall:** the arXiv HTML search rendition has different numbering (e.g. its Hamiltonian-flow equation is (98), while the downloaded published PDF has (88)). Cite the PDF numbering above consistently.

### Heslot: bibliography verified, full text not verified

André Heslot, *Quantum mechanics as a classical theory*, **Physical Review D 31**, 1341–1348 (1985), DOI https://doi.org/10.1103/PhysRevD.31.1341. Publisher page https://journals.aps.org/prd/abstract/10.1103/PhysRevD.31.1341 verifies author, title, date (15 March 1985), and abstract; access requires authorization. No accessible full-text copy was located in this bounded search. Ashtekar–Schilling p. 5 credits its finite-dimensional projective symplectic formulation. No Heslot page-level theorem or equation is claimed checked.

## Rigorous deductions, with explicit conventions

These formulas are derived here from the verified geometry; they are not verbatim numbered equations of the sources. Work in finite dimension, or on a regular common domain when differentiating unbounded generators. For a normalized representative ψ, horizontal vectors u, v satisfy orthogonality to ψ. Use Ashtekar–Schilling's scaled metric and symplectic form:

```latex
\begin{aligned}
g(u,v)&=2\hbar\,\operatorname{Re}\langle u|v\rangle,\\
\omega(u,v)&=2\hbar\,\operatorname{Im}\langle u|v\rangle,\\
Jv&=iv,\qquad \iota_{X_h}\omega=dh,\\
h([\psi])&=\langle\psi|\hat H|\psi\rangle,\\
X_h&=-J\operatorname{grad}_g h\\
&=-\frac{i}{\hbar}(\hat H-h)\psi .
\end{aligned}
```

For a smooth variation of initial rays, ξ is the pushed-forward initial tangent. The Levi-Civita connection is torsion-free and the Kähler complex structure is parallel. Define the Hessian endomorphism by raising one index with g:

```latex
\begin{aligned}
D_t\xi&=\nabla_\xi X_h\\
&=-J(\operatorname{Hess}_g h)^\sharp\xi,\\
(\operatorname{Hess}_g h)^\sharp\xi
&=\nabla_\xi\operatorname{grad}_g h .
\end{aligned}
```

This is a **first-order variational equation**, not the geodesic Jacobi equation; Schrödinger orbits need not be geodesics. It also holds for explicitly time-dependent h at fixed external time. In other sign/metric conventions its displayed prefactors change.

For expectation-value Hamiltonians the field is Killing, so the covariant linearization is skew-adjoint in g. Therefore:

```latex
\begin{aligned}
\frac{d}{dt}g(\xi,\xi)
&=2g(\nabla_\xi X_h,\xi)=0,\\
\|D\Phi_{t,s}\xi\|_{\rm FS}
&=\|\xi\|_{\rm FS},\\
\lambda_{\rm FS}(\xi)
&=\lim_{t\to\infty}\frac1t
\log\frac{\|D\Phi_{t,0}\xi\|_{\rm FS}}
{\|\xi\|_{\rm FS}}=0 .
\end{aligned}
```

Here the usual dimensionless FS metric is g/(2ħ); multiplying a metric by a fixed positive constant does not change this conclusion. Equivalently, the FS distance arccos of the absolute normalized overlap is exactly preserved by a common unitary. This holds also for time-dependent Hamiltonians and unitary Floquet maps. It does **not** assert absence of quantum-chaos signatures or constancy of sensitivity to a parameter in the Hamiltonian itself.

## Why coherent/variational dynamics can differ

Let the embedded coherent/variational manifold M have nondegenerate pulled-back symplectic form and restricted energy E. Its action-based Hamiltonian dynamics is exact quantum restriction only if the full quantum vector field is tangent to M everywhere along the relevant trajectories. For a Kähler submanifold the variational vector field is its metric tangent projection:

```latex
\begin{aligned}
X_{\rm var}&=P_{TM}X_h,\\
X_h|_M\in TM
&\quad\Longrightarrow\quad
\Phi_t|_M\ \text{preserves the induced metric}.
\end{aligned}
```

Consequently an *exactly invariant* coherent manifold cannot have positive exponents in its induced FS metric. If projection discards a normal component, the resulting Hamiltonian dynamics need not be Killing and may have stretching. On a Kähler M its covariant tangent equation still has the Hessian form above, now with E and the induced geometry, but its generator need not be skew-adjoint. On a general non-Kähler symplectic manifold, use the full covariant derivative of the Hamiltonian field; parallel J is unavailable.

Pitfalls for the main calculation:

1. A coordinate Hessian alone is not the covariant Hessian away from a stationary point; Christoffel terms matter. In a constant canonical Darboux matrix chart, ordinary coordinate linearization is the familiar symplectic matrix times ordinary Hessian, but its Euclidean norm is not automatically FS.
2. Symplecticity preserves phase-space volume, not lengths. Conversely, norm preservation of each state alone is insufficient; exact unitarity preserves **mutual** overlaps.
3. A fixed-width coherent manifold and the larger squeezed-Gaussian manifold are different. Quadratic dynamics may preserve the latter while leaving the former. For an inverted oscillator, exponentially separated centers with evolving covariance do not imply exponential FS separation under the same unitary.
4. In noncompact/infinite-dimensional settings, a center-coordinate Euclidean metric need not be uniformly equivalent to the induced FS metric along an escaping/squeezing orbit. Do not transfer Lyapunov exponents between them without bounds.
5. A positive exponent of an approximate classical/TDVP flow is not automatically an exponent of the exact quantum state flow. A controlled semiclassical comparison requires a time window, norm, error estimate and order of limits.

Scope: no simulations or Gaussian-covariance derivation performed here; full Heslot text and PDF visual equation rendering were not checked. The cited PDF equations were checked by text extraction, cross-reading surrounding statements, and independent sign-convention derivation.
