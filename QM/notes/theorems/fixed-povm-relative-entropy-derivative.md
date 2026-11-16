# Fixed and optimized measured relative entropy differentials

[Proved] Direct finite-sum differentiation, not a literature claim.

Assumptions: finite dimension, finite fixed POVM, natural logarithms, strictly positive outcome probabilities locally. Full rank is sufficient, not necessary. Physical tangent matrices are Hermitian and traceless.

For F(rho,sigma) = sum p_i log(p_i/q_i), the second-argument differential is -sum (p_i/q_i) Tr(E_i B). At rho=sigma it vanishes on trace-zero directions (ambient gradient -I).

For R=U rho U† and S_epsilon=U exp(-i epsilon G) rho exp(i epsilon G) U†, Delta=-i[UGU†,R]. At fixed time, F(R,S_epsilon) = epsilon²/2 sum Tr(E_i Delta)²/p_i + O(epsilon³). The first derivative is zero; the quadratic response is classical outcome Fisher information and may vanish for a measurement-invisible tangent.

Distinguish varying the second density matrix from varying the initial density matrix in both arguments. The latter requires both terms of the joint differential. Both derivatives vanish at coincident arguments on normalized states. Boundary probabilities need separate analysis; no uniform long-time expansion claimed.

Full derivation: `latex/fixed_povm_relative_entropy_derivative.tex`.

**Deliverable scope update (2026-11-15):** At the user's request, the LaTeX has been condensed to the requested fixed-POVM derivatives, equal-final versus equal-initial direction correction, substitutions before/after evolution, optimized derivatives, unitary time invariance, and the data-processing upper-bound proof. Repeated explanations, epsilon expansions, SLD/Fisher-information material, and the extended optimizer existence proof were removed from the deliverable. Older statements below describing those derivations as present in the LaTeX are historical; the mathematical research notes remain retained here. The upper bound is D(rho||V rho V†), finite for full-rank finite-dimensional states, and applies to any common record-producing channel. It does not automatically bound derivative graphs. No unrelated LaTeX file was changed.

User clarification: the requested variable is the initial rho in BOTH arguments, at fixed t and epsilon. The same initial direction A pushes forward to U_t A U_t† and U_t U_epsilon A U_epsilon† U_t†, not generally the same final direction. The LaTeX now explicitly emphasizes this derivative. With P_i=U_t† E_i U_t and Q_i=U_epsilon† P_i U_epsilon, its fixed-POVM gradient representative is sum [P_i log(p_i/q_i) - (p_i/q_i) Q_i]. For optimization, with Z_* the optimizer before common time evolution, the representative is log Z_* - U_epsilon† Z_* U_epsilon. Both act by trace pairing with traceless A; neither generally vanishes at fixed nonzero epsilon. These are not derivatives in epsilon.

## Fixed-POVM second differential

[Proved] For final-state directions (A_1,A_2) and (B_1,B_2), put a_i=Tr(E_i A_1), b_i=Tr(E_i A_2), c_i=Tr(E_i B_1), d_i=Tr(E_i B_2). The Hessian is sum [a_i c_i/p_i - (a_i d_i+b_i c_i)/q_i + p_i b_i d_i/q_i²]. Its diagonal is sum p_i(a_i/p_i-b_i/q_i)². Literal equal final directions (A,A) yield sum Tr(E_i A)²(p_i-q_i)²/(p_i q_i²). This is different from equal INITIAL directions, which require the two unitary pushforwards. Both vanish at epsilon=0 in the intended setup. No epsilon differentiation or curve acceleration is involved. Added to the same standalone LaTeX; compiler still unavailable.

## Optimized single-copy POVMs

[Proved] For full-rank states in finite dimension, the optimized divergence equals max over positive W of Tr(rho_1 log W) + 1 - Tr(rho_2 W). A self-contained derivation using classical variational relative entropy and operator Jensen is in the LaTeX. Strict concavity and interior compact level sets give a unique optimizer W_*. Its stationarity equation is Dlog(W_*)[rho_1] = rho_2. The negative-definite Hessian justifies smooth local optimizer dependence and the envelope differential.

On normalized variations (A,B), the joint differential is Tr(A log W_*) - Tr(B W_*). At equality W_*=I and the differential vanishes. The optimized local quadratic coefficient is SLD Fisher information. A uniform-in-POVM Taylor remainder justifies interchanging optimization and expansion; the SLD spectral measurement attains the coefficient. Common unitary evolution cancels exactly, so the optimized quantity has no time dependence in this setup.

Full rank matters: distinct pure states give infinite optimized measured relative entropy, even for arbitrarily small rotations. No unrestricted boundary differentiability is claimed.

## Optimized initial-state Hessian (explicit POVM probabilities)

[Proved] For V=U_epsilon fixed, p_i=Tr(M_i rho), q_i=Tr(M_i V rho V†), with M_i POVM effects, the initial-state differential is the usual joint formula evaluated at an optimal POVM, with a_i=Tr(M_i A) and b_i=Tr(M_i V A V†). It is NOT the supremum of directional derivatives.

[Proved] Let Z>0 be the unique optimizer of Tr(rho log Z)+1-Tr(V rho V† Z). Write R_s=(Z+sI)^(-1), T_Z(X)=integral R_s X R_s ds and L(H)=integral (R_s H R_s rho R_s+R_s rho R_s H R_s) ds. Stationarity is T_Z(rho)=V rho V†. L is positive definite and self-adjoint on Hermitian matrices. With C_A=T_Z(A)-V A V†, optimizer response H=DZ[A] solves L(H)=C_A. Thus DF[A]=Tr[A(log Z-V† Z V)] and D²F[A,A]=Tr[C_A L^(-1)(C_A)] >=0. This is an exact implicit formula, not a numerical computation. At V=I, Z=I and C_A=0. Common unitary evolution cancels from both differentials.

Spectral projections of Z form an optimal POVM, with its eigenvalues equal to the corresponding p_i/q_i. Z varies smoothly even where a spectral-projector choice is not smooth. The full derivation was appended to `latex/fixed_povm_relative_entropy_derivative.tex`; compilation remains unavailable.

## Deliverable correction and verification

[Proved; added draft check] Two further sections in the standalone LaTeX check the user's draft and explicitly substitute all probabilities in both the unevolved and time-evolved fixed-POVM cases. The common-final-direction Hessian is correct, but is not generally the initial-state Hessian, even at t=0. The unitary state map differentiates to A -> U A U†, not U rho U†. The +1 term drops precisely because Tr A=0. With F_i=U_epsilon† E_i U_epsilon (and time-dependent pulled-back analogues), the initial-state Hessian is sum Tr(E_i rho)[Tr(E_i A)/Tr(E_i rho)-Tr(F_i A)/Tr(F_i rho)]². At epsilon=0 both initial-state differentials vanish. No extra assumptions permit cancellation of the fixed-POVM time dependence.

The original LaTeX was only a section fragment despite the standalone-deliverable requirement. It has now been replaced, at the user's request, by a standalone article with amsmath/amssymb, using rho_1, rho_2, U_t and U_epsilon consistently. The right perturbation unitary is explicitly corrected to its adjoint. Both second-state and initial-state derivatives are included.

No TeX compiler was found on PATH (pdflatex, tectonic, latexmk, xelatex, lualatex checked); the replacement has NOT been compile-tested.

Local source screening: PDF page 1 of `papers/frechet_derivative.pdf`, `papers/1803.00665v2.pdf`, and `papers/2212.01585v3.pdf` was extracted with pypdf to identify topics. None is used as a source for the optimized derivative theorem here; the derivation is given explicitly instead. No new bibliographic claims or downloads.
