# Tasks
## General task
- Your goal is to derive a quantum version of Pesin theorem. 
- If not possible, focus on the study of these various entropies and its connection to chaos via OTOC, Quantum lyapunov exponents, POVM, density matrices, subspaces and others. 
- In the notes.md file in each subfolder you are currently working in, mateo will regularly write what you should focus on. If there is none, create one
- Reproduce all derivations step by step in LaTeX, without skipping any steps.
- Attempt to derive something new or extend the result — a sharper bound, a generalisation, a connection to another entropy measure, or an application to a specific class of quantum channels.
- Try to work with POVM, CP maps, OTOC, anything that is connected even remotely to this concept.
- Use small easy systems for fact checking - include the analytical matrices and the derivation in the latex output.
- For now we will focus on the CNT entropy. Create another folder in KS, name it CNT and focus on the basics of CNT entropy
- [x] CNT basics: folder created, Output.tex written (8 sections), Python verification complete (2026-06-01)
- [x] CNT deeper: proved h_AFL - h_CNT = log d exactly for all clustering states; verified numerically for XXX chain at 5 inverse temperatures; Section 10 added to Output.tex (2026-06-01)
- [x] CNT extension: proved h_CNT(alpha_t)=0 for finite systems; LR-Pesin bound h_CNT<=v_LR*log(d)*|t|; modular flow identified with imaginary-time evolution; Section 11 added to Output.tex (2026-06-01)
- [x] CNT vs quantum chaos: proved diagonal-U OPU triviality lemma; initial slope Delta_1 in [0,log d] as chaos diagnostic; KI dual-unitary achieves maximum Delta_1=log d; Section 12 added to Output.tex (2026-06-01)
- [x] CNT conjecture: disproved h_CNT >= Delta_1/|t|; proved Fekete bound h_CNT >= Delta_1/2; new finding: AFL orbit entropy under time evolution is SUBADDITIVE (not superadditive), unlike the shift; equality Delta_1 = h_CNT in linear-growth regime proved; Section 13 added (2026-06-01)
- [x] CNT summary: summary.md written; 4 new results, 5 open problems documented (2026-06-01)
- [ ] CNT open problem: prove rigorously that S(rho[Z^(n)]) is subadditive for the AFL orbit under time evolution (not superadditive), and show when equality S_2=2*S_1 holds; connect to operator entanglement and the dual-unitary property
