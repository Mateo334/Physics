# j=10 kicked top: fixed-Jy time-translation QFI after 50 complete kicks

[Numerically tested] Script `calculations/code/kicked_top_time_translation_qfi_j10.py`; run `OPENBLAS_NUM_THREADS=1; py -3.13 calculations/code/kicked_top_time_translation_qfi_j10.py`. Figure (visually inspected), 51-by-8 QFI CSV, compressed NPZ, and compact numerical JSON summary are `calculations/results/kicked_top_j10_Jy_time_translation_qfi_50kicks.{png,csv,npz,json}`. These results are reproducible and the script refuses to overwrite them.

Model: j=10, d=21, ħ=1; coherent pure initial state with θ=1.1, φ=0.7 made by `e^{-i φ Jz}e^{-i θ Jy}|j,j>` with the same Hermitian Jy used previously in `kicked_top_instantaneous_povm_pga.py`. All eight K=0.5,1,2,4,6,8,10,12 use this **same** state, and complete kicks evolve by `ψ_{n+1}=exp[-i K Jz²/(2j)]exp[-i π Jy/2]ψ_n`. Auxiliary translation generator H=Jy remains fixed, **not** the full periodic Hamiltonian; the computed QFI is `FQ(n,K)=4(⟨Jy²⟩-⟨Jy⟩²)` for n=0,...,50, **not** parameter-K QFI. Natural-number spin basis m=-10,...,10 and ℏ=1.

| K | FQ(0) | FQ(1) | FQ(10) | FQ(50) | max n<=50 |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 13.407459 | 20.835914 | 10.879534 | 19.272806 | 26.801501 |
| 1 | 13.407459 | 27.165972 | 33.516754 | 27.513072 | 49.193556 |
| 2 | 13.407459 | 27.103681 | 87.643502 | 93.922594 | 109.257611 |
| 4 | 13.407459 | 28.832398 | 163.223369 | 117.393393 | 224.188772 |
| 6 | 13.407459 | 102.576515 | 131.214336 | 159.777642 | 193.714655 |
| 8 | 13.407459 | 87.681532 | 124.607658 | 184.006428 | 198.464160 |
| 10 | 13.407459 | 117.551688 | 113.525479 | 105.498340 | 190.892247 |
| 12 | 13.407459 | 116.889490 | 165.303421 | 130.851120 | 194.959666 |

Checks: all eight FQ(0) match the coherent-state identity `2j[1-(sin θ sin φ)²]=13.40745939569193` within 1e-10; Floquet rotation unitary error below 1e-12; max state normalization error across all n,K ≤9.55e-15; all FQ within `[0,4j²]=[0,400]`, and moment imaginary errors ≤3.56e-15. No fitted growth or asymptotic chaos claim is inferred: the eight time series fluctuate and some exhibit revivals in finite dimension. The output CSV lists **every** computed numerical value, not just the selected kick numbers in this note.
