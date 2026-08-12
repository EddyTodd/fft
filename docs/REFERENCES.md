# Primary references

This bibliography is intentionally biased toward original papers and project-authored methodology rather than tertiary summaries.

1. J. W. Cooley and J. W. Tukey, **“An Algorithm for the Machine Calculation of Complex Fourier Series,”** *Mathematics of Computation* 19(90), 297–301, 1965. DOI: `10.1090/S0025-5718-1965-0178586-1`.
2. C. M. Rader, **“Discrete Fourier Transforms When the Number of Data Samples Is Prime,”** *Proceedings of the IEEE* 56(6), 1107–1108, 1968. DOI: `10.1109/PROC.1968.6477`.
3. L. I. Bluestein, **“A Linear Filtering Approach to the Computation of Discrete Fourier Transform,”** *IEEE Transactions on Audio and Electroacoustics* 18(4), 451–455, 1970. Commonly known as Bluestein’s/chirp-z algorithm.
4. P. Duhamel and H. Hollmann, **“‘Split Radix’ FFT Algorithm,”** *Electronics Letters* 20(1), 14–16, 1984. DOI: `10.1049/el:19840012`.
5. M. Frigo and S. G. Johnson, **“FFTW: An Adaptive Software Architecture for the FFT,”** *ICASSP*, 1998, vol. 3, 1381–1384.
6. M. Frigo and S. G. Johnson, **“The Design and Implementation of FFTW3,”** *Proceedings of the IEEE* 93(2), 216–231, 2005.
7. S. G. Johnson and M. Frigo, **“A Modified Split-Radix FFT With Fewer Arithmetic Operations,”** *IEEE Transactions on Signal Processing* 55(1), 111–119, 2007.
8. S. Haynal and H. Haynal, **“Generating and Searching Families of FFT Algorithms,”** 2011, arXiv:`1103.5740`.
9. J. Alman and K. Rao, **“Faster Walsh-Hadamard and Discrete Fourier Transforms From Matrix Non-Rigidity,”** 2022, arXiv:`2211.06459`.
10. M. Frigo and S. G. Johnson, **benchFFT speed methodology**, FFTW project, `https://www.fftw.org/speed/method.html`.
11. M. Frigo and S. G. Johnson, **benchFFT accuracy methodology**, FFTW project, `https://www.fftw.org/accuracy/method.html`.

## Why these matter here

- Cooley–Tukey is the baseline decomposition family.
- Rader and Bluestein provide distinct reductions for prime/arbitrary lengths.
- split-radix is important for arithmetic-count research; modified split-radix demonstrates that lower arithmetic count remains an active theoretical question rather than a settled proxy for wall-clock speed.
- FFTW is the canonical example of architecture-adaptive FFT engineering and supplies a mature precedent for separating setup/planning, execution speed, data format, and numerical accuracy.
- recent arithmetic-complexity work is tracked so the repository does not incorrectly imply that textbook split-radix is the theoretical endpoint.
