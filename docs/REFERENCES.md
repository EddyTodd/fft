# Primary references

This bibliography emphasizes original algorithm papers and specifications relevant to the permanent FFT mechanism catalog.

1. I. J. Good, **“The Interaction Algorithm and Practical Fourier Analysis,”** *Journal of the Royal Statistical Society, Series B* 20(2), 361–372, 1958. DOI: `10.1111/j.2517-6161.1958.tb00300.x`.
2. J. W. Cooley and J. W. Tukey, **“An Algorithm for the Machine Calculation of Complex Fourier Series,”** *Mathematics of Computation* 19(90), 297–301, 1965. DOI: `10.1090/S0025-5718-1965-0178586-1`.
3. C. M. Rader, **“Discrete Fourier Transforms When the Number of Data Samples Is Prime,”** *Proceedings of the IEEE* 56(6), 1107–1108, 1968. DOI: `10.1109/PROC.1968.6477`.
4. L. I. Bluestein, **“A Linear Filtering Approach to the Computation of Discrete Fourier Transform,”** *IEEE Transactions on Audio and Electroacoustics* 18(4), 451–455, 1970.
5. P. Duhamel and H. Hollmann, **“‘Split Radix’ FFT Algorithm,”** *Electronics Letters* 20(1), 14–16, 1984. DOI: `10.1049/el:19840012`.
6. S. G. Johnson and M. Frigo, **“A Modified Split-Radix FFT With Fewer Arithmetic Operations,”** *IEEE Transactions on Signal Processing* 55(1), 111–119, 2007.
7. M. Frigo and S. G. Johnson, **“FFTW: An Adaptive Software Architecture for the FFT,”** *ICASSP*, 1998, vol. 3, 1381–1384.
8. M. Frigo and S. G. Johnson, **“The Design and Implementation of FFTW3,”** *Proceedings of the IEEE* 93(2), 216–231, 2005.
9. S. Haynal and H. Haynal, **“Generating and Searching Families of FFT Algorithms,”** 2011, arXiv:`1103.5740`.
10. Intel, **Intel Intrinsics Guide**, primary ISA/intrinsic reference for the explicit x86 AVX2/AVX-512/FMA extension.

## Relevance to v1

- Good's interaction/Fourier work is the historical basis for the coprime-factor/Good-Thomas family represented by `GoodThomasPlan`.
- Cooley-Tukey supplies radix-2, radix-4, and mixed-radix decomposition structure.
- Duhamel-Hollmann and Johnson-Frigo distinguish classical and scaled modified split-radix mechanisms.
- Rader and Bluestein provide complementary prime/arbitrary-length reductions.
- FFTW remains useful historical context for reusable planning and codelet systems, but it is not a dependency of the v1 core.
