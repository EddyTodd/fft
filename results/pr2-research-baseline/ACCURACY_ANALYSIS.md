# Accuracy analysis

Ranking uses the worst forward L2 error across the deterministic signal families. Lower is better. This is descriptive of the checked-in long-double reference experiment, not a proof of a global stability ordering.

| N | Rank | Algorithm | Median forward L2 | Worst forward L2 | Worst backward L2 | Worst forward Linf |
|---:|---:|---|---:|---:|---:|---:|
| 64 | 1 | split-radix | 1.569e-16 | 2.160e-16 | 2.140e-16 | 3.197e-16 |
| 64 | 2 | radix4 | 1.615e-16 | 2.295e-16 | 2.287e-16 | 3.197e-16 |
| 64 | 3 | stockham-radix2 | 2.406e-16 | 2.792e-16 | 2.783e-16 | 4.814e-16 |
| 64 | 4 | mixed-radix | 5.179e-16 | 5.422e-16 | 5.422e-16 | 1.230e-15 |
| 64 | 5 | auto | 5.842e-16 | 7.500e-16 | 7.499e-16 | 1.341e-15 |
| 64 | 6 | radix2-iterative | 5.842e-16 | 7.500e-16 | 7.499e-16 | 1.341e-15 |
| 64 | 7 | radix2-recursive | 5.842e-16 | 7.500e-16 | 7.499e-16 | 1.341e-15 |
| 64 | 8 | bluestein | 1.309e-15 | 1.465e-15 | 1.465e-15 | 3.318e-15 |
| 64 | 9 | dft | 1.022e-14 | 1.101e-14 | 1.101e-14 | 2.328e-14 |
| 127 | 1 | auto | 4.310e-15 | 4.777e-15 | 4.776e-15 | 6.429e-15 |
| 127 | 2 | rader | 4.310e-15 | 4.777e-15 | 4.776e-15 | 6.429e-15 |
| 127 | 3 | bluestein | 4.289e-15 | 4.890e-15 | 4.890e-15 | 6.829e-15 |
| 127 | 4 | dft | 1.991e-14 | 2.624e-14 | 2.624e-14 | 3.657e-14 |
| 127 | 5 | mixed-radix | 1.991e-14 | 2.624e-14 | 2.624e-14 | 3.657e-14 |
| 256 | 1 | split-radix | 2.534e-16 | 2.970e-16 | 3.028e-16 | 3.987e-16 |
| 256 | 2 | radix4 | 2.714e-16 | 3.089e-16 | 3.118e-16 | 4.001e-16 |
| 256 | 3 | stockham-radix2 | 3.346e-16 | 3.520e-16 | 3.557e-16 | 7.008e-16 |
| 256 | 4 | mixed-radix | 6.453e-16 | 7.011e-16 | 6.980e-16 | 1.573e-15 |
| 256 | 5 | auto | 1.845e-15 | 2.556e-15 | 2.555e-15 | 4.562e-15 |
| 256 | 6 | radix2-iterative | 1.845e-15 | 2.556e-15 | 2.555e-15 | 4.562e-15 |
| 256 | 7 | radix2-recursive | 1.845e-15 | 2.556e-15 | 2.555e-15 | 4.562e-15 |
| 256 | 8 | bluestein | 5.771e-15 | 8.157e-15 | 8.158e-15 | 2.163e-14 |
| 256 | 9 | dft | 3.062e-14 | 6.236e-14 | 6.236e-14 | 1.115e-13 |
| 509 | 1 | bluestein | 2.372e-14 | 2.491e-14 | 2.491e-14 | 5.766e-14 |
| 509 | 2 | auto | 2.273e-14 | 2.495e-14 | 2.496e-14 | 5.205e-14 |
| 509 | 3 | rader | 2.273e-14 | 2.495e-14 | 2.496e-14 | 5.205e-14 |
| 509 | 4 | dft | 9.426e-14 | 9.691e-14 | 9.690e-14 | 2.399e-13 |
| 509 | 5 | mixed-radix | 9.426e-14 | 9.691e-14 | 9.690e-14 | 2.399e-13 |
