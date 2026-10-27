# Exact triangular-lattice finite lower-bound audit

All comparisons are exact rational comparisons. Ratios below are decimal displays only.
The all-n bound follows from the centering and homogeneity proof in the script; this finite table does not prove a universal comparison between lattice families.

| n | triangular / square-grid | triangular / prior envelope | envelope winner(s) |
|---:|---:|---:|:---|
| 2 | 1 | 1 | gaussian_average, regular_polygon, gaussian_integer_disk, triangular_integer_disk |
| 3 | 1.1851851852 | 1 | regular_polygon, triangular_integer_disk |
| 4 | 1.2671191692 | 0.75 | regular_polygon |
| 5 | 0.71428787541 | 0.54495840104 | regular_polygon |
| 6 | 1.1005475696 | 0.68961044718 | centered_polygon |
| 7 | 1.6162281065 | 1 | centered_polygon, triangular_integer_disk |
| 8 | 0.82958178248 | 0.38561041485 | centered_polygon |
| 9 | 0.67725412967 | 0.48449172288 | centered_polygon |
| 10 | 1.2969107454 | 0.49319628946 | old_explicit_center_improvement |
| 11 | 0.76300083297 | 0.44145378869 | old_explicit_center_improvement |
| 12 | 1.2078633172 | 0.57984832258 | old_explicit_center_improvement |
| 13 | 1.9087870623 | 1.9087870623 | triangular_integer_disk |
| 14 | 2.1416614590 | 2.1416614590 | triangular_integer_disk |
| 15 | 2.2799751846 | 2.2799751846 | triangular_integer_disk |
| 16 | 2.1491805479 | 2.1491805479 | triangular_integer_disk |
| 17 | 1.5389356275 | 1.5389356275 | triangular_integer_disk |
| 18 | 1.9609338541 | 1.9609338541 | triangular_integer_disk |
| 19 | 4.1337823858 | 4.1337823858 | triangular_integer_disk |
| 20 | 0.60662568242 | 0.60662568242 | gaussian_integer_disk |
| 21 | 0.069509992118 | 0.069509992118 | gaussian_integer_disk |
| 22 | 0.38751575458 | 0.38751575458 | gaussian_integer_disk |
| 23 | 1.4222155505 | 1.4222155505 | triangular_integer_disk |
| 24 | 5.3389828355 | 5.3389828355 | triangular_integer_disk |
| 25 | 5.7080848943 | 5.7080848943 | triangular_integer_disk |
| 26 | 3.9447396197 | 3.9447396197 | triangular_integer_disk |
| 27 | 1.3492747577 | 1.3492747577 | triangular_integer_disk |
| 28 | 0.64137612063 | 0.64137612063 | gaussian_integer_disk |
| 29 | 0.28731975942 | 0.28731975942 | gaussian_integer_disk |
| 30 | 2.8766835422 | 2.8766835422 | triangular_integer_disk |

Within this finite range, the triangular rule strictly improves the prior envelope at n = 13, 14, 15, 16, 17, 18, 19, 23, 24, 25, 26, 27, 30.

Special exact checks: T_2=4; T_3=27 with points (0,0),(-1,0),(-1,1); T_7=7^21/6^15. The square first shell gives its centered-square value at n=5. No claim of sharpness is made for n>=8.

The prior envelope includes the Gaussian average, unit regular polygon, centered regular polygon, Gaussian-integer disk rule, and the already-proved explicit center improvement for n>=10. Hermite is recorded in JSON but is dominated by the regular polygon for every n>=3.

The old n>=11 log improvement is converted to the rational lower bound C_n(1+x_m); it is a direct corollary of existing work, not a new perturbation theorem.

No assertion about the exact linear asymptotic coefficient or global crystallization is used.

Reproduce: python -B -S -O triangular_lattice.py --max-n 30
