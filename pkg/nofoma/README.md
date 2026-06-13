[![CI](https://github.com/gap-packages/nofoma/actions/workflows/CI.yml/badge.svg)](https://github.com/gap-packages/nofoma/actions/workflows/CI.yml)
[![Code Coverage](https://codecov.io/gh/gap-packages/nofoma/graph/badge.svg)](https://codecov.io/gh/gap-packages/nofoma)

# The GAP package nofoma

This package computes the Frobenius normal form and
the Jordan-Chevalley decomposition of a (square) matrix over any field
that is available in GAP, along with the Jordan normal form and primary decompositions
of matrices over finite fields.

To install the package, just unpack the tar file inside the pkg directory
of your GAP installation.

Then you can load the package into GAP by typing `LoadPackge("nofoma");`

## References

The algorithms in this package are based on, and a combination of:

- K. Bongartz, A direct approach to the rational normal form, preprint available at arXiv:1410.1683.
- M. Neunhoeffer and C. E. Praeger, Computing minimal polynomials of matrices, LMS J. Comput. Math. 11 (2008), 252-279;
- M. Geck, On Jacob's construction of the rational canonical form, Electron. J. Linear Algebra 36 (2020), 177-182.
- M. Geck, On the Jordan-Chevalley decomposition of a matrix, preprint.

## Feedback and support

If you have any bug reports, feature requests, or suggestions, then please
tell us via the
[issue tracker on GitHub](https://github.com/gap-packages/nofoma/issues).

## License

The nofoma package is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as published
by the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.
