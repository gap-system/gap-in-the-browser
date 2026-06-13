# nofoma - Normal forms of matrices
# Copyright (C) 2019-2022  Meinolf Geck <meinolf.geck@mathematik.uni-stuttgart.de>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, see <https://www.gnu.org/licenses/>.

#! @Chapter The nofoma package
#! @ChapterLabel The nofoma package
#! Let <M>K</M> be a field and <M>A</M> be an <M>n\times n</M>-matrix over 
#! <M>K</M>. This package provides functions for computing both the Frobenius normal form 
#! and the Jordan normal form of <M>A</M>. 
#! Furthermore, it also includes a functions for
#! the computation of a primary decomposition and the Jordan-Chevalley decomposition of 
#! <M>A</M>. 
#! 
#! @Section Installation of the &nofoma; package
#!
#! To install this package first unpack it inside some &GAP; root directory
#! in the subdirectory `pkg` (see <Ref Sect="Installing a GAP Package" BookName="ref"/>).
#! Then &nofoma; can already be loaded and used
#! (just type `LoadPackage("nofoma");`).
#!


DeclareInfoClass("Infonofoma");

#! @Chapter Normal forms of matrices
#! @Section The Frobenius normal form
#! Given a field <M>K</M> and an <M>n\times n</M>-matrix <M>A</M> 
#! over <M>K</M>, the <E>Frobenius normal form</E> of <M>A</M> is a block diagonal
#! matrix, where the diagonal blocks are companion matrices 
#! corresponding to the invariant factors of <M>A</M>. It reflects the minimal
#! decomposition of the vector space <M>K^n</M> into cyclic subspaces 
#! under the action of <M>A</M>.
#! The Frobenius normal form is also called the rational canonical form.

#! @Arguments A
#! @Description
#!  Returns the invariant factors of a matrix <A>A</A>
#!  and an invertible matrix <M>P</M> such that <M>PAP^{-1}</M> is the 
#!  Frobenius normal form of <A>A</A>. The algorithm first computes a maximal 
#!  vector and an <A>A</A>-invariant complement following Jacob's construction
#!  (as described in matrix language in <Cite Key ="Gec20"/>); then the 
#!  algorithm continues recursively. It works for matrices over any field 
#!  that is available in &GAP;. The output is a triple with
#!  * 1st component  = list of invariant factors; 
#!  * 2nd component = base change matrix <M>P</M>; and 
#!  * 3rd component = indices where the various blocks in the normal form 
#!       begin.
#! 
#! @BeginExampleSession
#! gap> A:=[ [  2,  2,  0,  1,  0,  2,  1 ],
#! >         [  0,  4,  0,  0,  0,  1,  0 ],
#! >         [  0,  1,  1,  0,  0,  1,  1 ],
#! >         [  0, -1,  0,  1,  0, -1,  0 ],
#! >         [  0, -7,  0,  0,  1, -5,  0 ],
#! >         [  0, -2,  0,  0,  0,  1,  0 ],
#! >         [  0, -1,  0,  0,  0, -1,  1 ] ];;
#! gap> f:=FrobeniusNormalForm(A);
#! [ [ x_1^4-7*x_1^3+17*x_1^2-17*x_1+6, x_1^2-3*x_1+2, x_1-1 ], 
#!   [ [    1,   -2,    1,    1,    0,    0,    1 ],
#!     [    2,   -7,    1,    2,    0,   -1,    3 ],
#!     [    4,  -26,    1,    4,    0,   -8,    6 ],
#!     [    8,  -89,    1,    8,    0,  -35,   11 ],
#!     [ -1/2,   -2,    0,  1/2,    0,   -2, -3/2 ],
#!     [   -1,   -4,    0,    0,    0,   -4,   -2 ],
#!     [    0,  9/4,    0,   -3,    1,  5/4,  1/4 ] ],
#!   [ 1, 5, 7 ]  ]                 
#! gap> PrintArray(f[2]*A*f[2]^-1);
#! [ [   0,   1,   0,   0,   0,   0,   0 ], 
#!   [   0,   0,   1,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   1,   0,   0,   0 ],
#!   [  -6,  17, -17,   7,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   1,   0 ],
#!   [   0,   0,   0,   0,  -2,   3,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   1 ] ]
#! @EndExampleSession
#! Note that the Frobenius normal form is unique up to the choice of the companion matrices
#! and the permutation of the blocks corresponding to the invariant factors. 
#! So while this function is significantly more efficient than the existing 'RationalCanonicalFormTransform',
#! the two functions yield slightly different results. 
#! In 'RationalCanonicalFormTransform', the companion matrices are consistent with the output of 'CompanionMat'. 
#! However, given an <M>n\times n</M> cyclic matrix <M>A</M>, along with a corresponding cyclic vector <M>v</M>, 
#! one can compute a change of basis matrix from A to a companion matrix of its minimal polynomial by computing 
#! <M>v</M> multiplied with powers of <M>A</M> (i.e., <M>v</M>, <M>vA</M>, ...., <M>vA^{n-1}</M>). This approach
#! follows GAP’s convention of right multiplication and yields a companion matrix in the form used by <Ref Func="FrobeniusNormalForm"/>.
#! Furthermore 'RationalCanonicalFormTransform' sorts the invariant factors in ascending order, while <Ref Func="FrobeniusNormalForm"/>
#! sorts them in descending order. 
#! @BeginExampleSession
#! gap> A := [ [ 0*Z(5), Z(5)^3, 0*Z(5), Z(5)^0, Z(5)^3 ], 
#! >   [ Z(5)^0, 0*Z(5), Z(5)^2, Z(5)^2, Z(5)^2 ], 
#! >   [ Z(5), Z(5), 0*Z(5), Z(5)^0, Z(5)^3 ], 
#! >   [ 0*Z(5), Z(5), Z(5)^2, Z(5)^2, Z(5) ], 
#! >   [ Z(5)^3, Z(5)^3, Z(5)^0, Z(5)^3, Z(5)^0 ] ];;
#! gap> T:=RationalCanonicalFormTransform(A);;
#! gap> S:=TransposedMat(FrobeniusNormalForm(TransposedMat(A))[2]);;
#! gap> Display(A^T);
#!  . . . . 3
#!  1 . . . 2
#!  . 1 . . 3
#!  . . 1 . 4
#!  . . . 1 .
#! gap> Display(A^S);
#!  . . . . 3
#!  1 . . . 2
#!  . 1 . . 3
#!  . . 1 . 4
#!  . . . 1 .
#! @EndExampleSession
#! 
#! Additionally, <C>RationalCanonicalFormTransform</C> sorts
#! the invariant factors in ascending order, whereas the 
#! <Ref Func="FrobeniusNormalForm"/> sorts them in 
#! descending order. Consequently, the outputs of the two functions 
#! agree up to a permutation of blocks and transposition.
#! To get a drop in replacement for 'RationalCanonicalFormTransform', see <Ref Func="FrobeniusNormalFormLikeRCFT"/>.
#! @BeginExampleSession
#! gap> aa:=[[  0, -8, 12, 40,-36,  4,  0, 59, 15, -9],
#! >         [ -2, -2, -2,  6,-11,  1, -1, 10,  1,  0],
#! >         [  1,  5,  0, -6, 12, -2,  0,-12, -4,  2],
#! >         [  0,  0,  0,  2,  0,  0,  0,  7,  0,  0],
#! >         [  0,  2, -3, -7,  8, -1,  0, -7, -3,  2],
#! >         [ -5, -4, -6, 18,-30,  2, -2, 35,  5, -1],
#! >         [ -1, -6,  6, 20,-28,  3,  0, 24, 10, -6],
#! >         [  0,  0,  0, -1,  0,  0,  0, -3,  0,  0],
#! >         [  0,  0, -1, -2, -2,  0, -1, -7,  0,  0],
#! >         [  0, -8,  9, 21,-36,  4, -2, 12, 12, -8]];;
#! gap> t:=RationalCanonicalFormTransform(aa);;
#! gap> Display(aa^t);
#! [ [   0,   0,   0,   1,   0,   0,   0,   0,   0,   0 ],
#!   [   1,   0,   0,   0,   0,   0,   0,   0,   0,   0 ],
#!   [   0,   1,   0,   0,   0,   0,   0,   0,   0,   0 ],
#!   [   0,   0,   1,   0,   0,   0,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   0,   0,   0,   1 ],
#!   [   0,   0,   0,   0,   1,   0,   0,   0,   0,   1 ],
#!   [   0,   0,   0,   0,   0,   1,   0,   0,   0,   1 ],
#!   [   0,   0,   0,   0,   0,   0,   1,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   0,   1,   0,  -1 ],
#!   [   0,   0,   0,   0,   0,   0,   0,   0,   1,  -1 ] ]
#! gap> res:=FrobeniusNormalForm(aa);;
#! gap> Display(aa^(res[2]^-1));
#! [ [   0,   1,   0,   0,   0,   0,   0,   0,   0,   0 ],
#!   [   0,   0,   1,   0,   0,   0,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   1,   0,   0,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   1,   0,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   1,   0,   0,   0,   0 ],
#!   [   1,   1,   1,   0,  -1,  -1,   0,   0,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   0,   1,   0,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   0,   0,   1,   0 ],
#!   [   0,   0,   0,   0,   0,   0,   0,   0,   0,   1 ],
#!   [   0,   0,   0,   0,   0,   0,   1,   0,   0,   0 ] ]
#! @EndExampleSession
DeclareGlobalFunction("FrobeniusNormalForm");

#! @Arguments A
#! @Description
#! This function returns the same result as <Ref Func="FrobeniusNormalForm"/>, except that the invariant factors 
#! are sorted in descending order and the companion matrices on the diagonal are transposed. Furthermore, if <M>P</M> is the 
#! computed base change matrix, the Frobenius normal form is obtained by <M>P^{-1}AP</M> (instead of <M>PAP^{-1}</M>).
#! This means that <M>P</M> can be used as a direct drop-in replacement for RationalCanonicalFormTransform. 
#!
#! Note that this function works by calling <Ref Func="FrobeniusNormalForm"/> and then modifying the computed 
#! transformation matrix and is thus potentially less efficient than using
#! the original function. 
#! @BeginExampleSession
#! gap> A := [ [ 0*Z(5), Z(5)^2, Z(5)^2, 0*Z(5), Z(5)^3, Z(5)^3, 0*Z(5), Z(5)^3, 
#! >       0*Z(5), Z(5) ], 
#! >   [ Z(5)^0, Z(5)^0, Z(5)^0, Z(5), Z(5)^3, Z(5), Z(5)^3, 0*Z(5), Z(5), 
#! >       0*Z(5) ], 
#! >   [ Z(5), 0*Z(5), 0*Z(5), Z(5)^3, Z(5)^2, Z(5)^0, 0*Z(5), Z(5)^0, 
#! >       Z(5), Z(5)^2 ], 
#! >   [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^2, 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^2, 
#! >       0*Z(5), 0*Z(5) ], 
#! >   [ 0*Z(5), Z(5)^2, Z(5)^2, Z(5)^0, Z(5)^0, Z(5)^3, 0*Z(5), Z(5)^0, 
#! >       Z(5)^2, Z(5)^2 ], 
#! >   [ 0*Z(5), Z(5), Z(5)^3, Z(5)^0, 0*Z(5), Z(5)^2, Z(5)^0, 0*Z(5), 
#! >       0*Z(5), Z(5)^3 ], 
#! >   [ Z(5)^3, Z(5)^3, Z(5), 0*Z(5), Z(5)^2, Z(5)^0, 0*Z(5), Z(5)^3, 
#!>       0*Z(5), Z(5)^3 ], 
#! >   [ 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^3, 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^2, 
#! >       0*Z(5), 0*Z(5) ], 
#! >   [ 0*Z(5), 0*Z(5), Z(5)^3, Z(5)^0, Z(5)^0, 0*Z(5), Z(5)^3, Z(5)^0, 
#! >       0*Z(5), 0*Z(5) ], 
#! >   [ 0*Z(5), Z(5)^2, Z(5)^3, Z(5), Z(5)^3, Z(5)^3, Z(5)^0, Z(5)^2, 
#! >       Z(5)^2, Z(5)^2 ] ];;
#! gap> frob := FrobeniusNormalFormLikeRCFT(A)[2];;
#! gap> rat := RationalCanonicalFormTransform(A);;
#! gap> Display(A^rat);
#!  . . . 1 . . . . . .
#!  1 . . . . . . . . .
#!  . 1 . . . . . . . .
#!  . . 1 . . . . . . .
#!  . . . . . . . . . 4
#!  . . . . 1 . . . . 2
#!  . . . . . 1 . . . 1
#!  . . . . . . 1 . . .
#!  . . . . . . . 1 . 1
#!  . . . . . . . . 1 3
#! gap> Display(A^frob);
#!  . . . 1 . . . . . .
#!  1 . . . . . . . . .
#!  . 1 . . . . . . . .
#!  . . 1 . . . . . . .
#!  . . . . . . . . . 4
#!  . . . . 1 . . . . 2
#!  . . . . . 1 . . . 1
#!  . . . . . . 1 . . .
#!  . . . . . . . 1 . 1
#!  . . . . . . . . 1 3
#! @EndExampleSession
DeclareGlobalFunction("FrobeniusNormalFormLikeRCFT");

#! @Arguments A
#! @Description
#!  Returns the invariant factors of the matrix <A>A</A>,
#!  i.e.,  the minimal polynomials of the  diagonal blocks in the Frobenius
#!  normal form  of <A>A</A>. Thus, 'InvariantFactorsMat' also specifies the
#!  rational canonical form of <A>A</A>, but without computing the base change.
#!
#! @BeginExampleSession
#! gap> A := [ [ 2,  2, 0, 1, 0,  2, 1 ],
#! >           [ 0,  4, 0, 0, 0,  1, 0 ],
#! >           [ 0,  1, 1, 0, 0,  1, 1 ],
#! >           [ 0, -1, 0, 1, 0, -1, 0 ],
#! >           [ 0, -7, 0, 0, 1, -5, 0 ],
#! >           [ 0, -2, 0, 0, 0,  1, 0 ],
#! >           [ 0, -1, 0, 0, 0, -1, 1 ] ];;
#! gap> InvariantFactorsMat(A);
#!   [ x_1^4-7*x_1^3+17*x_1^2-17*x_1+6, x_1^2-3*x_1+2, x_1-1 ]
#! @EndExampleSession
DeclareGlobalFunction("InvariantFactorsMat");

#! @Section The Jordan normal form

#! The Jordan normal form of a matrix <M>A</M> is a block diagonal matrix, 
#! where the diagonal blocks are Jordan blocks corresponding to the 
#! elementary divisors of <M>A</M>. It reflects the maximal decomposition of 
#! the vector space <M>K^n</M> into cyclic subspaces under the action of
#! <M>A</M>. For a more thorough definition of the Jordan normal form
#! and details about the algorithms used, see <Cite Key ="Bon26"/>. 

DeclareGlobalFunction("JordanNormalFormIrred");

#! @Arguments A
#! @Description
#!  Returns a list containing two entries. The first is a base change matrix
#!  <M>B</M> such that <M>B</M><A>A</A><M>B^{-1}</M> is the Jordan 
#!  normal form of <A>A</A>, i.e. a block diagonal matrix where the diagonal blocks
#!  are Jordan blocks corresponding to the elementary divisors of <A>A</A> in 
#!  descending order. The second entry is a list of the elementary divisors of <A>A</A>, also
#!  in descending order.
#!  The algorithm first computes a primary decomposition
#!  of <A>A</A> and then 
#!  computes a cyclic decomposition of the primary components. Finally it computes 
#!  Jordan block form for each of the cyclic components. The blocks are ordered in the 
#!  same order as in the list containing the elementary divisors. It works for matrices 
#!  over finite fields. 
#!
#!  Since all of the blocks on the resulting transformed matrix are cyclic, one 
#!  can retrieve their size by the degrees of the respective elementary divisor. 
#! 
#! @BeginExampleSession
#! gap> A := [ [ 0*Z(5), 0*Z(5), Z(5)^3, Z(5)^3, Z(5)^3, Z(5)^0 ], 
#! >    [ 0*Z(5), Z(5)^2, Z(5)^2, Z(5)^0, Z(5)^3, Z(5)^3 ], 
#! >    [ Z(5)^0, Z(5)^0, Z(5)^3, Z(5)^2, Z(5)^0, Z(5) ], 
#! >    [ 0*Z(5), Z(5)^3, Z(5), Z(5), 0*Z(5), Z(5)^2 ], 
#! >    [ Z(5)^2, Z(5)^0, Z(5)^0, 0*Z(5), Z(5), Z(5) ], 
#! >    [ 0*Z(5), Z(5)^0, Z(5)^2, Z(5), Z(5), Z(5) ] ];;
#! gap> B := JordanNormalForm(A);;
#! gap> Display(A^Inverse(B[1]));
#! 3 . . . . .
#! . 1 . . . .
#! . . . 1 . .
#! . . 2 . . .
#! . . . . . 1
#! . . . . 3 4
#! @EndExampleSession
#! This function computes the Jordan normal form of <M>A</M> 
#! significantly faster if <M>A</M> is either cyclic or has irreducible 
#! minimal polynomial. 
#! @BeginExampleSession
#! gap> B:= [ [ Z(5), Z(5)^3, 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
#! > [ Z(5)^0, Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
#! > [ 0*Z(5), 0*Z(5), Z(5), Z(5)^3, 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
#! > [ 0*Z(5), 0*Z(5), Z(5)^0, Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
#! > [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5), Z(5)^3, 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
#! > [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0, Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5) ], 
#! > [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5), Z(5)^3, 0*Z(5), 0*Z(5) ], 
#! > [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0, Z(5), 0*Z(5), 0*Z(5) ], 
#! > [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5), Z(5)^3 ], 
#! > [ 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), 0*Z(5), Z(5)^0, Z(5) ] ];;
#! gap> Factors(MinimalPolynomial(B));
#! [ x_1^2+x_1+Z(5)^0 ]
#! gap> Display(B^Inverse(JordanNormalForm(B)[1]));
#! . 1 . . . . . . . .
#! 4 4 . . . . . . . .
#! . . . 1 . . . . . .
#! . . 4 4 . . . . . .
#! . . . . . 1 . . . .
#! . . . . 4 4 . . . .
#! . . . . . . . 1 . .
#! . . . . . . 4 4 . .
#! . . . . . . . . . 1
#! . . . . . . . . 4 4
#! @EndExampleSession
DeclareGlobalFunction("JordanNormalForm");

#! @Chapter Other functionality

#! @Section Matrix decompositions

#! @Section The Jordan-Chevalley decomposition
#! @Arguments A,f
#! @Description
#!  Returns the unique pair of matrices <M>D</M>, 
#!  <M>N</M> such that the matrix <A>A</A> is written as <M>A=D+N</M>, where 
#!  <M>N</M> is a nilpotent matrix and <M>D</M> is a matrix that is 
#!  diagonalisable (over some extension field of the default field of 
#!  <A>A</A>), such that <M>D.N=N.D</M>; the argument <A>f</A> is a 
#!  polynomial such that <M>f(A)=0</M> (e.g., the minimal polynomial of 
#!  <A>A</A>). This is called the Jordan-Chevalley decomposition of <A>A</A>; 
#!  the algorithm is based on <Cite Key ="Gec22"/>. Note that this 
#!  algorithm does not require the knowledge of the eigenvalues of <A>A</A>; 
#!  it works over any perfect field that is available in &GAP;.
#!
#! @BeginExampleSession
#! gap> A:=[ [  6, -2,  6,  1,  1 ],
#! >         [  1, -1,  2,  1, -2 ],
#! >         [ -2,  0, -1,  0, -1 ],
#! >         [ -1,  0, -2,  2, -1 ],
#! >         [ -4,  4, -6, -2,  3 ] ];;
#! gap> jc:=JordanChevalleyDecMat(A,MinimalPolynomial(A));
#! [ [ [  4,  0,  4, -1,  1 ], 
#!     [  1,  0,  1,  1, -1 ], 
#!     [ -1, -1,  0,  1, -1 ], 
#!     [  0,  0, -2,  3,  0 ], 
#!     [ -3,  2, -4, -1,  2 ] ], 
#!   [ [  2, -2,  2,  2,  0 ], 
#!     [  0, -1,  1,  0, -1 ], 
#!     [ -1,  1, -1, -1,  0 ], 
#!     [ -1,  0,  0, -1, -1 ], 
#!     [ -1,  2, -2, -1,  1 ] ] ]
#! gap> MinimalPolynomial(jc[1]);
#! x_1^3-5*x_1^2+9*x_1-5
#! gap> Factors(last);
#! [ x_1-1, x_1^2-4*x_1+5 ]  
#! gap> MinimalPolynomial(jc[2]);
#! x_1^2                     
#! @EndExampleSession
#!  If the input matrix is very large, then <Ref Func="JordanChevalleyDecMatF"/>
#!  may be more efficient; this function first computes the Frobenius normal 
#!  form of <A>A</A> and then applies <C>JordanChevalleyDecMat</C> to each diagonal 
#!  block. (The result will be the same as that of 
#!  'JordanChevalleyDecMat(<A>A</A>);)'
DeclareGlobalFunction("JordanChevalleyDecMat");

#! @Arguments A
#! @Description
#!  First computes the Frobenius normal form and
#!  then applies <Ref Func="JordanChevalleyDecMat"/> to each diagonal block.
DeclareGlobalFunction("JordanChevalleyDecMatF");

#! @Section The primary decomposition
#! @Arguments A
#! @Description
#!  Returns a list containing three elements. The first element is
#!  a base change matrix <M>B</M> such that <M>B</M><A>A</A><M>B^{-1}</M> is a
#!  primary form of the matrix <A>A</A>, i.e., a block diagonal matrix where the minimal polynomials
#!  of the the diagonal blocks are precisely the powers of irreducible factors
#!  of the minimal polynomial of <A>A</A>, in descending order. The second element is a list containing
#!  the collected irreducible factors of the minimal polynomial of <A>A</A>, in the same order. The
#!  last element is a list containing the the size of each block. 
#!  The exact algorithm used in this function is described in <Cite Key ="Bon26"/>
#! 
#! @BeginExampleSession
#! gap> A := [ [ Z(5)^2, 0*Z(5), Z(5)^2, Z(5)^3, Z(5) ], 
#! >    [ 0*Z(5), 0*Z(5), Z(5)^3, Z(5), Z(5)^0 ],  
#! >    [ Z(5), Z(5)^0, 0*Z(5), Z(5)^0, 0*Z(5) ],
#! >    [ Z(5)^0, Z(5)^0, Z(5)^0, 0*Z(5), Z(5)^3 ],
#! >    [ Z(5), 0*Z(5), Z(5)^3, 0*Z(5), Z(5)^3 ] ];;
#! gap> B := PrimaryDecomposition(A);;
#! gap> Display(B[3]);
#! [ 1, 4 ]
#! gap> Factors(MinimalPolynomial(A));
#! [ x_1-Z(5)^0, x_1^4-x_1^3+Z(5)^3*x_1+Z(5)^3 ]
#! gap> PrimA := A^Inverse(B[1]);;
#! gap> MinimalPolynomial(PrimA{[1..1]}{[1..1]});
#! x_1-Z(5)^0
#! gap> MinimalPolynomial(PrimA{[2..5]}{[2..5]});
#! x_1^4-x_1^3+Z(5)^3*x_1+Z(5)^3
#! @EndExampleSession
DeclareAttribute("PrimaryDecomposition", IsMatrixOrMatrixObj);

#! @Chapter Auxiliary functions
#! @Section Vectors and matrices and their associated polynomials

#! @Arguments a,b
#! @Description
#! Computes a divisor <M>a_1</M> of the polynomial <A>a</A> and a
#! divisor <M>b_1</M> of the polynomial <A>b</A> such that <M>a_1</M> and <M>b_1</M> are coprime
#! and the lcm of <A>a</A>, <A>b</A> is <M>a_1*b_1</M>. This is based on Lemma 5 in <Cite Key ="Bon14"/>.
#! (See also Lemma 4.3 in <Cite Key ="Gec20"/>).
#!
#! (Note that it does not use the prime factorisation of polynomials but
#! only gcd computations.)
#!
#! @BeginExampleSession
#! gap> x:=X(Rationals);;
#! gap> a:=x^2*(x-1)^3*(x-2)*(x-3);
#! x_1^7-8*x_1^6+24*x_1^5-34*x_1^4+23*x_1^3-6*x_1^2
#! gap> b:=x^2*(x-1)^2*(x-2)^4*(x-4);
#! x_1^9-14*x_1^8+81*x_1^7-252*x_1^6+456*x_1^5-480*x_1^4+272*x_1^3-64*x_1^2
#! gap> GcdCoprimeSplit(a,b);
#! [ x_1^5-4*x_1^4+5*x_1^3-2*x_1^2, x_1^4-6*x_1^3+12*x_1^2-10*x_1+3, 
#!   x_1^7-12*x_1^6+56*x_1^5-128*x_1^4+144*x_1^3-64*x_1^2 ]
#! @EndExampleSession
DeclareGlobalFunction("GcdCoprimeSplit");

#! @Arguments A,pol,v
#! @Description
#! Returns the row vector obtained  by multiplying the row vector <A>v</A>
#! with the matrix <A>pol</A>(<A>A</A>), where p is a polynomial. The actual
#! computation is more efficient than this naive approach.
#!
#! @BeginExampleSession
#! gap> A:=[ [ 0, 1, 0, 1 ],
#! >         [ 0, 0, 0, 0 ],
#! >         [ 0, 1, 0, 1 ],
#! >         [ 1, 1, 1, 1 ] ];;
#! gap> x:=X(Rationals);;
#! gap> pol:=x^6-6*x^5+12*x^4-10*x^3+3*x^2;;
#! gap> v:=[ 1, 1, 1, 1 ];;
#! gap> PolynomialToMatVec(A,pol,v);
#! [ 8, -16, 8, -16 ]
#! @EndExampleSession
DeclareGlobalFunction("PolynomialToMatVec");

DeclareGlobalFunction("PolynomialToMat");

#! @Arguments A,v1,v2,pol1,pol2
#! @Description
#!  Returns,  given  a matrix  <A>A</A>,  vectors <A>v1</A>,
#!  <A>v2</A> with minimal polynomials <A>pol1</A>, <A>pol2</A>,  a new pair [<M>v</M>,<M>pol</M>],
#!  where <M>v</M> has minimal polynomial <M>pol</M>, and <M>pol</M> is the least common
#!  multiple of <A>pol1</A> and <A>pol2</A>.
#!  This crucially relies on <Ref Func="GcdCoprimeSplit"/> to avoid  factorisation of
#!  polynomials.
DeclareGlobalFunction("LcmMaximalVectorMat");

#! @Arguments A,v
#! @Description
#!  Computes  the  smallest subspace containing the vector
#!  <A>v</A> that is invariant under the matrix <A>A</A>. The  output is a
#!  quadruple, with
#!  * 1st component = basis of that subspace in row echelon form;
#!  * 2nd component = matrix  with  rows <M><A>v</A>, <A>v</A>.<A>A</A>, <A>v</A>.<A>A</A>^2,
#!     ..., <A>v</A>.<A>A</A>^{{d-1}}</M> (where <M>d</M> is the degree of the local
#!     minimal polynomial of <A>v</A>);
#!  * 3rd component = the coefficients <M>a_0</M>, <M>a_1</M>, ...,
#!     <M>a_d</M> of the local minimal polynomial; and
#!  * 4th component = the positions of the pivots of the first component.
#!
#! @BeginExampleSession
#! gap> A:=[ [   5,   2,  -4,   2 ],
#! >         [  -1,   0,   2,  -1 ],
#! >         [  -1,  -1,   3,  -1 ],
#! >         [ -13,  -7,  14,  -6 ] ];;
#! gap> SpinMatVector(A,[1,0,0,0]);
#! [ [ [ 1, 0, 0, 0 ], [ 0, 1, -2, 1 ] ],
#!   [ [ 1, 0, 0, 0 ], [ 5, 2, -4, 2 ] ],
#!   [ -1, 0, 1 ],
#!   [ 1, 2 ] ]
#! gap> SpinMatVector(A,[0,1,0,0]);
#! [ [ [ 0, 1, 0, 0 ], [ 1, 0, -2, 1 ], [ 0, 0, 1, -1/2 ] ],
#!   [ [ 0, 1, 0, 0 ], [ -1, 0, 2, -1 ], [ 6, 3, -4, 2 ] ],
#!   [ 1, -1, -1, 1 ],
#!   [ 2, 1, 3 ] ]
#! gap> SpinMatVector(A,[1,1,0,0]);
#! [ [ [ 1, 1, 0, 0 ], [ 0, 1, 1, -1/2 ] ],
#!   [ [ 1, 1, 0, 0 ], [ 4, 2, -2, 1 ] ],
#!   [ 1, -2, 1 ],
#!   [ 1, 2 ] ]
#! @EndExampleSession
DeclareGlobalFunction("SpinMatVector");

#! @Arguments A
#! @Description
#!  Repeatedly computes the smallest invariant subspaces containing different vectors
#!  to compute a chain of cyclic subspaces. The output is a triple
#!  <C>[B,C,svec]</C>  where  <M>C</M> is such that  <M>C<A>A</A>C^{-1}</M>  has a block
#!  triangular shape with companion matrices along the diagonal), <M>B</M> is the
#!  row echelon form of C and svec is the list of indices where the blocks
#!  begin.
#!
#! @BeginExampleSession
#! gap> A:=[ [ 0, 1, 0, 1 ],
#! >         [ 0, 0, 1, 0 ],
#! >         [ 0, 1, 0, 1 ],
#! >         [ 1, 1, 1, 1 ] ];;
#! gap> sp:=CyclicChainMat(A);
#! [ [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 1 ], [ 0, 0, 1, 0 ], [ 0, 0, 0, 1 ] ],
#!   [ [ 1, 0, 0, 0 ], [ 0, 1, 0, 1 ], [ 1, 1, 2, 1 ], [ 0, 0, 0, 1 ] ],
#!   [ 1, 4, 5 ] ]
#! gap> PrintArray(sp[2]*A*sp[2]^-1);
#! [ [    0,    1,    0,    0 ],
#!   [    0,    0,    1,    0 ],
#!   [    0,    3,    1,    0 ],
#!   [  1/2,  1/2,  1/2,    0 ] ]
#! @EndExampleSession
DeclareGlobalFunction("CyclicChainMat");

#! @Arguments A
#! @Description
#!  Returns the minimal polynomial and a maximal vector
#!  of the matrix <A>A</A>, that is, a vector whose local minimal polynomial
#!  is that of <A>A</A>. This is done by repeatedly spinning up vectors until
#!  a maximal one is found. The exact algorithm is a combination of
#!  * the minimal polynomial algorithm by Neunhoeffer-Praeger; see <Cite Key ="Neu08"/>; and
#!  * the method described by Bongartz
#!      (see <Cite Key ="Bon14"/>) for computing
#!      maximal vectors.
#!
#!  See also the article by Geck at <Cite Key ="Gec20"/>.
#!
#! @BeginExampleSession
#! gap> A:=[ [  2,  2,  0,  1,  0,  2,  1 ],
#! >         [  0,  4,  0,  0,  0,  1,  0 ],
#! >         [  0,  1,  1,  0,  0,  1,  1 ],
#! >         [  0, -1,  0,  1,  0, -1,  0 ],
#! >         [  0, -7,  0,  0,  1, -5,  0 ],
#! >         [  0, -2,  0,  0,  0,  1,  0 ],
#! >         [  0, -1,  0,  0,  0, -1,  1 ] ];;
#! gap> MaximalVectorMat(A);
#! [ [ 1, -2, 1, 1, 0, 0, 1 ], x_1^4-7*x_1^3+17*x_1^2-17*x_1+6 ]
#! gap> v:=last[1];;
#! gap> SpinMatVector(A,v)[3];
#! [ 6, -17, 17, -7, 1 ]
#! @EndExampleSession
#! In the following example, <M>M_2</M> is the (challenging) test matrix
#! from the paper by Neunhoeffer-Praeger:
#!
#! @BeginLogSession
#! gap> LoadPackage("AtlasRep");; g:=AtlasGroup("B",1); M2:=g.1+g.2+g.1*g.2;
#! <matrix group of size 4154781481226426191177580544000000 with 2 generators>
#! <an immutable 4370x4370 matrix over GF2>
#! gap> SetInfoLevel(Infonofoma,1);
#! gap> MaximalVectorMat(M2);;time;
#! #I Degree of minimal polynomial is 2097
#! 6725
#! gap> MinimalPolynomial(M2);;time;
#! 13415
#! gap> LoadPackage("cvec");
#! gap> MinimalPolynomial(CMat(M2));;time;
#! 9721
#! @EndLogSession
DeclareGlobalFunction("MaximalVectorMat");

#! @Arguments T,d
#! @Description
#! Modifies an already given  complementary subspace
#! to the  complementary subspace defined by  Jacob;  concretely, this is
#! realized by assuming that  <A>T</A>  is a matrix in block triangular shape,
#! where the upper left diagonal block is a companion matrix (as returned
#! by <Ref Func="RatFormStep1"/>; the variable <A>d</A> gives the size of that block.
#! (If <A>T</A> gives a  maximal cyclic subspace,  then  Jacob's complement is
#! also  <A>T</A>-invariant;  but even if not,  it appears  to be  very useful
#! because it produces many zeroes.)
DeclareGlobalFunction("JacobMatComplement");

#! @Arguments A,v
#! @Description
#! Spins up a vector  <A>v</A> under a  matrix  <A>A</A>,  computes
#! a complementary subspace  (using  Jacob's construction),  and performs
#! the base change. The output is a quadruple  <C>[A1,P,pol,str]</C> where <M>A1</M> is
#! the new matrix, <M>P</M> is the base change,  <M>pol</M> is  the minimal polynomial
#! and <M>str</M> is either 'split' or 'not', according to whether the extension
#! is split or not.
#!
#! @BeginExampleSession
#! gap> v:=[ 1, 1, 1, 1 ];;
#! gap> A:=[ [ 0, 1, 0, 1 ],
#! >         [ 0, 0, 1, 0 ],
#! >         [ 0, 1, 0, 1 ],
#! >         [ 1, 1, 1, 1 ] ];;
#! gap> PrintArray(RatFormStep1(A,v)[1]);
#! [ [  0,  1,  0,  0 ],
#!   [  0,  0,  1,  0 ],
#!   [  0,  3,  1,  0 ],
#!   [  1,  0,  0,  0 ] ]
#! @EndExampleSession
DeclareGlobalFunction("RatFormStep1");

DeclareGlobalFunction("RatFormStep1J");
DeclareGlobalFunction("SquareFreePol");
