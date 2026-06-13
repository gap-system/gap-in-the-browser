gap> START_TEST("representations");
gap> ReadPackage("nofoma", "tst/utils.g");
true

## Plain list over integers, and a matrix object over rationals
gap> rat := [ [ 1, 1, 0 ], [ 0, 1, 0 ], [ 0, 0, 2 ] ];;
gap> ratm := Matrix(Rationals, rat);;

## Plain list over a small extension field, with some entries in the prime field
gap> smallplain := [ [ Z(3^2), Z(3)^0, 0*Z(3) ],
> [ 0*Z(3), Z(3^2)^3, Z(3)^0 ],
> [ 0*Z(3), 0*Z(3), Z(3^2) ] ];;

## Compressed matrices over GF(2) and over a small field of order > 2
gap> gf2 := [ [ Z(2)^0, Z(2)^0, 0*Z(2) ],
> [ 0*Z(2), Z(2)^0, Z(2)^0 ],
> [ 0*Z(2), 0*Z(2), Z(2)^0 ] ];;
gap> ConvertToMatrixRep(gf2);;
gap> gf9 := StructuralCopy(smallplain);;
gap> ConvertToMatrixRep(gf9);;

## Plain list and matrix object over a finite field of order > 257
gap> bigplain := [ [ Z(5^5)^7, Z(5)^0, 0*Z(5) ],
> [ 0*Z(5), Z(5^5)^11, Z(5)^0 ],
> [ 0*Z(5), 0*Z(5), Z(5^5)^13 ] ];;
gap> bigm := Matrix(GF(5^5), bigplain);;

## Invariant factors should work for all these matrix families
gap> nfmCheckInvariantFactorsForMatrix(rat);
true
gap> nfmCheckInvariantFactorsForMatrix(ratm);
true
gap> nfmCheckInvariantFactorsForMatrix(smallplain);
true
gap> nfmCheckInvariantFactorsForMatrix(gf2);
true
gap> nfmCheckInvariantFactorsForMatrix(gf9);
true
gap> nfmCheckInvariantFactorsForMatrix(bigplain);
true
gap> nfmCheckInvariantFactorsForMatrix(bigm);
true

## Jordan-Chevalley decomposition should work for all these matrix families
gap> CheckJordanChev(rat, JordanChevalleyDecMat(rat, MinimalPolynomial(rat)));
[ 1, x_1^2 ]
gap> CheckJordanChev(ratm, JordanChevalleyDecMat(ratm, MinimalPolynomial(ratm)));
[ 1, x_1^2 ]
gap> CheckJordanChev(smallplain, JordanChevalleyDecMat(smallplain, MinimalPolynomial(smallplain)));
[ Z(3)^0, x_1^2 ]
gap> CheckJordanChev(gf2, JordanChevalleyDecMat(gf2, MinimalPolynomial(gf2)));
[ Z(2)^0, x_1^3 ]
gap> CheckJordanChev(gf9, JordanChevalleyDecMat(gf9, MinimalPolynomial(gf9)));
[ Z(3)^0, x_1^2 ]
gap> CheckJordanChev(bigplain, JordanChevalleyDecMat(bigplain, MinimalPolynomial(bigplain)));
[ Z(5)^0, x_1 ]
gap> CheckJordanChev(bigm, JordanChevalleyDecMat(bigm, MinimalPolynomial(bigm)));
[ Z(5)^0, x_1 ]

## Primary decomposition should work for all these matrix families
gap> nfmCheckPrimaryDecompositionForMatrix(rat);
true
gap> nfmCheckPrimaryDecompositionForMatrix(ratm);
true
gap> nfmCheckPrimaryDecompositionForMatrix(smallplain);
true
gap> nfmCheckPrimaryDecompositionForMatrix(gf2);
true
gap> nfmCheckPrimaryDecompositionForMatrix(gf9);
true
gap> nfmCheckPrimaryDecompositionForMatrix(bigplain);
true
gap> nfmCheckPrimaryDecompositionForMatrix(bigm);
true

## Jordan normal form should preserve the input matrix family
gap> nfmJordanNormalformFamilies(smallplain);
[ "plain", "plain" ]
gap> nfmJordanNormalformFamilies(gf2);
[ "gf2", "gf2" ]
gap> nfmJordanNormalformFamilies(gf9);
[ "8bit", "8bit" ]
gap> nfmJordanNormalformFamilies(bigplain);
[ "plain", "plain" ]
gap> nfmJordanNormalformFamilies(bigm);
[ "matobj", "matobj" ]

#
gap> STOP_TEST("representations");
