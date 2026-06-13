gap> START_TEST("nofoma");
gap> ReadPackage("nofoma", "tst/utils.g");
true

## Some test matrices; see also
##  A. Steel,  A new algorithm for the computation of  canonical forms of
##               matrices over fields, J. Symb. Comp. 24 (1997), 409-432.
## 
gap> bev:=TransposedMat(
> [[2,0,0,0,0,0,0], [2,4,1,-1,-7,-2,-1], [0,0,1,0,0,0,0], [1,0,0,1,0,0,0],
> [0,0,0,0,1,0,0], [2,1,1,-1,-5,1,-1], [1,0,1,0,0,0,1]]);;
gap> steel:=
> [[-23,19,-9,-75,34,9,56,15,-34,-9], [-2,2,-1,-6,3,1,4,2,-3,0], 
> [4,-4,3,10,-5,-1,-6,-4,5,1], [-2,2,-1,-5,3,1,3,2,-3,0], 
> [0,0,0,0,2,0,0,0,0,0], [12,-12,6,33,-18,-4,-18,-12,18,0], 
> [-1,-3,0,2,1,0,1,1,2,1], [-26,22,-10,-83,36,10,61,18,-39,-10], 
> [-1,-3,0,1,1,0,2,1,2,0], [8,-12,4,27,-12,-4,-12,-7,15,0]];;
gap> low:=Z(19)^0 *[[16,0,0,0,0,0,0,0,0,0],[5,2,0,0,0,0,0,0,0,0],
> [9,6,7,0,0,0,0,0,0,0],[15,9,1,2,0,0,0,0,0,0],[7,14,8,2,9,0,0,0,0,0],
> [12,3,17,5,1,12,0,0,0,0],[7,11,0,4,6,9,10,0,0,0],[8,3,15,16,17,18,18,12,0,0],
> [6,3,7,12,1,12,11,14,10,0],[18,14,7,17,16,15,13,13,3,8]];;
gap> ddd:=DiagonalMat([1,2,-1,0,4,3,3,4,5,6,7,-1,5,4,0,0,3,2,1]);;

## Testing Frobenius Normalform
gap> CheckFrobForm(bev, FrobeniusNormalForm(bev));
true
gap> CheckFrobForm(steel, FrobeniusNormalForm(steel));
true
gap> CheckFrobForm(low, FrobeniusNormalForm(low));
true
gap> CheckFrobForm(ddd, FrobeniusNormalForm(ddd));
true
gap> CheckFrobForm(Z(4)^0*steel, FrobeniusNormalForm(Z(4)^0*steel));
true
gap> CheckFrobForm(Z(4)^0*bev, FrobeniusNormalForm(Z(4)^0*bev));
true
gap> CheckFrobForm(Z(4)^0*ddd, FrobeniusNormalForm(Z(4)^0*ddd));
true
gap> CheckFrobForm(nfmmat1(steel), FrobeniusNormalForm(nfmmat1(steel)));
true
gap> CheckFrobForm(Z(29)*nfmmat1(steel), FrobeniusNormalForm(Z(29)*nfmmat1(steel)));
true

#
gap> bigfield := RandomInvertibleMat(10,GF(5^5));; #check non-list matrices
gap> bigfield := Matrix(GF(5^5), bigfield);;
gap> CheckFrobForm(bigfield, FrobeniusNormalForm(bigfield));
true

#
gap> issue24 := Matrix(GF(5^5), RandomInvertibleMat(5,GF(5^5)));;
gap> CheckFrobForm(issue24, FrobeniusNormalForm(issue24));
true

#
gap> issue68 := [ [ Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11^2)^30, 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11^2)^29, 0*Z(11), 0*Z(11), 0*Z(11), Z(11^2)^119, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11^2)^27, 0*Z(11), 0*Z(11), 0*Z(11), Z(11^2)^30, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11^2)^86, 0*Z(11), 0*Z(11), 0*Z(11), Z(11^2)^29, 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11), 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5, 0*Z(11) ],
> [ 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), 0*Z(11), Z(11)^5 ] ];;
gap> CheckFrobForm(issue68, FrobeniusNormalForm(issue68));
true

## Testing Jordan Chevalley Decomposition
gap> CheckJordanChev(bev,JordanChevalleyDecMatF(bev));
[ 1, x_1^2 ]
gap> CheckJordanChev(steel,JordanChevalleyDecMatF(steel));
[ 1, x_1^2 ]
gap> CheckJordanChev(low,JordanChevalleyDecMatF(low));
[ Z(19)^0, x_1^2 ]
gap> CheckJordanChev(ddd,JordanChevalleyDecMatF(ddd));
[ 1, x_1 ]
gap> CheckJordanChev(Z(4)^0*steel,JordanChevalleyDecMatF(Z(4)^0*steel));
[ Z(2)^0, x_1^4 ]
gap> CheckJordanChev(Z(4)^0*bev,JordanChevalleyDecMatF(Z(4)^0*bev));
[ Z(2)^0, x_1^2 ]
gap> CheckJordanChev(Z(4)^0*ddd,JordanChevalleyDecMatF(Z(4)^0*ddd));
[ Z(2)^0, x_1 ]
gap> CheckJordanChev(nfmmat1(steel),JordanChevalleyDecMatF(nfmmat1(steel)));
[ 1, x_1 ]
gap> CheckJordanChev(Z(29)*nfmmat1(steel),JordanChevalleyDecMatF(Z(29)*nfmmat1(steel)));
[ Z(29)^0, x_1 ]

## Stop test
gap> STOP_TEST("nofoma");
