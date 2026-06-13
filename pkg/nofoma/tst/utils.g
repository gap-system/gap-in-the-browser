# Returns the Frobenius normal form of a matrix with invariant factors
# <A>facs</A>. It works by building the block diagonal matrix with
# diagonal blocks given by the companion matrices according to the invariant
# factors.
CreateNormalForm := function(k, facs)
  local A,l,r,i,B;
  r:=Length(facs);
  l:=[1];
  for i in [1..r] do
    l[i+1]:=l[i]+Degree(facs[i]);
  od;
  A:=ZeroMatrix(k, l[r+1]-1, l[r+1]-1);
  for i in [1..r] do
    B := nfmCompanionMat1(k,CoefficientsOfUnivariatePolynomial(facs[i]));
    CopySubMatrix( B, A, [1..NrRows(B)], [l[i]..l[i+1]-1], [1..NrCols(B)], [l[i]..l[i+1]-1] );
  od;
  return A;
end;


# F=output of FrobeniusNormalForm
CheckFrobForm := function(A,F)
  local P,i,nf;
  nf:=CreateNormalForm(BaseDomain(A), F[1]);
  P:=F[2];
  A := Matrix(A, P);
  nf := Matrix(nf, P);
  if P*A*P^(-1)<>nf then
    Error("base change not ok!");
  fi;
  for i in [1..Length(F[1])-1] do
    if EuclideanRemainder(F[1][i],F[1][i+1])<>0*F[1][i] then
      Error("divisibility not ok!");
    fi;
  od;
  return true;
end;

CheckJordanChev := function(mat,jc)
  local m;
  m:=MinimalPolynomial(jc[1]);
  return [Gcd(m,Derivative(m)),MinimalPolynomial(jc[2])];
end;

nfmmat1 := function(mat)
  local a,a1,b,i;
  a1:=TransposedMat(Concatenation(mat,mat));
  a:=[];
  for i in [1..Length(a1)-1] do
    Add(a,a1[i]);
  od;
  Add(a,0*a1[1]);
  b:=TransposedMat(Concatenation(TransposedMat(mat),TransposedMat(mat)));
  return Concatenation(a,b);
end;

#TODO: now that primdecomp has sorted blocks this can be tested more efficiently
#check primary decomp function with field F and dimension n
#checks if the submatrices have the correct minimal polynomials
nfmCheckPrimaryDecomposition := function(F, n)
  local A,Prim,B,dim,k,minpolfacs,sub;
  A := Matrix(F,RandomInvertibleMat(n,F));
  Prim := PrimaryDecomposition(A);
  B := A^Inverse(Prim[1]);
  minpolfacs := Factors(MinimalPolynomial(F,A));
  k := 1;
  for dim in Prim[3] do
    sub := ExtractSubMatrix(B,[k..k+dim-1],[k..k+dim-1]);
    k := k+dim;
    if not Factors(MinimalPolynomial(F,sub))[1] in minpolfacs then
      return false;
    fi;
  od;
  return true;
end;

nfmCheckPrimaryDecompositionNonCyclic := function(F, n)
  local A,Prim,B,dim,k,minpolfacs,sub;
  A := Matrix(F,nfmGenerateNonCyclicMatrix(F,n));
  Prim := PrimaryDecomposition(A);
  B := A^Inverse(Prim[1]);
  minpolfacs := Factors(MinimalPolynomial(F,A));
  k := 1;
  for dim in Prim[3] do
    sub := ExtractSubMatrix(B,[k..k+dim-1],[k..k+dim-1]);
    k := k+dim;
    if not Factors(MinimalPolynomial(F,sub))[1] in minpolfacs then
      return false;
    fi;
  od;
  return true;
end;

nfmMatrixFamily := function(A)
  if IsGF2MatrixRep(A) then
    return "gf2";
  elif Is8BitMatrixRep(A) then
    return "8bit";
  elif IsMatrixObj(A) then
    return "matobj";
  elif IsMatrix(A) then
    return "plain";
  fi;
  Error("unknown matrix family");
end;

nfmCheckInvariantFactorsForMatrix := function(A)
  return InvariantFactorsMat(A) = FrobeniusNormalForm(A)[1];
end;

nfmCheckPrimaryDecompositionForMatrix := function(A)
  local Prim,B,dim,k,minpolfacs,sub,AA;
  Prim := PrimaryDecomposition(A);
  AA := Matrix(A, Prim[1]);
  B := AA^Inverse(Prim[1]);
  minpolfacs := Factors(MinimalPolynomial(AA));
  k := 1;
  for dim in Prim[3] do
    sub := ExtractSubMatrix(B,[k..k+dim-1],[k..k+dim-1]);
    k := k+dim;
    if not Factors(MinimalPolynomial(sub))[1] in minpolfacs then
      return false;
    fi;
  od;
  return true;
end;

nfmJordanNormalformFamilies := function(A)
  local J;
  J := JordanNormalForm(A);
  return [nfmMatrixFamily(A), nfmMatrixFamily(J[1])];
end;
