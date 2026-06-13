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

BindGlobal("nfmPolCoeffs", function(coeffs)
  return UnivariatePolynomialByCoefficients(FamilyObj(coeffs[1]),coeffs,1);
end);

# Computes divisor a_1 of the polynomial a and a divisor
# b_1 of the polynomial b such that a_1 and b_1 are coprime
# and the lcm of a,b is a_1*b_1.
# For details about this function, see nofoma.gd.
InstallGlobalFunction(GcdCoprimeSplit,function(a,b)
  local d,tb,bb;
  d:=Gcd(a,b);
  if IsZero(Degree(d)) then
    return [d,a,b];
  fi;
  if Degree(b)<=Degree(a) then
    if Degree(d)=Degree(b) then
      return [b,a,b^0];
    else
      tb:=Quotient(b,d);
      bb:=Gcd(tb^Degree(d),d);
      return [d,Quotient(a,bb),tb*bb];
    fi;
  else
    return GcdCoprimeSplit(b,a){[1,3,2]};
  fi;
end);

# Compute pol(A)*v, but more efficiently: the naive approach would compute
# first the matrix pol(A) and thus matrix powers A, A^2, A^3, etc. which is
# very expensive. The code here avoids this by computing a  linear combination
# of v, A*v, A*(A*v), ...,
InstallGlobalFunction(PolynomialToMatVec,function(A,pol,v)
  local n,v1,i,coeffs;
  coeffs := CoefficientsOfUnivariatePolynomial(pol);
  n:=Length(coeffs);
  v1:=ShallowCopy(coeffs[n]*v);
  for i in [n-1,n-2..1] do
    v1:=v1*A;
    if not IsZero(coeffs[i]) then
      AddRowVector(v1,v,coeffs[i]);
    fi;
  od;
  return v1;
end);

#TODO: take polynomial as input
# Applies polynomial pol to A.
InstallGlobalFunction(PolynomialToMat,function(A,pol)
  local A1,idm,n,i,coeffs;
  idm:= OneMutable(A);
  coeffs := CoefficientsOfUnivariatePolynomial(pol);
  n:=Length(coeffs);
  A1:=coeffs[n]*idm;
  for i in [n-1,n-2..1] do
    A1:=A1*A;
    if not IsZero(coeffs[i]) then
      A1:=A1+coeffs[i]*idm;
    fi;
  od;
  return A1;
end);

# Takes matrix A, vectors v1,v2 with minimal polynomials pol1,pol2 as inputs
# Returns [v,pol] such that v has minimal polynomial pol, pol is lcm
# of pol1 and pol2.
# For details about this function, see nofoma.gd.
InstallGlobalFunction(LcmMaximalVectorMat,function(A,v1,v2,pol1,pol2)
  local d,v;
  if EuclideanRemainder(pol1,pol2)=0*pol1 then
    return [v1,pol1];
  elif EuclideanRemainder(pol2,pol1)=0*pol2 then
    return [v2,pol2];
  else
    d:=GcdCoprimeSplit(pol1,pol2);
    if Degree(d[1])=0 then
      return [v1+v2,pol1*pol2];
    else
      v:=PolynomialToMatVec(A,Quotient(pol1,d[2]),v1)+
             PolynomialToMatVec(A,Quotient(pol2,d[3]),v2);
      return [v,d[2]*d[3]];
    fi;
  fi;
end);

# SpinMatVector takes a matrix A and a vector v as input.
# Returns quadruple where first component contains a basis in row
# echelon form, the second the matrix with rows v, Av, A^2v, ..., A^(d-1)v,
# the third one the coefficients of the minimal polynomial of v (of degree d)
# and the last one the positions of the pivots of the first component.
# For details about this function, see nofoma.gd.

# SpinMatVector1 does not compute the minimal polynomial, only the subspace;
# here, the last four arguments can be empty lists.
BindGlobal("SpinMatVector1",function(A,v,orbit1,orbit,pivot)
  local i,j,nv,coeff;

  #v := CompatibleVector(v, A);  # needs <https://github.com/gap-system/gap/issues/6302>
  nv := CompatibleVector(A);
  for i in [1..NrRows(A)] do
    nv[i] := v[i];
  od;
  v := nv;

  i:=PositionNonZero(v);
  if i>Length(v) then
    Error("# zero vector");
  fi;
  nv:=ShallowCopy(v);
  while i <= Length(nv) do
    Add(pivot,i);
    if not IsOne(nv[i]) then
      MultVector(nv,Inverse(nv[i]));
    fi;
    Add(orbit1,nv);
    Add(orbit,v);

    v:=v*A;
    nv:=ShallowCopy(v);
    for j in [1..Length(orbit1)] do
      coeff:=-nv[pivot[j]];
      if not IsZero(coeff) then
        AddRowVector(nv,orbit1[j],coeff);
      fi;
    od;
    i:=PositionNonZero(nv);
  od;
  return [orbit1,orbit,v,pivot];
end);

InstallGlobalFunction(SpinMatVector,function(mat,vec)
  local sp,minpol;
  sp:=SpinMatVector1(mat,vec,[],[],[]);
  minpol:=-SolutionMat(sp[2],sp[3]);
  Add(minpol,minpol[1]^0);
  return [sp[1],sp[2],minpol,sp[4]];
end);

# Returns triple [B,C,svec]  where C is such that  C*<A>*C^-1  has a block
# triangular shape with companian matrices along the diagonal), B is the
# row echelon form of C and svec is the list of indices where the blocks
# begin.
# For details about this function, see nofoma.gd.
InstallGlobalFunction(CyclicChainMat,function(A)
  local v,chain,j,idm,sp,svec,l;
  idm:=OneMutable(A);
  if IsLowerTriangularMat(A) then
    return [idm,idm,[1..NrRows(A)+1],[1..NrRows(A)]];
  fi;
  sp:=SpinMatVector1(A,idm[1],[],[],[]);
  svec:=[1,Length(sp[1])+1];
  j:=1;
  while Length(sp[1])<NrRows(A) do
    while j in Set(sp[4]) do
      j:=j+1;
    od;
    sp:=SpinMatVector1(A,idm[j],sp[1],sp[2],sp[4]);
    Add(svec,Length(sp[1])+1);
  od;
  return [sp[1],sp[2],svec];
end);

# Input: A block triangular matrix M with companion matrices along the diagonal
# and a list svec of containing the indices where the blocks begin
# Returns a list of the relative minimal polynomial of the block matrices on
# the diagonal, where each polynomial is given by its coefficients
BindGlobal("nfmRelMinPols",function(M,svec)
  local i,l,rpol,one;
  one:=OneOfBaseDomain(M);
  rpol:=[];
  for i in [1..Length(svec)-1] do
    l:=List(-M[svec[i+1]-1]{[svec[i]..svec[i+1]-1]});
    Add(l,one);
    Add(rpol,nfmPolCoeffs(l));
  od;
  return rpol;
end);

# This is OrdPoly from Neunhoeffer-Praeger
BindGlobal("nfmOrderPolM",function(M,svec,rpols,z,v)
  local i,f,v1,h,g,l;
  f:=[];
  v1:=ShallowCopy(v);
  for i in [z,z-1..1] do
    l:=v1{[svec[i]..svec[i+1]-1]};
    if not IsZero(l) then
      if IsOne(Degree(rpols[i])) then
        g:=rpols[i];
      else
        h:=nfmPolCoeffs(l);
        if IsZero(Degree(h)) then
          g:=rpols[i];
        else
          g:=Quotient(rpols[i],Gcd(rpols[i],h));
        fi;
      fi;
      Add(f,g);
      if i>1 then
        v1:=PolynomialToMatVec(M,g,v1);
      fi;
    fi;
  od;
  return Product(f);
end);

# Returns vector that has the same minimal polynomial as mat.
# For details about this function, see nofoma.gd.
InstallGlobalFunction(MaximalVectorMat,function(A)
  local M,k,idm,sp,l,np,i,v1,z,one,f,rpols,svec,lm,x;
  k:=DefaultFieldOfMatrix(A);
  one:=One(k);
  idm:=RowsOfMatrix(OneMutable(A));
  x:=Indeterminate(k);
  if IsDiagonalMat(A) then       # first deal with diagonal matrix
    if ForAll([2..NrRows(A)],i->A[i,i]=A[1,1]) then
      v1:=idm[1];
      np:=x-A[1,1];
    else
      v1:=CompatibleVector(A);
      for i in [1..Length(v1)] do v1[i] := one; od;
      l:=Set([1..NrRows(A)],i->A[i,i]);
      np:=Product(l, a -> x-a);
    fi;
    Info(Infonofoma,2,"#I Degree of minimal polynomial is ",Degree(np)," \n");
    return [v1,np];
  fi;
  sp:=CyclicChainMat(A);         # general case: transform to cyclic chain
  sp[2] := Matrix(sp[2], A);
  M:=sp[2]*A*sp[2]^-1;
  svec:=sp[3];
  rpols:=nfmRelMinPols(M,svec);
  lm:=[idm[1],rpols[1]];
  for z in [2..Length(svec)-1] do
    if Degree(lm[2])^3>Length(svec) or
         not IsZero(PolynomialToMatVec(M,lm[2],idm[svec[z]])) then
      f:=nfmOrderPolM(M,svec,rpols,z,idm[svec[z]]);
      if not IsZero(EuclideanRemainder(lm[2],f)) then
        lm:=LcmMaximalVectorMat(M,lm[1],idm[svec[z]],lm[2],f);
      fi;
    fi;
  od;
  #v:=lm[1]*sp[2];
  #if SpinMatVector(A,v)[3]<>CoefficientsOfUnivariatePolynomial(lm[2]) then
  #  Error("mist");
  #else
  #  Print("youpie ");
  #fi;
  Info(Infonofoma,2,"#I Degree of minimal polynomial is ",Degree(lm[2]),"\n");
  return [lm[1]*sp[2],lm[2]];
end);

# Modifies an already given complementary subspace to the complementary
# subspace defined by Jacob.
# For details about this function, see nofoma.gd.
InstallGlobalFunction(JacobMatComplement,function(T,d)
  local base,i,ii,j,F,tT;
  base:=OneMutable(T);
  tT:=TransposedMat(T);
  F:=[base[d]];
  for i in [2..d] do
    Add(F,F[i-1]*tT);
  od;
  #TriangulizeMat(F);
  for i in [1..d] do
    ii:=d+1-i;
    for j in [i+1..d] do
      if not IsZero(F[j][ii]) then
        AddRowVector(F[j],F[i],-F[j][ii]);
      fi;
    od;
  od;
  for i in [d+1..NrRows(T)] do
    for j in [1..d] do
      base[i,j]:=-F[d+1-j][i];
    od;
  od;
  return base;
end);

BindGlobal("BuildBlockDiagonalMat1",function(d,B)
  local new,n;
  n:=d+NrRows(B);
  new:=IdentityMatrix(n,B);
  CopySubMatrix(B, new, [1..NrRows(B)], [d+1..n], [1..NrCols(B)], [d+1..n]);
  return new;
end);

InstallGlobalFunction(RatFormStep1,function(A,v)
  local A1,sp,t,i,d,idm,minp;
  idm:=OneMutable(A);
  sp:=SpinMatVector1(A,v,[],[],[]);
  t:=Matrix(sp[2],A);
  d:=NrRows(t);
  Append(t,idm{Difference([1..NrRows(A)],sp[4])});
  A1:=t*A*t^(-1);
  minp:=-List(A1[d]{[1..d]});
  Add(minp,minp[1]^0);
  return [A1,t,minp];
end);

InstallGlobalFunction(RatFormStep1J,function(A,v)
  local A1,sp,t,i,j,d,idm,minp;
  idm:=OneMutable(A);
  sp:=SpinMatVector1(A,v,[],[],[]);
  t:=Matrix(sp[2],A);
  d:=NrRows(t);
  Append(t,idm{Difference([1..NrRows(A)],sp[4])});
  A1:=t*A*t^(-1);
  minp:=-List(A1[d]{[1..d]});
  Add(minp,minp[1]^0);
  j:=JacobMatComplement(A1,Length(minp)-1);
  #return [j*A1*j^-1,j*t,minp];
  return [ExtractSubMatrix(A1, [d+1..NrRows(A1)], [d+1..NrRows(A1)]),j*t,minp];
end);

#TODO: replace by CompanionMat
BindGlobal("nfmCompanionMat1", function(k,f)
  local n,i,one,mat;
  n:=Length(f)-1;
  one:=f[1]^0;
  mat:=ZeroMatrix(k,n,n);
  for i in [1..n-1] do
    mat[i,i+1]:=one;
    mat[n,i]:=-f[i];
  od;
  mat[n,n]:=-f[n];
  return mat;
end);

# Returns invariant factors of mat in descending order,
# invertible matrix P such that PAP^(-1) is in rational canonical form,
# and list of pivot indices.
# For details about this function, see nofoma.gd.
InstallGlobalFunction(FrobeniusNormalForm,function(A)
  local P,invf,i,mv,step1,rest,d,piv;
  if IsDiagonalMat(A) and ForAll([2..NrRows(A)],i->A[i,i]=A[1,1]) then
    mv:=nfmPolCoeffs([-A[1,1],A[1,1]^0]);
    return [ListWithIdenticalEntries(NrRows(A),mv),A^0,[1..NrRows(A)]];
  fi;
  mv:=MaximalVectorMat(A);
  invf:=[mv[2]];
  d:=Degree(mv[2]);
  if d=NrRows(A) then      # cyclic space: only one companion block
    step1:=[mv[1]];        # create base change
    for i in [2..d] do
      Add(step1,step1[i-1]*A);
    od;
    step1 := Matrix(step1, A);
    return [invf,step1,[1]];
  else
    step1:=RatFormStep1J(A,mv[1]);       # compute Jacob complement
    piv:=[1];
    rest:=FrobeniusNormalForm(step1[1]);
    P:=BuildBlockDiagonalMat1(d,rest[2]);
    Append(invf,rest[1]);
    for i in rest[3] do
      Add(piv,i+d);
    od;
    return [invf,P*step1[2],piv];
  fi;
end);

InstallGlobalFunction(FrobeniusNormalFormLikeRCFT, function(mat)
  local frob, blocks, n, Perm, Trans, ind, d, i, l;
  frob := FrobeniusNormalForm(TransposedMat(mat));
  blocks := frob[3];
  n := NrRows(mat);
  if Length(blocks) = 1 then
    return frob;
  fi;
  Perm := ZeroMatrix(n, n, mat);
  for i in [1..Length(blocks)-1] do
    ind := blocks[i];
    d   := blocks[i+1] - 1;
    l   := d - ind + 1;
    CopySubMatrix(IdentityMatrix(l, mat), Perm, [1..l], [ind..d], [1..l], [n-d+1..n-ind+1]);
  od;
  l := n - d;
  CopySubMatrix(IdentityMatrix(l, mat), Perm, [1..l], [Last(blocks)..n], [1..l], [1..n-Last(blocks)+1]);
  return [Reversed(frob[1]), TransposedMat(TransposedMat(Perm)*frob[2]), Reversed(frob[3])];
end);

# Returns the invariant factors of mat (i.e. the minimal polynomials of the
# diagonal blocks in the rational canonical form of mat).
# For details about this function, see nofoma.gd.
InstallGlobalFunction(InvariantFactorsMat,function(A)
  local A1,i,d,np,f,n,p;
  n:=NrRows(A);
  np:=MaximalVectorMat(A);

  d:=Degree(np[2]);
  f:=[np[2]];
  if d=n then
    return f;
  elif d=1 then
    return ListWithIdenticalEntries(n,np[2]);
  else
    A1:=RatFormStep1(A,np[1])[1];
    Append(f,InvariantFactorsMat(ExtractSubMatrix(A1, [d+1..n], [d+1..n])));
    return f;
  fi;
end);

## Now Jordan-Chevalley decomposition

# compute squarefree part sqf of f, that is: if f = \prod_{i=1}^k f_i^{n_i}
# where the f_i are pairwise coprime irreducible factor, then sqf is f_1
# \cdots f_k. Assumes that all coefficients of f are in the field K.
# Return a list [sqf,n] where n is is an integer such that f divides sqf^n
InstallGlobalFunction(SquareFreePol,function(K,f)
  local d,n,i,p,df,f1,g,g1,g2,cf,e;
  n:=Degree(f);
  if n=1 then
    return [f,1];
  fi;
  df:=Derivative(f);
  if not IsZero(df) then
    f1:=Gcd(f,df);

    # if f and its derivative are coprime, then f is squarefree
    if IsZero(Degree(f1)) then
      return [f,1];
    fi;

    g1:=SquareFreePol(K,f1);
    g2:=SquareFreePol(K,Quotient(f,f1));
    return [Lcm(g1[1],g2[1]),g1[2]+g2[2]];
  fi;

  # the derivative is zero, meaning all non-zero terms of f have an exponent
  # divisible by p. So compute a p-th root of the polynomial by dividing the
  # exponents by p, and taking p-th roots of each coefficient.
  cf:=CoefficientsOfUnivariatePolynomial(f);
  p:=Characteristic(K);
  e:=Size(K)/p;
  g:=SquareFreePol(K,nfmPolCoeffs(List([0..n/p],i->cf[p*i+1]^e)));
  return [g[1],p*g[2]];
end);

# Takes a matrix mat and a polynomial f such that f(mat) = 0 as input.
# Returns [D,N] such that D is diagonalisable over some extension field
# and N is a nilpotent matrix such that mat = D + N and DN=ND.
# For details about this function, see nofoma.gd.
InstallGlobalFunction(JordanChevalleyDecMat,function(A,f)
  local Ak,k0,gg,g,tg;
  gg:=SquareFreePol(DefaultField(CoefficientsOfUnivariatePolynomial(f)),f);
  g:=gg[1];
  tg:=GcdRepresentation(Derivative(gg[1]),gg[1])[1];
  k0:=0;
  Ak:=A;
  Info(Infonofoma,2,"Iterations (m=",gg[2],").");
  while 2^k0<gg[2] do
    Ak:=Ak-PolynomialToMat(Ak,g)*PolynomialToMat(Ak,tg);
    k0:=k0+1;
  od;
  return [Ak,A-Ak];
end);

InstallGlobalFunction(JordanChevalleyDecMatF,function(mat)
  local k,f,jc,p,N,D,src_range,dst_range;
  k:=DefaultFieldOfMatrix(mat);
  f:=FrobeniusNormalForm(mat);
  Info(Infonofoma,2,"Frobenius normal form complete\n");
  Add(f[3],NrRows(mat)+1);
  jc:=List(f[1],p->JordanChevalleyDecMat(nfmCompanionMat1(k,CoefficientsOfUnivariatePolynomial(p)),p));
  N:=0*f[2];
  D:=0*f[2];
  for p in [1..Length(f[1])] do
    src_range := [1..f[3][p+1]-f[3][p]];
    dst_range := [f[3][p]..f[3][p+1]-1];
    CopySubMatrix( jc[p][1], D, src_range, dst_range, src_range, dst_range );
    CopySubMatrix( jc[p][2], N, src_range, dst_range, src_range, dst_range );
  od;
  return [f[2]^-1*D*f[2],f[2]^-1*N*f[2]];
end);

##Jordan Normal form code from here

BindGlobal("nfmConvertVecToRowMat", function(vec)
    local mat, n, i;
    n := Length(vec);
    mat := ZeroMatrix(BaseDomain(vec), 1, n);
    for i in [1..n] do
        mat[1,i] := vec[i];
    od;
    return mat;
end);

#Avoids creating the vector space
BindGlobal("nfmGenerateRandomVector", function(d, A) #length d, fitting matrix A
    local bd;
    if IsPlistRep(A) then
      bd:=BaseDomain(A);
      return List([1..d], i -> Random(bd));
    fi;
    return Randomize(ZeroVector(d, A));
end);

BindGlobal("nfmGenerateNonZeroVector", function(d, A)
    local vec;
    repeat
        vec := nfmGenerateRandomVector(d, A);
    until not IsZero(vec);
    return vec;
end);

# TODO: this doesn't consistently generate non cyclic matrices, since glueing
# cyclic matrices together may result in a cyclic matrix again
BindGlobal("nfmGenerateNonCyclicMatrix", function(F,n) #Field F, dimension n
    local dim, A, num, subA, scr;
    A := ZeroMatrix(F,n,n);
    num := Random(1,n-1); #dont go until n because we dont want cyclic matrix
    subA := Matrix(F,RandomInvertibleMat(num,F));
    CopySubMatrix(subA, A, [1..num], [1..num], [1..num], [1..num]);
    dim := num;
    while not dim = n do
        num := Random(1,n-dim);
        subA := Matrix(F,RandomInvertibleMat(num,F));
        CopySubMatrix(subA, A, [1..num], [dim+1..dim+num], [1..num], [dim+1..dim+num]);
        dim := dim + num;
    od;
    scr := Matrix(F,RandomInvertibleMat(n,F)); #Conjugate matrix
    return A^scr;
end);

#Spinning algorithm
#Returns (vec, vec*A, ..., vec*A^(goal-1))
BindGlobal("nfmSpinUntil", function(vec, A, goal)
    local n, i, res, F;
    F := BaseDomain(A);
    if goal = 0 then
        return nfmConvertVecToRowMat(vec);
    fi;
    n := NrRows(A);
    res := ZeroMatrix(goal, n, A);
    CopySubMatrix(nfmConvertVecToRowMat(vec), res, [1..1], [1..1], [1..n], [1..n]);
    for i in [2..goal] do
        vec := vec * A;
        CopySubMatrix(nfmConvertVecToRowMat(vec), res, [1..1], [i..i], [1..n], [1..n]);
    od;
    return res;
end);

#Returns a vector v such that v is not in subspace spanned by gen
#Assumes that gen is already echelonised
BindGlobal("nfmFindVectorNotInSubspaceNC", function(gen) #assumes gen is already echelonised
    local w, i, n, F, r, zsf;
    r := NrRows(gen); #dimension of subspace
    F := BaseDomain(gen);
    n := Length(gen[1]);
    if (r = n) then
        return ZeroVector(n,gen);
    fi;
    zsf := ZeroMatrix(r+1,r+1,gen);
    CopySubMatrix(gen, zsf, [1..r], [1..r], [1..r], [1..r]);
    w := ZeroVector(n,gen);
    for i in [1..n] do
        if IsZero(zsf[i,i]) then
            w[i] := One(F);
            return w;
        fi;
    od;
    Error("Failed to find vector that is not in the generated subspace.");
    return fail;
end);

#Finds a cyclic vector for A
#Doesn't check if $A$ is cyclic
BindGlobal("nfmFindCyclicVectorNC", function(A) #Field F, Matrix A, n upper bound of loops;
    local checked,vec,gens,n,i,F;
    n := NrRows(A);
    F := BaseDomain(A);
    checked := []; #Avoid checking double
    vec := nfmGenerateRandomVector(n,A);
    for i in [1..40] do
        if not (vec in checked) then
            gens := nfmSpinUntil(vec, A, n); #potential basis
            if not (RankMat(gens) < n) then #check if linearly independent,
                return gens; #returns spinned basis, first element is the cyclic vector
            else
                Add(checked, vec);
            fi;
        fi;
        vec := nfmGenerateRandomVector(n,A);
    od;
    Info(Infonofoma,4,"Failed to find a cyclic vector!\n");
    return fail;
end);

BindGlobal("nfmRemoveZeroRows", function(mat)
    local n,i,matcopy;
    matcopy := MutableCopyMat(mat);
    n := NrRows(matcopy);
    for i in [n,n-1..1] do
        if IsZero(matcopy[i]) then
            Remove(matcopy, i);
        fi;
    od;
    return matcopy;
end);

#Evaluates v*p(A) VERY efficiently
#Takes (v,vA,...v^n-1A) and polynomial as input
BindGlobal("nfmPolyEvalFromSpan", function(span,pol)
    local coeffs, i, resu;
    coeffs := CoefficientsOfUnivariatePolynomial(pol);
    resu := coeffs[1]*span[1];
    for i in [2..Size(coeffs)] do
        resu := resu + coeffs[i]*span[i];
    od;
    return resu;
end);


BindGlobal("factoriseByKnownFactors", function(kFacs,pol)
  local fac,factup,i,resfacs,count,newpol,oldpol;
  resfacs := [];
  oldpol := pol;
  for factup in kFacs do
    fac := factup[1];
    count := 0;
    for i in [1..factup[2]] do
      newpol := Quotient(oldpol, fac);
      if newpol = fail then
        break;
      fi;
      count := count + 1;
      oldpol := newpol;
    od;
    if not count = 0 then
      Add(resfacs, [fac, count]);
    fi;
  od;
  return resfacs;
end);

#TODO: if size of minpolfacs is one we can return identity
#Primary Decomposition using a modified version of Allan Steel's algorithm for use in the Jordan normal form function
#Takes matrix along with its minimal polynomial and its factorised version as input
#Returns matrix B such that B*A*B^-1 is in primary decomp. form, dimensions of primary subspaces
#and factors of minimal polynomial in correct order
BindGlobal("nfmPrimaryDecompositionforJNF", function(A, minpol, minpolfacs)
    local rank,F,n,m,f,w,p,j,i,wspan,gens,facs,L_i,qi,k,v,
    COB,pot,gs,f2,toAdd,pos,dims,dim;
    rank := 0;
    n := NrRows(A);
    F := DefaultFieldOfMatrix(A);
    v := nfmGenerateNonZeroVector(n,A);
    gens := []; #Li_s as in Steel paper will go in here (but without separate U)
    gs := []; #distinct factors of minimal polynomial
    dims := [];
    while not rank = n do
        m := UnivariatePolynomial(F,SpinMatVector(A,v)[3]);
        p := One(PolynomialRing(F));
        for i in [1..Size(gens)] do
            L_i := gens[i];
            if not IsOne(Gcd(m,gs[i])) then
                f := m;
                pot := 0;
                for j in [1..n] do
                    f2 := Quotient(f,gs[i]);
                    if f2 = fail then
                        break;
                    fi;
                    pot := pot + 1;
                    f := f2;
                od;
                w := PolynomialToMatVec(A,f,v);
                p := p * gs[i]^pot;
                #minimal polynomial of w has degree smaller than or equal to minpol/f
                wspan := nfmSpinUntil(w,A,Degree(minpol)-Degree(f));
                toAdd := EcheloniseMat(Concatenation(wspan,gens[i]));
                if not IsMatrix(toAdd) then
                    toAdd := [toAdd];  # Convert vector to 1-row matrix
                    toAdd := Matrix(toAdd);
                fi;
                gens[i] := toAdd;
            fi;
        od;
        v := PolynomialToMatVec(A,p,v);
        m := Quotient(m,p);
        if not IsOne(m) then
            facs := factoriseByKnownFactors(minpolfacs,m);
            for i in [1..Size(facs)] do
                qi := (facs[i][1])^(facs[i][2]);
                w := PolynomialToMatVec(A,Quotient(m,qi),v);
                wspan := SpinMatVector1(A,w,[],[],[])[2];
                Add(gens, wspan);
                Add(gs, facs[i][1]);
            od;
        fi;
        # put everything into a single matrix
        COB := ZeroMatrix(F,n,n);
        k := 0;
        for L_i in gens do
            CopySubMatrix(L_i, COB, [1..NrRows(L_i)], [k+1..k+NrRows(L_i)], [1..n], [1..n]);
            k := k+NrRows(L_i);
        od;
        COB := nfmRemoveZeroRows(COB);
        rank := NrRows(COB);
        if not rank = n then
            COB := EcheloniseMat(COB);
            v := nfmFindVectorNotInSubspaceNC(COB);
        fi;
    od;
    COB := ZeroMatrix(F,n,n);
    #sort primary subspaces back into order and collect dims
    k := 0;
    for i in [1..Size(gens)] do
      pos := Position(gs, minpolfacs[i][1]);
      L_i := gens[pos];
      dim := NrRows(L_i); #?
      Add(dims, dim);
      CopySubMatrix(L_i, COB, [1..dim], [k+1..k+dim],[1..n], [1..n]);
      k := k + dim;
    od;
    return [COB, dims];
end);

#Primary Decomposition for cyclic matrices
#Jordan normal form will call this function if a cyclic matrix is detected
BindGlobal("nfmPrimaryDecompositionforJNFCyclic", function(A, minpol, minpolfacs)
    local vspan,n,w,mf,wspan,qi,k,COB,dims,elDivs;
    n := NrRows(A);
    vspan := nfmFindCyclicVectorNC(A);
    dims := [];
    COB := ZeroMutable(A);
    k := 0;
    elDivs := [];
    for mf in minpolfacs do
        qi := (mf[1])^(mf[2]);
        Add(elDivs, qi);
        w := nfmPolyEvalFromSpan(vspan,Quotient(minpol,qi));
        wspan := nfmSpinUntil(w,A,Degree(mf[1])*mf[2]);
        CopySubMatrix(wspan, COB, [1..NrRows(wspan)], [k+1..k+NrRows(wspan)], [1..n], [1..n]);
        Add(dims, NrRows(wspan));
        k := k + NrRows(wspan);
    od;
    return [COB, elDivs, dims];
end);

#TODO: recognise cyclic matrices
#Primary Decomposition using a modified version of Allan Steel's algorithm
#Standalone version
#Returns matrix B such that B*A*B^-1 is in primary decomposition form
#along with dimensions of primary subspaces
InstallMethod(PrimaryDecomposition,
    [ IsMatrixOrMatrixObj ],
function(A)
    local rank,F,n,m,f,w,p,j,i,wspan,gens,facs,L_i,qi,k,v,
    COB,pot,gs,f2,dims,toAdd,dim,minpol,collected;
    rank := 0;
    n := NrRows(A);
    F := DefaultFieldOfMatrix(A);
    v := nfmGenerateNonZeroVector(n,A);
    minpol := MinimalPolynomial(F,A);
    collected := Collected(Factors(minpol));
    if Degree(minpol) = n then
      return nfmPrimaryDecompositionforJNFCyclic(A,minpol,collected);
    fi;
    if IsZero(A) then
      return IdentityMat(n,F);
    fi;
    gens := []; #Li_s as in Steel paper will go in here
    gs := []; #distinct factors of minimal polynomial
    while not rank = n do
        m := UnivariatePolynomial(F,SpinMatVector(A,v)[3]);
        p := One(PolynomialRing(F));
        for i in [1..Size(gens)] do
            L_i := gens[i];
            if not IsOne(Gcd(m,gs[i])) then
                f := m;
                pot := 0;
                for j in [1..n] do
                    f2 := Quotient(f,gs[i]);
                    if f2 = fail then
                        break;
                    fi;
                    pot := pot + 1;
                    f := f2;
                od;
                w := PolynomialToMatVec(A,f,v);
                p := p * gs[i]^pot;
                #minimal polynomial of w has degree smaller than or equal to n - degree(f)
                wspan := nfmSpinUntil(w,A,n-Degree(f));
                toAdd := EcheloniseMat(Concatenation(wspan,gens[i]));
                if not IsMatrix(toAdd) then
                    toAdd := [toAdd];  # Convert vector to 1-row matrix
                    toAdd := Matrix(toAdd);
                fi;
                gens[i] := toAdd;
            fi;
        od;
        v := PolynomialToMatVec(A,p,v);
        m := Quotient(m,p);
        if not IsOne(m) then
            facs := Collected(Factors(m));
            for i in [1..Size(facs)] do
                qi := (facs[i][1])^(facs[i][2]);
                w := PolynomialToMatVec(A,Quotient(m,qi),v);
                wspan := SpinMatVector1(A,w,[],[],[])[2];
                Add(gens, wspan);
                Add(gs, facs[i][1]);
            od;
        fi;
        # put everything into a single matrix
        COB := ZeroMatrix(F,n,n);
        k := 0;
        for L_i in gens do
            CopySubMatrix(L_i, COB, [1..NrRows(L_i)], [k+1..k+NrRows(L_i)], [1..n], [1..n]);
            k := k+NrRows(L_i);
        od;
        COB := nfmRemoveZeroRows(COB);
        rank := NrRows(COB);
        if not rank = n then
            COB := EcheloniseMat(COB);
            v := nfmFindVectorNotInSubspaceNC(COB);
        fi;
    od;
    #sort blocks and count dims
    COB := ZeroMatrix(F,n,n);
    SortParallel(gs, gens);
    k := 0;
    dims := [];
    for L_i in gens do
      dim := NrRows(L_i);
      Add(dims, dim);
      CopySubMatrix(L_i, COB, [1..dim], [k+1..k+dim],[1..n], [1..n]);
      k := k + dim;
    od;
    return [COB, collected, dims];
end);

#Input: For matrix A with minimal polynomial p^m: m, vector v and p(A)
#Returns: v, length r of v, vp^(r-1)(A)
BindGlobal("GetMinPolPowerWithVec", function(m,v,Ainp)
    local j, veccopy, lastcopy;
    if IsZero(v) then
        return [v,0,v];
    fi;
    lastcopy := ShallowCopy(v);
    veccopy := ShallowCopy(v);
    for j in [1..m] do
        veccopy := veccopy * Ainp;
        if IsZero(veccopy) then #check if vector got nulled
            return [v, j, lastcopy];
        fi;
        lastcopy := veccopy;
    od;
    Info(Infonofoma,2,"Failed to find minimal polynomial of vector");
    return fail; #shouldn't happen
end);

#Returns linear dependence q_1,...,q_k as described in paper for cyclic decomposition
BindGlobal("FindLinearDependenceNC", function(vecs, A, d) #vecs, A, degree of p, returns coeffs of qis (ascending degree)
    local n,F,i,rel,tosolve,currdim,qis;
    F := DefaultFieldOfMatrix(A);
    n := NrRows(A);
    tosolve := ZeroMatrix(F,Length(vecs)*d,n);
    currdim := 1;
    for i in [1..Length(vecs)] do
        CopySubMatrix(nfmSpinUntil(vecs[i], A, d), tosolve, [1..d], [currdim..currdim +  d-1], [1..n], [1..n]);
        currdim := currdim +  d;
    od;
    rel := NullspaceMat(tosolve);
    if rel = fail or rel = [] then #shouldn't happen
        Error("Failed to find a linear dependence.");
    fi;
    rel := rel[1];
    #turn into usable polynomial coeffs
    qis := [];
    currdim := 1;
    for i in [1..Length(vecs)] do
        Add(qis, ExtractSubVector(rel, [currdim..currdim+d-1]));
        currdim := currdim + d;
    od;
    return qis;
end);

#Input: matrix A with minimalpolynomial p^m
#Returns matrix B such that A^Inverse(B) is in cyclic decomposition form, dimensions of cyclic subspaces
BindGlobal("CyclicDecompositionOfPrimarySubspace", function (A, p, m)
    local F,n,d,Ainp,ws,allspun,wspun,tomult,minpolpowers,
    dims,wtrip,i,conj,sumdim,qis,k,j,r,wstrich,currdim,w,vecs;
    F := DefaultFieldOfMatrix(A);
    n := NrRows(A);
    d := Degree(p);
    if m * d = n then #return if it's already cyclic
        return [One(A), [n]];
    fi;
    Ainp := p(A); #TODO: evaluate this using frobform? or maybe polyevalfromspan
    ws := [];
    w := nfmGenerateNonZeroVector(n,A);
    Add(ws,w);
    allspun := EcheloniseMat(nfmSpinUntil(w,A,n));
    while NrRows(allspun) < n do
        w := nfmFindVectorNotInSubspaceNC(allspun);
        wspun := nfmSpinUntil(w,A,n);
        Append(allspun, wspun);
        Add(ws, w);
        allspun := EcheloniseMat(allspun);
    od;
    minpolpowers := [];
    for i in [1..NrRows(ws)] do
        minpolpowers[i] := GetMinPolPowerWithVec(m,ws[i],Ainp);#for all ws: v, r such that p^r(A)(v)=0 and p^(r-1)(A)(v)
    od;
    dims := []; #dimensions of my cyclic subspaces
    sumdim := 0; #dimension of my direct sum
    conj := ZeroMutable(A);
    while not Sum(minpolpowers, v -> v[2]*d) = n do #generally bigger than n, working our way down
        SortBy(minpolpowers, v -> v[2]); #Sort ws by their A-length (ascending)
        vecs := List([1..Length(minpolpowers)], i -> minpolpowers[i][3]); #p(A)^(r-1)(w)
        qis := FindLinearDependenceNC(vecs,A,d); #coefficients of qis as described in theorem
        j := PositionNonZero(qis); #find j as described in theorem (works because we sorted this list beforehand)
        r := minpolpowers[j][2]; #r as described in theorem
        wstrich := ZeroVector(F,n);
        for i in [1..Length(qis)] do #calculate new reduced w_j
            if not IsZero(qis[i]) then #avoid unnecessary computations
                tomult :=  PolynomialToMatVec(A, nfmPolCoeffs(qis[i]), minpolpowers[i][1]);
                for k in [1..minpolpowers[i][2] - r] do
                    tomult := tomult * Ainp; #this too perhaps?
                od;
                wstrich := wstrich + tomult;
            fi;
        od;
        if IsZero(wstrich) then
            Remove(minpolpowers,j); #if it is zero vec we get no more information from it
        else
            minpolpowers[j] := GetMinPolPowerWithVec(m,wstrich,Ainp);
        fi;
    od;
    currdim := 1;
    SortBy(minpolpowers, v -> v[2]); #sort by A-length
    for i in [1..Length(minpolpowers)] do
        wtrip := minpolpowers[i];
        conj{[currdim..currdim+(wtrip[2]*d)-1]}{[1..n]} := nfmSpinUntil(wtrip[1],A,wtrip[2]*d);
        Add(dims, wtrip[2]*d);
        currdim := currdim + wtrip[2]*d;
    od;
    conj := Matrix(F,conj);
    return [conj,dims];
end);

#Input: Cyclic matrix A with minimal polynomial p^m
#Returns matrix B such that A^Inverse(B) is in Jordan block form
BindGlobal("JordanBlock", function(A, p, m) #For JordanNormalform
    local i,spun,n,basis,d,r;
    n := NrRows(A);
    d := Degree(p);
    spun := nfmFindCyclicVectorNC(A);
    basis := ZeroMutable(A);
    CopySubMatrix(spun,basis,[1..d],[1..d],[1..n],[1..n]);
    for r in [1..d] do
        for i in [1..m-1] do
            basis[d*i+r] := nfmPolyEvalFromSpan(ExtractSubMatrix(spun,[r..n],[1..n]),p^i);
        od;
    od;
    return basis;
end);

#Input: Matrix A with irreducible Minimal polynomial
#Returns matrix B such that A^Inverse(B) is in Jordan normal form
InstallGlobalFunction(JordanNormalFormIrred, function(A,minpol)
    local F,n,cobrank,COB,blockdim,spun,v,w,elDivs;
    n := NrRows(A);
    F := DefaultFieldOfMatrix(A); # get underlying field
    blockdim := Degree(minpol);
    if blockdim = n then
      return [nfmFindCyclicVectorNC(A), [minpol]];
    fi;
    v := nfmGenerateNonZeroVector(n,A);
    spun := nfmSpinUntil(v, A, blockdim);
    COB := ZeroMutable(A);
    CopySubMatrix(spun, COB, [1..blockdim],[1..blockdim],[1..n],[1..n]);
    cobrank := blockdim;
    elDivs := [minpol];  
    while not cobrank = n do 
        w := nfmFindVectorNotInSubspaceNC(
            EcheloniseMat(COB{[1..cobrank]}{[1..n]})
        );
        spun := nfmSpinUntil(w,A,blockdim);
        CopySubMatrix(spun, COB, [1..blockdim],[cobrank+1..cobrank+blockdim],[1..n],[1..n]);
        cobrank := cobrank + blockdim;
        Add(elDivs, minpol);
    od;
    return [COB,elDivs];
end);

#Input: Matrix A
#Returns matrix B such that A^Inverse(B) is in Jordan normal form
InstallGlobalFunction(JordanNormalForm, function(A)
    local n,F,pol,minpol,primary,primarydims,crhr,cyclicdims,elDivs,prim1,
    subsubCOB,COB,subA,cy,i,j,facOcc,subCOB,crcy,subsubA,prepreCOB,preCOB;
    F := DefaultFieldOfMatrix(A);
    n := NrRows(A);
    elDivs := [];
    if IsZero(A) then
      return [OneMutable(A), UnivariatePolynomial(F,One(F))];
    fi;
    minpol := MinimalPolynomial(F,A);
    facOcc := Collected(Factors(minpol));  #factors of minimalpolynomial and their multiplicity
    if Size(facOcc) = 1 and facOcc[1][2] = 1 then
        return JordanNormalFormIrred(A,minpol);
    fi;
    if Degree(minpol) = n then
        prim1 := nfmPrimaryDecompositionforJNFCyclic(A, minpol, facOcc);
        primary := [prim1[1],prim1[3]];
    else
        primary := nfmPrimaryDecompositionforJNF(A, minpol, facOcc);
    fi;
    primarydims := primary[2]; #dimensions of generalized eigenspaces
    COB := primary[1]; #Change of basis matrix, this will be the final COB from A to JNF
    A := COB*A*Inverse(COB); #A in hauptraumform
    crhr := 1; #current row (primary subspace)
    preCOB := ZeroMutable(A); #COB matrix from A in primary form to primary subspaces in cyclic form
    for i in [1..Size(primarydims)] do
        if primarydims[i] = 1 then
            preCOB[crhr,crhr] := One(F);
            crhr := crhr + 1;
            Add(elDivs, facOcc[i][1]);
            continue;
        fi;
        pol := facOcc[i][1];
        subA := ExtractSubMatrix(A,[crhr..crhr+primarydims[i]-1],[crhr..crhr+primarydims[i]-1]);
        if facOcc[i][2] = 1 then
            cy := JordanNormalFormIrred(subA,pol);
            elDivs := Concatenation(elDivs,cy[2]);
            CopySubMatrix(cy[1], preCOB, [1..primarydims[i]], [crhr..crhr+primarydims[i]-1], [1..primarydims[i]], [crhr..crhr+primarydims[i]-1]);
            crhr := crhr + primarydims[i];
            continue;
        fi;
        cy := CyclicDecompositionOfPrimarySubspace(subA, pol, facOcc[i][2]); #decompose primary spaces into cyclic ones
        cyclicdims := cy[2]; #dimensions of cyclic subspaces
        subCOB := cy[1]; #subCOB to be assembled
        subCOB := Matrix(F,subCOB); 
        subA := subCOB*subA*Inverse(subCOB); #subA in cyclic decomposition form
        crcy := 1; #current row (cyclic subspace)
        prepreCOB := ZeroMatrix(F,primarydims[i], primarydims[i]);
        for j in [1..Size(cyclicdims)] do
            if cyclicdims[j] = 1 then
                prepreCOB[crcy,crcy] := One(F);
                crcy := crcy + 1;
                Add(elDivs,pol);
                continue;
            fi;
            subsubA := ExtractSubMatrix(subA, [crcy..crcy+cyclicdims[j]-1], [crcy..crcy+cyclicdims[j]-1]);
            subsubCOB := JordanBlock(subsubA,pol,cyclicdims[j]/Degree(pol));
            CopySubMatrix(subsubCOB, prepreCOB, [1..cyclicdims[j]], [crcy..crcy+cyclicdims[j]-1], [1..cyclicdims[j]], [crcy..crcy+cyclicdims[j]-1]);
            crcy := crcy + cyclicdims[j];
            Add(elDivs,pol^(cyclicdims[j]/Degree(pol)));
        od;
        subCOB := prepreCOB * subCOB;
        CopySubMatrix(subCOB, preCOB, [1..primarydims[i]], [crhr..crhr+primarydims[i]-1], [1..primarydims[i]], [crhr..crhr+primarydims[i]-1]);
        crhr := crhr + primarydims[i];
    od;
    return [preCOB*COB, elDivs];
end);
