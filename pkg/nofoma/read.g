#
# nofoma: Normal forms of matrices
#
# Reading the implementation part of the package.
#
ReadPackage( "nofoma", "gap/nofoma.gi");


#
# for compatibility with MatrixObj in  GAP versions before 4.16, we need
# to add this "missing methods" (the downranking is to hopefully minimize
# conflicts with any other methods installed e.g. by packages)
#
if not CompareVersionNumbers(GAPInfo.Version, "4.16") then

  InstallOtherMethod( MinimalPolynomial, [ IsMatrixObj ], -1000,
      m -> MinimalPolynomial( BaseDomain( m ), m, 1 ) );
  InstallOtherMethod( MinimalPolynomial, [ IsRing, IsMatrixObj ], -1000,
      {R, m} -> MinimalPolynomial( R, m, 1 ) );
  InstallOtherMethod( MinimalPolynomial, [ IsRing, IsMatrixObj, IsPosInt ], -1000,
      {R, m, inum} -> MinimalPolynomial( R, Unpack( m ), inum ) );

fi;
