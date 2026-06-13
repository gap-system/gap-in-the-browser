#############################################################################
##
#W  qo-example.tst       GAP4 package IBNP       Gareth Evans & Chris Wensley
##

gap> START_TEST( "qo-example.tst" );
gap> ibnp_infolevel_saved := InfoLevel( InfoIBNP );; 
gap> SetInfoLevel( InfoIBNP, 0 );; 

gap> R := PolynomialRing( GF(2), ["x","y"] );;
gap> e := One(R);;  x:=R.1;;  y:=R.2;;
gap> ord := MonomialGrlexOrdering( [x,y] );;
gap> L2 := [ x^6*y^6+x^2*y^2, x^6+x ];;
gap> p := x^6*y^7 + x*y^7;;
gap> q := y^8;;
gap> CommutativeDivision := "Thomas";;
gap> ibasT := InvolutiveBasis( R, L2, ord );
rec( div := "Thomas", 
  mvars := [ [ 1 ], [ 2 ], [ 1 ], [ 2 ], [ 1 ], [ 2 ], [ 1 ], [ 2 ], [ 1 ], 
      [ 2 ], [ 1 ], [ 1, 2 ] ], 
  polys := [ x^6+x, x*y^6+x^2*y^2, x^6*y+x*y, x^2*y^6+x^3*y^2, x^6*y^2+x*y^2, 
      x^3*y^6+x^4*y^2, x^6*y^3+x*y^3, x^4*y^6+x^5*y^2, x^6*y^4+x*y^4, 
      x^5*y^6+x*y^2, x^6*y^5+x*y^5, x^6*y^6+x^2*y^2 ] )
gap> rp := LoggedIPolyReduceCP( R, p, ibasT, ord );
rec( logs := [ 0*Z(2), y, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 
      0*Z(2), 0*Z(2), 0*Z(2), y ], 
  polys := [ x^6+x, x*y^6+x^2*y^2, x^6*y+x*y, x^2*y^6+x^3*y^2, x^6*y^2+x*y^2, 
      x^3*y^6+x^4*y^2, x^6*y^3+x*y^3, x^4*y^6+x^5*y^2, x^6*y^4+x*y^4, 
      x^5*y^6+x*y^2, x^6*y^5+x*y^5, x^6*y^6+x^2*y^2 ], result := 0*Z(2) )
gap> rq := IPolyReduceCP( R, q, ibasT, ord );
y^8
gap> ## note that none of the powers of y reduce

gap> CommutativeDivision := "Janet";;
gap> ibasJ := InvolutiveBasis( R, L2, ord );
rec( div := "Janet", 
  mvars := [ [ 1 ], [ 2 ], [ 1 ], [ 2 ], [ 1 ], [ 2 ], [ 1 ], [ 2 ], [ 1 ], 
      [ 2 ], [ 1 ], [ 1, 2 ] ], 
  polys := [ x^6+x, x*y^6+x^2*y^2, x^6*y+x*y, x^2*y^6+x^3*y^2, x^6*y^2+x*y^2, 
      x^3*y^6+x^4*y^2, x^6*y^3+x*y^3, x^4*y^6+x^5*y^2, x^6*y^4+x*y^4, 
      x^5*y^6+x*y^2, x^6*y^5+x*y^5, x^6*y^6+x^2*y^2 ] )
gap> ibasT!.polys = ibasJ!.polys;
true
gap> ibasT!.mvars = ibasJ!.mvars;
true
gap> ##  so the same reductions apply using either division:
gap> rp := LoggedIPolyReduceCP( R, p, ibasJ, ord );
rec( logs := [ 0*Z(2), y, 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 0*Z(2), 
      0*Z(2), 0*Z(2), 0*Z(2), y ], 
  polys := [ x^6+x, x*y^6+x^2*y^2, x^6*y+x*y, x^2*y^6+x^3*y^2, x^6*y^2+x*y^2, 
      x^3*y^6+x^4*y^2, x^6*y^3+x*y^3, x^4*y^6+x^5*y^2, x^6*y^4+x*y^4, 
      x^5*y^6+x*y^2, x^6*y^5+x*y^5, x^6*y^6+x^2*y^2 ], result := 0*Z(2) )
gap> rq := IPolyReduceCP( R, q, ibasT, ord );
y^8

gap> CommutativeDivision := "Pommaret";;
gap> ibasP := InvolutiveBasisCP( R, L2, ord );
#I  reached the involutive abort limit 20
fail

gap> SetInfoLevel( InfoIBNP, ibnp_infolevel_saved );; 
gap> STOP_TEST( "qo-example.tst", 10000 );

#############################################################################
##
#E  qo-example.tst . . . . . . . . . . . . . . . . . . . . . . . . ends here
