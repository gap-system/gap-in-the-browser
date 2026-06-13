#############################################################################
##
#W  utils.tst            GAP 4 package AtlasRep                 Thomas Breuer
##
##  In order to run the tests, one starts GAP from the 'tst' subdirectory
##  of the 'pkg/atlasrep' directory, and calls 'Test( "utils.tst" );'.
##
gap> START_TEST( "utils.tst" );

# Load the necessary packages.
gap> LoadPackage( "atlasrep", false );
true
gap> LoadPackage( "ctbllib", false );
true

# Test the special 'ConjugacyClasses' and 'CharacterTable' methods.
gap> G:= AtlasGroup( "M11" );;
gap> ccl:= ConjugacyClasses( G );;
gap> List( ccl, Representative ) =
>    ResultOfStraightLineProgram( AtlasProgram( "M11", "classes" ).program,
>                                 GeneratorsOfGroup( G ) );
true
gap> if IsPackageLoaded( "CTblLib" ) then
>      t:= CharacterTable( G );;
>      if Identifier( t ) <> "M11" then
>        Error( "wrong OrdinaryCharacterTable method chosen" );
>      fi;
>    fi;

# Test that 'AGR_ChecksumFits' treats 'HexSHA256' values with leading zeros
# correctly.
gap> str:= "mu 2 1 3\nmu 3 3 4\nmu 2 4 3\nmu 3 2 4\nmu 3 4 2\n";;
gap> checksum:= HexSHA256( str );
gap> if checksum[1] = '0' then
>      Remove( checksum );
>    fi;
gap> Length( checksum ) < 64;
true
gap> AGR_ChecksumFits( str, checksum );
true

# Done.
gap> STOP_TEST( "utils.tst" );


#############################################################################
##
#E

