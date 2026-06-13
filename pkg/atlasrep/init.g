#############################################################################
##
#W  init.g               GAP 4 package AtlasRep                 Thomas Breuer
##

# Read the declaration part.
ReadPackage( "atlasrep", "gap/userpref.g"  );
ReadPackage( "atlasrep", "gap/bbox.gd"     );
ReadPackage( "atlasrep", "gap/access.gd"   );
if not IsBound( InfoCMeatAxe ) then
  # This file is also part of the C-MeaAxe package.
  ReadPackage( "atlasrep", "gap/scanmtx.gd" );
  ReadPackage( "atlasrep", "gap/scanmtx.gi" );
fi;
ReadPackage( "atlasrep", "gap/types.gd"    );
ReadPackage( "atlasrep", "gap/interfac.gd" );
ReadPackage( "atlasrep", "gap/mindeg.gd"   );
ReadPackage( "atlasrep", "gap/utils.gd"    );


#############################################################################
##
#E

