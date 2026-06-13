#############################################################################
##
#W  obsolete.gi          GAP 4 package AtlasRep                 Thomas Breuer
##
##  This file contains implementations of global variables
##  that had been documented in earlier versions of the AtlasRep package.
##


#############################################################################
##
#F  AtlasOfGroupRepresentationsTestClassScripts( ... )
#F  AtlasOfGroupRepresentationsTestCompatibleMaxes( ... )
#F  AtlasOfGroupRepresentationsTestFileHeaders( ... )
#F  AtlasOfGroupRepresentationsTestFiles( ... )
#F  AtlasOfGroupRepresentationsTestGroupOrders( ... )
#F  AtlasOfGroupRepresentationsTestStdCompatibility( ... )
#F  AtlasOfGroupRepresentationsTestSubgroupOrders( ... )
#F  AtlasOfGroupRepresentationsTestWords( ... )
##
##  These functions are deprecated since version 1.5 of the package.
##
BindGlobal( "AtlasOfGroupRepresentationsTestClassScripts",
    AGR.Test.ClassScripts );
BindGlobal( "AtlasOfGroupRepresentationsTestCompatibleMaxes",
    AGR.Test.CompatibleMaxes );
BindGlobal( "AtlasOfGroupRepresentationsTestFileHeaders",
    AGR.Test.FileHeaders );
BindGlobal( "AtlasOfGroupRepresentationsTestFiles",
    AGR.Test.Files );
BindGlobal( "AtlasOfGroupRepresentationsTestGroupOrders",
    AGR.Test.GroupOrders );
BindGlobal( "AtlasOfGroupRepresentationsTestStdCompatibility",
    AGR.Test.StdCompatibility );
BindGlobal( "AtlasOfGroupRepresentationsTestSubgroupOrders",
    AGR.Test.MaxesOrders );
BindGlobal( "AtlasOfGroupRepresentationsTestWords",
    AGR.Test.Words );


#############################################################################
##
#E

