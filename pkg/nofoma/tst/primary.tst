gap> START_TEST("PrimaryDecompositionosition");
gap> ReadPackage("nofoma", "tst/utils.g");
true

##Test Primary decomposition
gap> nfmCheckPrimaryDecomposition(GF(5),50);
true
gap> nfmCheckPrimaryDecomposition(GF(25),25);
true
gap> nfmCheckPrimaryDecomposition(GF(7),150);
true
gap> nfmCheckPrimaryDecompositionNonCyclic(GF(5),50);
true
gap> nfmCheckPrimaryDecompositionNonCyclic(GF(5^5),10);
true

## Stop test
gap> STOP_TEST("PrimaryDecompositionosition");
