#
# nofoma: Normal forms of matrices
#
# This file contains package meta data. For additional information on
# the meaning and correct usage of these fields, please consult the
# manual of the "Example" package as well as the comments in its
# PackageInfo.g file.
#
SetPackageInfo( rec(

PackageName := "nofoma",
Subtitle := "Normal forms of matrices",
Version := "1.0.1",
Date := "17/05/2026", # dd/mm/yyyy format
License := "GPL-2.0-or-later",

Persons := [
  rec(
    FirstNames := "Meinolf",
    LastName := "Geck",
    WWWHome := "https://pnp.mathematik.uni-stuttgart.de/idsr/idsr1/geckmf/",
    Email := "meinolf.geck@mathematik.uni-stuttgart.de",
    IsAuthor := true,
    IsMaintainer := true,
    PostalAddress := "Fachbereich Mathematik, Pfaffenwaldring 57, 70569 Stuttgart, Germany",
    Place := "Stuttgart",
    Institution := "University of Stuttgart",
  ),
  rec(
    FirstNames := "Alia",
    LastName := "Bonnet",
    Email := "alia.bonnet@rwth.aachen.de",
    IsAuthor := true,
    IsMaintainer := true,
    Place := "Aachen",
    Institution := "RWTH Aachen"
  )
],

SourceRepository := rec(
    Type := "git",
    URL := Concatenation( "https://github.com/gap-packages/", ~.PackageName ),
),
IssueTrackerURL := Concatenation( ~.SourceRepository.URL, "/issues" ),
PackageWWWHome  := Concatenation( "https://gap-packages.github.io/", ~.PackageName ),
README_URL      := Concatenation( ~.PackageWWWHome, "/README.md" ),
PackageInfoURL  := Concatenation( ~.PackageWWWHome, "/PackageInfo.g" ),
ArchiveURL      := Concatenation( ~.SourceRepository.URL,
                                 "/releases/download/v", ~.Version,
                                 "/", ~.PackageName, "-", ~.Version ),
ArchiveFormats := ".tar.gz .tar.bz2",

##  Status information. Currently the following cases are recognized:
##    "accepted"      for successfully refereed packages
##    "submitted"     for packages submitted for the refereeing
##    "deposited"     for packages for which the GAP developers agreed
##                    to distribute them with the core GAP system
##    "dev"           for development versions of packages
##    "other"         for all other packages
##
Status := "dev",

AbstractHTML   :=
  "This package computes the Frobenius normal form and\
  the Jordan-Chevalley decomposition of a (square) matrix over any field\
  that is available in GAP. It also computes the Jordan normal form\
  of matrices over finite fields.",

BannerString := Concatenation(
"──────────────────────────────────────────────────────────────────────────\n",
"Loading  nofoma ", ~.Version, " (Normal forms of matrices), \n",
"by Meinolf Geck (https://pnp.mathematik.uni-stuttgart.de/idsr/idsr1/geckmf/)\n",
"and Alia Bonnet (https://github.com/AliaBonnet).\n",
"Help about the main functions in this package is obtained by typing:\n",
"    ?FrobeniusNormalForm     ?JordanNormalform     ?JordanChevalleyDecMat\n",
"──────────────────────────────────────────────────────────────────────────\n"),

PackageDoc := rec(
  BookName  := "nofoma",
  ArchiveURLSubset := ["doc"],
  HTMLStart := "doc/chap0_mj.html",
  PDFFile   := "doc/manual.pdf",
  SixFile   := "doc/manual.six",
  LongTitle := "Normal forms of matrices",
),

Dependencies := rec(
  GAP := ">= 4.15",
  NeededOtherPackages := [
    ["AutPGrp", ">= 1.5"],
  ],
  SuggestedOtherPackages := [ ],
  ExternalConditions := [ ],
),

AvailabilityTest := ReturnTrue,

TestFile := "tst/testall.g",

Keywords := [
    "maximal vectors",
    "Frobenius normal form",
    "Jordan normal form",
    "Jordan-Chevalley decomposition",
],

AutoDoc := rec(
    TitlePage := rec(
        # We reuse the AbstractHTML for the GAPDoc title page abstract -- we can
        # only do that as long as we are careful and keep it to pure ASCII, as any
        # use of HTML won't work inside of GAPDoc.
        Abstract := Concatenation(~.AbstractHTML,
            """
            <P/>
            Bug reports, suggestions and comments are, welcome.
            Please submit them to our issue tracker
            """,
            "<URL>", ~.IssueTrackerURL, "</URL>."
        ),

        Copyright := """
            &copyright; 2026 by Meinolf Geck.
            The &nofoma; package is free software; you can redistribute it and/or
            modify it under the terms of the GNU General Public License as published
            by the Free Software Foundation; either version 2 of the License, or (at
            your option) any later version.
        """,
        Acknowledgements := """
            The original code by Meinolf Geck was cleaned up and made ready for inclusion
            in the &GAP; package distribution Alia Bonnet in 2025, with some assistance
            by Max Horn.
        """,
    ),
),

));
