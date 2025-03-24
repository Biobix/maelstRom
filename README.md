
<!-- README.md is generated from README.Rmd. Please edit that file -->

# maelstRom

<!-- badges: start -->
<!-- badges: end -->

Welcome to maelstRom’s readme!

```

                                                                                        ,,ggAAggg,
                                                                                      gAACAAAAAAGAA
                                                                                    /U¨       `*GUCU
                           _         _      _____                          ,gAUG¨               ¨
                          | ||      | ||   | ___ \\                       AGUG'   ,gggGGGGGgg
  _ __ ___    __ _   ___  | || ___  | ||_  | |_/ //  ___   _ __ ___      ]GA'   gAAACAAAAUCCCAA
 | '_ ` _ \\ / _` ||/ _ \\| ||/ __||| ___|||  _ <<  / _ \\| '_ ` _ \\     C'  gAAAAU¨'      ¨]C[
 | ||| ||| || ((| || ___//| ||\__ \\| ||_  | | \ \\| (() || ||| ||| ||       CCAGG            *
 |_|||_|||_||\__,_||\___|||_|||___//\____|||_|  \_\\\___//|_|||_|||_||       UCAC    ,aggg,
                                                                              ]AC   AAAACCAAa
                                                                                    *'-¨*ACAA
                                                                                         ]CAC
         G C C                           U G G                           G C A           CAA
     A U | | | U U       | | |       C C | | | U A       | | |       A G | | | A A     _]U`
 | C | | | | | | | A | | | | | | | G | | | | | | | A | | | | | | | U | | | | | | | Co-´
 C       | | |       C G | | | U A       | | |       C G | | | C A       | | |       
                          A C A                           A C A

 Modeller of allele-specific transcriptomics maelstRom, version 1.1.11
    
```

maelstRom is the Modeller of Allele-Specific Transcriptomics, an R package providing
extensive functions for various RNAseq-based allelic analyses. This
ranges from basic tasks such as (solely) RNAseq-based genotyping, to the
analysis of more complex population-level phenomena such as
(differential) allelic bias, allellic divergence and (loss of)
imprinting analyses. More information about what these are and how to
model them using maelstRom can be found in the package vignette.

## System requirements

### Hardware requirements

maelstRom requires a standard computer (server) with sufficient RAM for
in-memory operations.

### OS requirements

This package is supported for Windows, Linux and MacOS. The package has
been tested on the following systems:

-   Linux: Ubuntu 18.04 (bionic), Ubuntu 16.04 (Xenial)
-   Windows: Windows 10

### Other software requirements

maelstRom is an R software package with additional C/C++ code under-the-hood,
and as such relies on (versions listed are those used for maelstRom’s latest
test, though more recent versions should work fine as well):

-   CMake (\>= v3.0.0)
-   R (\>= v4.0.2)
-   GNU multiple precision arithmetic library (gmp; \>= v6.1.2)
-   Boost Multiple Precision Floating-Point Reliable Library (mpfr;
    \>= v4.0.1)

Furthermore, maelstRom relies on other R packages for its operations, which
are listed in the tutorial:
<https://biobix.github.io/maelstRom/articles/maelstRom_Allelic_Dispersion_tutorial.html#session-info>

## Installation

If you want to install maelstRom in a conda environment including an
R-installation, the following bash commands sets up the previously
listed software dependencies, then opens the environment, and launches
R. **An important** exception to this is CMake, which needs to be
installed on your machine.

``` bash
conda create -n <envname> r-essentials r-base r-devtools gmp mpfr
conda activate <envname>
R
```

Setting up this environment takes about 9 minutes on a linux-based
server (Ubuntu 18.04; bionic)

maelstRom source code can be found on its associated Github page:
<https://github.com/Biobix/maelstRom>. As such, maelstRom can be installed using
the `install_github` function from R’s `devtools` package:

``` r
library(devtools)
install_github("BioBix/maelstRom")
```

From a clean R installation (e.g. as in the conda environment set up
previously), this installation takes about 2 minutes.

For a local installation, rtools will be required. In case the code
chunk above throws the error `ERROR: loading failed for 'i386'`, running
the following instead potentially solves the issue:

``` r
library(devtools)
install_github("BioBix/maelstRom", INSTALL_opts=c("--no-multiarch"))
```

## Getting started

maelstRom contains a vignette going over its entire dAD anaylis pipeline, which
can be found here:
<https://biobix.github.io/maelstRom/articles/maelstRom_Allelic_Dispersion_tutorial.html>


Besides this vignette, there's also a tutorial on further exploration of maelstRom's dAD results:
<https://biobix.github.io/maelstRom/articles/maelstRom_results_exploration.html>
as well as a tutorial on (differential) imprinting analyses using maelstRom:
<https://biobix.github.io/maelstRom/articles/maelstRom_imprinting_tutorial.html>.

Running the entire vignette takes about an hour on a standard computer
(local installation; Windows 10) if parallellization is not used
(single-core). Enabling parallellization when running maelstRom on
large-scale datasets on server is advised.

maelstRom has an associated paper pending publication. Coming soon…


## License

This repo's main license, AGPL-3.0, applies to the entire repo,
except for subfolders that have their own license file. 
In such cases, the license file in the subfolder takes precedence.

## Contact

For theoretical/technical questions or issues while using the maelstRom
package, please contact me at <cedric.stroobandt@ugent.be>; for
everything else, contact Tim De Meyer at <tim.demeyer@ugent.be>.
