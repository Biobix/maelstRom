
<!-- README.md is generated from README.Rmd. Please edit that file -->

# maelstRom

<!-- badges: start -->
<!-- badges: end -->

Welcome to maelstRom’s readme!

maelstRom is the Modeller of Allelic Gene Expression, an R package providing
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

-   R (v4.0.2)
-   GNU multiple precision arithmetic library (gmp; v6.1.2)
-   Boost Multiple Precision Floating-Point Reliable Library (mpfr;
    v4.0.1)

Furthermore, maelstRom relies on other R packages for its operations, which
are listed in the tutorial:
<https://biobix.github.io/maelstRom/articles/maelstRom_tutorial.html#session-info>

## Installation

If you want to install maelstRom in a conda environment including an
R-installation, the following bash commands sets up the previously
listed software dependencies, then opens the environment, and launches
R:

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

From a clean R installation (e.g. as in the conda environment set up
previously), this installation takes about 2 minutes.

For a local installation, rtools will be required. In case the code
chunk above throws the error `ERROR: loading failed for 'i386'`, running
the following instead potentially solves the issue:

``` r
library(devtools)
install_github("BioBix/maelstRom", INSTALL_opts=c("--no-multiarch"))
```

## Getting started

maelstRom contains a vignette going over its entire anaylis pipeline, which
can be found here:
<https://biobix.github.io/maelstRom/articles/maelstRom_tutorial.html>

Besides this “regular” vignette an “expanded” one is included as well
(<https://biobix.github.io/maelstRom/articles/maelstRom_expanded_tutorial.html>),
containing more in-depth code using maelstRom’s base functions while the
regular one uses wrapper functions that handle many analyses and
intermediary steps under-the hood. As such, the regular vignette is
recommended for first-time users or users that just want a plug-and-play
pipeline, while the expanded vignette provides more insight into maelstRom’s
analyses and is recommended when setting out to create your own
custom/specialized analysis pipeline.

Running the entire vignette takes about 2.5 hours on a standard computer
(local installation; Windows 10) if parallellization is not used
(single-core). Enabling parallellization when running maelstRom on
large-scale datasets on server is advised.

maelstRom has an associated paper pending publication. Coming soon…

## Contact

For theoretical/technical questions or issues while using the maelstRom
package, please contact me at <cedric.stroobandt@ugent.be>; for
everything else, contact Tim De Meyer at <tim.demeyer@ugent.be>.
