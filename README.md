
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MAGE

<!-- badges: start -->
<!-- badges: end -->

Welcome to MAGE’s readme!

MAGE is the Modeller of Allelic Gene Expression, an R package providing
extensive functions for various RNAseq-based allelic analyses. This
ranges from basic tasks such as (solely) RNAseq-base genotyping, to the
analysis of more complex population-level phenomena such as
(differential) allelic bias, allellic divergence and (loss of)
imprinting analyses. More information about what these are and how to
model them using MAGE can be found in the package vignette.

## Installation

MAGE can be installed using the `install_github` function from the
`devtools` package:

``` r
library(devtools)
install_github("BioBix/MAGE")
```

## Getting started

MAGE contains a vignette going over its entire anaylis pipeline, which
can be found under “articles” on <https://biobix.github.io/MAGE>.

Both a “regular” and “expanded” vignette are included, the latter
writing out all analyses much more in-depth using MAGE’s base functions
while the former makes use of wrapper functions that handle many
analyses and intermediary steps under-the hood. As such, the regular
vignette is recommended for first-time users or users that just want a
plug-and-play pipeline, while the expanded vignette provides more
insight into MAGE’s analyses and is recommended when setting out to
create your own custom/specialized analysis pipeline.

MAGE has an associated paper pending publication. Coming soon…
