
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MAGE

<!-- badges: start -->
<!-- badges: end -->

MAGE: Modeller of Allelic Gene Expression

## Installation

Put installation instructions here:

``` r
install.packages("MAGE")
```

## Examples

Link to vignette…

MAGE does stuff:

``` r
library(MAGE)
data("MAGE", package = "MAGE")

ChromPlot <- MAGE::MAGE_ADChromplot(AD_Data, DE_Data, Meth_Data, CNAgain_Data, CNAloss_Data,
                                    pvalSIG = 0.05, roll_median = 15)

ChromPlot[["ADDE_plot"]] / ChromPlot[["MethCNA_plot"]] / (ChromPlot[["LEG1"]] + 
  ChromPlot[["LEG2"]] + ChromPlot[["LEG3"]]) + 
  patchwork::plot_layout(heights = c(2,1,0.5), widths = c(1,1,1))
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" style="display: block; margin: auto;" />

Render `README.Rmd` regularly, to keep `README.md` up-to-date.
`devtools::build_readme()` is handy for this. You could also use GitHub
Actions to re-render `README.Rmd` every time you push. An example
workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots. In that case, don’t forget to commit and push
the resulting figure files, so they display on GitHub and CRAN.
