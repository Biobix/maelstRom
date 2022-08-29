#' Plot MAGE's imprinting detection results.
#'
#' \code{MAGE_imprintplot} plots the results of MAGE's imprinting detection analysis. More specifically, it plots a histogram of the observed reference allele fraction,
#' als well as the PMF of the fitted imprinted model, which is a binomial mixture model using a heterozygous peak with no variant bias (p-parameter = 0.5) but that is
#' split up according to a fitted degree of imprinting (see \code{imprinting_est} for this fitting procedure and \code{pmf_impr} for the imprinted pmf). Optionally,
#' the non-imprinted fit (so just a binomial mixture model without variant bias) is plotted as well, if \code{plot_NoImpr} is TRUE.
#'
#' @param ref_counts Numeric vector. Reference counts of the locus.
#' @param var_counts Numeric vector. Variant counts of the locus.
#' @param allelefreq Number. Allele frequency.
#' @param impr Number. Degree of imprinting.
#' @param SE Number. Sequencing error rate.
#' @param inbr Number. Degree of inbreeding (default = 0).
#' @param PlotCov Number. Even though the observed number of reads can vary per sample, we have to assume a single number of observed reads to plot imprinting PMFs.
#' One can opt to pick a set value for this (the default is 50) or e.g. use the locus' median or mean coverage here, so as to best correspond to the observed data.
#' @param Plot_NoImpr Logical. If TRUE, also plot the non-imprintend PMF (a binomial mixture model without variant bias i.e. binomial p-parameter of heterozygotes = 0.5)
#' @param SplitPeaks Logical.  If TRUE, the PMFs of both homozygotes and the heterozygotes are plotted seperately and using different colors. If FALSE, one total imprinted PMF is plotted.
#' @param ImrCols String (vector). Colors for the plotted PMFs of the imprinted fit. Three values are required if \code{SplitPeaks} is TRUE, otherwise just one value is required.
#' Default colors are used in case none are given as input.
#' @param NoImprCol String. Color for the non-imprinted PMF (also no allelic bias; heterozygous binomial p-parameter = 0.5). There's no option to split the peaks of this PMF as it would make
#' the final figure way too clustered, so only one color value is required. A default color is used in case non is given as input.
#' @param wd_res String. Working directory where plots are saved; if non is given, the plot itself (as a \code{ggplot2} object) is returned instead.
#' @param position Number. Position of the locus, to be used in naming the outpud file if a \code{wd_res} is given.
#' @param chr Number. Chromosome of the locus, to be used in naming the outpud file if a \code{wd_res} is given.
#' @param gene String. Gene the locus is part of, to be used in naming the outpud file if a \code{wd_res} is given; optional.
#' @export

MAGE_imprintplot <- function (ref_counts, var_counts, allelefreq, impr, SE, inbr = 0, PlotCov = 50, plot_NoImpr = FALSE, SplitPeaks = FALSE, 
                              ImprCols = NULL, NoImprCol = NULL, wd_res = NULL, chr = "", position = "", gene = "") {
  
  if(SplitPeaks){
    if(is.null(ImprCols) | length(ImprCols)!= 3){
      ImprCols <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF")
    }
  }else{
    if(is.null(ImprCols) | length(ImprCols)!= 1){
      ImprCols <- "#FFB300"
    }
  }
  
  if(is.null(NoImprCol) | length(NoImprCol)!=1){
    NoImprCol <- "purple"
  }
  
  
  PlotCov <- ceiling(PlotCov)
  
  pr <- allelefreq ^ 2 + (inbr * allelefreq * (1 - allelefreq))
  pv <- (1 - allelefreq) ^ 2 + (inbr * allelefreq * (1 - allelefreq))
  prv <- allelefreq * (1 - allelefreq) * (1 - inbr)

  TheoryDat_impr <- prv*(dbinom(0:PlotCov, PlotCov, 0.5-0.5*impr)+dbinom(0:PlotCov, PlotCov, 0.5+0.5*impr)) * length(ref_counts)
  TheoryDat_impr_RR <- pr*(dbinom(0:PlotCov, PlotCov, 1-SE)) * length(ref_counts)
  TheoryDat_impr_VV <- pv*(dbinom(0:PlotCov, PlotCov, SE)) * length(ref_counts)
  
  Ratios <- ref_counts / (ref_counts + var_counts)
  
  #Ratios <- Ratios[Ratios < 0.1 | Ratios > 0.9]
  
  ImprDat <- data.frame("Ratio" = Ratios)
  
  ImprDat_Theory <- data.frame("DataRV" = TheoryDat_impr, "Pos" = (0:PlotCov)/PlotCov, "DataRR" = TheoryDat_impr_RR, "DataVV" = TheoryDat_impr_VV, 
                               "DataRRVV" = TheoryDat_impr_RR+TheoryDat_impr_VV, "DataALL" = TheoryDat_impr+TheoryDat_impr_RR+TheoryDat_impr_VV)
  

  PosV <- (0:PlotCov)/PlotCov
  AdvDF <- data.frame("Val" = c(TheoryDat_impr, TheoryDat_impr_RR, TheoryDat_impr_VV),
                      "Pos" = c(PosV, PosV, PosV),
                      "Geno" = as.factor(c(rep("RV", length(TheoryDat_impr)), rep("RR", length(TheoryDat_impr)), rep("VV", length(TheoryDat_impr)))),
                      "Hom" = as.factor(c(rep(0, length(TheoryDat_impr)), rep(1, 2*(length(TheoryDat_impr)))) ))
  
  
  if(!plot_NoImpr){
    
    if(SplitPeaks){
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = ImprDat, ggplot2::aes(Ratio), bins = PlotCov) +
        ggplot2::geom_line(data = AdvDF[AdvDF$Geno=="RR",], ggplot2::aes(x=Pos,y=Val, color= Geno, group = Geno), size = 1) +
        ggplot2::geom_line(data = AdvDF[AdvDF$Geno=="VV",], ggplot2::aes(x=Pos,y=Val, color= Geno, group = Geno), size = 1) +
        ggplot2::geom_line(data = AdvDF[AdvDF$Geno=="RV",], ggplot2::aes(x=Pos,y=Val, color= Geno, group = Geno), size = 1) +
        ggplot2::labs(x="Reference allele fraction", y="Frequency") +
        #ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left") +
        #ggplot2::theme(legend.margin=margin(2,12,2,2), legend.key.size = unit(0.7, 'cm'), legend.text = element_text(size=11)) +
        #ggplot2::scale_linetype_manual(name = "Genotype", values = c(1,2), guide = ggplot2::guide_legend(reverse = TRUE), label = c("Heterozygous fit, imprinted", "Homozygous fit")) + 
        ggplot2::scale_color_manual(guide = "none", values = ImprCols[c(1,3,2)]) +
        ggplot2::geom_area(data = AdvDF[AdvDF$Geno=="RR",], ggplot2::aes(x=Pos,y=Val), color = NA, fill = ImprCols[1], alpha = 0.45, size = 0, show.legend = FALSE) +
        ggplot2::geom_area(data = AdvDF[AdvDF$Geno=="VV",], ggplot2::aes(x=Pos,y=Val), color = NA, fill = ImprCols[2], alpha = 0.45, size = 0, show.legend = FALSE) +
        ggplot2::geom_area(data = AdvDF[AdvDF$Geno=="RV",], ggplot2::aes(x=Pos,y=Val), color = ImprCols[3], fill = ImprCols[3], alpha = 0.45, size = 0, show.legend = FALSE)
    }else{
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = ImprDat, ggplot2::aes(Ratio), bins = PlotCov) +
        ggplot2::geom_line(data = ImprDat_Theory, ggplot2::aes(x=Pos,y=DataALL, color= ImprCols), size = 1) +
        ggplot2::labs(x="Reference allele fraction", y="Frequency") +
        #ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left") +
        #ggplot2::theme(legend.margin=margin(2,12,2,2), legend.key.size = unit(0.7, 'cm'), legend.text = element_text(size=11)) +
        #ggplot2::scale_linetype_manual(name = "Genotype", values = c(1,2), guide = ggplot2::guide_legend(reverse = TRUE), label = c("Heterozygous fit, imprinted", "Homozygous fit")) + 
        ggplot2::scale_color_manual(guide = "none", values = ImprCols) +
        ggplot2::geom_area(data = ImprDat_Theory, ggplot2::aes(x=Pos,y=DataALL), color = NA, fill = ImprCols, alpha = 0.45, size = 0, show.legend = FALSE) 
    }
    
  }else{
    
    pr_H0 <- allelefreq ^ 2
    pv_H0 <- (1 - allelefreq) ^ 2
    prv_H0 <- allelefreq * (1 - allelefreq) * 2
    
    TheoryDat_impr_H0 <- prv_H0*dbinom(0:PlotCov, PlotCov, 0.5) * length(ref_counts)
    TheoryDat_impr_RR_H0 <- pr_H0*(dbinom(0:PlotCov, PlotCov, 1-SE)) * length(ref_counts)
    TheoryDat_impr_VV_H0 <- pv_H0*(dbinom(0:PlotCov, PlotCov, SE)) * length(ref_counts)
    
    ImprDat_Theory_H0 <- data.frame("Pos" = (0:PlotCov)/PlotCov, "DataALL" = TheoryDat_impr_H0+TheoryDat_impr_RR_H0+TheoryDat_impr_VV_H0)
    
    AdvDF_H0 <- data.frame("Val" = c(ImprDat_Theory_H0$DataALL, ImprDat_Theory_H0$DataALL),
                        "Pos" = c(PosV, PosV),
                        "Hypothesis" = c(rep("Imprinting", nrow(ImprDat_Theory_H0)), rep("No Imprinting", nrow(ImprDat_Theory_H0))))
    
    if(SplitPeaks){
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = ImprDat, ggplot2::aes(Ratio), bins = PlotCov) +
        ggplot2::geom_line(data = AdvDF[AdvDF$Geno=="RR",], ggplot2::aes(x=Pos,y=Val, linetype = "Imprinting"), color= ImprCols[1], size = 1) +
        ggplot2::geom_line(data = AdvDF[AdvDF$Geno=="VV",], ggplot2::aes(x=Pos,y=Val, linetype = "Imprinting"), color= ImprCols[2], size = 1) +
        ggplot2::geom_line(data = AdvDF[AdvDF$Geno=="RV",], ggplot2::aes(x=Pos,y=Val, linetype = "Imprinting"), color= ImprCols[3], size = 1) +
        ggplot2::labs(x="Reference allele fraction", y="Frequency") +
        #ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left") +
        #ggplot2::theme(legend.margin=margin(2,12,2,2), legend.key.size = unit(0.7, 'cm'), legend.text = element_text(size=11)) +
        #ggplot2::scale_linetype_manual(name = "Genotype", values = c(1,2), guide = ggplot2::guide_legend(reverse = TRUE), label = c("Heterozygous fit, imprinted", "Homozygous fit")) + 
        ggplot2::geom_area(data = AdvDF[AdvDF$Geno=="RR",], ggplot2::aes(x=Pos,y=Val), color = ImprCols[1], fill = ImprCols[1], alpha = 0.45, size = 0, show.legend = FALSE) +
        ggplot2::geom_area(data = AdvDF[AdvDF$Geno=="VV",], ggplot2::aes(x=Pos,y=Val), color = ImprCols[2], fill = ImprCols[2], alpha = 0.45, size = 0, show.legend = FALSE) +
        ggplot2::geom_area(data = AdvDF[AdvDF$Geno=="RV",], ggplot2::aes(x=Pos,y=Val), color = ImprCols[3], fill = ImprCols[3], alpha = 0.45, size = 0, show.legend = FALSE) +
        ggplot2::geom_line(data = ImprDat_Theory_H0, ggplot2::aes(x=Pos,y=DataALL, linetype = "No Imprinting", color= NoImprCol), size = 1) +
        ggplot2::scale_linetype_manual(name = "Hypothesis", values = c("solid", "dashed"), labels = c("Imprinting", "No Imprinting")) +
        ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left") +
        ggplot2::scale_color_identity(guide = "none")
    }else{
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = ImprDat, ggplot2::aes(Ratio), bins = PlotCov) +
        ggplot2::geom_line(data = ImprDat_Theory, ggplot2::aes(x=Pos,y=DataALL, linetype = "Imprinting"), color= ImprCols, size = 1) +
        ggplot2::labs(x="Reference allele fraction", y="Frequency") +
        #ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left") +
        #ggplot2::theme(legend.margin=margin(2,12,2,2), legend.key.size = unit(0.7, 'cm'), legend.text = element_text(size=11)) +
        #ggplot2::scale_linetype_manual(name = "Genotype", values = c(1,2), guide = ggplot2::guide_legend(reverse = TRUE), label = c("Heterozygous fit, imprinted", "Homozygous fit")) + 
        ggplot2::geom_area(data = ImprDat_Theory, ggplot2::aes(x=Pos,y=DataALL), color = ImprCols[1], fill = ImprCols, alpha = 0.45, size = 0, show.legend = FALSE) +
        ggplot2::geom_line(data = ImprDat_Theory_H0, ggplot2::aes(x=Pos,y=DataALL, linetype = "No Imprinting", color= NoImprCol), size = 1) +
        ggplot2::scale_linetype_manual(name = "Hypothesis", values = c("solid", "dashed"), labels = c("Imprinting", "No Imprinting")) +
        ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left") +
        ggplot2::scale_color_identity(guide = "none")
    }
    
  }
  
  if(!is.null(wd_res)){
    ggplot2::ggsave(filename = paste(wd_res, "chr", chr, "_", position, "_", gene, ".png", sep = ""), plot = pl, width = 10, height = 8, units = "in", dpi = 300)
  } else{
    return(pl)
  }
}
