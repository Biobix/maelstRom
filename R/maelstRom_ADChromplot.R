#' Make a plot summarizing allelic divergence (-related phenomena) across a chromosome
#'
#' \code{maelstRom_ADChromplot} returns various ggplot-objects visualizing maelstRom's differential allelic divergence results across a chromosome,
#' together with differential expression logFC which can be determined using the same RNAseq data, as well as (optionally) methylation- and copy-number-alteration-data if available.
#' All of this data should be supplied as separate dataframes, which enables plotting even when some of this data is missing for some of the loci or genes.
#'
#' @param AD_Data Dataframe. Should at least contain chromosomal position ("pos"), overdispersion in controls and cases expressed as rho which lies between
#' 0 and 1 ("rho_control", "rho_case") and the p-value testing for differential allelic divergence (difference in overdispersion parameter between controls and cases; "pval")
#' @param DE_Data Dataframe. Should at least contain chromosomal position ("pos") and the log(2)FC of expression between cases and controls ("logFC")
#' @param Meth_Data Dataframe. Optional; should at least contain chromosomal position ("pos") and the p-value testing for case hypermethylation ("logFC")
#' @param CNAgain_Data Dataframe. Optional; should at least contain chromosomal position ("pos") and the average gain in cases ("CNAgain")
#' @param CNAloss_Data Dataframe. Optional; should at least contain chromosomal position ("pos") and the average loss in cases ("CNAloss")
#' @param allelefreq Number. Allele frequency.
#' @param pvalSIG Number. Significance cutoff of provided p-values (these should be FDR-corrected in advance)
#' @param roll_median integer. All measures are plotted as rolling medians across \code{roll_median} data points.
#' @export
#' @return A list containing the following components:
#' \item{ADDE_plot}{A plot vizualizing Allelic Divergence and Differential Expression results.}
#' \item{MethCNA_plot}{A plot vizualizing hypermethylation and copy-number-alteration results.}
#' \item{LEG1, LEG2, LEG3}{Various legend items, to be combined with \code{ADDE_plot} and \code{MethCNA_plot} if the user so desires, e.g. using the \code{patchwork} package.}

maelstRom_ADChromplot <- function(AD_Data, DE_Data, Meth_Data=NULL, CNAgain_Data=NULL, CNAloss_Data=NULL, pvalSIG = 0.05, roll_median = 15){
  
  WO1 <- 0
  WO2 <- 0
  WOX <- 0
  
  if(is.null(Meth_Data)){
    WOX <- 1
    #Meth_Data <- data.frame("gene"=NA, "pos"=NA, "logp_meth"=NA, "MOI_RM"=NA)
  }
  if(is.null(CNAgain_Data)){
    WO1 <- 1
    #CNAgain_Data <- data.frame("gene"=NA, "pos"=0, "CNAgain"=-100, "MOI_RM"=-100)
  }
  if(is.null(CNAloss_Data)){
    WO2 <- 1
    #CNAloss_Data <- data.frame("gene"=NA, "pos"=0, "CNAloss"=-100, "MOI_RM"=-100)
  }
  
  ColLim <- c(0,13)
  
  LS <- 0.8 # LineSize
  
  AD_Data$RollMedRhoT <- zoo::rollmedian(AD_Data[,colnames(AD_Data)%in%c("pos","rho_case")],roll_median,fill=NA)[,2]
  AD_Data$RollMedRhoC <- zoo::rollmedian(AD_Data[,colnames(AD_Data)%in%c("pos","rho_control")],roll_median,fill=NA)[,2]
  
  DE_Data$MOI_RM <- zoo::rollmedian(DE_Data[, c("pos","logFC")][,c(2,1)],roll_median,fill=NA)[,1]
  Meth_Data$MOI_RM <- zoo::rollmedian(Meth_Data[, c("pos","logp_meth")][,c(2,1)],roll_median,fill=NA)[,1]
  CNAgain_Data$MOI_RM <- zoo::rollmedian(CNAgain_Data[, c("pos","CNAgain")][,c(2,1)],roll_median,fill=NA)[,1]
  CNAloss_Data$MOI_RM <- zoo::rollmedian(CNAloss_Data[, c("pos","CNAloss")][,c(2,1)],roll_median,fill=NA)[,1]
  
  
  AD_Data$ST <- NA
  AD_Data$EN <- NA
  AD_Data$FILL_meanminlogpval <- NA
  AD_Data$FILL_Sig <- NA
  MH <- (roll_median - 1)/2
  for(i in (MH+1):(nrow(AD_Data)-MH)){
    AD_Data$ST[i] <- mean(c(AD_Data$pos[i-1], AD_Data$pos[i]))
    AD_Data$EN[i] <- mean(c(AD_Data$pos[i+1], AD_Data$pos[i]))
    AD_Data$FILL_meanminlogpval[i] <- -mean(log(AD_Data$pval[(i-MH):(i+MH)]))
  }
  
  AD_Data$FILL_meanminlogpval_clone <- AD_Data$FILL_meanminlogpval
  AD_Data$FILL_meanminlogpval_clone[AD_Data$FILL_meanminlogpval > ColLim[2]] <- ColLim[2]
  AD_Data$FILL_meanminlogpval_clone[AD_Data$FILL_meanminlogpval < ColLim[1]] <- ColLim[1]
  
  MidP <- -log(pvalSIG)
  
  
  # PLOTS
  
  if(WOX==1){
    Meth_Data <- data.frame("gene"=NA, "pos"=0, "logp_meth"=-100, "MOI_RM"=-100)
  }
  if(WO1==1){
    CNAgain_Data <- data.frame("gene"=NA, "pos"=0, "CNAgain"=-100, "MOI_RM"=-100)
  }
  if(WO2==1){
    CNAloss_Data <- data.frame("gene"=NA, "pos"=0, "CNAloss"=-100, "MOI_RM"=-100)
  }
  
  P2 <- ggplot2::ggplot() + ggplot2::geom_rect(data = AD_Data, ggplot2::aes(xmin = ST, xmax = EN, ymax=RollMedRhoT, ymin=0, fill = FILL_meanminlogpval_clone)) +
    ggplot2::scale_fill_gradientn(colors=c("#0033FF","#D8E0FF","#FFD8D8", "#FF0000"),
                                  values=scales::rescale(c(ColLim[1], MidP-0.00001, MidP+0.00001,ColLim[2])),
                                  limits=ColLim, guide = ggplot2::guide_colorbar(ticks.colour = NA),
                                  breaks = (ColLim[1]+ColLim[2])/2, labels = c("diseased")) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_area(data = AD_Data, ggplot2::aes(x = pos, y = RollMedRhoC, fill = FILL_meanminlogpval_clone)) +
    ggplot2::scale_fill_gradientn(colors=c("#000000","#000000","#000000", "#000000"),
                                  values=scales::rescale(c(ColLim[1], MidP-0.00001, MidP+0.00001,ColLim[2])),
                                  limits=ColLim, guide = ggplot2::guide_colorbar(ticks.colour = NA),
                                  breaks = (ColLim[1]+ColLim[2])/2, labels = c("control")) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                   axis.text.x=ggplot2::element_blank()) +
    ggplot2::ylab("Overdispersion (Rho)") +
    ggplot2::coord_cartesian(ylim=c(0,0.4), xlim = c(0, max(AD_Data$EN))) +
    ggplot2::theme(legend.key.height = ggplot2::unit(1.2, "mm"), legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::theme(legend.position = c(.75, .955), legend.justification = c("left", "top"), legend.box.just = "left",
                   legend.title = ggplot2::element_blank(), legend.margin = ggplot2::margin(5,5,5,5), legend.box.background = ggplot2::element_rect(colour = "white") ) +
    ggplot2::theme(legend.spacing.y = ggplot2::unit(-0.045, "mm"))
  
  P6 <- ggplot2::ggplot() + ggplot2::geom_line(data = CNAgain_Data, ggplot2::aes(x = pos, y = MOI_RM, linetype = "gain"), size = LS) +
    ggplot2::geom_line(data = CNAloss_Data, ggplot2::aes(x = pos, y = MOI_RM, linetype = "loss"), size = LS) +
    ggplot2::ylab("Nr gain/loss") + ggplot2::scale_y_continuous(breaks = c(0,0.5,1)) + ggplot2::coord_cartesian(ylim=c(0,1), xlim = c(0, max(AD_Data$EN))) +
    ggplot2::xlab("Chromosomal Location") +
    ggplot2::scale_linetype_manual(name = NULL, values = c("solid", "dashed")) +
    ggplot2::theme(legend.position = c(.75, .955), legend.justification = c("left", "top"), legend.box.just = "left") +
    ggplot2::theme(legend.margin=ggplot2::margin(-1,5,5,5), legend.key.width = ggplot2::unit(0.8, 'cm'),
                   legend.key.height = ggplot2::unit(0.5, "cm"), legend.text = ggplot2::element_text(size=12)) +
    ggplot2::theme(legend.key = ggplot2::element_rect(colour = '#bdbdbd', size = 0.6)) +
    ggplot2::guides(linetype = ggplot2::guide_legend(nrow = 1))
  
  P_A <- P2 + ggplot2::geom_line(data = DE_Data, ggplot2::aes(x = pos, y = (MOI_RM+1.5)*(0.4/3)), color = "#8B5723", size = LS) +
    ggplot2::geom_hline(yintercept=(0+1.5)*(0.4/3), linetype = "dotted", color = "#8B5723", size = LS) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~(.*3/0.4)-1.5, name="Log2(FC)", breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5))) +
    ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "#8B5723"),
                   axis.title.y.right = ggplot2::element_text(color = "#8B5723"),
                   axis.ticks.y.right = ggplot2::element_line(color = "#8B5723"),
                   axis.line.y.right = ggplot2::element_line(color = "#8B5723")) +
    ggplot2::theme(legend.position = "none")
  
  P_B <- P6 + ggplot2::geom_line(data = Meth_Data, ggplot2::aes(x = pos, y = (MOI_RM/log10(exp(1))/(12))), color = "#E64B35FF", size = LS) +
    ggplot2::geom_hline(yintercept=(-log(0.05))/(12), linetype = "dotted", color = "#E64B35FF", size = LS) +
    ggplot2::scale_y_continuous(sec.axis = ggplot2::sec_axis(~(.*(12)), name="-Log(p_meth)", breaks = c(0,3,6,9,12))) +
    ggplot2::theme(axis.text.y.right = ggplot2::element_text(color = "#E64B35FF"),
                   axis.title.y.right = ggplot2::element_text(color = "#E64B35FF"),
                   axis.ticks.y.right = ggplot2::element_line(color = "#E64B35FF"),
                   axis.line.y.right = ggplot2::element_line(color = "#E64B35FF")) +
    ggplot2::theme(legend.position = "none")
  
  
  # LEGENDS
  
  LEG1 <- ggplot2::ggplot() + ggplot2::geom_rect(data = AD_Data, ggplot2::aes(xmin = ST, xmax = EN, ymax=RollMedRhoT, ymin=0, fill = FILL_meanminlogpval_clone)) +
    ggplot2::scale_fill_gradientn(colors=c("#0033FF","#D8E0FF","#FFD8D8", "#FF0000"),
                         values=scales::rescale(c(ColLim[1], MidP-0.00001, MidP+0.00001,ColLim[2])),
                         limits=ColLim, labels = c("0", "3", "6", "9", ">12"),
                         breaks = c(0,3,6,9,12)) +
    ggplot2::geom_area(data = AD_Data, ggplot2::aes(x = pos, y = RollMedRhoC), fill = "black") +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank()) +
    ggplot2::ylab("Overdispersion (Rho)") +
    ggplot2::coord_cartesian(ylim=c(0,0.4)) +
    ggplot2::theme(legend.position="bottom", legend.box = "horizontal", legend.key.width = ggplot2::unit(10, "mm"), legend.text = ggplot2::element_text(size = 12),
          legend.key.height = ggplot2::unit(7, "mm"), legend.title = ggplot2::element_text()) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(title.position = "top", title.hjust = 0.5)) + ggplot2::labs(fill = "Mean -Log(p_dAD)") +
    ggplot2::theme(legend.position=c(0.5,0.5), legend.direction = "horizontal")
  
  LEG2 <- ggplot2::ggplot() + ggplot2::geom_rect(data = AD_Data, ggplot2::aes(xmin = ST, xmax = EN, ymax=RollMedRhoT, ymin=0, fill = FILL_meanminlogpval_clone)) +
    ggplot2::scale_fill_gradientn(colors=c("#0033FF","#D8E0FF","#FFD8D8", "#FF0000"),
                         values=scales::rescale(c(ColLim[1], MidP-0.00001, MidP+0.00001,ColLim[2])),
                         limits=ColLim, guide = ggplot2::guide_colorbar(ticks.colour = NA, label.position = "right"),
                         breaks = (ColLim[1]+ColLim[2])/2, labels = expression(overdispersion[case])) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_area(data = AD_Data, ggplot2::aes(x = pos, y = RollMedRhoC, fill = FILL_meanminlogpval_clone)) +
    ggplot2::scale_fill_gradientn(colors=c("#000000","#000000","#000000", "#000000"),
                         values=scales::rescale(c(ColLim[1], MidP-0.00001, MidP+0.00001,ColLim[2])),
                         limits=ColLim, guide = ggplot2::guide_colorbar(ticks.colour = NA),
                         breaks = (ColLim[1]+ColLim[2])/2, labels = expression(overdispersion[control])) +
    ggplot2::theme(axis.title.x=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank()) +
    ggplot2::ylab("Overdispersion (Rho)") +
    ggplot2::coord_cartesian(ylim=c(0,0.4)) +
    ggplot2::theme(legend.key.height = ggplot2::unit(1.2, "mm"), legend.text = ggplot2::element_text(size = 12)) +
    ggplot2::theme(legend.box.just = "left",
          legend.title = ggplot2::element_blank(), legend.margin = ggplot2::margin(5,5,5,5), legend.box.background = ggplot2::element_rect(colour = "white") ) +
    ggplot2::theme(legend.spacing.y = ggplot2::unit(-0.045, "mm"))

  LEG2 <- LEG2 + ggplot2::theme(legend.position=c(0.5,0.5))

  LEG3 <- ggplot2::ggplot() + ggplot2::geom_line(data = NULL, ggplot2::aes(x = 1, y = -100, linetype = "mean gain", color = "mean gain"), size = 1) +
    ggplot2::geom_line(data = NULL, ggplot2::aes(x = 1, y = -100, linetype = "mean loss", color = "mean loss"), size = 1) +
    ggplot2::geom_line(data = NULL, ggplot2::aes(x = 1, y = -100, linetype = "bli", color = "bli"), size = 1) +
    ggplot2::geom_line(data = NULL, ggplot2::aes(x = 1, y = -100, linetype = "bla", color = "bla"), size = 1) +
    ggplot2::ylab("Nr gain/loss") + ggplot2::scale_y_continuous(breaks = c(0,0.5,1)) + ggplot2::coord_cartesian(ylim=c(0,1)) +
    ggplot2::xlab("Chromosomal Location") +
    ggplot2::scale_linetype_manual(name = NULL, values = c("solid", "solid", "solid", "dashed"), labels = expression(log2(FC[case/control]), -log(pval[methylation]), "mean gain", "mean loss")) +
    ggplot2::scale_color_manual(name = NULL, values = c("#8B5723", "#E64B35FF", "black", "black"), labels = expression(log2(FC[case/control]), -log(pval[methylation]), "mean gain", "mean loss")) +
    ggplot2::theme(legend.key.width = ggplot2::unit(0.8, 'cm'),
          legend.key.height = ggplot2::unit(0.5, "cm"), legend.text = ggplot2::element_text(size=12)) +
    ggplot2::theme(legend.key = ggplot2::element_rect(colour = '#bdbdbd', size = 0.6)) +
    ggplot2::guides(linetype = ggplot2::guide_legend(nrow = 2), color = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::theme(legend.text.align = 0) +
    ggplot2::theme(legend.position=c(0.5,0.5))
  
  LEG1 <- LEG1 + ggplot2::theme(aspect.ratio = 10^-16, legend.margin = ggplot2::margin(5,0,0,0), axis.line=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), panel.background=ggplot2::element_blank(), panel.border=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), plot.background=ggplot2::element_blank())
  
  LEG2 <- LEG2 + ggplot2::theme(aspect.ratio = 10^-16, legend.margin = ggplot2::margin(5,5,5,5), axis.line=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), panel.background=ggplot2::element_blank(), panel.border=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), plot.background=ggplot2::element_blank())
  
  LEG3 <- LEG3 + ggplot2::theme(aspect.ratio = 10^-16, line = ggplot2::element_blank(), legend.margin = ggplot2::margin(5,0,0,0), axis.line=ggplot2::element_blank(), axis.text.x=ggplot2::element_blank(), axis.text.y=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank(), axis.title.x=ggplot2::element_blank(), axis.title.y=ggplot2::element_blank(), panel.background=ggplot2::element_blank(), panel.border=ggplot2::element_blank(), panel.grid.major=ggplot2::element_blank(), panel.grid.minor=ggplot2::element_blank(), plot.background=ggplot2::element_blank())

  return(list("ADDE_plot"=P_A, "MethCNA_plot"=P_B, "LEG1"=LEG1, "LEG2"=LEG2, "LEG3"=LEG3))
}
