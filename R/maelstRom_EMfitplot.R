#' Plot maelstRom's EM-fit results.
#'
#' \code{maelstRom_EMfitplot} plots the results of maelstRom's beta-binomial EM-fit. More specifically, it plots both the observed reference allele fractions as a histogram (optionally scaled
#' to a single total count parameter, see \code{ScaleHist} and \code{ScaleCount}) and the parameters yielded by the EM-fit as beta-binomial PMFs (assuming a total read counts \code{ScaleCount}
#' and compressed to the zero-to-one interval of the reference allele fraction instead of integer reference counts). Optionally, the unshifted fit assuming no allelic bias (heterozygous
#' pi-parameter = 0.5) is plotted as well, if \code{plot_NoShift} is TRUE and its distributional parameters are provided as input.
#'
#' @param ref_counts Numeric vector. Reference counts.
#' @param var_counts Numeric vector. Variant counts.
#' @param pr Number. Reference homozygote genotype probability of the locus.
#' @param pv Number. Variant homozygote genotype probability of the locus.
#' @param prv Number. Heterozygote genotype probability of the locus.
#' @param theta_hom Number. The dispersion parameter of the homozygous peaks.
#' @param theta_het Number. The dispersion parameter of the heterozygous peak.
#' @param pr_NoShift Number. Reference homozygote genotype probability of the locus under the null hypothesis fit (no allelic bias; heterozygous pi-parameter = 0.5).
#' @param pv_NoShift Number. Variant homozygote genotype probability of the locus under the null hypothesis fit (no allelic bias; heterozygous pi-parameter = 0.5).
#' @param prv_NoShift Number. Heterozygote genotype probability of the locus under the null hypothesis fit (no allelic bias; heterozygous pi-parameter = 0.5).
#' @param theta_hom_NoShift Number. The dispersion parameter of the homozygous peaks under the null hypothesis fit (no allelic bias; heterozygous pi-parameter = 0.5).
#' @param theta_het_NoShift Number. The dispersion parameter of the heterozygous peak under the null hypothesis fit (no allelic bias; heterozygous pi-parameter = 0.5).
#' @param probshift Number. The reference allele fraction in heterozygotes, indicating allelic bias when deviating from 0.5
#' @param SE Number. Sequencing error rate.
#' @param MinCount Number. A minimal count filter for plotting only; samples with a lower count will not be included in the plotted histogram.
#' @param ScaleCount Number. Even though the observed number of reads can vary per sample, we have to assume a single number of observed reads to plot PMFs. One can opt to pick
#' a set value for this (the default is 50) or e.g. use the locus' median or mean coverage here, so as to best correspond to the observed data.
#' @param ScaleHist Logical. See the remark given at \code{ScaleCount}; if TRUE, tranforms the input data to all have the same total allele counts \code{ScaleCount}, with the actual number
#' of reference reads being determined as the quantile function assuming \code{ScaleCount} total counts with the same probability as the observed data (this probability being calculated
#' with the actually observed number of counts). The idea here is to better reflect how well the observed data corresponds to the fitted PMF, which is plotted for a single number of observed reads.
#' @param nbins Number. Number of bins for the histogram depicting the observed reference allele fractions.
#' @param plot_NoShift Logical. If TRUE, also plots the PMF of the unshifted fit (no allelic bias; heterozygous pi-parameter = 0.5)
#' @param SplitPeaks Logical. If TRUE, the PMFs of both homozygotes and the heterozygotes are plotted seperately and using different colors. If FALSE, one total PMF of the mixture model is plotted.
#' @param ShiftCols String (vector). Colors for the plotted PMFs of the fit incorporating allelic bias. Three values are required if \code{SplitPeaks} is TRUE, otherwise just one value is required.
#' Default colors are used in case none are given as input.
#' @param NoShiftCol String. Color for the unshifted PMF (no allelic bias; heterozygous pi-parameter = 0.5). There's no option to split the peaks of this unshifted fit as it would make
#' the final figure way too clustered, so only one color value is required. A default color is used in case non is given as input.
#' @param wd_res String. Working directory where plots are saved; if non is given, the plot itself (as a \code{ggplot2} object) is returned instead.
#' @param position Number. Position of the locus, to be used in naming the output file if a \code{wd_res} is given.
#' @param chr Number. Chromosome of the locus, to be used in naming the output file if a \code{wd_res} is given.
#' @param gene String. Gene the locus is part of, to be used in naming the output file if a \code{wd_res} is given; optional.
#' @param DataList_out Dataframe. Alternatively, \code{maelstRom_EMfitplot} accepts a dataframes containing per-sample data of the locus of interest. From this, it can
#' infer its \code{ref_counts}, \code{var_counts}, \code{pr_NoShift}, \code{prv_NoShift} and \code{pv_NoShift} arguments. This is possible
#' as long as \code{EMfit_betabinom_robust} was run on this dataframe with its \code{fitH0} argument set to TRUE.
#' @param Geno_AB_res Dataframe with 1 row. Alternatively, \code{maelstRom_EMfitplot} accepts a dataframe containing beta-binomial mixture model EM-fit results from calling \code{EMfit_betabinom_robust}'s
#' wrapper function \code{BetaBinomGenotyping}. As such, expected columns are listed on \code{BetaBinomGenotyping}'s help page, more specifically
#' its produced \code{Geno_AB_res} output. The one row provided must correspond to the locus of interest. From this, \code{maelstRom_EMfitplot} can
#' infer its \code{pr}, \code{prv}, \code{pv}, \code{theta_hom}, \code{theta_het}, \code{theta_hom_NoShift}, \code{theta_het_NoShift} and \code{probshift} arguments.
#' @param dAD_res Dataframe with 1 row. Alternatively, \code{maelstRom_EMfitplot} accepts a dataframe containing beta-binomial mixture model EM-fit results from calling \code{EMfit_betabinom_popcomb}'s
#' wrapper function \code{dAD_analysis}. As such, expected columns are listed on \code{dAD_analysis}'s help page, more specifically
#' its produced \code{dAD_res} output. The one row provided must correspond to the locus of interest. From this, \code{maelstRom_EMfitplot} can
#' infer its \code{pr}, \code{prv}, \code{pv}, \code{theta_hom}, \code{theta_het}, and \code{probshift} arguments. The "NoShift" case can not be plotted in this scenario,
#' as dAD-detection happens without AB-detection. If \code{dAD_res} is not NULL, \code{Geno_AB_res} has to be NULL and vice-versa. Also, when using \code{dAD_res},
#' you need to specify whether to plot controls or cases using the \code{PlotWhich} argument.
#' @param PlotWhich String. Either "control" or "case", to know which data to fetch from \code{dAD_res} for plotting.
#' The dataframe given as the \code{DataList_out} argument should of course be of the corresponding population.
#' @export

maelstRom_EMfitplot <- function(ref_counts = NULL, var_counts = NULL, pr = NULL, prv = NULL, pv = NULL, theta_hom = NULL, theta_het = NULL, pr_NoShift = NULL, prv_NoShift = NULL, pv_NoShift = NULL, 
                                theta_hom_NoShift = NULL, theta_het_NoShift = NULL, probshift, SE, MinCount = 0, ScaleCount = 50, ScaleHist = FALSE, nbins = 100, 
                                plot_NoShift = FALSE, SplitPeaks = TRUE, ShiftCols = NULL, NoShiftCol = NULL, wd_res = NULL, chr = "", position = "", gene = "",
                                DataList_out = NULL, Geno_AB_res = NULL, dAD_res = NULL, PlotWhich = "choose", BG = "#F2F2F2") {
  
  if(!is.null(DataList_out)){
    ref_counts<-DataList_out$ref_count
    var_counts<-DataList_out$var_count
    pr_NoShift<-mean(DataList_out$prr_H0)
    prv_NoShift<-mean(DataList_out$prv_H0)
    pv_NoShift<-mean(DataList_out$pvv_H0)
  }
  if(!is.null(Geno_AB_res)){
    pr<-Geno_AB_res$rho_rr
    prv<-Geno_AB_res$rho_rv
    pv<-Geno_AB_res$rho_vv
    theta_hom<-Geno_AB_res$theta_hom
    theta_het<-Geno_AB_res$theta_het
    theta_hom_NoShift<-Geno_AB_res$theta_hom_NoShift
    theta_het_NoShift<-Geno_AB_res$theta_het_NoShift
    probshift<-as.numeric(Geno_AB_res$probshift)
  }
  if(!is.null(dAD_res)){
    pr<-dAD_res$pr
    prv<-dAD_res$prv
    pv<-dAD_res$pv
    theta_hom<-dAD_res$ThetaHomH0
    if(PlotWhich == "control"){
      theta_het<-dAD_res$ThetaHetCTRL
    }else if(PlotWhich == "case"){
      theta_het<-dAD_res$ThetaHetCASE
    }else{
      stop("PlotWhich has to be either \"control\" or \"case\"")
    }
    probshift<-as.numeric(dAD_res$PiFitH1)
    plot_NoShift <- FALSE
  }
  
  if(SplitPeaks){
    if(is.null(ShiftCols) | length(ShiftCols)!= 3){
      ShiftCols <- c("#3C6CC2", "#FFC000", "#00B050")
    }
  }else{
    if(is.null(ShiftCols) | length(ShiftCols)!= 1){
      ShiftCols <- "#FFB300"
    }
  }
  
  if(is.null(NoShiftCol) | length(NoShiftCol)!=1){
    NoShiftCol <- "purple"
  }
  
  ScaleCount <- ceiling(ScaleCount)
  # First: filtering!
  RC <- ref_counts[ref_counts+var_counts >= MinCount]
  VC <- var_counts[ref_counts+var_counts >= MinCount]
  
  if(length(RC)==0 | length(VC)==0){
    stop("MinCount too high; no data remains")
  }
  
  if(ScaleHist){
    # Next: calculating weights!
    spr <- pr * maelstRom::dBetaBinom(RC, RC + VC, pi = 1-SE, theta = theta_hom, LOG = FALSE)
    spv <- pv * maelstRom::dBetaBinom(VC, RC + VC, pi = 1-SE, theta = theta_hom, LOG = FALSE)
    sprv <- prv * maelstRom::dBetaBinom(RC, RC + VC, pi = probshift, theta = theta_het, LOG = FALSE)
    pdata <- rowSums(cbind(spv, sprv, spr)) 
    if (any(pdata==0)){ 
      ProblemCases <- which(pdata==0)
      for(case in ProblemCases){
        var_part<-VC[case]
        ref_part<-RC[case]
        spr_part<-maelstRom::dBetaBinom(ref_part, ref_part+var_part, pi = 1-SE, theta = theta_hom, LOG = TRUE)
        spv_part<-maelstRom::dBetaBinom(var_part, ref_part+var_part, pi = 1-SE, theta = theta_hom, LOG = TRUE)
        sprv_part<-maelstRom::dBetaBinom(ref_part, ref_part+var_part, pi = probshift, theta = theta_het, LOG = TRUE)
        spvec<-c(spr_part, sprv_part, spv_part)
        if(spr_part==max(spvec)){
          spr[case]<-1
          pdata[case]<-1
        } else if(sprv_part==max(spvec)){
          sprv[case]<-1
          pdata[case]<-1
        } else{
          spv[case]<-1
          pdata[case]<-1
        }
      }
    }
    spv <- spv/pdata
    sprv <- sprv/pdata
    spr <- spr/pdata
    # Now: for each observation, get the observation most closely corresponding to this observation,
    # should its total count have been ScaleCount
    Helper <- maelstRom::pBetaBinom(ms = RC, ns = RC+VC, pi = 1-SE, theta = theta_hom)
    HomRR_ref <- maelstRom::qBetaBinom(Helper, ns = rep(ScaleCount, length(RC)), pi = 1-SE, theta = theta_hom)
    HomRR_var <- ScaleCount-HomRR_ref
    Helper <- maelstRom::pBetaBinom(ms = RC, ns = RC+VC, pi = probshift, theta = theta_het)
    Het_ref <- maelstRom::qBetaBinom(Helper, ns = rep(ScaleCount, length(RC)), pi = probshift, theta = theta_het)
    Het_var <- ScaleCount-Het_ref
    Helper <- maelstRom::pBetaBinom(ms = VC, ns = RC+VC, pi = 1-SE, theta = theta_hom)
    HomVV_var <- maelstRom::qBetaBinom(Helper, ns = rep(ScaleCount, length(RC)), pi = 1-SE, theta = theta_hom)
    HomVV_ref <- ScaleCount-HomVV_var
    
    
    
    
    Dat_TheoryData <- (pr * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE) +
                         pv * maelstRom::dBetaBinom(ScaleCount:0, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE) +
                         prv * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = probshift, theta = theta_het, LOG = FALSE)) * length(RC)
    Dat_RRData <- (pr * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE)) * length(RC)
    Dat_VVData <- (pv * maelstRom::dBetaBinom(ScaleCount:0, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE)) * length(RC)
    Dat_RVData <- (prv * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = probshift, theta = theta_het, LOG = FALSE)) * length(RC)
    
    MDX <- data.frame("Count" = 0:ScaleCount, "Data" = Dat_TheoryData, "DataRR" = Dat_RRData, "DataVV" = Dat_VVData, "DataRV" = Dat_RVData)
    MDX$St <- 0:ScaleCount/ScaleCount
    MDX$St <- MDX$St - (1/ScaleCount)/2
    MDX$St[1] <- 0
    MDX$En <- 0:ScaleCount/ScaleCount
    MDX$En <- MDX$En + (1/ScaleCount)/2
    MDX$En[nrow(MDX)] <- 1
    
    ScaleCountV <- (0:ScaleCount)/ScaleCount
    HV <- seq((1/nbins)/2,1,(1/nbins))
    HV2 <- seq(0,1,(1/nbins))
    HV2X <- HV2[HV2 != 1]
    HV2Y <- HV2[HV2 != 0]
    
    MD <- data.frame("Binstart" = HV2X, "Binend" = HV2Y, "Data" = 0, "Pos" = HV, "DataRR" = 0, "DataVV" = 0, "DataRV" = 0)
    MD$Pos[1]<- 0
    MD$Pos[nrow(MD)]<- 1
    
    for(i2 in 1:nrow(MD)){
      
      Alls <- which(MDX$St < MD$Binend[i2] & !(MDX$St < MD$Binstart[i2] & MDX$En <= MD$Binstart[i2]))
      
      for(AA in Alls){
        
        over <- min(MD$Binend[i2], MDX$En[AA]) - max(MD$Binstart[i2], MDX$St[AA])
        over <- over/(MDX$En[AA]-MDX$St[AA])
        
        
        MD$Data[i2] <- MD$Data[i2] + over*Dat_TheoryData[AA]
        MD$DataRR[i2] <- MD$DataRR[i2] + over*Dat_RRData[AA]
        MD$DataVV[i2] <- MD$DataVV[i2] + over*Dat_VVData[AA]
        MD$DataRV[i2] <- MD$DataRV[i2] + over*Dat_RVData[AA]
        
      }
      
    }
    
    Dat_Theory <- MD
    
    if(plot_NoShift){
      Dat_TheoryData_H0 <- (pr_NoShift * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom_NoShift, LOG = FALSE) +
                              pv_NoShift * maelstRom::dBetaBinom(ScaleCount:0, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom_NoShift, LOG = FALSE) +
                              prv_NoShift * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 0.5, theta = theta_het_NoShift, LOG = FALSE)) * length(RC)
      
      Dat_Theory_Ho <- data.frame("Data" = Dat_TheoryData_H0, "Pos" = (0:ScaleCount)/ScaleCount) 
    }
    
    
    
    
    
    
    
    Dat <- data.frame("Ratios" = c(HomRR_ref/(HomRR_ref+HomRR_var), Het_ref/(Het_ref+Het_var),
                                   HomVV_ref/(HomVV_ref+HomVV_var)), "Weights" = c(spr, sprv, spv))
    if(SplitPeaks){
      
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = Dat, ggplot2::aes(Ratios, weight = Weights), breaks = seq(from = 0, to = 1, by = (1/(nbins))),
                                                        fill = "#7F7F7F") + 
        ggplot2::labs(x="reference allele fraction", y="frequency") +
        ggplot2::geom_vline(xintercept=probshift, linetype="dotted", col = "#00B050", size = 1) +
        ggplot2::geom_vline(xintercept=0.5, linetype="dotted", col = "black", size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRR), color= ShiftCols[1], fill = ShiftCols[1], alpha = 0.45, size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataVV), color= ShiftCols[2], fill = ShiftCols[2], alpha = 0.45, size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRV), color= ShiftCols[3], fill = ShiftCols[3], alpha = 0.45, size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRR, linetype = "shifted"), color= ShiftCols[1], size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataVV, linetype = "shifted"), color= ShiftCols[2], size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRV, linetype = "shifted"), color= ShiftCols[3], size = 1) +
        ggplot2::theme(legend.position="none")
      if(plot_NoShift){
        pl <- pl + ggplot2::geom_line(data = Dat_Theory_H0, ggplot2::aes(x = Pos, y = Data, linetype = "unshifted", color = NoShiftCol), size = 1) +
          ggplot2::scale_linetype_manual(name = "Fit type", values = c("solid", "dashed"), labels = c("shifted", "unshifted")) +
          ggplot2::scale_color_identity(guide = "none") +
          ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left")
      }
    }else{
      
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = Dat, ggplot2::aes(Ratios, weight = Weights), breaks = seq(from = 0, to = 1, by = (1/(nbins))),
                                                        fill = "#7F7F7F") + 
        ggplot2::labs(x="reference allele fraction", y="frequency") +
        ggplot2::geom_vline(xintercept=probshift, linetype="dotted", col = "#00B050", size = 1) +
        ggplot2::geom_vline(xintercept=0.5, linetype="dotted", col = "black", size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=Data), color= ShiftCols, fill = ShiftCols, alpha = 0.45, size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=Data, linetype = "shifted"), color= ShiftCols, size = 1) +
        ggplot2::theme(legend.position="none")
      if(plot_NoShift){
        pl <- pl + ggplot2::geom_line(data = Dat_Theory_H0, ggplot2::aes(x = Pos, y = Data, linetype = "unshifted", color = NoShiftCol), size = 1) +
          ggplot2::scale_linetype_manual(name = "Fit type", values = c("solid", "dashed"), labels = c("shifted", "unshifted")) +
          ggplot2::scale_color_identity(guide = "none") +
          ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left")
      }
    }
    
  }else{
    
    Dat_TheoryData <- (pr * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE) +
                         pv * maelstRom::dBetaBinom(ScaleCount:0, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE) +
                         prv * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = probshift, theta = theta_het, LOG = FALSE)) * length(RC)
    Dat_RRData <- (pr * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE)) * length(RC)
    Dat_VVData <- (pv * maelstRom::dBetaBinom(ScaleCount:0, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom, LOG = FALSE)) * length(RC)
    Dat_RVData <- (prv * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = probshift, theta = theta_het, LOG = FALSE)) * length(RC)
    
    MDX <- data.frame("Count" = 0:ScaleCount, "Data" = Dat_TheoryData, "DataRR" = Dat_RRData, "DataVV" = Dat_VVData, "DataRV" = Dat_RVData)
    MDX$St <- 0:ScaleCount/ScaleCount
    MDX$St <- MDX$St - (1/ScaleCount)/2
    MDX$St[1] <- 0
    MDX$En <- 0:ScaleCount/ScaleCount
    MDX$En <- MDX$En + (1/ScaleCount)/2
    MDX$En[nrow(MDX)] <- 1
    
    ScaleCountV <- (0:ScaleCount)/ScaleCount
    HV <- seq((1/nbins)/2,1,(1/nbins))
    HV2 <- seq(0,1,(1/nbins))
    HV2X <- HV2[HV2 != 1]
    HV2Y <- HV2[HV2 != 0]
    
    MD <- data.frame("Binstart" = HV2X, "Binend" = HV2Y, "Data" = 0, "Pos" = HV, "DataRR" = 0, "DataVV" = 0, "DataRV" = 0)
    MD$Pos[1]<- 0
    MD$Pos[nrow(MD)]<- 1
    
    for(i2 in 1:nrow(MD)){
      
      Alls <- which(MDX$St < MD$Binend[i2] & !(MDX$St < MD$Binstart[i2] & MDX$En <= MD$Binstart[i2]))
      
      for(AA in Alls){
        
        over <- min(MD$Binend[i2], MDX$En[AA]) - max(MD$Binstart[i2], MDX$St[AA])
        over <- over/(MDX$En[AA]-MDX$St[AA])
        
        
        MD$Data[i2] <- MD$Data[i2] + over*Dat_TheoryData[AA]
        MD$DataRR[i2] <- MD$DataRR[i2] + over*Dat_RRData[AA]
        MD$DataVV[i2] <- MD$DataVV[i2] + over*Dat_VVData[AA]
        MD$DataRV[i2] <- MD$DataRV[i2] + over*Dat_RVData[AA]
        
      }
      
    }
    
    Dat_Theory <- MD
    
    if(plot_NoShift){
      Dat_TheoryData_H0 <- (pr_NoShift * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom_NoShift, LOG = FALSE) +
                              pv_NoShift * maelstRom::dBetaBinom(ScaleCount:0, rep(ScaleCount, ScaleCount+1), pi = 1-SE, theta = theta_hom_NoShift, LOG = FALSE) +
                              prv_NoShift * maelstRom::dBetaBinom(0:ScaleCount, rep(ScaleCount, ScaleCount+1), pi = 0.5, theta = theta_het_NoShift, LOG = FALSE)) * length(RC)
      
      Dat_Theory_Ho <- data.frame("Data" = Dat_TheoryData_H0, "Pos" = (0:ScaleCount)/ScaleCount) 
    }
    
    Dat <- data.frame("Ratios" = RC / (RC+VC))
    
    if(SplitPeaks){
      
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = Dat, ggplot2::aes(Ratios), breaks = seq(from = 0, to = 1, by = (1/(nbins))),
                                                        fill = "#7F7F7F") + 
        ggplot2::labs(x="reference allele fraction", y="frequency") +
        ggplot2::geom_vline(xintercept=probshift, linetype="dotted", col = "#00B050", size = 1) +
        ggplot2::geom_vline(xintercept=0.5, linetype="dotted", col = "black", size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRR), color= ShiftCols[1], fill = ShiftCols[1], alpha = 0.45, size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataVV), color= ShiftCols[2], fill = ShiftCols[2], alpha = 0.45, size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRV), color= ShiftCols[3], fill = ShiftCols[3], alpha = 0.45, size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRR, linetype = "shifted"), color= ShiftCols[1], size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataVV, linetype = "shifted"), color= ShiftCols[2], size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=DataRV, linetype = "shifted"), color= ShiftCols[3], size = 1) +
        
        ggplot2::theme(legend.position="none")
      
      if(plot_NoShift){
        pl <- pl + ggplot2::geom_line(data = Dat_Theory_H0, ggplot2::aes(x = Pos, y = Data, linetype = "unshifted", color = NoShiftCol), size = 1) +
          ggplot2::scale_linetype_manual(name = "Fit type", values = c("solid", "dashed"), labels = c("shifted", "unshifted")) +
          ggplot2::scale_color_identity(guide = "none") +
          ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left")
        
      }
      
    }else{
      pl <- ggplot2::ggplot() + ggplot2::geom_histogram(data = Dat, ggplot2::aes(Ratios), breaks = seq(from = 0, to = 1, by = (1/(nbins))),
                                                        fill = "#7F7F7F") + 
        ggplot2::labs(x="reference allele fraction", y="frequency") +
        ggplot2::geom_vline(xintercept=probshift, linetype="dotted", col = "#00B050", size = 1) +
        ggplot2::geom_vline(xintercept=0.5, linetype="dotted", col = "black", size = 1) +
        ggplot2::geom_area(data = Dat_Theory, ggplot2::aes(x=Pos,y=Data), color= ShiftCols, fill = ShiftCols, alpha = 0.45, size = 1) +
        ggplot2::geom_line(data = Dat_Theory, ggplot2::aes(x=Pos,y=Data, linetype = "shifted"), color= ShiftCols, size = 1) +
        ggplot2::theme(legend.position="none")
      if(plot_NoShift){
        pl <- pl + ggplot2::geom_line(data = Dat_Theory_H0, ggplot2::aes(x = Pos, y = Data, linetype = "unshifted", color = NoShiftCol), size = 1) +
          ggplot2::scale_linetype_manual(name = "Fit type", values = c("solid", "dashed"), labels = c("shifted", "unshifted")) +
          ggplot2::scale_color_identity(guide = "none") +
          ggplot2::theme(legend.position = c(.1, .955), legend.justification = c("left", "top"), legend.box.just = "left")
      }
    }
  }
  
  pl <- pl + ggplot2::theme(plot.background = ggplot2::element_rect(fill =BG, color = BG)) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = BG, colour = BG))
  
  if(!is.null(wd_res)){
    ggplot2::ggsave(filename = paste(wd_res, "eQTL_chr", chr, "_", position, "_", gene, "BetaBinom.png", sep = ""), plot = pl, width = 10, height = 8, units = "in", dpi = 300)
  } else{
    return(pl)
  }
  
}
