# TumPurHelpFun_Final <- function(RC, TC, TumReads_oi, a, b, q, SCP_oi, n = 1000){
#   
#   
#   logSums_MaxMethod <- function(logvec){
#     # EXPECTED: a vector of logs
#     logvec <- logvec[!is.na(logvec)]
#     X1 <- max(logvec)
#     H <- logvec - X1
#     return(X1 + log(sum(exp(H))))
#   }
#   
#   LogTrapezoidalInt <- function(f, lower, upper, n, ...){
#     Xs <- seq(from = lower, to = upper, length.out = n)
#     h <- (upper-lower)/(n-1)
#     Ys <- f(Xs, ...)
#     Ys[1] <- Ys[1] + log(5) - log(12)
#     Ys[n] <- Ys[n] + log(5) - log(12)
#     Ys[2] <- Ys[2] + log(13) - log(12)
#     Ys[n-1] <- Ys[n-1] + log(13) - log(12)
#     return(log(h)+logSums_MaxMethod(Ys))
#   }
#   
#   
#   
#   RC_mat <- matrix(rep(0:RC,length(TumReads_oi)), ncol=length(TumReads_oi), byrow=FALSE)
#   TumReads_mat <- matrix(rep(TumReads_oi, RC+1), ncol=length(TumReads_oi), byrow=TRUE)
#   RC_mat_C <- RC_mat; TumReads_mat_C <- TumReads_mat
#   
#   RC_mat[RC_mat_C > TumReads_mat_C] <- NA; RC_mat[RC-RC_mat_C > TC-TumReads_mat_C] <- NA
#   TumReads_mat[RC_mat_C > TumReads_mat_C] <- NA; TumReads_mat[RC-RC_mat_C > TC-TumReads_mat_C] <- NA
#   
#   
#   
#   ToInt <- function(Ppar_X, a,b,q,TC,RC,curTR,curRC){
#     
#     Tumprob <- qbeta(Ppar_X, q*a, q*b)
#     Hprob <- qbeta(Ppar_X, a, b)
#     
#     Res <- dbinom(x = curRC, size = curTR, prob = Tumprob, log = TRUE) +
#       dbinom(x = RC-curRC, size = TC-curTR, prob = Hprob, log = TRUE)
#     
#     return(Res)
#     
#   }
#   ToInt <- Vectorize(ToInt, vectorize.args = c("Ppar_X"))
#   
#   
#   
#   teller <- 0
#   probsvec <- rep(0, length(TumReads_oi))
#   for(curTR in TumReads_oi){
#     teller <- teller+1
#     RCcol <- RC_mat[,teller]
#     RCcol <- RCcol[!is.na(RCcol)]
#     
#     probsvecELS <- c()
#     for(RCcolEL in RCcol){
#       probsvecELS <- c(probsvecELS, LogTrapezoidalInt(ToInt, lower=0, upper=1, n = n, 
#                                                       a=a,b=b,q=q,TC=TC,RC=RC,curTR=curTR,curRC=RCcolEL))
#     }
#     
#     probsvec[teller] <- logSums_MaxMethod(probsvecELS)
#       
#   }
#   
#   
#   return(logSums_MaxMethod(log(SCP_oi)+probsvec))
#   
#   
# }
