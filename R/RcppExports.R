# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Exact beta-binomial density using a multiprecision library.
#'
#' \code{dBetaBinom_MP} calculates the beta-binomial density while
#' avoiding numerical mistakes (catastrophic cancellations) due to 
#' extreme parameter values. This function is called by \code{dBetaBinom}
#' if necessary, and should not be called outside of this.
#'
#' @param ms Numeric vector. Vector of number of successes
#' @param ns Numeric vector. Vector of number of trials
#' @param piX Number. Probability of success; \code{0 >= piX >= 1}
#' @param thetaX Number. Dispersion parameter; \code{0 >= thetaX > +Inf}
#' @param LOG Logical. if TRUE, return log-densities
#' @param NecPres Number. Necessary Precision, i.e. number of bits, for an accurate density calculatation, as determined by the function \code{dBetaBinom}
#' @return A numeric vector of the same length as ms and ns, containing (log-)beta-binomial densities
#' @export
dBetaBinom_MP <- function(ms, ns, piX, thetaX, LOG, NecPres) {
    .Call(`_maelstRom_dBetaBinom_MP`, ms, ns, piX, thetaX, LOG, NecPres)
}

#' Exact gradient of the beta-binomial log-likelihood function for pi using a multiprecision library.
#'
#' \code{grad_pi_MP} calculates the value of the gradient of the beta-binomial log-likelihood function to pi
#' at given data points, while avoiding numerical mistakes (catastrophic cancellations) due to 
#' extreme parameter values. This function is called by \code{grad_pi}
#' if necessary, and should not be called outside of this.
#'
#' @param ms Numeric vector. Vector of number of successes
#' @param ns Numeric vector. Vector of number of trials
#' @param piX Number. Probability of success; \code{0 >= piX >= 1}
#' @param thetaX Number. Dispersion parameter; \code{0 >= thetaX > +Inf}
#' @param NecPres Number. Necessary Precision, i.e. number of bits, for an accurate gradient calculatation, as determined by the function \code{grad_pi}
#' @return A numeric vector of the same length as ms and ns, containing the gradient to pi in the give data points
#' @export
grad_pi_MP <- function(ms, ns, piX, thetaX, NecPres) {
    .Call(`_maelstRom_grad_pi_MP`, ms, ns, piX, thetaX, NecPres)
}

#' Exact gradient of the beta-binomial log-likelihood function for theta using a multiprecision library.
#'
#' \code{grad_theta_MP} calculates the value of the gradient of the beta-binomial log-likelihood function to theta
#' at given data points, while avoiding numerical mistakes (catastrophic cancellations) due to 
#' extreme parameter values. This function is called by \code{grad_theta}
#' if necessary, and should not be called outside of this.
#'
#' @param ms Numeric vector. Vector of number of successes
#' @param ns Numeric vector. Vector of number of trials
#' @param piX Number. Probability of success; \code{0 >= piX >= 1}
#' @param thetaX Number. Dispersion parameter; \code{0 >= thetaX > +Inf}
#' @param NecPres Number. Necessary Precision, i.e. number of bits, for an accurate gradient calculatation, as determined by the function \code{grad_pi}
#' @return A numeric vector of the same length as ms and ns, containing the gradient to theta in the give data points
#' @export
grad_theta_MP <- function(ms, ns, piX, thetaX, NecPres) {
    .Call(`_maelstRom_grad_theta_MP`, ms, ns, piX, thetaX, NecPres)
}

#' Exact beta-binomial density using sums.
#'
#' \code{dBetaBinom_cpp_old} calculates the beta-binomial density via
#' a number of sums, which is slow for high-value data
#' but fast for low-value data.
#' This function is called by \code{dBetaBinom}
#' if necessary, and should not be called outside of this.
#'
#' @param ms Numeric vector. Vector of number of successes
#' @param ns Numeric vector. Vector of number of trials
#' @param pi Number. Probability of success; \code{0 >= pi >= 1}
#' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
#' @param LOG Logical. if TRUE, return log-densities
#' @return A numeric vector of the same length as ms and ns, containing (log-)beta-binomial densities
#' @export
dBetaBinom_cpp_old <- function(ms, ns, pi, theta, LOG) {
    .Call(`_maelstRom_dBetaBinom_cpp_old`, ms, ns, pi, theta, LOG)
}

#' Exact gradient of the beta-binomial log-likelihood function for pi using sums.
#'
#' \code{grad_pi_old} calculates the value of the gradient of the beta-binomial log-likelihood function to pi
#' at given data points via a number of sums, which is slow for high-value data
#' but fast for low-value data.
#' This function is called by \code{grad_pi} if necessary, and should not be called outside of this.
#'
#' @param ms Numeric vector. Vector of number of successes
#' @param ns Numeric vector. Vector of number of trials
#' @param pi Number. Probability of success; \code{0 >= pi >= 1}
#' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
#' @return A numeric vector of the same length as ms and ns, containing the gradient to pi in the give data points
#' @export
grad_pi_old <- function(ms, ns, pi, theta) {
    .Call(`_maelstRom_grad_pi_old`, ms, ns, pi, theta)
}

#' Exact gradient of the beta-binomial log-likelihood function for theta using sums.
#'
#' \code{grad_theta_old} calculates the value of the gradient of the beta-binomial log-likelihood function to theta
#' at given data points via a number of sums, which is slow for high-value data
#' but fast for low-value data.
#' This function is called by \code{grad_theta}
#' if necessary, and should not be called outside of this.
#'
#' @param ms Numeric vector. Vector of number of successes
#' @param ns Numeric vector. Vector of number of trials
#' @param pi Number. Probability of success; \code{0 >= pi >= 1}
#' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
#' @return A numeric vector of the same length as ms and ns, containing the gradient to theta in the give data points
#' @export
grad_theta_old <- function(ms, ns, pi, theta) {
    .Call(`_maelstRom_grad_theta_old`, ms, ns, pi, theta)
}

#' UNUSED test function for C++ driven numeric optimization
#' @export
MyOptTest <- function(parX) {
    .Call(`_maelstRom_MyOptTest`, parX)
}

#' dBetaBinom implementation for internal use by C++, accounting for numerical precision problems via increased memory usage (number of bits per number)
#' @export
dBetaBin_cppi_MP <- function(M, N, PI, THETA, NecPres) {
    .Call(`_maelstRom_dBetaBin_cppi_MP`, M, N, PI, THETA, NecPres)
}

#' dBetaBinom implementation for internal use by C++
#' @export
dBetaBin_cppi <- function(M, N, PI, THETA, LOG, MemLim, Xtra) {
    .Call(`_maelstRom_dBetaBin_cppi`, M, N, PI, THETA, LOG, MemLim, Xtra)
}

#' @rdname dpqrBetaBinom
#' @export
dBetaBinom <- function(ms, ns, pi, theta, LOG = FALSE, MemLim = 2048L, Xtra = 7L) {
    .Call(`_maelstRom_dBetaBinom`, ms, ns, pi, theta, LOG, MemLim, Xtra)
}

#' Beta-binomial gradient with respect to its pi parameter; for internal use by C++, accounting for numerical precision problems via increased memory usage (number of bits per number)
#' @export
GradPi_cppi_MP <- function(M, N, PI, THETA, NecPres) {
    .Call(`_maelstRom_GradPi_cppi_MP`, M, N, PI, THETA, NecPres)
}

#' Beta-binomial gradient with respect to its pi parameter; for internal use by C++
#' @export
GradPi_cppi <- function(M, N, PI, THETA, MemLim, Xtra) {
    .Call(`_maelstRom_GradPi_cppi`, M, N, PI, THETA, MemLim, Xtra)
}

#' Beta-binomial gradient with respect to its pi parameter
#' @export
grad_pi <- function(ms, ns, pi, theta, MemLim = 2048L, Xtra = 7L) {
    .Call(`_maelstRom_grad_pi`, ms, ns, pi, theta, MemLim, Xtra)
}

#' Beta-binomial gradient with respect to its theta parameter; for internal use by C++, accounting for numerical precision problems via increased memory usage (number of bits per number)
#' @export
GradTheta_cppi_MP <- function(M, N, PI, THETA, NecPres) {
    .Call(`_maelstRom_GradTheta_cppi_MP`, M, N, PI, THETA, NecPres)
}

#' Beta-binomial gradient with respect to its theta parameter; for internal use by C++
#' @export
GradTheta_cppi <- function(M, N, PI, THETA, MemLim, Xtra) {
    .Call(`_maelstRom_GradTheta_cppi`, M, N, PI, THETA, MemLim, Xtra)
}

#' Beta-binomial gradient with respect to its theta parameter
#' @export
grad_theta <- function(ms, ns, pi, theta, MemLim = 2048L, Xtra = 7L) {
    .Call(`_maelstRom_grad_theta`, ms, ns, pi, theta, MemLim, Xtra)
}

#' Optimizer for pi, theta_control, and theta_case, for heterozygous samples in the differential control-case beta-binomial mixture fit using expectation-maximization
#' @export
CppCnT_Optim <- function(StartVals, ref_counts, var_counts, isCase, sprv, MemLim = 2048L, Xtra = 7L, step_size = 0.01, tol = 0.1, epsabs = 1e-3) {
    .Call(`_maelstRom_CppCnT_Optim`, StartVals, ref_counts, var_counts, isCase, sprv, MemLim, Xtra, step_size, tol, epsabs)
}

#' Optimizer for theta for homozygous samples in the beta-binomial mixture fit using expectation-maximization
#' @export
CppHom_Optim <- function(ThetaHomStart, SE, ref_counts, var_counts, spr, spv, MemLim = 2048L, Xtra = 7L, step_size = 0.01, tol = 0.1, epsabs = 1e-3) {
    .Call(`_maelstRom_CppHom_Optim`, ThetaHomStart, SE, ref_counts, var_counts, spr, spv, MemLim, Xtra, step_size, tol, epsabs)
}

#' Optimizer for pi, theta for heterozygous samples in the beta-binomial mixture fit using expectation-maximization
#' @export
CppHet_Optim <- function(StartVals, ref_counts, var_counts, sprv, MemLim = 2048L, Xtra = 7L, step_size = 0.01, tol = 0.1, epsabs = 1e-3) {
    .Call(`_maelstRom_CppHet_Optim`, StartVals, ref_counts, var_counts, sprv, MemLim, Xtra, step_size, tol, epsabs)
}

#' Optimizer for theta for heterozygous samples assuming a fixed pi parameter in the beta-binomial mixture fit using expectation-maximization
#' @export
CppHetH0_Optim <- function(ThetaHetStart, probshift, ref_counts, var_counts, sprv, MemLim = 2048L, Xtra = 7L, step_size = 0.01, tol = 0.1, epsabs = 1e-3) {
    .Call(`_maelstRom_CppHetH0_Optim`, ThetaHetStart, probshift, ref_counts, var_counts, sprv, MemLim, Xtra, step_size, tol, epsabs)
}

#' Beta-binomial second-order derivative with respect to its pi parameter (twice); for internal use by C++, accounting for numerical precision problems via increased memory usage (number of bits per number)
#' @export
GradPiPi_cppi_MP <- function(M, N, PI, THETA, NecPres) {
    .Call(`_maelstRom_GradPiPi_cppi_MP`, M, N, PI, THETA, NecPres)
}

#' Beta-binomial second-order derivative with respect to its pi parameter (twice); for internal use by C++
#' @export
GradPiPi_cppi <- function(M, N, PI, THETA, MemLim, Xtra) {
    .Call(`_maelstRom_GradPiPi_cppi`, M, N, PI, THETA, MemLim, Xtra)
}

#' Beta-binomial second-order derivative with respect to its pi parameter (twice)
#' @export
grad_pi_pi <- function(ms, ns, pi, theta, MemLim = 2048L, Xtra = 7L) {
    .Call(`_maelstRom_grad_pi_pi`, ms, ns, pi, theta, MemLim, Xtra)
}

#' Beta-binomial second-order derivative with respect to its pi parameter (once) and theta-parameter (once); for internal use by C++, accounting for numerical precision problems via increased memory usage (number of bits per number)
#' @export
GradPiTheta_cppi_MP <- function(M, N, PI, THETA, NecPres) {
    .Call(`_maelstRom_GradPiTheta_cppi_MP`, M, N, PI, THETA, NecPres)
}

#' Beta-binomial second-order derivative with respect to its pi parameter (once) and theta-parameter (once); for internal use by C++
#' @export
GradPiTheta_cppi <- function(M, N, PI, THETA, MemLim, Xtra) {
    .Call(`_maelstRom_GradPiTheta_cppi`, M, N, PI, THETA, MemLim, Xtra)
}

#' Beta-binomial second-order derivative with respect to its pi parameter (once) and theta-parameter (once)
#' @export
grad_pi_theta <- function(ms, ns, pi, theta, MemLim = 2048L, Xtra = 7L) {
    .Call(`_maelstRom_grad_pi_theta`, ms, ns, pi, theta, MemLim, Xtra)
}

#' Helper function for GradThetaTheta_cppi
#' @export
MPhelper_GradThetaTheta <- function(M, N, PI, NecPres) {
    .Call(`_maelstRom_MPhelper_GradThetaTheta`, M, N, PI, NecPres)
}

#' Beta-binomial second-order derivative with respect to its theta-parameter (twice); for internal use by C++, accounting for numerical precision problems via increased memory usage (number of bits per number)
#' @export
GradThetaTheta_cppi_MP <- function(M, N, PI, THETA, NecPres) {
    .Call(`_maelstRom_GradThetaTheta_cppi_MP`, M, N, PI, THETA, NecPres)
}

#' Beta-binomial second-order derivative with respect to its theta-parameter (twice); for internal use by C++
#' @export
GradThetaTheta_cppi <- function(M, N, PI, THETA, MemLim, Xtra) {
    .Call(`_maelstRom_GradThetaTheta_cppi`, M, N, PI, THETA, MemLim, Xtra)
}

#' Beta-binomial second-order derivative with respect to its theta-parameter (twice)
#' @export
grad_theta_theta <- function(ms, ns, pi, theta, MemLim = 2048L, Xtra = 7L) {
    .Call(`_maelstRom_grad_theta_theta`, ms, ns, pi, theta, MemLim, Xtra)
}

#' Returns an (approximate) log of a sum based on individual logs of the terms being summed
#'
#' \code{logSums_MaxMethod_CPP} is an internal C++ implementation of logSum_MaxMethod, and
#' calculates the (approximate) log of a sum of terms, based on the individual logs of these terms.
#' Simply taking the exponent of these individual terms, then summing them, then taking the log again, is usually not possible if taking the exponent of the
#' supplied values results in a numeric overflow. The rationale here is:
#' 
#' log(x + y + z) =(suppose x is the MAXIMUM amongst {x,y,z}) log(((x+y+z)/x) * x) = log(x) + log((x+y+z)/x) = log(x) + log(x/x + y/x + z/x) =
#' log(x) + log(exp(log(x/x)) + exp(log(y/x)) + exp(log(z/x))) = log(x) + log(exp(log(x)-log(x)) + exp(log(y)-log(x)) + exp(log(z)-log(x)))
#' 
#' Now, log(x)-log(x) cancels out, and the remaining (log(y)-log(x)) and (log(z)-log(x)) are negative (x is the maximum amontst {x,y,z}) meaning that,
#' worst case scenario, x is that much larger than y and z that these terms are high and negative, resulting in a calculated exp-value of zero;
#' this isn't too bad though, as the first term (log(x)) is still there, and if x is so much higher than the other terms, then the log of sums will be
#' mainly dictated by x anyway (when remaining on R's default numerical precision, and given a limited number of terms being summed; once again, this value is approximate).
#'
#' @param logvec Numeric vector. separate log values, of which we want to calculate the log of the sum of their (non-log) values
#' @export
#' @return The log of the sum of the values for which the separate logs were supplied as input.
logSums_MaxMethod_CPP <- function(logvec) {
    .Call(`_maelstRom_logSums_MaxMethod_CPP`, logvec)
}

#' Returns an (approximate) log of a sum based on individual logs of the terms being summed
#'
#' \code{logSums_MaxMethod_CPP_R} is a C++ implementation of logSum_MaxMethod callable from R, and
#' calculates the (approximate) log of a sum of terms, based on the individual logs of these terms.
#' Simply taking the exponent of these individual terms, then summing them, then taking the log again, is usually not possible if taking the exponent of the
#' supplied values results in a numeric overflow. The rationale here is:
#' 
#' log(x + y + z) =(suppose x is the MAXIMUM amongst {x,y,z}) log(((x+y+z)/x) * x) = log(x) + log((x+y+z)/x) = log(x) + log(x/x + y/x + z/x) =
#' log(x) + log(exp(log(x/x)) + exp(log(y/x)) + exp(log(z/x))) = log(x) + log(exp(log(x)-log(x)) + exp(log(y)-log(x)) + exp(log(z)-log(x)))
#' 
#' Now, log(x)-log(x) cancels out, and the remaining (log(y)-log(x)) and (log(z)-log(x)) are negative (x is the maximum amontst {x,y,z}) meaning that,
#' worst case scenario, x is that much larger than y and z that these terms are high and negative, resulting in a calculated exp-value of zero;
#' this isn't too bad though, as the first term (log(x)) is still there, and if x is so much higher than the other terms, then the log of sums will be
#' mainly dictated by x anyway (when remaining on R's default numerical precision, and given a limited number of terms being summed; once again, this value is approximate).
#'
#' @param logvec Numeric vector. separate log values, of which we want to calculate the log of the sum of their (non-log) values
#' @export
#' @return The log of the sum of the values for which the separate logs were supplied as input.
logSums_MaxMethod_CPP_R <- function(logvec) {
    .Call(`_maelstRom_logSums_MaxMethod_CPP_R`, logvec)
}

#' Returns an (approximate) log of a sum based on individual logs of the terms being summed, but can include negative terms (see further)
#'
#' \code{logSums_MaxMethod_CPP} is an internal C++ implementation of logSum_MaxMethod, and
#' calculates the (approximate) log of a sum of terms, based on the individual logs of these terms.
#' Simply taking the exponent of these individual terms, then summing them, then taking the log again, is usually not possible if taking the exponent of the
#' supplied values results in a numeric overflow. The rationale here is:
#' 
#' log(x + y + z) =(suppose x is the MAXIMUM amongst {x,y,z}) log(((x+y+z)/x) * x) = log(x) + log((x+y+z)/x) = log(x) + log(x/x + y/x + z/x) =
#' log(x) + log(exp(log(x/x)) + exp(log(y/x)) + exp(log(z/x))) = log(x) + log(exp(log(x)-log(x)) + exp(log(y)-log(x)) + exp(log(z)-log(x)))
#' 
#' Now, log(x)-log(x) cancels out, and the remaining (log(y)-log(x)) and (log(z)-log(x)) are negative (x is the maximum amontst {x,y,z}) meaning that,
#' worst case scenario, x is that much larger than y and z that these terms are high and negative, resulting in a calculated exp-value of zero;
#' this isn't too bad though, as the first term (log(x)) is still there, and if x is so much higher than the other terms, then the log of sums will be
#' mainly dictated by x anyway (when remaining on R's default numerical precision, and given a limited number of terms being summed; once again, this value is approximate).
#' 
#' In this implementation, negative terms are actually allowed, but the sign of terms should be supplied as a separate input "signvec"
#' If negative, a term's separate contribution to the final big log-term in the formula above (e.g. exp(log(y)-log(x)) for term y)
#' is subtracted instead of added. This is okay as long as the final result of which the log is eventually taken still ends up being positive.
#' Accordingly, this method should only be used if this is actually expected.
#' For example, this function is used by maelstRom's implementation for calculating the probability distrubution of a sum
#' of two completely correlated beta-binomials using the  Newton-Cotes numerical integration method,
#' for which we expect positive outcomes for maelstRom's applications (integration of probability distributions). 
#'
#' @param logvec Numeric vector. separate log values, of which we want to calculate the log of the sum of their (non-log) values
#' @param signvec Numeric vector. separate signs of the values being summed/subtracted (as +1 or -1), as their log-values provided in logvec can,
#' of course, not contain this information.
#' @export
#' @return The log of the sum of the values for which the separate logs were supplied as input.
logSums_MaxMethodSigned_CPP <- function(logvec, signvec) {
    .Call(`_maelstRom_logSums_MaxMethodSigned_CPP`, logvec, signvec)
}

#' Internal helper function for TumPur_LogLik_CPP when NumIntMethod == "Gregory"
#' @export
LogTrapezoidalInt_CPP <- function(lower, upper, n, a, b, q, TC, RC, curTR, curRC) {
    .Call(`_maelstRom_LogTrapezoidalInt_CPP`, lower, upper, n, a, b, q, TC, RC, curTR, curRC)
}

#' Internal helper function for TumPur_LogLik_CPP when NumIntMethod == "NewtonCotes"
#' @export
LogNewtonCotes_CPP <- function(lower, upper, a, b, q, TC, RC, curTR, curRC, Wvec) {
    .Call(`_maelstRom_LogNewtonCotes_CPP`, lower, upper, a, b, q, TC, RC, curTR, curRC, Wvec)
}

#' Internal helper function for TumPur_LogLik_CPP when NumIntMethod == "GaussianQuad"
#' @export
LogGaussianQuad_CPP <- function(lower, upper, a, b, q, TC, RC, curTR, curRC, Wvec, Nvec) {
    .Call(`_maelstRom_LogGaussianQuad_CPP`, lower, upper, a, b, q, TC, RC, curTR, curRC, Wvec, Nvec)
}

#' Internal helper function for TumPur_LogLik_CPP when NumIntMethod == "TanhSinhQuad"
#' @export
LogTanhSinhQuad_CPP <- function(lower, upper, n, a, b, q, TC, RC, curTR, curRC, prec) {
    .Call(`_maelstRom_LogTanhSinhQuad_CPP`, lower, upper, n, a, b, q, TC, RC, curTR, curRC, prec)
}

#' Internal helper function for TumPur_LogLik_CPP
#' @export
TumPurHelpFun_CPP <- function(RCcol, TC, TumReads_oi, a, b, q, SCP_oi, n, NumIntMethod, prec, Wvec, Nvec) {
    .Call(`_maelstRom_TumPurHelpFun_CPP`, RCcol, TC, TumReads_oi, a, b, q, SCP_oi, n, NumIntMethod, prec, Wvec, Nvec)
}

#' Internal helper function for TumPur_LogLik_CPP
#' @export
TumPurHelpFun_CPP_R <- function(RCcol, TC, TumReads_oi, a, b, q, SCP_oi, n, NumIntMethod, prec, Wvec, Nvec) {
    .Call(`_maelstRom_TumPurHelpFun_CPP_R`, RCcol, TC, TumReads_oi, a, b, q, SCP_oi, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' UNUSED alternative implementation of the qbeta funtion
#' @export
qbeta_C <- function(qs, a, b) {
    .Call(`_maelstRom_qbeta_C`, qs, a, b)
}

#' UNUSED alternative implementation of the qbeta funtion
#' @export
qbeta_C3 <- function(qs, a, b) {
    .Call(`_maelstRom_qbeta_C3`, qs, a, b)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP2 <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP2`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP2X <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP2X`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP3 <- function(q, a, b, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP3`, q, a, b, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP3X <- function(q, a, b, qlim, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP3X`, q, a, b, qlim, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP4 <- function(q, a, b, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP4`, q, a, b, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP_DB1 <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP_DB1`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP_DB2 <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP_DB2`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP_DB3 <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP_DB3`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' UNUSED
#' @export
BrolDB <- function(TC, TP) {
    .Call(`_maelstRom_BrolDB`, TC, TP)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP2_10 <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP2_10`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP2_100 <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP2_100`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP2_1000 <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP2_1000`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

#' @rdname TumPur_LogLik_CPPs
#' @export
TumPur_LogLik_CPP2Y <- function(optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n = 0L, NumIntMethod = "Gregory", prec = 0.0001, Wvec = 0L, Nvec = 0L) {
    .Call(`_maelstRom_TumPur_LogLik_CPP2Y`, optpars, ref_counts, var_counts, tumpur, weights, SCPthreshold, n, NumIntMethod, prec, Wvec, Nvec)
}

