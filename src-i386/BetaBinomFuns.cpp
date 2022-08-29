
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <gmp.h>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/digamma.hpp>
using namespace Rcpp;

//' Exact beta-binomial density using a multiprecision library.
//'
//' \code{dBetaBinom_MP} calculates the beta-binomial density while
//' avoiding numerical mistakes (catastrophic cancellations) due to 
//' extreme parameter values. This function is called by \code{dBetaBinom}
//' if necessary, and should not be called outside of this.
//'
//' @param ms Numeric vector. Vector of number of successes
//' @param ns Numeric vector. Vector of number of trials
//' @param piX Number. Probability of success; \code{0 >= piX >= 1}
//' @param thetaX Number. Dispersion parameter; \code{0 >= thetaX > +Inf}
//' @param LOG Logical. if TRUE, return log-densities
//' @param NecPres Number. Necessary Precision, i.e. number of bits, for an accurate density calculatation, as determined by the function \code{dBetaBinom}
//' @return A numeric vector of the same length as ms and ns, containing (log-)beta-binomial densities
//' @export
// [[Rcpp::export]]
NumericVector dBetaBinom_MP(NumericVector ms, NumericVector ns, double piX, double thetaX, bool LOG, double NecPres) {
  
  boost::multiprecision::mpfr_float::default_precision(NecPres);
  
  int l = ms.size();
  NumericVector out(l);
  
  for(int i = 0; i < l; ++i) {
    boost::multiprecision::mpfr_float m = ms[i];
    boost::multiprecision::mpfr_float n = ns[i];
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PV3 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = piX;
    boost::multiprecision::mpfr_float theta = thetaX;
    
    if(m>0){
      PV1 = -lgamma(m)-lgamma(pi/theta)+lgamma(m+pi/theta)-log(m);
    }
    if((n-m)>0){
      PV2 = -lgamma(n-m)-lgamma((1-pi)/theta)+lgamma((n-m)+((1-pi)/theta))-log(n-m);
    }
    if(n>1){
      PV3 = lgamma(n)+lgamma(1/theta)-lgamma(n+(1/theta))+log(n);
    }
    PVtot = PV1+PV2+PV3;
    out[i] = PVtot.convert_to<double>();
    
  }
  if (LOG) {
    return out;
  } else {
    return exp(out);
  }
  
}


//' Exact gradient of the beta-binomial log-likelihood function for pi using a multiprecision library.
//'
//' \code{grad_pi_MP} calculates the value of the gradient of the beta-binomial log-likelihood function to pi
//' at given data points, while avoiding numerical mistakes (catastrophic cancellations) due to 
//' extreme parameter values. This function is called by \code{grad_pi}
//' if necessary, and should not be called outside of this.
//'
//' @param ms Numeric vector. Vector of number of successes
//' @param ns Numeric vector. Vector of number of trials
//' @param piX Number. Probability of success; \code{0 >= piX >= 1}
//' @param thetaX Number. Dispersion parameter; \code{0 >= thetaX > +Inf}
//' @param NecPres Number. Necessary Precision, i.e. number of bits, for an accurate gradient calculatation, as determined by the function \code{grad_pi}
//' @return A numeric vector of the same length as ms and ns, containing the gradient to pi in the give data points
//' @export
// [[Rcpp::export]]
NumericVector grad_pi_MP(NumericVector ms, NumericVector ns, double piX, double thetaX, double NecPres) {
  
  boost::multiprecision::mpfr_float::default_precision(NecPres);
  
  int l = ms.size();
  NumericVector out(l);
  
  for(int i = 0; i < l; ++i) {
    boost::multiprecision::mpfr_float m = ms[i];
    boost::multiprecision::mpfr_float n = ns[i];
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = piX;
    boost::multiprecision::mpfr_float theta = thetaX;
    
    if(m>0){
      PV1 = (1/theta) * (boost::math::digamma(m+(pi/theta)) - boost::math::digamma(pi/theta));
    }
    if((n-m)>0){
      PV2 = (1/theta) * (boost::math::digamma((n-m)+((1-pi)/theta)) - boost::math::digamma((1-pi)/theta));
    }
    PVtot = PV1-PV2;
    out[i] = PVtot.convert_to<double>();
  }
  return out;
}


//' Exact gradient of the beta-binomial log-likelihood function for theta using a multiprecision library.
//'
//' \code{grad_theta_MP} calculates the value of the gradient of the beta-binomial log-likelihood function to theta
//' at given data points, while avoiding numerical mistakes (catastrophic cancellations) due to 
//' extreme parameter values. This function is called by \code{grad_theta}
//' if necessary, and should not be called outside of this.
//'
//' @param ms Numeric vector. Vector of number of successes
//' @param ns Numeric vector. Vector of number of trials
//' @param piX Number. Probability of success; \code{0 >= piX >= 1}
//' @param thetaX Number. Dispersion parameter; \code{0 >= thetaX > +Inf}
//' @param NecPres Number. Necessary Precision, i.e. number of bits, for an accurate gradient calculatation, as determined by the function \code{grad_pi}
//' @return A numeric vector of the same length as ms and ns, containing the gradient to theta in the give data points
//' @export
// [[Rcpp::export]]
NumericVector grad_theta_MP(NumericVector ms, NumericVector ns, double piX, double thetaX, double NecPres) {
  
  boost::multiprecision::mpfr_float::default_precision(NecPres);
  
  int l = ms.size();
  NumericVector out(l);
  
  for(int i = 0; i < l; ++i) {
    boost::multiprecision::mpfr_float m = ms[i];
    boost::multiprecision::mpfr_float n = ns[i];
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PV3 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = piX;
    boost::multiprecision::mpfr_float theta = thetaX;
    
    if(m>0){
      PV1 = -(pi/theta/theta)*(boost::math::digamma(m+(pi/theta)) - boost::math::digamma(pi/theta));
    }
    if((n-m)>0){
      PV2 = -((1-pi)/theta/theta)*(boost::math::digamma((n-m)+((1-pi)/theta)) - boost::math::digamma((1-pi)/theta));
    }
    if(n>1){
      PV3 = -(1/theta/theta)*(boost::math::digamma(n+(1/theta)) - boost::math::digamma(1/theta));
    }
    PVtot = PV1+PV2-PV3;
    out[i] = PVtot.convert_to<double>();
  }
  return out;
}


//' Exact beta-binomial density using sums.
//'
//' \code{dBetaBinom_cpp_old} calculates the beta-binomial density via
//' a number of sums, which is slow for high-value data
//' but fast for low-value data.
//' This function is called by \code{dBetaBinom}
//' if necessary, and should not be called outside of this.
//'
//' @param ms Numeric vector. Vector of number of successes
//' @param ns Numeric vector. Vector of number of trials
//' @param pi Number. Probability of success; \code{0 >= pi >= 1}
//' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
//' @param LOG Logical. if TRUE, return log-densities
//' @return A numeric vector of the same length as ms and ns, containing (log-)beta-binomial densities
//' @export
// [[Rcpp::export]]
NumericVector dBetaBinom_cpp_old(NumericVector ms, NumericVector ns, double pi, double theta, bool LOG) {
  int l = ms.size();
  NumericVector out(l);
  
  for(int i = 0; i < l; ++i) {
    
    int m = ms[i];
    int n = ns[i];
    double PV1 = 0;
    double PV2 = 0;
    double PV3 = 0;
    
    if(m>0){
      for(double k = 0; k < (m); ++k) {
        PV1 += log(pi + k*theta) + log(1/(k+1));
      }
    }
    if((n-m) > 0){
      for(double k = 0; k < (n-m); ++k) {
        PV2 += log(1 - pi + k*theta) + log(1/(k+1));
      }  
    }
    if((n) > 1){
      for(double k = 1; k < (n); ++k) {
        PV3 += log(1/(1 + k*theta)) + log(k+1);
      }  
    }
    out[i] = PV1 + PV2 + PV3;
  }
  
  if (LOG) {
    return out;
  } else {
    return exp(out);
  }
  
}


//' Exact gradient of the beta-binomial log-likelihood function for pi using sums.
//'
//' \code{grad_pi_old} calculates the value of the gradient of the beta-binomial log-likelihood function to pi
//' at given data points via a number of sums, which is slow for high-value data
//' but fast for low-value data.
//' This function is called by \code{grad_pi} if necessary, and should not be called outside of this.
//'
//' @param ms Numeric vector. Vector of number of successes
//' @param ns Numeric vector. Vector of number of trials
//' @param pi Number. Probability of success; \code{0 >= pi >= 1}
//' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
//' @return A numeric vector of the same length as ms and ns, containing the gradient to pi in the give data points
//' @export
// [[Rcpp::export]]
NumericVector grad_pi_old(NumericVector ms, NumericVector ns, double pi, double theta) {
  int l = ms.size();
  NumericVector out(l);
  
  for(int i = 0; i < l; ++i) {
    
    int m = ms[i];
    int n = ns[i];
    double PV1 = 0;
    double PV2 = 0;
    
    if(m>0){
      for(double k = 0; k < (m); ++k) {
        PV1 += 1/(pi + k*theta);
      }
    }
    if((n-m) > 0){
      for(double k = 0; k < (n-m); ++k) {
        PV2 += 1/(1 - pi + k*theta);
      }  
    }
    out[i] = PV1 - PV2;
  }
  
  return out;
  
}


//' Exact gradient of the beta-binomial log-likelihood function for theta using sums.
//'
//' \code{grad_theta_old} calculates the value of the gradient of the beta-binomial log-likelihood function to theta
//' at given data points via a number of sums, which is slow for high-value data
//' but fast for low-value data.
//' This function is called by \code{grad_theta}
//' if necessary, and should not be called outside of this.
//'
//' @param ms Numeric vector. Vector of number of successes
//' @param ns Numeric vector. Vector of number of trials
//' @param pi Number. Probability of success; \code{0 >= pi >= 1}
//' @param theta Number. Dispersion parameter; \code{0 >= theta > +Inf}
//' @return A numeric vector of the same length as ms and ns, containing the gradient to theta in the give data points
//' @export
// [[Rcpp::export]]
NumericVector grad_theta_old(NumericVector ms, NumericVector ns, double pi, double theta) {
  int l = ms.size();
  NumericVector out(l);
  
  for(int i = 0; i < l; ++i) {
    
    int m = ms[i];
    int n = ns[i];
    double PV1 = 0;
    double PV2 = 0;
    double PV3 = 0;
    
    if(m>0){
      for(double k = 0; k < (m); ++k) {
        PV1 += k/(pi + k*theta);
      }
    }
    if((n-m) > 0){
      for(double k = 0; k < (n-m); ++k) {
        PV2 += k/(1 - pi + k*theta);
      }  
    }
    if((n) > 1){
      for(double k = 1; k < (n); ++k) {
        PV3 += k/(1 + k*theta);
      }  
    }
    out[i] = PV1 + PV2 - PV3;
  }
  
  return out;
  
}
