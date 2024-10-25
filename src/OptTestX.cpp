
// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>
#include <RcppGSL.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_interp.h>

#include <gmp.h>
#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>

using namespace Rcpp;

struct myCnT_P {NumericVector ref_counts; NumericVector var_counts; NumericVector isCase; NumericVector sprv; int MemLim; int Xtra;};

struct myHom_P {double SE; NumericVector ref_counts; NumericVector var_counts; NumericVector spr; NumericVector spv; int MemLim; int Xtra;};

struct myHet_P {NumericVector ref_counts; NumericVector var_counts; NumericVector sprv; int MemLim; int Xtra;};

struct myHetH0_P {double probshift; NumericVector ref_counts; NumericVector var_counts; NumericVector sprv; int MemLim; int Xtra;};



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  dBetaBin_cppi_MP(double M, double N, double PI, double THETA, int NecPres)
  {
    
    double OUT = 0;

    boost::multiprecision::mpfr_float::default_precision(NecPres);
    
    boost::multiprecision::mpfr_float m = M;
    boost::multiprecision::mpfr_float n = N;
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PV3 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = PI;
    boost::multiprecision::mpfr_float theta = THETA;
    
    if(m>0){
      PV1 = -lgamma(m)-lgamma(pi/theta)+lgamma(m+pi/theta)-log(m);
    }
    if((n-m)>0){
      PV2 = -lgamma(n-m)-lgamma((1-pi)/theta)+lgamma((n-m)+((1-pi)/theta))-log(n-m);
    }
    if(n>0){
      PV3 = lgamma(n)+lgamma(1/theta)-lgamma(n+(1/theta))+log(n);
    }
    PVtot = PV1+PV2+PV3;
    
    OUT = PVtot.convert_to<double>();
    
    return OUT;
      
  }




//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  dBetaBin_cppi(double M, double N, double PI, double THETA, bool LOG, int MemLim, int Xtra)
  {
    
    double OUT = 0;
    
    /*if(N < 50){
     
     if((M) > 0){
     for(double k = 0; k < (M); ++k) {
     OUT += log(PI + k*THETA) + log(1/(k+1));
     }
     }
     if((N-M) > 0){
     for(double k = 0; k < (N-M); ++k) {
     OUT += log(1 - PI + k*THETA) + log(1/(k+1));
     }  
     }
     if((N) > 1){
     for(double k = 1; k < (N); ++k) {
     OUT += log(1/(1 + k*THETA)) + log(k+1);
     }  
     }
     
    } else { */
    
    if(THETA==0){
      
      OUT = lgamma(N+1)-lgamma(M+1)-lgamma(N-M+1)+M*log(PI)+(N-M)*log(1-PI);
      
    } else {
      
      int Hlp2 = ceil(log10(fabs(R::lbeta(N, 1/THETA))));
      int NecBits = ceil((Hlp2+Xtra)/log10(2));
      
      if(Hlp2 == -2147483648 | Hlp2 == 2147483647){
        
        OUT = lgamma(N+1)-lgamma(M+1)-lgamma(N-M+1)+M*log(PI)+(N-M)*log(1-PI);
        
      } else if(NecBits <= 53){
        
        if((M) > 0){
          OUT += -R::lbeta(M, PI/THETA) - log(M);
        }
        if((N-M) > 0){
          OUT += -R::lbeta(N-M, ((1-PI)/THETA)) - log(N-M);
        }
        if((N) > 0){
          OUT += R::lbeta(N, 1/THETA) + log(N);
        }
        
      } else{
        
        int Hlp = ceil(log10(lgamma(N+(1/THETA))));
        int NecBitsTRUE = ceil((Hlp+Xtra)/log10(2));
        
        if(Hlp == -2147483648 | Hlp == 2147483647){
          
          OUT = lgamma(N+1)-lgamma(M+1)-lgamma(N-M+1)+M*log(PI)+(N-M)*log(1-PI);
          
        } else if ((NecBitsTRUE > MemLim)){ /* One might think I could combine this with the previous if-statement, but this one could give a TRUE in the case of NAN's */
          
          OUT = lgamma(N+1)-lgamma(M+1)-lgamma(N-M+1)+M*log(PI)+(N-M)*log(1-PI);
          
        } else {
          
          OUT = dBetaBin_cppi_MP(M, N, PI, THETA, Hlp + Xtra);
          
        }
        
      } 
      
    }
    
    /*} */
    
    if(LOG){
      return OUT;
    } else {
      return(exp(OUT));
    }
    
  }



//' @rdname dpqrBetaBinom
//' @export
// [[Rcpp::export]]
NumericVector
  dBetaBinom(NumericVector ms, NumericVector ns, double pi, double theta, bool LOG = false, int MemLim = 2048, int Xtra = 7)
  {
    
    int l = ms.size();
    NumericVector OUT(l);
    
    for(int i = 0; i < l; ++i) {
      
      OUT[i] = dBetaBin_cppi(ms[i], ns[i], pi, theta, LOG, MemLim, Xtra);
      
    }
    
    return OUT;
    
  }



double
  LogLikComp_het_CnTcpp (const gsl_vector *v, void *params)
  {
    double PiHet, ThetaHetControl, ThetaHetCase;
    
    struct myCnT_P *p = (struct myCnT_P *)params;
    
    PiHet = gsl_vector_get(v, 0);
    ThetaHetControl = gsl_vector_get(v, 1);
    ThetaHetCase = gsl_vector_get(v, 2);
    
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts); 
    NumericVector isCase = (p->isCase); 
    NumericVector sprv = (p->sprv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    PiHet = std::max(std::min((tanh(PiHet/2)+1)/2, 1-(1e-16)), 1e-16);
    ThetaHetControl = exp(ThetaHetControl);
    ThetaHetCase = exp(ThetaHetCase);
    
    double out = 0;
    
    for(int i = 0; i < sprv.size(); ++i) {
      
      if(isCase[i]==0){
        out += sprv[i] * dBetaBin_cppi(ref_counts[i], ref_counts[i]+var_counts[i], PiHet, ThetaHetControl, true, MemLim, Xtra);
      } else {
        out += sprv[i] * dBetaBin_cppi(ref_counts[i],  ref_counts[i]+var_counts[i], PiHet, ThetaHetCase, true, MemLim, Xtra);
      }
      
    }
    
    return -(out);
    
  }




/*double
  my_fY (const gsl_vector *v, void *params)
  {
    double x, y;
    double *p = (double *)params;
    
    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);
    
    return p[2] * (x - p[0]) * (x - p[0]) +
      p[3] * (y - p[1]) * (y - p[1]) + p[4];
  }
*/

/* The gradient of f, df = (df/dx, df/dy). */
/*void
  my_dfY (const gsl_vector *v, void *params,
         gsl_vector *df)
  {
    double x, y;
    double *p = (double *)params;
    
    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);
    
    gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
    gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
  }
*/







//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradPi_cppi_MP(double M, double N, double PI, double THETA, int NecPres)
  {
    
    double OUT = 0;
    
    boost::multiprecision::mpfr_float::default_precision(NecPres);
    
    boost::multiprecision::mpfr_float m = M;
    boost::multiprecision::mpfr_float n = N;
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = PI;
    boost::multiprecision::mpfr_float theta = THETA;
    
    if(m>0){
      PV1 = (1/theta) * (boost::math::digamma(m+(pi/theta)) - boost::math::digamma(pi/theta));
    }
    if((n-m)>0){
      PV2 = (1/theta) * (boost::math::digamma((n-m)+((1-pi)/theta)) - boost::math::digamma((1-pi)/theta));
    }
    PVtot = PV1-PV2;
    OUT = PVtot.convert_to<double>();
    
    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradPi_cppi(double M, double N, double PI, double THETA, int MemLim, int Xtra)
  {
    
    double OUT = 0;
    
    /*if(N < 1){
      
      if((M) > 0){
        for(double k = 0; k < (M); ++k) {
          OUT += 1/(PI + k*THETA);
        }
      }
      if((N-M) > 0){
        for(double k = 0; k < (N-M); ++k) {
          OUT -= 1/(1 - PI + k*THETA);
        }  
      }
      
    } else { */
      
      if(THETA == 0){
        
        OUT = M/PI - (N-M)/(1-PI);
        
      } else {
        
        int Hlp = 0;
        
        if(((N-M)+((1-PI)/THETA))>(M+(PI/THETA))){
          /* Asymptotic behaviour of digamma function, see wikipedia */
          Hlp = ceil(log10((1/THETA)*(log(((N-M)+((1-PI)/THETA)))-(1/(2*((N-M)+((1-PI)/THETA)))))));
        } else {
          Hlp = ceil(log10((1/THETA)*(log((M+(PI/THETA)))-(1/(2*(M+(PI/THETA)))))));
        }
        
        int NecBitsTRUE = ceil((Hlp+Xtra)/log10(2));
        
        if(Hlp == 2147483647 | Hlp == -2147483648){
          
          OUT = M/PI - (N-M)/(1-PI);
          
        } else if((PI/THETA) < 1e-100 & ((1-PI)/THETA) < 1e-100){
          
          OUT = (1/PI) - (1/(1-PI));
        
        } else if(NecBitsTRUE <= 53){
          
          if((M) > 0){
            OUT += (1/THETA) * (boost::math::digamma(M+(PI/THETA)) - boost::math::digamma(PI/THETA));
          }
          if((N-M) > 0){
            OUT -= (1/THETA) * (boost::math::digamma((N-M)+((1-PI)/THETA)) - boost::math::digamma((1-PI)/THETA));
          }
          
        } else if( (M-1)*THETA/PI < 0.01 & (N-M-1)*THETA/(1-PI) < 0.01 &
          (M < 1e12 | PI*N < 1e12) & /* Only if these are both of a similar and large order of magnitude can catastrophic cancellation occur. */
          fabs((pow(THETA, 4)*pow(M-1, 5)))/pow(PI, 4) < 1e-7 & fabs((pow(THETA, 4)*pow(N-M-1, 5)))/pow((1-PI), 4) < 1e-7 ){
          
          /* See the Theta gradient function... */
          
          OUT += (M-PI*N)/(PI*(1-PI)) +
                 THETA/2*( -(M-1)*M/pow(PI, 2) + (N-M-1)*(N-M)/pow((1-PI), 2) ) +
                 pow(THETA, 2)/6*( (M-1)*M*(2*M-1)/pow(PI, 3) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 3) ) +
                 pow(THETA, 3)/4*( -M*M*(M-1)*(M-1)/pow(PI, 4) + (N-M)*(N-M)*(N-M-1)*(N-M-1)/pow((1-PI), 4) );
          
        } else if ((NecBitsTRUE > MemLim)){
          
          OUT = M/PI - (N-M)/(1-PI);
          
        } else {
          
          OUT = GradPi_cppi_MP(M, N, PI, THETA, Hlp+Xtra);
          
        }
        
      }
      
    /*} */

    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  grad_pi(NumericVector ms, NumericVector ns, double pi, double theta, int MemLim = 2048, int Xtra = 7)
  {
    
    int l = ms.size();
    NumericVector OUT(l);
    
    for(int i = 0; i < l; ++i) {
      
      OUT[i] = GradPi_cppi(ms[i], ns[i], pi, theta, MemLim, Xtra);
      
    }
    
    return OUT;
    
  }




//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradTheta_cppi_MP(double M, double N, double PI, double THETA, int NecPres)
  {
    
    double OUT = 0;
    
    boost::multiprecision::mpfr_float::default_precision(NecPres);
    
    boost::multiprecision::mpfr_float m = M;
    boost::multiprecision::mpfr_float n = N;
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PV3 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = PI;
    boost::multiprecision::mpfr_float theta = THETA;
    
    if(m>0){
      PV1 = -(pi/theta/theta)*(boost::math::digamma(m+(pi/theta)) - boost::math::digamma(pi/theta));
    }
    if((n-m)>0){
      PV2 = -((1-pi)/theta/theta)*(boost::math::digamma((n-m)+((1-pi)/theta)) - boost::math::digamma((1-pi)/theta));
    }
    if(n>0){
      PV3 = -(1/theta/theta)*(boost::math::digamma(n+(1/theta)) - boost::math::digamma(1/theta));
    }
    PVtot = PV1+PV2-PV3;
    OUT = PVtot.convert_to<double>();
    
    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradTheta_cppi(double M, double N, double PI, double THETA, int MemLim, int Xtra)
  {
    
    double OUT = 0;
    
    /*if(N < 1){
      
      if((M) > 0){
        for(double k = 0; k < (M); ++k) {
          OUT += k/(PI + k*THETA);
        }
      }
      if((N-M) > 0){
        for(double k = 0; k < (N-M); ++k) {
          OUT += k/(1 - PI + k*THETA);
        }  
      }
      if((N) > 1){
        for(double k = 1; k < (N); ++k) {
          OUT -= k/(1 + k*THETA);
        } 
      }
      
    } else { */
      
      if(THETA == 0){
        
        if((M) > 0){
          OUT += (M*(M-1)/2)/PI;
        }
        if((N-M) > 0){
          OUT += ((N-M)*((N-M)-1)/2)/(1-PI);
        }
        if((N) > 0){
          OUT -= (N*(N-1)/2);
        }
        
      } else {
        
        int Hlp = ceil(log10((1/THETA/THETA)*(log((N+(1/THETA)))-(1/(2*(N+(1/THETA)))))));
        int NecBitsTRUE = ceil((Hlp+Xtra)/log10(2));
        
        if(Hlp == 2147483647 | Hlp == -2147483648){ /* Because of an Inf - Inf, I presume... so VERY small THETA */
          
          if((M) > 0){
            OUT += (M*(M-1)/2)/PI;
          }
          if((N-M) > 0){
            OUT += ((N-M)*((N-M)-1)/2)/(1-PI);
          }
          if((N) > 0){
            OUT -= (N*(N-1)/2);
          }
          
          
        } else if((PI/THETA) < 1e-100 & ((1-PI)/THETA) < 1e-100){
          
          OUT = -1/THETA;
          
        } else if (NecBitsTRUE <= 53){
          
          if((M) > 0){
            OUT += -(PI/THETA/THETA)*(boost::math::digamma(M+(PI/THETA)) - boost::math::digamma(PI/THETA));
          }
          if((N-M) > 0){
            OUT += -((1-PI)/THETA/THETA)*(boost::math::digamma((N-M)+((1-PI)/THETA)) - boost::math::digamma((1-PI)/THETA));
          }
          if((N) > 0){
            OUT -= -(1/THETA/THETA)*(boost::math::digamma(N+(1/THETA)) - boost::math::digamma(1/THETA));
          }
          
          
        } else if( (N-1)*THETA < 0.01 & (M-1)*THETA/PI < 0.01 & (N-M-1)*THETA/(1-PI) < 0.01 &
          N*(N-1)/2 < 1e12 & 
          fabs((pow(THETA, 3)*pow(N-1, 5))) < 1e-7 & fabs((pow(THETA, 3)*pow(M-1, 5)))/pow(PI, 4) + fabs((pow(THETA, 3)*pow(N-M-1, 5)))/pow((1-PI), 4) < 1e-7 ){
          
          /* Small x required for accurate 1/(1+x) Taylor expansion.
           *  I mean... actually, we'll calculate the error at a later step,
           *  but it doesn't hurt to check for a small x to begin with, otherwise it's not really worth the trouble. */
          
          /* So far so good. Next up, we don't want catastrophic cancellation to occur among the 1st order Taylor terms
           * using default double precision. For this, we check to order of magnitude of the negative term (containing N and no PI).
           * If this is large, one of the positive terms will - probably - also be large, and catastrophic cancelation can occur.
           * If it isn't large, but one of the positive terms is (this can happen due to PI), this isn't necessarily a problem.
           * Sure, the results won't necessarily be precise until a certain number of decimal digits, BUT it will be very large while the negative term isn't,
           * and the very problem of catastrophic cancellation is terms of similar magnitude but opposite sign yielding a small result... which won't happen here.
           * suppose we want precision up until 5 decimal digits, then this order of magnitude shouldn't be higher than 1e10, as a conservative measure (extreme scenario) */
          
          /* Finally, I'll include Taylor Terms up until the third order, meaning the error is dictated by O(x^4)...
           * see notes / articles / PhD / whatever for the specifics, including the actual code inside this if-loop (Faulhaber's formula etc.). */
          
          OUT += 0.5*( -(N-1)*N + (M-1)*M/PI + (N-M-1)*(N-M)/(1-PI) ) +
            THETA/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) ) +
            pow(THETA, 2)/4*( -N*N*(N-1)*(N-1) + M*M*(M-1)*(M-1)/pow(PI, 3) + (N-M)*(N-M)*(N-M-1)*(N-M-1)/pow((1-PI), 3) );
          
        } else if ((NecBitsTRUE > MemLim)){
          
          if((M) > 0){
            OUT += (M*(M-1)/2)/PI;
          }
          if((N-M) > 0){
            OUT += ((N-M)*((N-M)-1)/2)/(1-PI);
          }
          if((N) > 0){
            OUT -= (N*(N-1)/2);
          }
          
        } else {
          
          OUT = GradTheta_cppi_MP(M, N, PI, THETA, Hlp+Xtra);
          
        }
        
      }
      
      
    /*} */
    
    return OUT;
    
  }



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  grad_theta(NumericVector ms, NumericVector ns, double pi, double theta, int MemLim = 2048, int Xtra = 7)
  {
    
    int l = ms.size();
    NumericVector OUT(l);
    
    for(int i = 0; i < l; ++i) {
      
      OUT[i] = GradTheta_cppi(ms[i], ns[i], pi, theta, MemLim, Xtra);
      
    }
    
    return OUT;
    
  }




void
  LogLikComp_het_CnTcpp_Df (const gsl_vector *v, void *params, gsl_vector *df)
  {
    double PiHetFactor, PiHet, ThetaHetControl, ThetaHetCase;
    
    struct myCnT_P *p = (struct myCnT_P *)params;
    
    PiHet = gsl_vector_get(v, 0);
    ThetaHetControl = gsl_vector_get(v, 1);
    ThetaHetCase = gsl_vector_get(v, 2);
    
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts); 
    NumericVector isCase = (p->isCase); 
    NumericVector sprv = (p->sprv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    PiHetFactor = 1/(2+2*cosh(PiHet));
    PiHet = std::max(std::min((tanh(PiHet/2)+1)/2, 1-(1e-16)), 1e-16);
    ThetaHetControl = exp(ThetaHetControl);
    ThetaHetCase = exp(ThetaHetCase);
    
    double GradPi = 0;
    double GradThetaControl = 0;
    double GradThetaCase = 0;
    
    for(int i = 0; i < sprv.size(); ++i) {
      
      if(isCase[i]==0){
        GradPi += sprv[i] * PiHetFactor * GradPi_cppi(ref_counts[i], ref_counts[i]+var_counts[i], PiHet, ThetaHetControl, MemLim, Xtra);
        GradThetaControl += sprv[i] * ThetaHetControl * GradTheta_cppi(ref_counts[i], ref_counts[i]+var_counts[i], PiHet, ThetaHetControl, MemLim, Xtra);
      } else {
        GradPi += sprv[i] * PiHetFactor * GradPi_cppi(ref_counts[i],  ref_counts[i]+var_counts[i], PiHet, ThetaHetCase, MemLim, Xtra);
        GradThetaCase += sprv[i] * ThetaHetCase * GradTheta_cppi(ref_counts[i], ref_counts[i]+var_counts[i], PiHet, ThetaHetCase, MemLim, Xtra);
      }
      
    }
    
    gsl_vector_set(df, 0, -GradPi);
    gsl_vector_set(df, 1, -GradThetaControl);
    gsl_vector_set(df, 2, -GradThetaCase);
    
  }





/* Compute both f and df together. */
void
  LogLikComp_het_CnTcpp_fdfY (const gsl_vector *x, void *params,
          double *f, gsl_vector *df)
  {
    *f = LogLikComp_het_CnTcpp(x, params);
    LogLikComp_het_CnTcpp_Df(x, params, df);
  }




/* struct myCnT_P {NumericVector ref_counts; NumericVector var_counts; NumericVector isCase; NumericVector sprv;}; */

/* OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))), 
 fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, method = "BFGS",
 ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)},
 error = function(e) NULL) */
 


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  CppCnT_Optim (NumericVector StartVals, NumericVector ref_counts, NumericVector var_counts, NumericVector isCase, NumericVector sprv, int MemLim = 2048, int Xtra = 7,
                double step_size = 0.01, double tol = 0.1, double epsabs = 1e-3) /* set epsabs to 1e-1 for permutation procedures!!! */
  {
    
    double PiStart = StartVals[0];
    double ThetaControlStart = StartVals[1];
    double ThetaCaseStart = StartVals[2];
    
    NumericVector R(4);   // allocate a return vector
    
    size_t iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    /* Position of the minimum (1,2), scale factors
     10,20, height 30. */
    struct myCnT_P params = {ref_counts, var_counts, isCase, sprv, MemLim, Xtra};
    
    gsl_vector *x;
    gsl_multimin_function_fdf my_func;
    
    my_func.n = 3;
    my_func.f = LogLikComp_het_CnTcpp;
    my_func.df = LogLikComp_het_CnTcpp_Df;
    my_func.fdf = LogLikComp_het_CnTcpp_fdfY;
    my_func.params = &params;
    
    /* Starting point, x = (5,7) */
    x = gsl_vector_alloc (3);
    gsl_vector_set (x, 0, PiStart);
    gsl_vector_set (x, 1, ThetaControlStart);
    gsl_vector_set (x, 2, ThetaCaseStart);
    
    double m = 1000;
    
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, 3);
    
    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tol);
    
    do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
      
      if (status)
        break;
      
      status = gsl_multimin_test_gradient (s->gradient, epsabs);
      
      /*
       if (status == GSL_SUCCESS)
       printf ("Minimum found at:\n");
       
       printf ("%5d %.5f %.5f %10.5f\n", iter,
       gsl_vector_get (s->x, 0),
       gsl_vector_get (s->x, 1),
       s->f);
       */
      
    }
    while (status == GSL_CONTINUE && iter < 100);
    
    m = gsl_multimin_fdfminimizer_minimum(s);
    
    R(0) = m;
    R(1) = gsl_vector_get (s->x, 0);
    R(2) = gsl_vector_get (s->x, 1);
    R(3) = gsl_vector_get (s->x, 2);
    
    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
    
    return wrap(R);
  }










double
  LogLikComp_homcpp (const gsl_vector *v, void *params)
  {
    double ThetaHom;
    
    struct myHom_P *p = (struct myHom_P *)params;
    
    ThetaHom = gsl_vector_get(v, 0);
    
    double SE = (p->SE); 
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts); 
    NumericVector spr = (p->spr);
    NumericVector spv = (p->spv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    ThetaHom = exp(ThetaHom);

    double out = 0;
    
    for(int i = 0; i < spr.size(); ++i) {
      
      out += spr[i] * dBetaBin_cppi(ref_counts[i], ref_counts[i]+var_counts[i], 1-SE, ThetaHom, true, MemLim, Xtra) + 
        spv[i] * dBetaBin_cppi(var_counts[i],  ref_counts[i]+var_counts[i], 1-SE, ThetaHom, true, MemLim, Xtra);

    }
    
    return -(out);
    
  }



void
  LogLikComp_homcpp_Df (const gsl_vector *v, void *params, gsl_vector *df)
  {
    double ThetaHom;
    
    struct myHom_P *p = (struct myHom_P *)params;
    
    ThetaHom = gsl_vector_get(v, 0);

    double SE = (p->SE); 
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts); 
    NumericVector spr = (p->spr);
    NumericVector spv = (p->spv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    ThetaHom = exp(ThetaHom);

    double Grad = 0;
    
    for(int i = 0; i < spr.size(); ++i) {
      
      Grad += spr[i] * ThetaHom * GradTheta_cppi(ref_counts[i], ref_counts[i]+var_counts[i], 1-SE, ThetaHom, MemLim, Xtra) + 
        spv[i] * ThetaHom * GradTheta_cppi(var_counts[i], ref_counts[i]+var_counts[i], 1-SE, ThetaHom, MemLim, Xtra);
      
    }
    
    gsl_vector_set(df, 0, -Grad);
    
  }




void
  LogLikComp_homcpp_fdfY (const gsl_vector *x, void *params,
                              double *f, gsl_vector *df)
  {
    *f = LogLikComp_homcpp(x, params);
    LogLikComp_homcpp_Df(x, params, df);
  }



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  CppHom_Optim (double ThetaHomStart, double SE, NumericVector ref_counts, NumericVector var_counts, NumericVector spr, NumericVector spv,
                int MemLim = 2048, int Xtra = 7, double step_size = 0.01, double tol = 0.1, double epsabs = 1e-3) /* set epsabs to 1e-1 for permutation procedures!!! */
  {
    
    NumericVector R(2);   // allocate a return vector
    
    size_t iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    struct myHom_P params = {SE, ref_counts, var_counts, spr, spv, MemLim, Xtra};
    
    gsl_vector *x;
    gsl_multimin_function_fdf my_func;
    
    my_func.n = 1;
    my_func.f = LogLikComp_homcpp;
    my_func.df = LogLikComp_homcpp_Df;
    my_func.fdf = LogLikComp_homcpp_fdfY;
    my_func.params = &params;
    
    x = gsl_vector_alloc (1);
    gsl_vector_set (x, 0, ThetaHomStart);
    
    double m = 1000;
    
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, 1);
    
    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tol);
    
    do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
      
      if (status)
        break;
      
      status = gsl_multimin_test_gradient (s->gradient, epsabs);
      
      
    }
    while (status == GSL_CONTINUE && iter < 100);
    
    m = gsl_multimin_fdfminimizer_minimum(s);
    
    R(0) = m;
    R(1) = gsl_vector_get (s->x, 0);
    
    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
    
    return wrap(R);
  }







double
  LogLikComp_hetcpp (const gsl_vector *v, void *params)
  {
    double PiHet, ThetaHet;
    
    struct myHet_P *p = (struct myHet_P *)params;
    
    PiHet = gsl_vector_get(v, 0);
    ThetaHet = gsl_vector_get(v, 1);
    
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts);
    NumericVector sprv = (p->sprv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    PiHet = std::max(std::min((tanh(PiHet/2)+1)/2, 1-(1e-16)), 1e-16);
    ThetaHet = exp(ThetaHet);
    
    double out = 0;
    
    for(int i = 0; i < sprv.size(); ++i) {
      
      out += sprv[i] * dBetaBin_cppi(ref_counts[i], ref_counts[i]+var_counts[i], PiHet, ThetaHet, true, MemLim, Xtra);

    }
    
    return -(out);
    
  }


void
  LogLikComp_hetcpp_Df (const gsl_vector *v, void *params, gsl_vector *df)
  {
    double PiHetFactor, PiHet, ThetaHet;
    
    struct myHet_P *p = (struct myHet_P *)params;
    
    PiHet = gsl_vector_get(v, 0);
    ThetaHet = gsl_vector_get(v, 1);
    
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts);
    NumericVector sprv = (p->sprv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    PiHetFactor = 1/(2+2*cosh(PiHet));
    PiHet = std::max(std::min((tanh(PiHet/2)+1)/2, 1-(1e-16)), 1e-16);
    ThetaHet = exp(ThetaHet);
    
    double GradPi = 0;
    double GradTheta = 0;
    
    for(int i = 0; i < sprv.size(); ++i) {
      
      GradPi += sprv[i] * PiHetFactor * GradPi_cppi(ref_counts[i], ref_counts[i]+var_counts[i], PiHet, ThetaHet, MemLim, Xtra);
      GradTheta += sprv[i] * ThetaHet * GradTheta_cppi(ref_counts[i], ref_counts[i]+var_counts[i], PiHet, ThetaHet, MemLim, Xtra);
      
    }
    
    gsl_vector_set(df, 0, -GradPi);
    gsl_vector_set(df, 1, -GradTheta);
    
  }





/* Compute both f and df together. */
void
  LogLikComp_hetcpp_fdfY (const gsl_vector *x, void *params,
                              double *f, gsl_vector *df)
  {
    *f = LogLikComp_hetcpp(x, params);
    LogLikComp_hetcpp_Df(x, params, df);
  }




/* struct myCnT_P {NumericVector ref_counts; NumericVector var_counts; NumericVector isCase; NumericVector sprv;}; */

/* OptObj <- tryCatch( {optim(par = c(gtools::logit(probshift), log(min(max(theta_het_control, ResetThetaMin), ResetThetaMax)), log(min(max(theta_het_tumor, ResetThetaMin), ResetThetaMax))), 
 fn = LogLikComp_het_CnT, gr = GradComp_het_CnT, method = "BFGS",
 ref_counts = data_counts$ref_count, var_counts=data_counts$var_count, sprv = sprv, isCase = data_counts$isCase)},
 error = function(e) NULL) */



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  CppHet_Optim (NumericVector StartVals, NumericVector ref_counts, NumericVector var_counts, NumericVector sprv, int MemLim = 2048, int Xtra = 7,
                double step_size = 0.01, double tol = 0.1, double epsabs = 1e-3) /* set epsabs to 1e-1 for permutation procedures!!! */
  {
    
    double PiStart = StartVals[0];
    double ThetaStart = StartVals[1];
    
    NumericVector R(3);   // allocate a return vector
    
    size_t iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    /* Position of the minimum (1,2), scale factors
     10,20, height 30. */
    struct myHet_P params = {ref_counts, var_counts, sprv, MemLim, Xtra};
    
    gsl_vector *x;
    gsl_multimin_function_fdf my_func;
    
    my_func.n = 2;
    my_func.f = LogLikComp_hetcpp;
    my_func.df = LogLikComp_hetcpp_Df;
    my_func.fdf = LogLikComp_hetcpp_fdfY;
    my_func.params = &params;
    
    /* Starting point, x = (5,7) */
    x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, PiStart);
    gsl_vector_set (x, 1, ThetaStart);
    
    double m = 1000;
    
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, 2);
    
    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tol);
    
    do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
      
      if (status)
        break;
      
      status = gsl_multimin_test_gradient (s->gradient, epsabs);
      
      /*
       if (status == GSL_SUCCESS)
       printf ("Minimum found at:\n");
       
       printf ("%5d %.5f %.5f %10.5f\n", iter,
       gsl_vector_get (s->x, 0),
       gsl_vector_get (s->x, 1),
       s->f);
       */
      
    }
    while (status == GSL_CONTINUE && iter < 100);
    
    m = gsl_multimin_fdfminimizer_minimum(s);
    
    R(0) = m;
    R(1) = gsl_vector_get (s->x, 0);
    R(2) = gsl_vector_get (s->x, 1);
    
    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
    
    return wrap(R);
  }









double
  LogLikComp_hetH0cpp (const gsl_vector *v, void *params)
  {
    double ThetaHet;
    
    struct myHetH0_P *p = (struct myHetH0_P *)params;
    
    ThetaHet = gsl_vector_get(v, 0);
    
    double probshift = (p->probshift); 
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts); 
    NumericVector sprv = (p->sprv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    ThetaHet = exp(ThetaHet);
    
    double out = 0;
    
    for(int i = 0; i < sprv.size(); ++i) {
      
      out += sprv[i] * dBetaBin_cppi(ref_counts[i], ref_counts[i]+var_counts[i], probshift, ThetaHet, true, MemLim, Xtra);
      
    }
    
    return -(out);
    
  }


void
  LogLikComp_hetH0cpp_Df (const gsl_vector *v, void *params, gsl_vector *df)
  {
    double ThetaHet;
    
    struct myHetH0_P *p = (struct myHetH0_P *)params;
    
    ThetaHet = gsl_vector_get(v, 0);
    
    double probshift = (p->probshift); 
    NumericVector ref_counts = (p->ref_counts); 
    NumericVector var_counts = (p->var_counts); 
    NumericVector sprv = (p->sprv);
    int MemLim = (p->MemLim);
    int Xtra = (p->Xtra);
    
    ThetaHet = exp(ThetaHet);
    
    double Grad = 0;
    
    for(int i = 0; i < sprv.size(); ++i) {
      
      Grad += sprv[i] * ThetaHet * GradTheta_cppi(ref_counts[i], ref_counts[i]+var_counts[i], probshift, ThetaHet, MemLim, Xtra);
      
    }
    
    gsl_vector_set(df, 0, -Grad);
    
  }




void
  LogLikComp_hetH0cpp_fdfY (const gsl_vector *x, void *params,
                          double *f, gsl_vector *df)
  {
    *f = LogLikComp_hetH0cpp(x, params);
    LogLikComp_hetH0cpp_Df(x, params, df);
  }



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  CppHetH0_Optim (double ThetaHetStart, double probshift, NumericVector ref_counts, NumericVector var_counts, NumericVector sprv,
                int MemLim = 2048, int Xtra = 7, double step_size = 0.01, double tol = 0.1, double epsabs = 1e-3) /* set epsabs to 1e-1 for permutation procedures!!! */
  {
    
    NumericVector R(2);   // allocate a return vector
    
    size_t iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    struct myHom_P params = {probshift, ref_counts, var_counts, sprv, MemLim, Xtra};
    
    gsl_vector *x;
    gsl_multimin_function_fdf my_func;
    
    my_func.n = 1;
    my_func.f = LogLikComp_homcpp;
    my_func.df = LogLikComp_homcpp_Df;
    my_func.fdf = LogLikComp_homcpp_fdfY;
    my_func.params = &params;
    
    x = gsl_vector_alloc (1);
    gsl_vector_set (x, 0, ThetaHetStart);
    
    double m = 1000;
    
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, 1);
    
    gsl_multimin_fdfminimizer_set (s, &my_func, x, step_size, tol);
    
    do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
      
      if (status)
        break;
      
      status = gsl_multimin_test_gradient (s->gradient, epsabs);
      
      
    }
    while (status == GSL_CONTINUE && iter < 100);
    
    m = gsl_multimin_fdfminimizer_minimum(s);
    
    R(0) = m;
    R(1) = gsl_vector_get (s->x, 0);
    
    gsl_multimin_fdfminimizer_free (s);
    gsl_vector_free (x);
    
    return wrap(R);
  }







//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradPiPi_cppi_MP(double M, double N, double PI, double THETA, int NecPres)
  {
    
    double OUT = 0;
    
    boost::multiprecision::mpfr_float::default_precision(NecPres);
    
    boost::multiprecision::mpfr_float m = M;
    boost::multiprecision::mpfr_float n = N;
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = PI;
    boost::multiprecision::mpfr_float theta = THETA;
    
    if(m > 0){
      PV1 = (1/pow(theta, 2)) * (boost::math::trigamma(m+(pi/theta)) - boost::math::trigamma(pi/theta));
    }
    if((n-m) > 0){
      PV2 = (1/pow(theta, 2)) * (boost::math::trigamma((n-m)+((1-pi)/theta)) - boost::math::trigamma((1-pi)/theta));
    }
    
    PVtot = PV1+PV2;
    OUT = PVtot.convert_to<double>();
    
    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradPiPi_cppi(double M, double N, double PI, double THETA, int MemLim, int Xtra)
  {
    
    double OUT = 0;
    
    if(THETA == 0){
      
      OUT = -M/pow(PI, 2) - (N-M)/pow((1-PI), 2);
      
    } else {
      
      int Hlp = 0;
      
      if((((1-PI)/THETA))<((PI/THETA))){ /* Watch out!!! The trigamma is decreasing for x>0... so look for the lowest possible trigamma-arguments... */
        Hlp = ceil(log10((1/pow(THETA, 2))*boost::math::trigamma(((1-PI)/THETA))));
      } else {
        Hlp = ceil(log10((1/pow(THETA, 2))*boost::math::trigamma(((PI)/THETA))));
      }
      
      /*
      if(((N-M)+((1-PI)/THETA))>(M+(PI/THETA))){
      Hlp = ceil(log10((1/pow(THETA, 2))* (1/(2*pow(((N-M)+((1-PI)/THETA)),2)) + 1/((N-M)+((1-PI)/THETA)))));
      } else {
      Hlp = ceil(log10((1/pow(THETA, 2))* (1/(2*pow(((M)+((PI)/THETA)),2)) + 1/((M)+((PI)/THETA)))));
      }
      */
      
      int NecBitsTRUE = ceil((Hlp+Xtra)/log10(2));
      
      if(Hlp == 2147483647 | Hlp == -2147483648){
        
        OUT = -M/pow(PI, 2) - (N-M)/pow((1-PI), 2);
        
      } else if((PI/THETA) < 1e-100 & ((1-PI)/THETA) < 1e-100){
        
        OUT = -(1/pow(PI, 2)) - (1/pow((1-PI), 2));
        
      } else if(NecBitsTRUE <= 53){
        
        if((M) > 0){
          OUT += (1/pow(THETA, 2)) * (boost::math::trigamma(M+(PI/THETA)) - boost::math::trigamma(PI/THETA));
        }
        if((N-M) > 0){
          OUT += (1/pow(THETA, 2)) * (boost::math::trigamma((N-M)+((1-PI)/THETA)) - boost::math::trigamma((1-PI)/THETA));
        }
        
      } else if( (M-1)*THETA/PI < 0.01 & (N-M-1)*THETA/(1-PI) < 0.01 &
        (pow(PI, 2)*N+M < 1e12 | 2*PI*M < 1e12) & /* Only if these are both of a similar and large order of magnitude can catastrophic cancellation occur. */
        fabs((pow(THETA, 4)*pow(M-1, 5)))/pow(PI, 4) < 1e-7 & fabs((pow(THETA, 4)*pow(N-M-1, 5)))/pow((1-PI), 4) < 1e-7 ){
        
        /* See the Theta gradient function... */
        
        OUT += (-N*pow(PI, 2) + M*(-1 + 2*PI))/(pow((PI-1), 2)*pow(PI, 2)) -
        (((M-N)*(1+M-N)*THETA)/pow((-1+PI),3)) + ((-1+M)*M*THETA)/pow(PI,3) +
        (1/6)*((3*(1+2*M-2*N)*(M-N)*(1+M-N))/pow((-1+PI),4)-(3*(-1+M)*M*(-1+2*M))/pow(PI,4))*pow(THETA,2) +
        (1/4)*(-((4*pow((M-N),2)*pow((1+M-N),2))/pow((-1+PI),5))+(4*pow((-1+M),2)*pow(M,2))/pow(PI,5))*pow(THETA,3);
        
      } else if ((NecBitsTRUE > MemLim)){
        
        OUT = -M/pow(PI, 2) - (N-M)/pow((1-PI), 2);
        
      } else {
        
        OUT = GradPiPi_cppi_MP(M, N, PI, THETA, Hlp+Xtra);
        
      }
      
    }
    

    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  grad_pi_pi(NumericVector ms, NumericVector ns, double pi, double theta, int MemLim = 2048, int Xtra = 7)
  {
       
    int l = ms.size();
    NumericVector OUT(l);
       
    for(int i = 0; i < l; ++i) {
         
      OUT[i] = GradPiPi_cppi(ms[i], ns[i], pi, theta, MemLim, Xtra);
         
    }
       
    return OUT;
       
  }











//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradPiTheta_cppi_MP(double M, double N, double PI, double THETA, int NecPres)
  {
    
    double OUT = 0;
    
    boost::multiprecision::mpfr_float::default_precision(NecPres);
    
    boost::multiprecision::mpfr_float m = M;
    boost::multiprecision::mpfr_float n = N;
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = PI;
    boost::multiprecision::mpfr_float theta = THETA;
    
    if(m > 0){
      PV1 = (-theta*boost::math::digamma(m+pi/theta) + theta*boost::math::digamma(pi/theta) + pi*(-boost::math::trigamma(m+pi/theta) + boost::math::trigamma(pi/theta)))/pow(theta,3);
    }
    if(n-m > 0){
      PV2 = (theta*boost::math::digamma(-m+n+(1-pi)/theta)-theta*boost::math::digamma((1-pi)/theta)-(-1+pi)*(boost::math::trigamma(-m+n+(1-pi)/theta)-boost::math::trigamma((1-pi)/theta)))/pow(theta,3);
    }
    
    PVtot = PV1+PV2;
    OUT = PVtot.convert_to<double>();
    
    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradPiTheta_cppi(double M, double N, double PI, double THETA, int MemLim, int Xtra)
  {
    
    double OUT = 0;
    
    if(THETA == 0){
      
      OUT = 0.5*(((M-N)*(1+M-N))/pow((-1+PI),2)-((-1+M)*M)/pow(PI,2));
      
    } else {
      

      int Hlp = 0;
      
      /*
      if(((N-M)+((1-PI)/THETA))>(M+(PI/THETA))){
        Hlp = ceil(log10((1/THETA)*(log(((N-M)+((1-PI)/THETA)))-(1/(2*((N-M)+((1-PI)/THETA)))))));
      } else {
        Hlp = ceil(log10((1/THETA)*(log((M+(PI/THETA)))-(1/(2*(M+(PI/THETA)))))));
      }
      */
      
      int Hlp1 = ceil(log10((1/pow(THETA,2))*boost::math::digamma(M+PI/THETA)));
      int Hlp2 = ceil(log10((1/pow(THETA,2))*boost::math::digamma(-M+N+(1-PI)/THETA)));
      int Hlp3 = ceil(log10((1/pow(THETA,3))*boost::math::trigamma(PI/THETA)*PI));
      int Hlp4 = ceil(log10((1/pow(THETA,3))*boost::math::trigamma((1-PI)/THETA)*(1-PI)));
      
      Hlp = std::max(Hlp1, std::max(Hlp2, std::max(Hlp3, Hlp4)));
      
      int NecBitsTRUE = ceil((Hlp+Xtra)/log10(2));
      
      if((Hlp1 == 2147483647) | (Hlp1 == -2147483648) | (Hlp2 == 2147483647) | (Hlp2 == -2147483648) | (Hlp3 == 2147483647) | (Hlp3 == -2147483648) | (Hlp4 == 2147483647) | (Hlp4 == -2147483648)){
        
        OUT = 0.5*(((M-N)*(1+M-N))/pow((-1+PI),2)-((-1+M)*M)/pow(PI,2));
        
      } else if((PI/THETA) < 1e-100 & ((1-PI)/THETA) < 1e-100){
        
        OUT = 0; /* Maybe take a look at the Taylor expansion of 1/(1+x) for large x to double-check this...? */
        
      } else if(NecBitsTRUE <= 53){
      
        if((M) > 0){
          OUT += (-THETA*boost::math::digamma(M+PI/THETA) + THETA*boost::math::digamma(PI/THETA) + PI*(-boost::math::trigamma(M+PI/THETA) + boost::math::trigamma(PI/THETA)))/pow(THETA,3);
        }
        if((N-M) > 0){
          OUT += (THETA*boost::math::digamma(-M+N+(1-PI)/THETA)-THETA*boost::math::digamma((1-PI)/THETA)-(-1+PI)*(boost::math::trigamma(-M+N+(1-PI)/THETA)-boost::math::trigamma((1-PI)/THETA)))/pow(THETA,3);
        }
        
      } else if( ((M-1)*THETA/PI < 0.01) & ((N-M-1)*THETA/(1-PI) < 0.01) &
        ((0.5*pow(N,2)/pow((PI-1),2) < 1e12) | (std::max(N*M/pow((PI-1),2), 0.5*pow(M,2)/pow((PI-1),2)) < 1e12)) & /* Only if these are both of a similar and large order of magnitude can catastrophic cancellation occur. */
        (fabs((pow(THETA, 4)*pow(M-1, 5)))/pow(PI, 4) < 1e-7) & (fabs((pow(THETA, 4)*pow(N-M-1, 5)))/pow((1-PI), 4) < 1e-7) ){
        
        /* See the Theta gradient function... */
        
        OUT += 0.5*(((M-N)*(1+M-N))/pow((-1+PI),2)-((-1+M)*M)/pow(PI,2)) +
        (1/3)*(-(((1+2*M-2*N)*(M-N)*(1+M-N))/pow((-1+PI),3))+((-1+M)*M*(-1+2*M))/pow(PI,3))*THETA +
        (3/4)*((pow((M-N),2)*pow((1+M-N),2))/pow((-1+PI),4)-(pow((-1+M),2)*pow(M,2))/pow(PI,4))*pow(THETA,2);
        
      } else if ((NecBitsTRUE > MemLim)){
        
        OUT = 0.5*(((M-N)*(1+M-N))/pow((-1+PI),2)-((-1+M)*M)/pow(PI,2));
        
      } else {
        
        OUT = GradPiTheta_cppi_MP(M, N, PI, THETA, Hlp+Xtra);
        
      }
      
    }
    
    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  grad_pi_theta(NumericVector ms, NumericVector ns, double pi, double theta, int MemLim = 2048, int Xtra = 7)
  {
    int l = ms.size();
    NumericVector OUT(l);
     
    for(int i = 0; i < l; ++i) {
     
     OUT[i] = GradPiTheta_cppi(ms[i], ns[i], pi, theta, MemLim, Xtra);
     
    }
    
    return OUT;
     
  }







//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  MPhelper_GradThetaTheta(double M, double N, double PI, int NecPres)
  {
    
    double OUT = 0;
    
    boost::multiprecision::mpfr_float::default_precision(NecPres);
    
    boost::multiprecision::mpfr_float m = M;
    boost::multiprecision::mpfr_float n = N;
    boost::multiprecision::mpfr_float OUTMP = 0;
    boost::multiprecision::mpfr_float pi = PI;

    OUTMP = 1/6*( (n-1)*n*(2*n-1) - (m-1)*m*(2*m-1)/pow(pi, 2) - (n-m-1)*(n-m)*(2*(n-m)-1)/pow((1-pi), 2) );
    
    OUT = OUTMP.convert_to<double>();
    
    return OUT;
    
  }



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradThetaTheta_cppi_MP(double M, double N, double PI, double THETA, int NecPres)
  {
    
    double OUT = 0;
    
    boost::multiprecision::mpfr_float::default_precision(NecPres);
    
    boost::multiprecision::mpfr_float m = M;
    boost::multiprecision::mpfr_float n = N;
    boost::multiprecision::mpfr_float PV1 = 0;
    boost::multiprecision::mpfr_float PV2 = 0;
    boost::multiprecision::mpfr_float PV3 = 0;
    boost::multiprecision::mpfr_float PVtot = 0;
    boost::multiprecision::mpfr_float pi = PI;
    boost::multiprecision::mpfr_float theta = THETA;
    
    if((m) > 0){
      PV1 = (pi*(2*theta*boost::math::digamma(m+pi/theta)-2*theta*boost::math::digamma(pi/theta)+pi*(boost::math::trigamma(m+pi/theta)-boost::math::trigamma(pi/theta))))/pow(theta,4);
    }
    if((n-m) > 0){
      PV2 = ((-1+pi)*(-2*theta*boost::math::digamma(-m+n+(1-pi)/theta)+2*theta*boost::math::digamma((1-pi)/theta)+(-1+pi)*(boost::math::trigamma(-m+n+(1-pi)/theta)-boost::math::trigamma((1-pi)/theta))))/pow(theta,4);
    }
    if((n) > 0){
      PV3 = (-2*theta*boost::math::digamma(n+1/theta)+2*theta*boost::math::digamma(1/theta)-boost::math::trigamma(n+1/theta)+boost::math::trigamma(1/theta))/pow(theta,4);
    }
    
    PVtot = PV1+PV2+PV3;
    OUT = PVtot.convert_to<double>();
    
    return OUT;
    
  }


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double
  GradThetaTheta_cppi(double M, double N, double PI, double THETA, int MemLim, int Xtra)
  {
    
    double OUT = 0;
    
    if(THETA == 0){
      
      int HlpBU = ceil(log10(1/6*((N-1)*N*(2*N-1))));
      int NecBitsBU = ceil((HlpBU + 3)/log10(2));
      
      if(NecBitsBU <= 53){
        OUT = 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) );
      } else if (NecBitsBU <= MemLim){
        OUT = MPhelper_GradThetaTheta(M, N, PI, HlpBU + 3);
      } else {
        OUT = 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) );
      }
      
    } else {
      
      int Hlp = 0;
      
      int Hlp1 = ceil(log10((1/pow(THETA,3))*boost::math::digamma(M+PI/THETA)*2*PI));
      int Hlp2 = ceil(log10((1/pow(THETA,3))*boost::math::digamma(-M+N+(1-PI)/THETA)*2*(1-PI)));
      int Hlp3 = ceil(log10((1/pow(THETA,3))*boost::math::digamma(N+1/THETA)*2));
      int Hlp4 = ceil(log10((1/pow(THETA,4))*boost::math::trigamma(PI/THETA)*(PI)));
      int Hlp5 = ceil(log10((1/pow(THETA,4))*boost::math::trigamma((1-PI)/THETA)*(1-PI)));
      int Hlp6 = ceil(log10((1/pow(THETA,4))*boost::math::trigamma(1/THETA)));
      
      Hlp = std::max(Hlp1, std::max(Hlp2, std::max(Hlp3, std::max(Hlp4, std::max(Hlp5, Hlp6)))));
      
      int NecBitsTRUE = ceil((Hlp+Xtra)/log10(2));
      
      if((Hlp1 == 2147483647) | (Hlp1 == -2147483648) | (Hlp2 == 2147483647) | (Hlp2 == -2147483648) | (Hlp3 == 2147483647) | (Hlp3 == -2147483648) | 
         (Hlp4 == 2147483647) | (Hlp4 == -2147483648) | (Hlp5 == 2147483647) | (Hlp5 == -2147483648) | (Hlp6 == 2147483647) | (Hlp6 == -2147483648)){
    
        int HlpBU = ceil(log10(1/6*((N-1)*N*(2*N-1))));
        int NecBitsBU = ceil((HlpBU + 3)/log10(2));
        
        if(NecBitsBU <= 53){
          OUT = 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) );
        } else if (NecBitsBU <= MemLim){
          OUT = MPhelper_GradThetaTheta(M, N, PI, (HlpBU + 3));
        } else {
          OUT = 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) );
        }
    
      } else if (((PI/THETA) < 1e-100) & (((1-PI)/THETA) < 1e-100)){
        
        OUT = 1/pow(THETA,2);
        
      } else if (NecBitsTRUE <= 53){
        
        if((M) > 0){
          OUT += (PI*(2*THETA*boost::math::digamma(M+PI/THETA)-2*THETA*boost::math::digamma(PI/THETA)+PI*(boost::math::trigamma(M+PI/THETA)-boost::math::trigamma(PI/THETA))))/pow(THETA,4);
        }
        if((N-M) > 0){
          OUT += ((-1+PI)*(-2*THETA*boost::math::digamma(-M+N+(1-PI)/THETA)+2*THETA*boost::math::digamma((1-PI)/THETA)+(-1+PI)*(boost::math::trigamma(-M+N+(1-PI)/THETA)-boost::math::trigamma((1-PI)/THETA))))/pow(THETA,4);
        }
        if((N) > 0){
          OUT += (-2*THETA*boost::math::digamma(N+1/THETA)+2*THETA*boost::math::digamma(1/THETA)-boost::math::trigamma(N+1/THETA)+boost::math::trigamma(1/THETA))/pow(THETA,4);
        }
        
      } else if( ((N-1)*THETA < 0.01) & ((M-1)*THETA/PI < 0.01) & ((N-M-1)*THETA/(1-PI) < 0.01) &
        (N*(N-1)*(2*N-1)/6 < 1e12) &
        (fabs((pow(THETA, 3)*pow(N-1, 5))) < 1e-7) & (fabs((pow(THETA, 3)*pow(M-1, 5)))/pow(PI, 4) + fabs((pow(THETA, 3)*pow(N-M-1, 5)))/pow((1-PI), 4) < 1e-7) ){
        
        int HlpBU = ceil(log10(1/6*((N-1)*N*(2*N-1))));
        int NecBitsBU = ceil((HlpBU + 3)/log10(2));
        
        if(NecBitsBU <= 53){
          OUT += 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) ) +
            0.5*THETA*( -N*N*(N-1)*(N-1) + M*M*(M-1)*(M-1)/pow(PI, 3) + (N-M)*(N-M)*(N-M-1)*(N-M-1)/pow((1-PI), 3) );
        } else if (NecBitsBU <= MemLim){
          OUT += MPhelper_GradThetaTheta(M, N, PI, (HlpBU + 3)) +
            0.5*THETA*( -N*N*(N-1)*(N-1) + M*M*(M-1)*(M-1)/pow(PI, 3) + (N-M)*(N-M)*(N-M-1)*(N-M-1)/pow((1-PI), 3) );
        } else {
          OUT += 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) ) +
            0.5*THETA*( -N*N*(N-1)*(N-1) + M*M*(M-1)*(M-1)/pow(PI, 3) + (N-M)*(N-M)*(N-M-1)*(N-M-1)/pow((1-PI), 3) );
        }
        
        
      } else if ((NecBitsTRUE > MemLim)){
        
        int HlpBU = ceil(log10(1/6*((N-1)*N*(2*N-1))));
        int NecBitsBU = ceil((HlpBU + 3)/log10(2));
        
        if(NecBitsBU <= 53){
          OUT = 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) );
        } else if (NecBitsBU <= MemLim){
          OUT = MPhelper_GradThetaTheta(M, N, PI, (HlpBU + 3));
        } else {
          OUT = 1/6*( (N-1)*N*(2*N-1) - (M-1)*M*(2*M-1)/pow(PI, 2) - (N-M-1)*(N-M)*(2*(N-M)-1)/pow((1-PI), 2) );
        }
        
      } else {
        
        OUT = GradThetaTheta_cppi_MP(M, N, PI, THETA, Hlp+Xtra);
        
      }
      
    }
    
    
    /*} */
    
    return OUT;
    
  }



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  grad_theta_theta(NumericVector ms, NumericVector ns, double pi, double theta, int MemLim = 2048, int Xtra = 7)
  {
    int l = ms.size();
    NumericVector OUT(l);
       
    for(int i = 0; i < l; ++i) {
         
    OUT[i] = GradThetaTheta_cppi(ms[i], ns[i], pi, theta, MemLim, Xtra);
         
    }
       
    return OUT;
       
  }



