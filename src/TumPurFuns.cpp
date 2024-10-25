
// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <boost/math/special_functions/beta.hpp>
using namespace Rcpp;


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double logSums_MaxMethod_CPP(std::vector<double> logvec) {
  
  std::vector<double> r(logvec.size());
  int k=0;
  for (int i = 0; i < logvec.size(); ++i) {
    if (! std::isnan(logvec[i])) {
      r[k] = logvec[i];
      k++;
    }
  }
  r.resize(k);
  
  double X1 = *std::max_element(r.begin(), r.end());
  
  std::vector<double> H(r.size());
  for (int i = 0; i < r.size(); ++i) {
    H[i] = r[i]-X1;
  }
  
  double out = 0;
  for (int i = 0; i < H.size(); ++i) {
    out = out+exp(H[i]);
  }
  out = log(out) + X1;
  return(out);
  
}



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double logSums_MaxMethod_CPP_R(NumericVector logvec) {
  
  std::vector<double> r(logvec.size());
  int k=0;
  for (int i = 0; i < logvec.size(); ++i) {
    if (! std::isnan(logvec[i])) {
      r[k] = logvec[i];
      k++;
    }
  }
  r.resize(k);
  
  double X1 = *std::max_element(r.begin(), r.end());
  
  std::vector<double> H(r.size());
  for (int i = 0; i < r.size(); ++i) {
    H[i] = r[i]-X1;
  }
  
  double out = 0;
  for (int i = 0; i < H.size(); ++i) {
    out = out+exp(H[i]);
  }
  out = log(out) + X1;
  return(out);
  
}



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double logSums_MaxMethodSigned_CPP(std::vector<double> logvec, std::vector<int> signvec) {
  
  std::vector<double> r(logvec.size());
  std::vector<double> rs(logvec.size());
  int k=0;
  for (int i = 0; i < logvec.size(); ++i) {
    if (! std::isnan(logvec[i])) {
      r[k] = logvec[i];
      rs[k] = signvec[i];
      k++;
    }
  }
  r.resize(k);
  rs.resize(k);
  
  double X1 = *std::max_element(r.begin(), r.end());
  
  std::vector<double> H(r.size());
  for (int i = 0; i < r.size(); ++i) {
    H[i] = r[i]-X1;
  }
  
  double out = 0;
  for (int i = 0; i < H.size(); ++i) {
    if(signvec[i] == 1){
      out = out+exp(H[i]);
    }else if(signvec[i] == -1){
      out = out-exp(H[i]);
    }
    
  }
  out = log(out) + X1;
  return(out);
  
}


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double LogTrapezoidalInt_CPP(double lower, double upper, int n, double a, double b, double q, int TC, int RC, int curTR, int curRC) {
  
  double step = (upper-lower)/(n-1);
  std::vector<double> Xs(n);
  Xs[0] = lower;
  Xs[(Xs.size()-1)] = upper;
  double filler = 0;
  for (int i = 1; i < (Xs.size()-1); ++i) {
    Xs[i] = Xs[i-1]+step;
  }
  
  double pH = 0;
  double pT = 0;
  
  std::vector<double> Ys(n);
  for (int i = 0; i < Ys.size(); ++i) {
    
    /*pH = boost::math::ibeta_inv(a, b, Xs[i]);
    pT = boost::math::ibeta_inv(q*a, q*b, Xs[i]); */
    
    pH = R::qbeta(Xs[i], a, b, 1, 0);
    pT = R::qbeta(Xs[i], q*a, q*b, 1, 0);
    
    Ys[i] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+curRC*log(pT)+(curTR-curRC)*log(1-pT) +
      lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(RC-curRC)*log(pH)+((TC-curTR)-(RC-curRC))*log(1-pH);
    
    if(std::isnan(Ys[i])){ /* Somewhere, a zero times infinity occured */
      
      Ys[i] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+((curRC==0)?(0):(curRC*log(pT)))+(((curTR-curRC)==0)?(0):((curTR-curRC)*log(1-pT))) +
        lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(((RC-curRC)==0)?(0):((RC-curRC)*log(pH)))+((((TC-curTR)-(RC-curRC))==0)?(0):(((TC-curTR)-(RC-curRC))*log(1-pH)));
      
    }
    
  }
  
  /*Ys[0] = Ys[0] + log(5) - log(12);
  Ys[n-1] = Ys[n-1] + log(5) - log(12);
  Ys[1] = Ys[1] + log(13) - log(12);
  Ys[n-2] = Ys[n-2] + log(13) - log(12);*/
  
  Ys[0] = Ys[0] + log(1070017) - log(3628800);
  Ys[n-1] = Ys[n-1] + log(1070017) - log(3628800);
  
  Ys[1] = Ys[1] + log(5537111) - log(3628800);
  Ys[n-2] = Ys[n-2] + log(5537111) - log(3628800);
  
  Ys[2] = Ys[2] + log(103613) - log(403200);
  Ys[n-3] = Ys[n-3] + log(103613) - log(403200);
  
  Ys[3] = Ys[3] + log(261115) - log(145152);
  Ys[n-4] = Ys[n-4] + log(261115) - log(145152);
  
  Ys[4] = Ys[4] + log(298951) - log(725760);
  Ys[n-5] = Ys[n-5] + log(298951) - log(725760);
  
  Ys[5] = Ys[5] + log(515677) - log(403200);
  Ys[n-6] = Ys[n-6] + log(515677) - log(403200);
  
  Ys[6] = Ys[6] + log(3349879) - log(3628800);
  Ys[n-7] = Ys[n-7] + log(3349879) - log(3628800);
  
  Ys[7] = Ys[7] + log(3662753) - log(3628800);
  Ys[n-8] = Ys[n-8] + log(3662753) - log(3628800);
  
  return(log(step)+logSums_MaxMethod_CPP(Ys));
  
}




//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double LogNewtonCotes_CPP(double lower, double upper, double a, double b, double q, int TC, int RC, int curTR, int curRC, NumericVector Wvec) {
  
  int n = Wvec.size();
  double h = (upper-lower)/(n+1);
  
  std::vector<double> Xs(n);
  Xs[0] = lower+h;
  Xs[(Xs.size()-1)] = upper-h;
  for (int i = 1; i < (Xs.size()-1); ++i) {
    Xs[i] = Xs[i-1]+h;
  }
  
  std::vector<double> TWvec(n);
  std::vector<int> Svec(n);
  for (int i = 0; i < TWvec.size(); ++i) {
    TWvec[i] = fabs(Wvec[i]);
    if(Wvec[i]<0){
      Svec[i] = -1;
    }else{
      Svec[i]=1;
    }
  }
  
  double pH = 0;
  double pT = 0;
  
  std::vector<double> Ys(n);
  for (int i = 0; i < Ys.size(); ++i) {
    
    pH = R::qbeta(Xs[i], a, b, 1, 0);
    pT = R::qbeta(Xs[i], q*a, q*b, 1, 0);
    
    Ys[i] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+curRC*log(pT)+(curTR-curRC)*log(1-pT) +
      lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(RC-curRC)*log(pH)+((TC-curTR)-(RC-curRC))*log(1-pH) +
      log(TWvec[i]);
    
    if(std::isnan(Ys[i])){ /* Somewhere, a zero times infinity occured */
  
      Ys[i] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+((curRC==0)?(0):(curRC*log(pT)))+(((curTR-curRC)==0)?(0):((curTR-curRC)*log(1-pT))) +
        lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(((RC-curRC)==0)?(0):((RC-curRC)*log(pH)))+((((TC-curTR)-(RC-curRC))==0)?(0):(((TC-curTR)-(RC-curRC))*log(1-pH))) +
        log(TWvec[i]);
      
    }
    
  }
  
  return(log(h)+logSums_MaxMethodSigned_CPP(Ys, Svec));
  
}




//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double LogGaussianQuad_CPP(double lower, double upper, double a, double b, double q, int TC, int RC, int curTR, int curRC, NumericVector Wvec, NumericVector Nvec) {
  
  double bMa = (upper-lower)/2;
  double bPa = (upper+lower)/2;
  
  double pH = 0;
  double pT = 0;
  
  std::vector<double> Xs(Nvec.size());
  for (int i = 0; i < (Xs.size()); ++i) {
    Xs[i] = bPa + bMa*Nvec[i];
  }
  
  std::vector<double> Ys(Xs.size());
  for (int i = 0; i < Ys.size(); ++i) {
    
    pH = R::qbeta(Xs[i], a, b, 1, 0);
    pT = R::qbeta(Xs[i], q*a, q*b, 1, 0);
    
    Ys[i] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+curRC*log(pT)+(curTR-curRC)*log(1-pT) +
      lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(RC-curRC)*log(pH)+((TC-curTR)-(RC-curRC))*log(1-pH) +
      log(Wvec[i]);
    
    if(std::isnan(Ys[i])){ /* Somewhere, a zero times infinity occured */
  
      Ys[i] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+((curRC==0)?(0):(curRC*log(pT)))+(((curTR-curRC)==0)?(0):((curTR-curRC)*log(1-pT))) +
        lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(((RC-curRC)==0)?(0):((RC-curRC)*log(pH)))+((((TC-curTR)-(RC-curRC))==0)?(0):(((TC-curTR)-(RC-curRC))*log(1-pH))) +
        log(Wvec[i]);
      
    }
    
  }
  
  return(log(bMa)+logSums_MaxMethod_CPP(Ys));
  
}




//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double LogTanhSinhQuad_CPP(double lower, double upper, int n, double a, double b, double q, int TC, int RC, int curTR, int curRC, double prec) {
  
  const double pi = boost::math::constants::pi<double>();
  
  double bMa = (upper-lower)/2;
  double bPa = (upper+lower)/2;
  
  double pH = 0;
  double pT = 0;
  
  double h = 1;
  
  std::vector<double> Xs((2*n)+1);
  std::vector<double> XTs(Xs.size());
  std::vector<double> Ws(Xs.size());
  for (int i = 0; i < Xs.size(); ++i) {
    Xs[i] = -n + i*h;
    XTs[i] = bPa + bMa * tanh(0.5*pi*sinh(Xs[i]));
    if(Xs[i] < 0){
      Ws[i] = 0.5*pi*cosh(Xs[i])/pow(cosh(0.5*pi*sinh(Xs[i])), 2);
      Ws[Xs.size()-1-i] = Ws[i];
    }
    if(Xs[i] < 0){
      Ws[i] = pi/2;
    }
  }
  
  
  std::vector<double> Ys(XTs.size());
  for (int i2 = 0; i2 < Ys.size(); ++i2) {
    
    pH = R::qbeta(XTs[i2], a, b, 1, 0);
    pT = R::qbeta(XTs[i2], q*a, q*b, 1, 0);
    
    Ys[i2] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+curRC*log(pT)+(curTR-curRC)*log(1-pT) +
      lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(RC-curRC)*log(pH)+((TC-curTR)-(RC-curRC))*log(1-pH) +
      log(Ws[i2]);
    
  }
  
  double Intgrl = log(bPa*h) + logSums_MaxMethod_CPP(Ys);
  double Intgrl_old = Intgrl + 1000*prec;
  
  
  while(fabs(Intgrl - Intgrl_old) > prec){
    
    Intgrl_old = Intgrl;
    h = h/2;
    std::vector<double> Xs_new(Xs.size()-1);
    std::vector<double> XTs_new(Xs.size()-1);
    std::vector<double> Ws_new(Xs.size()-1);
    for (int i_n = 0; i_n < Xs_new.size(); ++i_n) {
      Xs_new[i_n] = Xs[i_n] + h;
      XTs_new[i_n] = bPa + bMa * tanh(0.5*pi*sinh(Xs_new[i_n]));
      if(Xs_new[i_n] < 0){
        Ws_new[i_n] = 0.5*pi*cosh(Xs_new[i_n])/pow(cosh(0.5*pi*sinh(Xs_new[i_n])), 2);
        Ws_new[Xs_new.size()-1-i_n] = Ws_new[i_n];
      }
    }
    
    std::vector<double> Ys_new(XTs_new.size());
    for (int i2_n = 0; i2_n < Ys_new.size(); ++i2_n) {
      
      pH = R::qbeta(XTs_new[i2_n], a, b, 1, 0);
      pT = R::qbeta(XTs_new[i2_n], q*a, q*b, 1, 0);
      
      Ys_new[i2_n] = lgamma(curTR+1)-lgamma(curRC+1)-lgamma(curTR-curRC+1)+curRC*log(pT)+(curTR-curRC)*log(1-pT) +
        lgamma((TC-curTR)+1)-lgamma((RC-curRC)+1)-lgamma((TC-curTR)-(RC-curRC)+1)+(RC-curRC)*log(pH)+((TC-curTR)-(RC-curRC))*log(1-pH) +
        log(Ws_new[i2_n]);
      
    }
    
    Ys.insert( Ys.end(), Ys_new.begin(), Ys_new.end() );
    Xs.insert( Xs.end(), Xs_new.begin(), Xs_new.end() );
    sort(Xs.begin(), Xs.end());
    
    Intgrl = log(bMa*h) + logSums_MaxMethod_CPP(Ys);
    
  }
  
  return(Intgrl);
  
}




//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPurHelpFun_CPP(std::vector<int> RCcol, int TC, std::vector<int> TumReads_oi, double a, double b, double q, std::vector<double> SCP_oi, 
                         int n, std::basic_string<char> NumIntMethod, double prec, NumericVector Wvec, NumericVector Nvec) {
  
  double lower = 0;
  double upper = 1;
  
  std::vector<double> probsvec(TumReads_oi.size());
  for (int i = 0; i < TumReads_oi.size(); ++i) {
    
    double CTROI = TumReads_oi[i];
    
    std::vector<double> probsvecELS(RCcol.size());
    int k=0;
    
    if(NumIntMethod.compare("Gregory") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
      
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogTrapezoidalInt_CPP(lower, upper, n, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol);
            k++;
          }
        }
      }
    }
    
    if(NumIntMethod.compare("TanhSinhQuad") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
        
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogTanhSinhQuad_CPP(lower, upper, n, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol, prec);
            k++;
          }
        }
      }
    }
    
    if(NumIntMethod.compare("NewtonCotes") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
        
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogNewtonCotes_CPP(lower, upper, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol, Wvec);
            k++;
          }
        }
      }
    }
    
    if(NumIntMethod.compare("GaussianQuad") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
        
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogGaussianQuad_CPP(lower, upper, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol, Wvec, Nvec);
            k++;
          }
        }
      }
    }
    
    probsvecELS.resize(k);
    probsvec[i] = logSums_MaxMethod_CPP(probsvecELS) + log(SCP_oi[i]);
    
  }
  
  return(logSums_MaxMethod_CPP(probsvec));
  
}






//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPurHelpFun_CPP_R(NumericVector RCcol, int TC, NumericVector TumReads_oi, double a, double b, double q, NumericVector SCP_oi, 
                         int n, std::basic_string<char> NumIntMethod, double prec, NumericVector Wvec, NumericVector Nvec) {
  
  double lower = 0;
  double upper = 1;
  
  std::vector<double> probsvec(TumReads_oi.size());
  for (int i = 0; i < TumReads_oi.size(); ++i) {
    
    double CTROI = TumReads_oi[i];
    
    std::vector<double> probsvecELS(RCcol.size());
    int k=0;
    
    if(NumIntMethod.compare("Gregory") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
        
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogTrapezoidalInt_CPP(lower, upper, n, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol);
            k++;
          }
        }
      }
    }
    
    if(NumIntMethod.compare("TanhSinhQuad") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
        
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogTanhSinhQuad_CPP(lower, upper, n, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol, prec);
            k++;
          }
        }
      }
    }
    
    if(NumIntMethod.compare("NewtonCotes") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
        
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogNewtonCotes_CPP(lower, upper, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol, Wvec);
            k++;
          }
        }
      }
    }
    
    if(NumIntMethod.compare("GaussianQuad") == 0){
      for (int j = 0; j < RCcol.size(); ++j) {
        double CRCcol = RCcol[j];
        
        if(CRCcol<=CTROI){
          if((RCcol[RCcol.size()-1]-CRCcol)<=(TC-CTROI)){
            probsvecELS[k] = LogGaussianQuad_CPP(lower, upper, a, b, q, TC, RCcol[RCcol.size()-1], CTROI, CRCcol, Wvec, Nvec);
            k++;
          }
        }
      }
    }
    
    probsvecELS.resize(k);
    probsvec[i] = logSums_MaxMethod_CPP(probsvecELS) + log(SCP_oi[i]);
    
  }
  
  return(logSums_MaxMethod_CPP(probsvec));
  
}





//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = exp(optpars[0]) + 1.01;
  double b = exp(optpars[1]) + 1.01;
  double q = exp(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output;
  
}




//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector qbeta_C(NumericVector qs, double a, double b) {
  
  int l = qs.size();
  NumericVector out(l);
  double q = 0;
  double p = 0;
  
  for(int i = 0; i < l; ++i) {
    
    q = qs[i];
    p = boost::math::ibeta_inv(a, b, q);
    out[i] = p;
  }

  return out;

}


//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector qbeta_C3(NumericVector qs, double a, double b) {
  
  int l = qs.size();
  NumericVector out(l);
  double q = 0;
  double p = 0;
  
  for(int i = 0; i < l; ++i) {
    
    q = qs[i];
    p = R::qbeta(q, a, b, 1, 0);
    out[i] = p;
  }
  
  return out;
  
}









//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP2(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = fabs(optpars[0]) + 1.01;
  double b = fabs(optpars[1]) + 1.01;
  double q = fabs(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output;
  
}









//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP2X(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = fabs(optpars[0]) + 1.01;
  double b = fabs(optpars[1]) + 1.01;
  double q = exp(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output;
  
}









//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP3(double q, double a, double b, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  q = fabs(q) ;
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output;
  
}





//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP3X(double q, double a, double b, double qlim, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  q = fabs(q) + qlim;
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output;
  
}







//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP4(double q, double a, double b, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  q = exp(q);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output;
  
}









//' DO SOMETHING
//' @export
// [[Rcpp::export]]
std::vector<double> TumPur_LogLik_CPP_DB1(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = exp(optpars[0]) + 1.01;
  double b = exp(optpars[1]) + 1.01;
  double q = exp(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return SampLiks;
  
}










//' DO SOMETHING
//' @export
// [[Rcpp::export]]
std::list<std::vector<double> > TumPur_LogLik_CPP_DB2(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = exp(optpars[0]) + 1.01;
  double b = exp(optpars[1]) + 1.01;
  double q = exp(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  std::list<std::vector<double> > listOfVectors;
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
    
    listOfVectors.push_back(SCP_oi_S);
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return listOfVectors;
  
}






//' DO SOMETHING
//' @export
// [[Rcpp::export]]
std::list<std::vector<double> > TumPur_LogLik_CPP_DB3(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = exp(optpars[0]) + 1.01;
  double b = exp(optpars[1]) + 1.01;
  double q = exp(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  std::list<std::vector<double> > listOfVectors;
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
    
    listOfVectors.push_back(SCP_oi);
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return listOfVectors;
  
}



//' DO SOMETHING
//' @export
// [[Rcpp::export]]
std::vector<double> BrolDB(int TC, double TP){
  
  std::vector<double> SCP_oi(TC+1);
  
  for (int j = 0; j < (TC+1); ++j) {
    SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
  }
  
  return SCP_oi;
  
}





//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP2_10(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = fabs(optpars[0]) + 1.01;
  double b = fabs(optpars[1]) + 1.01;
  double q = fabs(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output*10;
  
}






//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP2_100(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = fabs(optpars[0]) + 1.01;
  double b = fabs(optpars[1]) + 1.01;
  double q = fabs(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output*100;
  
}






//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP2_1000(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = fabs(optpars[0]) + 1.01;
  double b = fabs(optpars[1]) + 1.01;
  double q = fabs(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
      
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output*1000;
  
}









//' DO SOMETHING
//' @export
// [[Rcpp::export]]
double TumPur_LogLik_CPP2Y(NumericVector optpars, NumericVector ref_counts, NumericVector var_counts, NumericVector tumpur, NumericVector weights, double SCPthreshold, int n = 0, std::basic_string<char> NumIntMethod = "Gregory", double prec = 0.0001, NumericVector Wvec = 0, NumericVector Nvec = 0) {
  
  double a = exp(optpars[0]) + 1.01;
  double b = exp(optpars[1]) + 1.01;
  double q = fabs(optpars[2]);
  
  std::vector<double> SampLiks(ref_counts.size());
  
  for (int i = 0; i < SampLiks.size(); ++i) {
    
    int RC = ref_counts[i];
    int VC = var_counts[i];
    int TC = RC+VC;
    
    double TP = tumpur[i];
    
    std::vector<int> TumReads_oi(TC+1);
    std::vector<double> SCP_oi(TC+1);
    
    for (int j = 0; j < (TC+1); ++j) {
      TumReads_oi[j] = j;
      SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + j*log(TP) + (TC-j)*log(1-TP));
      
      if(std::isnan(SCP_oi[j])){ /* Somewhere, a zero times infinity occured */
        SCP_oi[j] = exp(lgamma(TC+1)-lgamma(j+1)-lgamma(TC-j+1) + ((j==0)?(0):(j*log(TP))) + (((TC-j)==0)?(0):((TC-j)*log(1-TP))) );
      }
    }
    
    std::vector<int> index(SCP_oi.size(), 0);
    for (int i = 0 ; i != index.size() ; i++) {
      index[i] = i;
    }
    sort(index.begin(), index.end(),
         [&](const int& a, const int& b) {
           return (SCP_oi[a] > SCP_oi[b]);
         }
    );
    
    int k = 0;
    double cumsum = 0;
    std::vector<int> TumReads_oi_S(TC+1);
    std::vector<double> SCP_oi_S(TC+1);
    while (cumsum <= SCPthreshold && cumsum != 1 && k != (TC+1)) {
      SCP_oi_S[k] = SCP_oi[index[k]];
      TumReads_oi_S[k] = TumReads_oi[index[k]];
      cumsum = cumsum + SCP_oi_S[k];
      k++;
    }
    
    SCP_oi_S.resize(k);
    TumReads_oi_S.resize(k);
    
    for (int j2 = 0; j2 < SCP_oi_S.size(); ++j2) {
      SCP_oi_S[j2] = (SCP_oi_S[j2]/cumsum);
    }
    
    
    std::vector<int> RCcol(RC+1);
    std::iota(RCcol.begin(), RCcol.end(), 0);
    SampLiks[i] = TumPurHelpFun_CPP(RCcol, TC, TumReads_oi_S, a, b, q, SCP_oi_S, n, NumIntMethod, prec, Wvec, Nvec);
    
  }
  
  double Output = 0;
  for (int OV = 0; OV < weights.size(); ++OV) {
    Output = Output - weights[OV]*SampLiks[OV];
  }
  return Output;
  
}



