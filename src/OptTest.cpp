
// [[Rcpp::depends(RcppGSL)]]

#include <Rcpp.h>
#include <RcppGSL.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_interp.h>

using namespace Rcpp;
 

double
  my_f (const gsl_vector *v, void *params)
  {
    double x, y;
    double *p = (double *)params;
    
    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);
    
    return p[2] * (x - p[0]) * (x - p[0]) +
      p[3] * (y - p[1]) * (y - p[1]) + p[4];
  }

/* The gradient of f, df = (df/dx, df/dy). */
void
  my_df (const gsl_vector *v, void *params,
         gsl_vector *df)
  {
    double x, y;
    double *p = (double *)params;
    
    x = gsl_vector_get(v, 0);
    y = gsl_vector_get(v, 1);
    
    gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
    gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
  }

/* Compute both f and df together. */
void
  my_fdf (const gsl_vector *x, void *params,
          double *f, gsl_vector *df)
  {
    *f = my_f(x, params);
    my_df(x, params, df);
  }






//' DO SOMETHING
//' @export
// [[Rcpp::export]]
NumericVector
  MyOptTest (NumericVector parX)
  {
    
    NumericVector R(3);   // allocate a return vector
    
    size_t iter = 0;
    int status;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    
    /* Position of the minimum (1,2), scale factors
     10,20, height 30. */
    double par[5] = { parX(0), parX(1), parX(2), parX(3), parX(4) };
    
    gsl_vector *x;
    gsl_multimin_function_fdf my_func;
    
    my_func.n = 2;
    my_func.f = my_f;
    my_func.df = my_df;
    my_func.fdf = my_fdf;
    my_func.params = par;
    
    /* Starting point, x = (5,7) */
    x = gsl_vector_alloc (2);
    gsl_vector_set (x, 0, 5.0);
    gsl_vector_set (x, 1, 7.0);
    
    double m = 1000;
    
    T = gsl_multimin_fdfminimizer_vector_bfgs2;
    s = gsl_multimin_fdfminimizer_alloc (T, 2);
    
    gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);
    
    do
    {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);
      
      if (status)
        break;
      
      status = gsl_multimin_test_gradient (s->gradient, 1e-3);
      
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
