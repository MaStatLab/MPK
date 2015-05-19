#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;




const double log2pi = std::log(2.0 * M_PI);

double log_exp_x_plus_exp_y(double x, double y) 
{
  double result;
  if ( x - y >= 100 ) result = x;
  else if ( x - y <= -100 ) result = y;
  else {
    if (x > y) {
      result = y + log( 1 + exp(x-y) );
    }
    else result = x + log( 1 + exp(y-x) );
  }
  return result;
}


arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) 
{
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

arma::mat mvrnormScaling( arma::mat Y, arma::vec mu, arma::mat sigma) 
{
  int n = Y.n_rows;
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

arma::mat mvrnormScaling2( arma::mat Y, arma::vec mu, arma::mat sigma, double c = 1.0 ) 
{
  int n = Y.n_rows;
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma) * sqrt(c);
}

arma::mat rWishartArma(arma::mat Sigma, int df)
{
  int p = Sigma.n_rows;
  vec m(p); m.zeros();
  mat X = mvrnormArma(df, m, Sigma );
  return X.t()*X;
}

arma::mat WishartScaling(arma::mat Y, arma::mat Sigma)
{
  // int df = Y.n_rows;
  int p = Sigma.n_rows;
  vec m(p); m.zeros();
  mat X = mvrnormScaling(Y, m, Sigma); 
  return X.t()*X;
}

arma::vec dmvnrm_arma_precision(  arma::mat x,  
                                  arma::rowvec mean,  
                                  arma::mat omega, 
                                  bool logd = true) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = trimatu(arma::chol(omega));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}



double rgammaBayes(double shape, double rate) 
{
  return rgamma(1, shape, 1.0/rate )(0);
}


double rgammaBayesTruncated(double shape, double rate, double left = 0, double right = -1) 
{
  // right = -1 means right = Infinity
  double u = 0;
  if(right == -1)
    u = runif(1, R::pgamma(left, shape, 1/rate, 1, 0))(0);
  else
    u = runif(1, R::pgamma(left, shape, 1/rate, 1, 0), R::pgamma(left, shape, 1/rate, 1, 0) )(0);
  double output = R::qgamma(u, shape, 1/rate, 1, 0);
  if(R::pgamma(left, shape, 1/rate, 1, 0) > 1 - 10E-10)
    output = left;

  return  output;
}



double beta_fun(arma::vec alpha, bool logB = true)
{
  double output = - lgamma(sum(alpha));
  int J = alpha.n_elem;
  
  for(int j=0; j<J; j++)
    output += lgamma(alpha(j));
  if(logB)
    return output;
  else
    return exp(output);
}



double marginalLikeDirichlet(arma::vec data, arma::vec alpha, bool logM = true)
{
  double output = beta_fun( alpha + data ) - beta_fun( alpha );
  if(logM)
    return output;
  else
    return exp(output);
}


arma::vec rDirichlet(arma::vec alpha, bool logR = true)
{
  int n = alpha.n_elem;
  vec output(n);
  for(int i=0; i<n; i++)
    output(i) = rgamma(1,alpha(i),1)(0);
  if(logR)  
    return ( log( output) - log(sum(output) ) );
  else
    return exp( log( output) - log(sum(output) ) );
}


double dGeneralizedBeta(double x, double a, double b, arma::vec extremes, bool logD = true  )
{
  vec alpha(2); alpha(0) = a; alpha(1) = b;
  double output = - beta_fun(alpha) - (sum(alpha) - 1)*log( extremes(1) - extremes(0) ) + 
    (alpha(0)-1)*log( x-extremes(0) ) + (alpha(1)-1)*log( extremes(1) - x );
    
  if(logD)
    return output;
  else
    return exp(output);
}


int sampling(vec probs)
{  
  int Kf = probs.n_elem;
  int v;
  Rcpp::NumericVector probsum(Kf);
  double x = R::runif(0.0, sum(probs));
  probsum(0)=probs(0);
  for(int k = 1; k < Kf; k++)
  {  probsum(k)=probs(k)+probsum(k-1);
  }
  if(x < probsum(0)){ v=0;}
  for (int k = 0; k < (Kf-1); k++)
  {  if(x > probsum(k)){ if(x< probsum(k+1)){v=k+1;} }
  }
  if(x > probsum(Kf-1)){v=Kf-1;}
    return v;
}




double KL(  arma::vec mu_1, 
            arma::vec mu_2, 
            arma::mat Sigma_1, 
            arma::mat Omega_1, 
            arma::mat Sigma_2, 
            arma::mat Omega_2)
{
  double output;
  int p = mu_1.n_elem;
  double val_1, val_2;
  double sign_1, sign_2;

  log_det(val_1, sign_1, Sigma_1);
  log_det(val_2, sign_2, Sigma_2);
  mat mu(p,1); 
  mu.col(0) = (mu_1 - mu_2);
  output = trace( Omega_1 * Sigma_2) + as_scalar( mu.t() * Omega_1 * mu ) - p + val_1 - val_2;
  return (output/2.0);
    
}


double dInvGamma( double x, double alpha, double beta, bool log_like = true)
{
  double output = 0;
  if(log_like)
    output = alpha * log(beta) - (alpha + 1) * log(x) - lgamma(alpha) - beta/x;
  else
    output = exp(output);
  return output;  
}



// log Multivariate gamma function (used for Wishart distn)
double multiGamma(double x, const int p = 1, bool logScale = true) {
  double ans = log(M_PI)*p*(p - 1)/4.0;
  for (int j = 1; j <= p; j++) 
    ans += lgamma(x - (1-j)/2.0);
  if(logScale == false)
    ans = exp(ans);
  return ans;
}




