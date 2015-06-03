#ifndef GIBBS_H
#define GIBBS_H

#include "RcppArmadillo.h"
#include <omp.h>

using namespace Rcpp;
using namespace arma;


struct kernel_coeffs_type
{
  arma::mat mu;
  arma::vec mu_0;
  arma::mat Omega;
  arma::mat Sigma;
  double epsilon;
};


class MCMC
{
  private:
  mat Y;            // the data (n,J)
  vec C;            // the group membership (n,1)
  int K;            // number of mixture components 
  int J;            // the number of groups
  int n;            // number of observations
  int p;            // observation dimension
  int num_iter, num_burnin, num_thin, num_display;   // number of iteration, burnin and thinning
  int seed;         // initial random seed

  /* --- hyperparameters --- */

  double nu_2;        // nu_2 = p + 2
  int nu_1;   // nu_1 = p + 2
  // Hyperparameter of the Inverse Wishart on Psi_1
  mat Psi_2; //  = eye<mat>(p,p);  
  // mean of the Normal prior on m_1
  vec m_2; // (p);  m_2.zeros();    
  // Covariance of the Normal prior on m_1
  mat inv_S_2; // =  eye<mat>(p,p); S_2 = S_2/1000;   
  // k_0 prior parameters
  vec tau_k0;  //  tau.fill(4);
  // alpha parameters
  vec tau_alpha;  // (2); tau_alpha.fill(1);
  double tau_epsilon;
  double trunc_epsilon;
  vec tau_epsilon0;
  double epsilon0_par;
  bool merge_step;  // 
  double merge_par;
  // latent indicator initial values
  uvec Z_input; 
  
  int length_chain; 
  
  double saveAlphaHM, saveEpsilon0HM;
  vec saveK0, saveEpsilon0, saveAlpha;
  mat saveOdds, saveEpsilon, saveM1, saveShifts;
  cube saveW, saveMu, saveMu0, saveOmega, saveOmega1;
  umat saveZ;
  
  
  void main_loop(Rcpp::List state);

  Rcpp::List GenerateZetas( mat loglike,
                            arma::mat logW );

  Rcpp::List UpdateZetas(   arma::cube mu, 
                            arma::cube Omega, 
                            arma::mat logW );
  
  double UpdateAlpha(double alpha, arma::mat N, arma::uvec R, double alpha_par);

  
  
  arma::mat UpdateLogWs(   arma::mat N, 
                               arma::uvec R,
                               double rho,
                               double alpha  );
                            
  kernel_coeffs_type UpdateMuSigmaEpsilon(  arma::uvec Z,
                              int k,  
                              arma::vec mu_0,
                              arma::mat Sigma_1, 
                              arma::mat Omega_1, 
                              double k_0, 
                              double epsilon,
                              double epsilon0,
                              arma::vec m_1,
                              arma::mat mean_std,
                              arma::mat cov_std  );      
                              
  double UpdateK0(  arma::cube Omega, 
                    arma::mat mu_0,
                    arma::vec m_1  );        
                    
  arma::mat UpdateSigma1(arma::cube Omega);
  
  arma::vec UpdateM1(  double k_0, 
                        arma::mat mu_0, 
                        arma::cube Omega );
  
                  
  double UpdateEpsilon0(  double epsilon0, 
                          arma::vec epsilon, 
                          double tau_epsilon, 
                          arma::vec tau_epsilon0,
                          double epsilon0_par  );      
                          
  double updateEta( arma::uvec R);
  
  double updateRho( arma::mat N, arma::uvec R);


  
  kernel_coeffs_type PriorMuSigmaEpsilon( arma::mat Sigma_1, 
                                          arma::mat Omega_1, 
                                          double k_0, 
                                          double epsilon0,
                                          arma::vec m_1 ); 
                                          
  arma::uvec UpdateR( arma::uvec R, arma::mat N, double alpha, double rho, double eta  );
                                        
                                          
  double merge_step_new( arma::uvec Z,
                             int k_1,
                             int k_2,
                             arma::vec mu_01,
                             arma::vec mu_02,
                             arma::mat Sigma_1,
                             double epsilon
                            );                                        

  
  public:
  
  // constructor 
  MCMC( mat Y, 
        vec C, 
        Rcpp::List prior,
        Rcpp::List mcmc,
        Rcpp::List state );
        
  Rcpp::List get_chain();
  Rcpp::List get_diff();
      
  
};



#endif


