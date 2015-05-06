#include <RcppArmadillo.h>
#include <omp.h>
#include "helpers.h"
#include "gibbs.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


MCMC::MCMC( mat Y, 
            vec C, 
            Rcpp::List prior,
            Rcpp::List mcmc,
            Rcpp::List state ) :
            Y(Y), 
            C(C)
{
       
      p = Y.n_cols;  
      n = Y.n_rows;

      // Number of data sets
      J = C.max() + 1;

      K = Rcpp::as<int>(prior["K"]);
      num_iter = Rcpp::as<int>(mcmc["nskip"]) * Rcpp::as<int>(mcmc["nsave"]);
      num_burnin = Rcpp::as<int>(mcmc["nburn"]);
      num_thin = Rcpp::as<int>(mcmc["nskip"]);    
      num_display = Rcpp::as<int>(mcmc["ndisplay"]);    

      // set Armadillo seed (Rcpp seed inherited from R in mpg.R)
      // NOTE: Rstudio seems to touch Armadillo's underlying RNG, so
      // differing results may be seen when run from Rstudio, see comments
      // from Dirk Eddelbuettel here:
      //  https://github.com/RcppCore/RcppArmadillo/issues/11
      seed = Rcpp::as<int>(mcmc["seed"]);
      arma::arma_rng::set_seed(seed);

      m_2 = Rcpp::as<vec>(prior["m_2"]);
      nu_2 = Rcpp::as<double>(prior["nu_2"]);    
      nu_1 = Rcpp::as<double>(prior["nu_1"]);    
      Psi_2 = Rcpp::as<mat>(prior["Psi_2"]);
      inv_S_2 = inv(Rcpp::as<mat>(prior["S_2"]));
      tau_k0 = Rcpp::as<vec>(prior["tau_k0"]);
      tau_alpha = Rcpp::as<vec>(prior["tau_alpha"]);
      tau_epsilon = Rcpp::as<double>(prior["tau_epsilon"]);
      tau_epsilon0 = Rcpp::as<vec>(prior["tau_epsilon0"]);
      merge_step = Rcpp::as<bool>(prior["merge_step"]);
      merge_par = Rcpp::as<double>(prior["merge_par"]);
      Z_input = Rcpp::as<uvec>(state["Z"]);      
      trunc_epsilon = Rcpp::as<double>(prior["trunc_epsilon"]);
      
      length_chain =  num_iter/num_thin;
      saveAlpha.set_size(length_chain);
      saveW.set_size(J,K,length_chain);
      saveK0.set_size(length_chain);
      saveEpsilon0.set_size(length_chain);
      saveEpsilon.set_size(length_chain,K);
      saveZ.set_size(length_chain,n);
      saveMu.set_size(J,K*p,length_chain);
      saveMu0.set_size(p,K,length_chain);
      saveOmega.set_size(p,K*p,length_chain);   
      saveOmega1.set_size(p,p,length_chain);  
      saveM1.set_size(length_chain,p);   
      saveOdds.set_size(length_chain,n);
      saveShifts.set_size(length_chain,n);
      
      main_loop(state);
                                  
}

void MCMC::main_loop(Rcpp::List state)
{    
  
  /* --- support variables --- */
  
  // counter
  int km = 0;
   // number of observations per group and per component
  mat N(J,K);
  // used in the swap step
  mat temp;
  vec indices;
  // latent assignments
  uvec Z = Z_input;
  
  /* --- parameters --- */
  
  // link between mean and covariance
  double k_0 = Rcpp::as<double>(state["k_0"]);
  // perturbation parameter for the mean
  double epsilon0  = Rcpp::as<double>(state["epsilon_0"]);
  double epsilon0_old = epsilon0;
  double epsilon0_par = Rcpp::as<double>(state["epsilon0MH"]);
  vec epsilon = Rcpp::as<vec>(state["epsilon"]);
  int epsilon0_count = 0; 
  int epsilon0_tot = 100;
  // mass parameter for the dirichlet prior on the mixture weights
  double alpha = Rcpp::as<double>(state["alpha"]);
  double alpha_old = alpha;
  double alpha_par = Rcpp::as<double>(state["alphaMH"]);
  double alpha_count = 0; 
  int alpha_tot = 100; 
  // mixture weights
  mat logW = log( Rcpp::as<mat>(state["w"]) );   // J x K matrix
  // mean \mu_{j,k}
  mat mus = Rcpp::as<mat>(state["mus"]);
  cube mu(J,p, K);
  for(int k = 0; k < K; k++)
    mu.slice(k) = mus.cols( p*k , p*(k+1)-1);
  // centering of mean locations \mu_k
  mat mu_0 = Rcpp::as<mat>(state["mu_0"]);     // p x K matrix
  // covariance locations
  cube Sigma(p,p,K);
  // precision locations
  mat Omegas = Rcpp::as<mat>(state["Omegas"]); 
  cube Omega(p,p,K);
  for(int k = 0; k < K; k++)
  {
    Omega.slice(k) = Omegas.cols(p*k, p*(k+1)-1);
    Sigma.slice(k) = inv_sympd(Omega.slice(k));
  }
    
  // centering for the Wishart prior 
  mat Omega_1 = Rcpp::as<mat>(state["Omega_1"]); 
  mat Sigma_1 = inv_sympd(Omega_1);  
  // mean of the mean
  vec m_1 = Rcpp::as<vec>(state["m_1"]); 
    

  kernel_coeffs_type kc;  
  
  // assign each observation to a mixture component
  
  List tempZ = UpdateZetas(mu, Omega, logW);   
  Z = Rcpp::as<uvec>(tempZ["Z"]);  
  N = Rcpp::as<mat>(tempZ["N"]);
  
  int max_threads = omp_get_max_threads();
  cout << "Maximum threads available: " << max_threads << endl;
      
  /* --- let the chain run --- */

  for(int it=0; it<(num_iter + num_burnin); it++)
  {          
    if((it+1)%num_display==0)
      cout << "Iteration: " << it + 1 << " of " << num_iter + num_burnin << endl;
            
    // MERGE STEP
    if(merge_step)
    {
      for( int k=0; k< K - 1 ; k++ )
      {
        if( sum(N.col(k)) > 0 )
        {
          for(int kk=k+1; kk < K ; kk++)
          {
            if( sum(N.col(kk)) > 0  )
            {              
              double kl_div = KL( mu_0.col(k), 
                                  mu_0.col(kk), 
                                  Sigma.slice(k), 
                                  Omega.slice(k), 
                                  Sigma.slice(kk), 
                                  Omega.slice(kk) ) / ( 0.5*epsilon(k) +  0.5*epsilon(kk));
              if( kl_div < R::qchisq(merge_par, (double)p, 1, 0) )
              {
                
                mu_0.col(k) = mu_0.col(k)*sum(N.col(k)) + mu_0.col(kk)*sum(N.col(kk));
                for(int j = 0; j<J; j++)
                  mu.slice(k).row(j) =  mu.slice(k).row(j)*( (N(j,k)+0.5 )/(N(j,k)+N(j,kk)+1.0) ) + mu.slice(k).row(j)*( (N(j,kk)+0.5 )/(N(j,k)+N(j,kk)+1.0) );
                Sigma.slice(k) = Sigma.slice(k)*sum(N.col(k)) + Sigma.slice(kk)*sum(N.col(kk));
                Omega.slice(k) = inv_sympd(Sigma.slice(k));
                epsilon(k) = epsilon(k)*sum(N.col(k)) + epsilon(kk)*sum(N.col(kk));
                N.col(k) = N.col(k) + N.col(kk);
                N.col(kk) = zeros<vec>(J);                  
                kc = PriorMuSigmaEpsilon(  Sigma_1, 
                                           Omega_1, 
                                           k_0, 
                                           epsilon0,
                                           m_1 );
                mu.slice(kk) = kc.mu;
                mu_0.col(kk) = kc.mu_0;
                Omega.slice(kk) = kc.Omega;   
                Sigma.slice(kk) = kc.Sigma;       
                epsilon(kk) = kc.epsilon;
                        
              }
  
            }
  
          }
        }
  
      }
    }
    
                
    alpha_old = alpha;  
    alpha = UpdateAlpha( alpha, N, alpha_par );
    if( it <= num_burnin )
    {
      if( alpha != alpha_old )
        alpha_count++;
           
      if( (it+1)  % alpha_tot == 0)
      {
        if( alpha_count < 30 )
          alpha_par *= 1.1;
        if( alpha_count > 50 )  
          alpha_par *= 0.9;
        alpha_count = 0;        
      }      
    }  
    else
    {
        if( alpha != alpha_old )
          alpha_count++;
    }        
    
    logW = UpdateLogWs( N, alpha );
    
    List tempZ = UpdateZetas(mu, Omega, logW);   
    Z = Rcpp::as<uvec>(tempZ["Z"]);  
    N = Rcpp::as<mat>(tempZ["N"]);  
    
    mat norm_draws = randn<mat>( K*(J+1), p);
    mat norm_draws_cov = randn<mat>( nu_1*K + n, p);
    vec dfs = zeros<vec>( K + 1 );
    dfs.rows(1, K)  =  cumsum( sum(N, 0).t() + nu_1 );
    
    // #pragma omp parallel for private(tempSMuSigma)
    for(int k=0; k < K; k++)
    {      
      kc = UpdateMuSigmaEpsilon(  Z,
                                  k,
                                  mu_0.col(k),
                                  Sigma_1, 
                                  Omega_1,
                                  k_0, 
                                  epsilon(k),
                                  epsilon0,
                                  m_1,
                                  norm_draws.rows( (J+1)*k, (J+1)*(k+1) - 1 ),
                                  norm_draws_cov.rows(dfs(k), dfs(k+1) - 1 )
                                ); 
      mu.slice(k) = kc.mu; 
      mu_0.col(k) = kc.mu_0;
      Omega.slice(k) = kc.Omega;   
      Sigma.slice(k) = kc.Sigma;  
      epsilon(k) = kc.epsilon;

    }  
    
    k_0 =  UpdateK0(Omega, mu_0, m_1);
    
    Sigma_1 = UpdateSigma1(Omega);
    Omega_1 = inv_sympd(Sigma_1);
    
    m_1 = UpdateM1( k_0, mu_0, Omega );
    
    
    epsilon0_old = epsilon0;
    epsilon0 = UpdateEpsilon0(  epsilon0_old, 
                                epsilon, 
                                tau_epsilon, 
                                tau_epsilon0,
                                epsilon0_par  );     
    
    if( it <= num_burnin )
    {
      if( epsilon0 != epsilon0_old )
        epsilon0_count++;
        
      if( (it+1)  % epsilon0_tot == 0)
      {
        if( epsilon0_count < 30 )
          epsilon0_par *= 1.1;
        if( epsilon0_count > 50 )  
          epsilon0_par *= 0.9;
        epsilon0_count = 0;
      }      
    }  
    else
    {
      if( epsilon0 != epsilon0_old )
        epsilon0_count++;
    }
    
  
    
    if( (it+1 > num_burnin) && ((it+1) % num_thin == 0))
    {  
      // save chain
      saveK0(km) = k_0;
      saveEpsilon0(km) = epsilon0;
      saveEpsilon.row(km) = epsilon.t();
      saveM1.row(km) = m_1.t();
      saveAlpha(km) = alpha;
      saveW.slice(km) = exp(logW);   
      saveOmega.slice(km) = reshape( mat(Omega.memptr(), Omega.n_elem, 1, false), p, K*p); 
      saveOmega1.slice(km) = Omega_1;
      saveMu.slice(km) = reshape( mat(mu.memptr(), mu.n_elem, 1, false), J, K*p);   
      saveMu0.slice(km) = mu_0;
      saveZ.row(km) = Z.t();  
      
      mat W = exp(logW);
      for(int i = 0; i < n; i++)
        saveOdds(km, i) = log(W(C(i),Z(i)) + DBL_MIN) 
          - log(1 + DBL_MIN - W(C(i),Z(i)) )
          - log( (sum( W.col(Z(i)) ) - W(C(i),Z(i)))/(J-1) + DBL_MIN  )
          + log( 1 + DBL_MIN - (sum( W.col(Z(i)) ) - W(C(i),Z(i)))/(J-1) );
          
      for(int i = 0; i < n; i++)
        saveShifts(km, i) = epsilon(Z(i));
      
      km++;        
    }
      
      
  }
  
  saveAlphaHM = alpha_par;
  saveEpsilon0HM = epsilon0_par;
  
  cout << endl << "MH acceptance rate " << endl;
  cout << "epsilon0: " << (double)epsilon0_count/num_iter << endl;
  cout << "alpha: " << alpha_count / (double)num_iter << endl;
  
}     




Rcpp::List MCMC::GenerateZetas( arma::mat log_like,
                                arma::mat logW  )
{
  // J is the number of data sets,
  // K is the number of mixture components
  // So, N has a row for each data set, a column for each component
  mat N(J, K);
  N.fill(0);

  // Zeta vector assigning each data point to a component
  uvec Z(n);
  Z.fill(0);

  // generate a new random uniform distribution for this iteration
  NumericVector U = runif(n);

  // log likelihood
  double tot_log_like = 0.0;
  vec prob;
  vec probsum;
  double x;
  bool not_assigned;

  int i;
  int k;
  #pragma omp parallel for private(k, prob, probsum, x, not_assigned)
  for(i = 0; i < n; i++)
  {
    prob = exp(log_like.row(i).t() + logW.row(C(i)).t());
    probsum = cumsum(prob);
    x = U(i) * sum(prob);
    not_assigned = true;
    for (k = 0; (k < K) && not_assigned; k++)
    {
      if(x <= probsum(k))
      {
        Z(i) = k;
        not_assigned = false;
      }
    }
  }

  for(i = 0; i < n; i++) {
    N(C(i), Z(i))++;
    tot_log_like += log_like(i,Z(i));
  }

  return Rcpp::List::create(  Rcpp::Named( "Z" ) = Z,
                              Rcpp::Named( "N" ) = N,
                              Rcpp::Named( "tot_log_like" ) = tot_log_like ) ;
}


Rcpp::List MCMC::UpdateZetas(   arma::cube mu,
                                arma::cube Omega,
                                arma::mat logW )
{
  mat log_like(n, K);
  uvec C_j;
  int j;  // used as private index of for loop inside omp below
  
  #pragma omp parallel for private(j, C_j)
  for(int k = 0; k < K; k++)
  {
    uvec index(1);
    index(0) = k;
    for(j=0; j < J; j++)
    {
      C_j = arma::find(C==j);
      log_like.submat(C_j,  index) = dmvnrm_arma_precision(
        Y.rows(C_j),
        mu.slice(k).row(j),
        Omega.slice(k)  );
    }
  }  
  
  Rcpp::List zetas_output = GenerateZetas(log_like, logW);
  
  return zetas_output;
}


double MCMC::UpdateAlpha(double alpha, arma::mat N, double alpha_par)
{
  double output = alpha;    
  double log_ratio = 0;
  double temp = rgammaBayes(  pow( alpha, 2 ) * alpha_par, 
                        alpha * alpha_par );                      
  
  log_ratio += R::dgamma(alpha, pow(temp,2)* alpha_par, 1/temp/alpha_par, 1);                          
  log_ratio -= R::dgamma(temp, pow(alpha,2)* alpha_par, 1/alpha/alpha_par, 1);  
  log_ratio += R::dgamma(temp, tau_alpha(0), 1/tau_alpha(1), 1);
  log_ratio -= R::dgamma(alpha, tau_alpha(0), 1/tau_alpha(1), 1);
  
  for(int j = 0; j < J; j++)
  {
      log_ratio += marginalLikeDirichlet( N.row(j).t(), (temp/K)*ones<vec>(K)  );
      log_ratio -= marginalLikeDirichlet( N.row(j).t(), (alpha/K)*ones<vec>(K)  );
  }
  if( exp(log_ratio) > R::runif(0,1) )
      output = temp;
      
  return output;
}




arma::mat MCMC::UpdateLogWs(   arma::mat N, 
                               double alpha  )
{
  
  mat logW(J,K);  

  for(int j=0; j<J; j++)
    logW.row(j) = rDirichlet( N.row(j).t() +  alpha * ones<vec>(K) / K  ).t();
  
  return logW ;      
}





kernel_coeffs_type MCMC::UpdateMuSigmaEpsilon(    arma::uvec Z,
                                                  int k, 
                                                  arma::vec mu_0,
                                                  arma::mat Sigma_1, 
                                                  arma::mat Omega_1, 
                                                  double k_0, 
                                                  double epsilon,
                                                  double epsilon0,
                                                  arma::vec m_1,
                                                  arma::mat mean_std,
                                                  arma::mat cov_std  ) 
{ 
  kernel_coeffs_type output;
  uvec Z_k = arma::find(Z==k);  
  mat data_group = Y.rows(Z_k);
  vec C_k = C(Z_k);
  int p = data_group.n_cols;
  int N_k = data_group.n_rows;   
  mat Omega(p,p);
  mat Sigma(p,p);
  mat mu(p,J);
  mat mu_0new(p,1);  
  double epsilon_new;
  
  if(N_k == 0)
  {
    Omega = WishartScaling(cov_std, Omega_1);
    Sigma = inv_sympd( Omega );       
    mu_0new = trans( mvrnormScaling(mean_std.row(0).cols(0,p-1), m_1, Sigma/k_0) );    
    
    for(int j=0; j<J; j++)
      mu.col(j) = trans( mvrnormScaling(mean_std.row(j+1).cols(0,p-1), mu_0new,  Sigma*epsilon ));      
      
    epsilon_new = 1/ rgammaBayesTruncated(tau_epsilon + 1, tau_epsilon*epsilon0, 1.0/trunc_epsilon, -1 );  
    
  }
  else  // N_k > 0
  {
        
    mat  Psi_1(p,p);
    double extra_piece_var_1 = 0;
    vec  m1_1(p);
    
    vec n_jk(J);
    mat mean_jk(p,J); mean_jk.fill(0);
    vec mean_k = mean(data_group,0).t();
    mat SS_jk(p,p), ss_jk_1(p,p);
    
    // marginal likelihood under model 0
    mat mean_k_rep = repmat( mean_k.t(), N_k, 1);
    mat SS_k = ( data_group - mean_k_rep ).t() * ( data_group - mean_k_rep );
    mat ss_k = N_k * ( mean_k - mu_0 ) * ( mean_k - mu_0 ).t();        
    
    // marginal likelihood under model 1 
    extra_piece_var_1 = k_0;
    m1_1 = k_0 * m_1;     
    SS_jk.fill(0);
    ss_jk_1.fill(0);

    for(int j=0; j<J; j++)
    {
      uvec indices = find(C_k==j);
      n_jk(j) = indices.n_elem;     
      
      if (n_jk(j) > 0)
      {
          mean_jk.col(j) = mean(data_group.rows(indices),0).t();
          mat mean_jk_rep = repmat( trans(mean_jk.col(j)),(int)n_jk(j), 1);
          SS_jk = SS_jk + (data_group.rows(indices) - mean_jk_rep).t() * ( data_group.rows(indices) - mean_jk_rep );
          ss_jk_1 = ss_jk_1 + (mean_jk.col(j) - mu_0) * (mean_jk.col(j) - mu_0).t() / (epsilon + 1.0/n_jk(j));
          extra_piece_var_1 +=  n_jk(j) / (epsilon * n_jk(j) + 1.0);
          m1_1 = m1_1 +  n_jk(j) / (epsilon * n_jk(j) + 1.0) * mean_jk.col(j);
         
      }
    }          
    Psi_1 = inv_sympd( Sigma_1 + SS_jk + ss_jk_1 );     

    Omega = WishartScaling(cov_std, Psi_1);
    Sigma = inv_sympd( Omega ); 
    m1_1 = m1_1 / extra_piece_var_1;
    mu_0new = trans( mvrnormScaling(mean_std.row(0).cols(0,p-1), m1_1, Sigma/extra_piece_var_1));
    
    double temp_ss = 0;
    for(int j=0; j<J; j++)
    {
      if( n_jk(j) > 0 )
      {
        mu.col(j) = trans( mvrnormScaling(mean_std.row(j+1).cols(0,p-1), 
          (n_jk(j)*mean_jk.col(j) + 1.0/epsilon*mu_0new)/(n_jk(j) + 1.0/epsilon), 
          Sigma/(n_jk(j) + 1.0/epsilon)));
      }
      else
        mu.col(j) = trans( mvrnormScaling(mean_std.row(j+1).cols(0,p-1), mu_0new,  Sigma*epsilon) );           
        
      temp_ss +=  as_scalar( ( mu.col(j).t() - mu_0new.t()) * Omega * ( mu.col(j) - mu_0new ) ); 
    }  

    epsilon_new = 1 / rgammaBayesTruncated(tau_epsilon + 1 + (double)p*J/2, 
                                            tau_epsilon*epsilon0 + temp_ss/2, 1.0/trunc_epsilon, -1);    
  }
  
  output.mu_0 = mu_0new;
  output.mu = mu.t();
  output.Sigma = Sigma;
  output.Omega = Omega;
  output.epsilon = epsilon_new;
  return output;    
           

};





double MCMC::UpdateK0(  arma::cube Omega, 
                        arma::mat mu_0,
                        arma::vec m_1  )
{
  double tau_2_tot = tau_k0(1);
  double tau_1_tot = tau_k0(0) + p*K;
  for(int k=0; k < K; k++)
      tau_2_tot += as_scalar( (mu_0.col(k) - m_1).t() * Omega.slice(k) * (mu_0.col(k) - m_1));  

  return rgammaBayes(tau_1_tot/2, tau_2_tot/2);
};






arma::mat MCMC::UpdateSigma1(arma::cube Omega)
{
  mat psi_2_tot = Psi_2;
  for(int k=0; k< K; k++)
      psi_2_tot += Omega.slice(k);

  return( rWishartArma(inv_sympd(psi_2_tot), K*nu_1 + nu_2) );
  
};


arma::vec MCMC::UpdateM1(   double k_0, 
                            arma::mat mu_0, 
                            arma::cube Omega )
{
  mat precision = inv_S_2;
  vec meanM = inv_S_2*m_2;
  for(int k=0; k< K; k++)
  {
      precision += k_0*Omega.slice(k);
      meanM += k_0 * ( Omega.slice(k)*mu_0.col(k) );
  }
  mat variance = inv_sympd(precision);
  mat output = mvrnormArma(1, variance*meanM, variance);
  return( output.row(0).t() );
};



double MCMC::UpdateEpsilon0(  double epsilon0, 
                              arma::vec epsilon, 
                              double tau_epsilon, 
                              arma::vec tau_epsilon0,
                              double epsilon0_par  )                                      
{
  double output = epsilon0;
  // proposal 
  double e0_new = as<double>(rbeta(1, epsilon0 * epsilon0_par, epsilon0_par * (1 - epsilon0) ));
  
  double log_acc = R::dbeta(e0_new, tau_epsilon0(0), tau_epsilon0(1), 1 );
  log_acc -= R::dbeta(epsilon0, tau_epsilon0(0), tau_epsilon0(1), 1 );    
  log_acc -= R::dbeta(e0_new, epsilon0 * epsilon0_par, epsilon0_par * (1 - epsilon0), 1 );
  log_acc += R::dbeta(epsilon0, e0_new * epsilon0_par, epsilon0_par * (1 - e0_new), 1 );
  
//  for(int k=0; k<K; k++)
//   log_acc += dInvGamma(epsilon(k), tau_epsilon + 1, tau_epsilon*e0_new) - 
//    dInvGamma(epsilon(k), tau_epsilon + 1, tau_epsilon*epsilon0);
    
  for(int k=0; k<K; k++)
  {
    log_acc += R::dgamma(1/epsilon(k), tau_epsilon + 1, 1/(tau_epsilon*e0_new), 1) - 
      R::dgamma(1/epsilon(k), tau_epsilon + 1, 1/(tau_epsilon*epsilon0),1);
    log_acc += R::pgamma( 1.0/2, tau_epsilon + 1, 1/(tau_epsilon*epsilon0), 0, 1  )-  
      R::pgamma( 1.0/2, tau_epsilon + 1, 1/(tau_epsilon*e0_new), 0, 1  );
  }

    
  if( exp(log_acc) > R::runif(0,1) )
    output = e0_new;
    
  return output;  
}






Rcpp::List MCMC::get_chain()
{
  return Rcpp::List::create(  
    Rcpp::Named( "alpha" ) = saveAlpha,
    Rcpp::Named( "epsilon0" ) = saveEpsilon0,
    Rcpp::Named( "epsilon" ) = saveEpsilon,
    Rcpp::Named( "k_0" ) = saveK0,
    Rcpp::Named( "m_1" ) = saveM1,
    Rcpp::Named( "mus" ) = saveMu,
    Rcpp::Named( "mu_0" ) = saveMu0,
    Rcpp::Named( "Omega_1" ) = saveOmega1,
    Rcpp::Named( "Omegas" ) = saveOmega,
    Rcpp::Named( "Z" ) = saveZ,
    Rcpp::Named( "w" ) = saveW,
    Rcpp::Named( "epsilon0HM" ) = saveEpsilon0HM,
    Rcpp::Named( "alphaHM" ) = saveAlphaHM
  );
};


Rcpp::List MCMC::get_diff()
{
  return Rcpp::List::create(  
    Rcpp::Named( "logOdds" ) = saveOdds,
    Rcpp::Named( "shifts" ) = saveShifts
    );
};



kernel_coeffs_type MCMC::PriorMuSigmaEpsilon(   arma::mat Sigma_1, 
                                                arma::mat Omega_1, 
                                                double k_0, 
                                                double epsilon0,
                                                arma::vec m_1 ) 
{ 
  kernel_coeffs_type output;
  mat Omega(p,p);
  mat Sigma(p,p);
  mat mu(p,J);
  mat mu_0new(p,1);  
  double epsilon;

  Omega = rWishartArma(Omega_1, nu_1);
  Sigma = inv_sympd( Omega );   
  
//  while(epsilon > trunc_epsilon)
//    epsilon = 1 / rgammaBayes(tau_epsilon + 1, tau_epsilon*epsilon0);
  epsilon = 1 / rgammaBayesTruncated(tau_epsilon + 1, tau_epsilon*epsilon0, 1.0/trunc_epsilon, -1);  
  
  mu_0new = trans(mvrnormArma(1, m_1, Sigma/k_0));    
  for(int j=0; j<J; j++)
    mu.col(j) = trans( mvrnormArma(1, mu_0new,  Sigma*epsilon ));      
  
  output.mu_0 = mu_0new;
  output.mu = mu.t();
  output.Sigma = Sigma;
  output.Omega = Omega;
  output.epsilon = epsilon;
  return output;    
};

