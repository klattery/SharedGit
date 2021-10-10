// Kevin Lattery June 2020
// Conjoint Model in Stan
// Wishart prior with Barlett Decomposition (needs revision)
// CovBlock defines blocks of parameters
// Hinge function for constraints
// Parallel threads using reduce_sum

functions{
  matrix logistic_hinge(matrix x, matrix delta) {
    return delta .* log1p_exp(x ./ delta);
    // https://statmodeling.stat.columbia.edu/2017/05/19/continuous-hinge-function-bayesian-modeling/
  }
  
  real MNL_LL_par(int[] array_slice,
                  int a_beg, int a_end, // Stan determines values of these for specific core
                  matrix beta_ind,  matrix X,  vector dep,
                  int[] start, int[] end,
                  int[] task_individual
  ) {
    real ll=0; // Log Likelihood to accumulate
    real util_alt; // Utility of alternative (was ab)
    real denom; // sum of e^util_alt
    real util_alty; // Utility of alterative * dep, cumulated
   
    for (t in a_beg:a_end){ 
      denom=0;
      util_alty=0;
      for (n in start[t]:end[t]){
          util_alt = X[n] * col(beta_ind,task_individual[t]);
          denom += exp(util_alt);
         
          if(dep[n]>0){
            util_alty += util_alt * dep[n];
          }
      }
      ll+=(util_alty-log(denom));
    }
    return ll;
  }
  
  matrix make_chol(real[] d, real[] tri_val, int[ , ] tri_pos){
  // d is diagonal,
  // tri_val is vector of n elements for lower tri, tri_pos is nx2 array of positions
   int K = size(d);
   matrix[K, K] A = diag_matrix(to_vector(d)); 
   int count = 1;
   if (num_elements(tri_val) > 0){
     for (i in 1:size(tri_val)) {
       A[tri_pos[i,1], tri_pos[i,2]] = tri_val[count];
       count += 1;
     }
   }
   return A;
  }

  int tri_sum(matrix sym){
  // get count of lower triangular elements  
   int K = cols(sym);
   int count = 0;
   for (i in 2:K){
     for (j in 1:(i-1)){
      if (fabs(sym[i,j]) > .00000001) count += 1;
     }
   }
   return count;
  }
  
  int vec_not0(vector vec){
    int count = 0;
    for (i in 1:num_elements(vec)){
      if (fabs(vec[i]) > .00000001) count += 1;
    }
    return count;
  }
  
  int[,] get_pos(matrix sym, int tri_n){
    int K = cols(sym);
    int pos = 1;
    int result [tri_n, 2];
    for (i in 2:K){
      for (j in 1:(i-1)){
        if (fabs(sym[i,j]) > .00000001){
          result[pos,1] = i;
          result[pos,2] = j;
          pos += 1;
        }
      }
    }
    return(result);
  }
} // End Functions

data {
  // Sizing Constants, Number of:
  int N; // Rows in data file
  int P; // Parameters = ncol(coded independent data)
  int I; // Individuals
  int T; // Total unique tasks (individuals x mean tasks per individual)  
  
  // Main data
  matrix[N, P] ind; // Independent Data (coded)  
  vector<lower = 0, upper = 1>[N] dep; // Dep variable
  
  // Upper level priors
  cov_matrix[P] prior_cov;  // Prior covariance matrix, recommend Lenk adjustment
  int<lower = 0> df; // df over P, recommend 5
  vector[P] prior_alpha;
  real<lower = 0> a_sig;  // alpha ~ N(prior_alpha, a_sig)
 
  // Control Options: constraints, weights, blocking
  vector[P] con_sign; // Sign of constraints -1, +1 or 0
  vector[N] wts; // weights of each row
  matrix[P,P] cov_block; // Specifies blocks of covariance items
  
  // Ragged array matching, For each task in 1:T, Specify:
  int<lower = 1, upper = I> task_individual[T]; // which i in 1:I does task
  int<lower = 1, upper = N> start[T]; // which row in [1:N] task begins
  int<lower = 1, upper = N> end[T];   // which row in [1:N] task ends
  
  int<lower = 0> splitsize; // grainsize for parallel processing
}

transformed data{
  vector[N] dep_wt = dep .* wts;
  matrix[P,P] L = cholesky_decompose(prior_cov/(P + df)); // For Wishart
  real df_chi[P];
  int tri_n = tri_sum(cov_block); // # of lower tri elements 
  int tri_pos[tri_n,2] = get_pos(cov_block, tri_n); // position of elements

  int con_n = vec_not0(con_sign); // Number of parameters constrained
  int con_p[con_n];               // Array of which parameters are constrained
  matrix[con_n, I] con_delta;     // Con function scale for parameter and respondent
  
  int array_slice[T]; // Parallel threads slice across tasks 1:T  
  int count = 1;  
  for (i in 1:P){
    df_chi[i] = P + df - i + 1;
    if (fabs(con_sign[i]) > .00000001){
      con_p[count] = i;
      con_delta[count,1:I] = rep_row_vector(con_sign[i], I);
      count += 1;
    }
  }
  for (i in 1:T){
    array_slice[i] = i;
  }
}

parameters {
  vector[P] alpha; // upper MVN: mean vector of ind betas
  real<lower=0> bart_c [P]; // Bartlett diag^2
  real bart_z [tri_n]; // Bartlett lower tri
  matrix[P, I] z; // individual z scores N(0,1)
}

transformed parameters {
  matrix[P, I] beta_ind;
  {
    matrix[P,P] cov_chol;
    if (tri_n == 0){ // No off-diagonal
      cov_chol = diag_post_multiply(L, to_vector(sqrt(bart_c))); // L * diag_matrix()
    } else {
      cov_chol = L * make_chol(sqrt(bart_c),bart_z,tri_pos); 
    }
    beta_ind = rep_matrix(alpha, I) + cov_chol * z; // Unconstrained
    if (con_n >0){
      beta_ind[con_p,1:I] = con_delta .* log1p_exp(beta_ind[con_p,1:I] ./ con_delta);
    } 
  }
}

model {
  // priors on the parameters
  alpha ~ normal(prior_alpha, a_sig); // PriorAlpha can be vector of 0's or AggModel
  target+= chi_square_lpdf(bart_c|df_chi); // for (i in 1:P) bart_c[i] ~ chi_square(P + df - i + 1);
  if (tri_n > 0) bart_z ~ std_normal();
  to_vector(z) ~ std_normal(); // log probabilities of each choice in the dataset
  target += reduce_sum(MNL_LL_par, array_slice, splitsize, 
                       beta_ind, ind, dep_wt, start, end, task_individual);
} // End Model
