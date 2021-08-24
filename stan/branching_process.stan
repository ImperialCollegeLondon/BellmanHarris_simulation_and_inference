functions{
  matrix revmat(matrix X){
    int n = rows(X);
    int n2 = cols(X);
    matrix[n,n2] s;
    for(i in 1:n2){
      s[1:n,i] = X[1:n,n2-i+1];
    }
    return s;
  }
}
data {
  int<lower=1> N; 
  int<lower=1> T; 
  int<lower=1> B; 
  int prev[N]; 
  matrix[N,N] G; 
  real g[N]; 
  matrix[N,N] C; // ifr
  real x[N];
  matrix[N,B] Basis; // ifr
  
}
transformed data {
  matrix[N,N]  one_mat = rep_matrix(0,N,N);
  matrix[N,N]  zero_mat = rep_matrix(0,N,N);
  vector[N] zero_vec = rep_vector(0,N);
  vector[N] one_vec = rep_vector(1, N);
  vector[N] delta_vec = rep_vector(1e-8, N);

}
parameters {
  vector[B] eta;
}

transformed parameters {
  vector[N] E_prev=one_vec;
  vector[N] R = one_vec;
  {
    matrix[N,N] A = zero_mat;
    matrix[N,N] B1 = zero_mat;
    R[1:N] = exp(Basis*eta);
    for(i in 1:N){
      for(j in 1:i){
        A[i,j] = R[i-j+1];
      }
    }
    B1[1:N,1] = G[1:N,1];
    B1[2:N,2] = G[2:N,2]+(C[2:N,1] .* A[2:N,1] .* B1[2:N,1]);
    for(i in 3:T){	
      B1[i:N,i] = G[i:N,i] + (C[i:N,(i-1):1] .* A[i:N,1:(i-1)] .* B1[i:N,1:(i-1)]) * one_vec[1:(i - 1)];
    }
    for(i in (T+1):N){	
      B1[i:N,i] = G[i:N,i] + (C[i:N,T:1] .* A[i:N,(i-T):(i-1)] .* B1[i:N,(i-T):(i-1)]) * one_vec[1:T];
    }
    E_prev = diagonal(B1)+delta_vec;    
  }
}
model {
  eta ~ normal(0,1);
  prev ~ poisson(E_prev);   
} 
generated quantities {
  vector[N] output1;
  vector[N] output2;

  
  for(i in 1:N){
    output1[i]=E_prev[i];
    output2[i]=R[i];
  }
}

