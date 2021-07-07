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
  int<lower=1> B; 
  int prev[N]; 
  int inc[N]; 
  real G[N]; 
  real g[N]; 
  matrix[N,N] C; // ifr
  matrix[N,B] Basis; // ifr
  real x[N];
}
transformed data {
  real delta = 1e-6;
}
parameters {
  vector[B] eta;
}

transformed parameters {
  matrix[N,N] A = rep_matrix(0,N,N);
  matrix[N,N] B1 = rep_matrix(0,N,N);
  matrix[N,N] B2 = rep_matrix(0,N,N);
  vector[N] E_prev=rep_vector(0,N);
  vector[N] E_inc=rep_vector(0,N);
  vector[N] R = rep_vector(1e-5,N);
  vector[N] convolution1=rep_vector(0,N);
  vector[N] convolution2=rep_vector(0,N);
  {
     R[1:N] = exp(Basis*eta);
  //  for(i in 1:N) {
  //    R[i] = 1.3 + sin(0.25*x[i]);
  //  }
    for(i in 1:N){
      for(j in 1:i){
        A[i,j] = R[i-j+1];
      }
    }
    
    B1[1:N,1] = rep_vector(G[1],N);
    B1[1:N,2] = rep_vector(G[2],N)+(C[1:N,1].*A[1:N,1].*B1[1:N,1]);
    B2[1:N,1] = rep_vector(g[1],N).*R[1:N];
    B2[1:N,2] = rep_vector(g[2],N).*R[1:N]+(C[1:N,1].*A[1:N,1].*B2[1:N,1]);
    
    for(i in 3:N){	
      convolution1 = (revmat(C[1:N,1:(i-1)]).*A[1:N,1:(i-1)].*B1[1:N,1:(i-1)])*rep_vector(1.0,(i-1));
      B1[1:N,i] = rep_vector(G[i],N) + convolution1;
      convolution2 = (revmat(C[1:N,1:(i-1)]).*A[1:N,1:(i-1)].*B2[1:N,1:(i-1)])*rep_vector(1.0,(i-1));
      B2[1:N,i] = rep_vector(g[i],N).*R[1:N] + convolution2;
    }
    E_prev = diagonal(B1);
    E_inc = diagonal(B2);
    
  }
}
model {
  eta ~ normal(0,1);
  prev ~ poisson(E_prev);   
  inc ~ poisson(E_inc);   
} 
generated quantities {
  vector[N] output1;
  vector[N] output2;
  vector[N] output3;
  
  for(i in 1:N){
    output1[i]=E_prev[i];
    output2[i]=E_inc[i];
    output3[i]=R[i];
  }
}
