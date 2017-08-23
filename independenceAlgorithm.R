library(fda)

D=21
BASIS = create.fourier.basis(rangeval = c(0,1), nbasis=D)

rarh1 = function(n , Sigma, k){
   #burnin
   l_burnin = n/2
   data = matrix(0, D, (n+l_burnin))
   for (i in 1:(n+l_burnin)){
      data[,i] = rnorm(D, 0, Sigma)
   }
   # Psi operator
   Psi = matrix(0, D, D)
   for (i in 1:D){
      for (j in 1:D){
         Psi[i, j] = rnorm(1, 0, t(Sigma[i])%*%Sigma[j])
      }
   }
   #adjust the norm
   Psi = k*Psi/(norm(Psi))
   coef = matrix(0, D, (n+l_burnin))
   coef[ ,1] = data[,1]
   # recursion
   for (i in 2:(n+l_burnin)){
      coef[,i] = Psi %*% coef[ ,i-1] + data[,i]
   }
   dat = coef[, (l_burnin+1):(n+l_burnin)]
   fd(dat, BASIS)
}

Cov.h = function(X, h){
   n = ncol(X)
   D = nrow(X)
   A = matrix(0, D, D)
   for (i in 1:(n-h)){
      A = A + X[,i] %*% t(X[,i+h])
   }
   A/n
}

Mult = function(A, B){
   d = dim(A)[1]
   val = sapply(1:d, function(i) sapply(1:d, function(j) A[i,j]*B[i,j]))
   sum(val)
}

Q_N = function(SC, H){
   SC = as.matrix(SC)
   n = nrow(SC)
   d = ncol(SC)
   C_0 = Cov.h(t(SC), 0)
   q_n = 0
   for (h in 1:H){
      C_h = Cov.h(t(SC), h)
      A = solve(C_0)%*%C_h
      B = C_h%*%solve(C_0)
      q_n = q_n + Mult(A,B)
   }
   stat = q_n*n
   df = d^2*H
   c = qchisq(0.975, df=df)
   if (stat>c){1} else {0}
}

fdata = rarh1(n=100, Sigma=3^-(1:D), k=0.5)
c.data = center.fd(fdata)
fPCA = pca.fd(c.data, nharm = D, centerfns = TRUE)


sapply(1:15, function(d) Q_N(fPCA$scores[,1:d], H=5))


