
library(fda)

D=21
BASIS = create.fourier.basis(rangeval = c(0,1), nbasis=D)

#---------------------------------- FAR(1) process -----------------------------------------------

far1 = function(n , Sigma){
   
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
   
   Psi =0.8* Psi/norm(Psi)
   
   
   coef = matrix(0, D, (n+l_burnin))
   coef[ ,1] = data[,1]
   
   # recursion
   for (i in 2:(n+l_burnin)){
      coef[,i] = Psi %*% coef[ ,i-1]  + data[,i]
   }
   dat = coef[, (l_burnin+1):(n+l_burnin)]
   fd(dat, BASIS)
}


#---------------------------------- FAR(1) process -----------------------------------------------

fma1 = function(n , Sigma){
   
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
   
   Psi =0.8* Psi/norm(Psi)
   
   
   coef = matrix(0, D, (n+l_burnin))
   coef[ ,1] = data[,1]
   
   # recursion
   for (i in 2:(n+l_burnin)){
      coef[,i] = Psi %*% data[ ,i-1]  + data[,i]
   }
   dat = coef[, (l_burnin+1):(n+l_burnin)]
   fd(dat, BASIS)
}




