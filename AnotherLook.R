library(fda)
library(corpcor) # for pseudo inverse
library(vars)

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

trace = function(res, p){
   res = t(t(res))
   d = length(res[1,])
   n = length(res[,1])
   if(d==1){
      out=(p*d+n)/(n-p*d)*var(res)
   }
   else{
      out=(p*d+n)/(n-p*d)*sum(diag(cov(res)))
   }
   out
}

fFPE = function(fdata, Pmax){
   n = ncol(fdata$coef)
   D = nrow(fdata$coef)
   
   #center the data
   data = center.fd(fdata)
   
   #fPCA
   fpca = pca.fd(data, nharm=D, centerfns=TRUE)
   
   #total variance
   vartot = sum(fpca$values)
   
   #create a matrix containing different values 
   #of fFPE for different d and p values
   values = matrix(0, D, (Pmax+1))
   for (d in 1:D){
      scores = fpca$scores[, 1:d]
      var.explain = sum(fpca$values[1:d])
      for(p in 0:Pmax)
      {
         if(d==1)
         {
            res = arima(scores, order = c(p, 0, 0), method="ML")$residuals
         }
         else
         {
            if(p==0)
            {
               mean = t(matrix(rep(colMeans(scores),n),d))
               res = scores - mean
            }
            else
            {
               # we have to give colnames to the scores, 
               # otherwise we get warnings from vars resid function
               colnames(scores) = as.character(seq(1:d))
               res = resid(VAR(scores, p=p, type="const")) + vartot - var.explain
            }
         }
         values[d, p+1] = trace(res = res, p = p) 
      }
   } 
   # compute the estimates hat.p and hat.d for optimal order p and dimension d
   hat.p = (which.min(values)-1)%/%D
   hat.d = which.min(values)%%D
   if(hat.d==0) hat.d = hat.d + D
   return(c(hat.d, hat.p))
}


algorithm2 <- function(fundata){
     n = ncol(fundata$coefs)
     D = nrow(fundata$coefs)

     fdata = center.fd(fundata)
     mu = mean.fd(fundata)
     fpca = pca.fd(fdata, nharm=D, centerfns=TRUE)
     ffpe = fFPE(fdata, Pmax=4)
     d.hat = ffpe[1]
     p.hat = ffpe[2]

     scores = fpca$scores[, 1:d.hat]
     scores.F = fpca$scores[, (d.hat+1):D]
         if(d.hat==1)
         {
            VAR.pre = predict(arima(scores, order=c(p.hat,0,0)), n.ahead=1)$pred[1]
            yhat = VAR.pre
         }else
         {
            if(p.hat==0)
            {
               yhat = colMeans(scores)
            }else
            {
               # we need to give colnames to the scores 
               # to avoid warnings from vars predict function below
               colnames(scores) <- as.character(seq(1:d.hat))
               VAR.pre= predict(VAR(scores, p.hat), n.ahead=1, type="const")$fcst
               yhat=c()
               for(i in 1:d.hat)
               {
                  yhat=c(yhat,VAR.pre[[i]][1])
               }
            }
         }
      FAR.pre = fpca$harmonics[1]*0
      for(i in 1:d.hat)
      {
         FAR.pre = FAR.pre + yhat[i]*fpca$harmonics[i]
      }

    kk = 0.9
    w = sapply(1:100, function(t) kk*(1-kk)^(n-t))
    X.w = c()
    for (j in 1:(D-d.hat)){
        X.w[j] = weighted.mean(scores.F[,j], weights=w)
    }

   X.F1 = fpca$harmonics[1]*0
   X.F2 = fpca$harmonics[1]*0
   X.F3 = fpca$harmonics[1]*0
   X.bar = colMeans(scores.F)
   if (p.hat == 1){
      X_n = scores.F[n,]
    } else {
        X_n = colMeans(scores.F[(n-p.hat+1):n,])
    }

    for(i in (d.hat+1):D)
        {
         X.F1 = X.F1 + X.bar[i-d.hat]*fpca$harmonics[i]
         X.F2 = X.F2 + X_n[i-d.hat]*fpca$harmonics[i]
         X.F3 = X.F3 + X.w[i-d.hat]*fpca$harmonics[i]
        }
  X.hat1 = FAR.pre + X.F1 + mu
  X.hat2 = FAR.pre + X.F2 + mu
  X.hat3 = FAR.pre + X.F3 + mu
  X.ffpe = FAR.pre + mu

  list(allmean = X.hat1, last = X.hat2, wmean = X.hat3, ffpe = X.ffpe)
}


simulation <- function(n){
   fundata = rarh1(n=n, Sigma=2^-(1:D), k=0.8)
   m = n*0.9
   E1 = E2 = E3 = E4 = c()
   for (s in m:(n-1)){
      ffdata = fundata[1:j]
      res = algorithm2(ffdata)
      pred1 = res$allmean; E1[s-m+1] = inprod(fundata[s+1]-pred1, fundata[s+1]-pred1); 
      pred2 = res$last; E2[s-m+1] = inprod(fundata[s+1]-pred2, fundata[s+1]-pred2)
      pred3 = res$wmean; E3[s-m+1] = inprod(fundata[s+1]-pred3, fundata[s+1]-pred3)
      pred4 = res$ffpe; E4[s-m+1] = inprod(fundata[s+1]-pred4, fundata[s+1]-pred4)
   }
   c(mean(E1), mean(E2), mean(E3), mean(E4))
}


R = sapply(1:100, function(j) simulation(n=200))
 




