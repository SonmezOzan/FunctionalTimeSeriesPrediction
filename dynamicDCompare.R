library(fda)
library(vars)

D=21
BASIS = create.fourier.basis(rangeval = c(0,1), nbasis=D)

#---------------------------------- FAR(1) process -----------------------------------------------

rarh1 = function(n , Sigma){
   
   #burnin
   l_burnin = n/2
   
   data = matrix(0, D, (n+l_burnin))
   for (i in 1:(n+l_burnin)){
      data[,i] = rnorm(D, 0, Sigma)
   }
   
   # Psi operator
   v1 = rnorm(D, 0, Sigma)
   v2 = rnorm(D, 0, Sigma)
   Psi = v1 %*% t(v2)
   
   Psi =0.8* Psi/(base::norm(Psi, type="2"))
   
   
   coef = matrix(0, D, (n+l_burnin))
   coef[ ,1] = data[,1]
   
   # recursion
   for (i in 2:(n+l_burnin)){
      coef[,i] = Psi %*% coef[ ,i-1]  + data[,i]
   }
   dat = coef[, (l_burnin+1):(n+l_burnin)]
   fd(dat, BASIS)
}


reduction = function(data, d){
   dd = length(data)
   if (d==0){
      out = rep(0, dd)
   } else {
      out = c(data[1:d], rep(0, dd-d))
   }
   out
}

g = function(x, y){exp(-x*y)}


# --------------------------------------------------------------------------
# --------------------------- fFPE method ----------------------------------
# --------------------------------------------------------------------------

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

fFPE = function(fdata, Pmax, Dmax=10){
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
   values = matrix(0,Dmax,(Pmax+1))
   for (d in 1:Dmax){
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
               res = resid(VAR(scores, p=p, type="const"))
            }
         }
         values[d, p+1] = trace(res = res, p = p) + vartot - var.explain
      }
   } 
   # compute the estimates hat.p and hat.d for optimal order p and dimension d
   hat.p = (which.min(values)-1)%/%Dmax
   hat.d = which.min(values)%%Dmax
   if(hat.d==0) hat.d = hat.d + Dmax
   return(c(hat.d, hat.p))
}



# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# --------------------------------------------------------------------------


simulation = function(n, Pmax, Dmax=10, sigma){
   
   ffdata = rarh1(n=n, Sigma=sigma)
   fdata = center.fd(ffdata)
   
   n.p = n/10
   
   ffpe = fFPE(fdata[1:(n-n.p)], Pmax)
   d.ffpe = ffpe[1]
   p.ffpe = ffpe[2]
   
   #fPCA
   fpca = pca.fd(fdata, nharm=D, centerfns=TRUE)
   v.hat = fpca$harmonics
   
   Qmax = 5
   alpha = seq(0, Qmax, 0.5)
   ell = seq(0 , 1, length=Pmax+1)
   values = c()
   
   d = d.ffpe
      full.scores = as.matrix(fpca$scores[,1:d]) # full projected data
      scores = full.scores[1:(n-n.p),] # training set for estimation
      beta.hat = list()
      if(d==1)
      {
         ar = arima(scores, order = c(Pmax, 0, 0), method="ML")
         betas = ar$coef
         for (j in 1:Pmax){
            beta.hat[[j]] = betas[((j-1)*d+1):(j*d)]
         }
      } else {
         # give colnames to the scores to avoid warnings 
         colnames(scores) = as.character(seq(1:d))
         var.p = VAR(scores, p=Pmax, type="const")
         betas = t(sapply(1:d, function(k) coef(var.p)[[k]][,1]))
         for (j in 1:Pmax){
            beta.hat[[j]] = betas[, ((j-1)*d+1):(j*d)]
         }
      }
      for (l in 1:length(alpha)){
         d.prime = round(g(alpha[l], ell)*d)
         
         # accumulation
         error = 0
         for (m in (n-n.p):(n-1)){
            X.mv = matrix(0, d, 1)
            for (i in 1:Pmax){
               V = beta.hat[[i]] %*% reduction(full.scores[m-i+1,],d.prime[i])
               X.mv = X.mv + V
            }
            X.f = 0*v.hat[1]
            for (el in 1:d){
               X.f = X.f + X.mv[el] * v.hat[el]
            }
            eps = inprod(fdata[m+1]-X.f, fdata[m+1]-X.f)
            error = error + eps 
         }
         values[l] = error
      }

   d.hat = d
   alpha.hat = min(which(values==min(values)))
   di.s = round(g(alpha[alpha.hat], ell[-1])*d.hat)
   error.new = min(values)/n.p
   
   
   # Prediction error based on fFPE criteria
   beta.hat.ffpe = list()
   scores.ffpe = as.matrix(fpca$scores[,1:d.ffpe])
   if(p.ffpe==0)
   {
      err.ffpe = sum(sapply((n-n.p):(n-1), function(m) inprod(fdata[m+1]-mean(fdata[1:m]),fdata[m+1]-mean(fdata[1:m]))))/n.p
   } else {
      
      if(d.ffpe==1)
      {
         ar.ffpe = arima(scores.ffpe, order = c(p.ffpe, 0, 0), method="ML")
         betas.f = ar.ffpe$coef
         for (j in 1:p.ffpe){
            beta.hat.ffpe[[j]] = betas.f[((j-1)*d.ffpe+1):(j*d.ffpe)]
         }
      } else {
         # give colnames to the scores to avoid warnings 
         colnames(scores.ffpe) = as.character(seq(1:d.ffpe))
         var.p.f = VAR(scores.ffpe, p=p.ffpe, type="const")
         betas.f = t(sapply(1:d.ffpe, function(k) coef(var.p.f)[[k]][,1]))
         for (j in 1:p.ffpe){
            beta.hat.ffpe[[j]] = betas.f[, ((j-1)*d.ffpe+1):(j*d.ffpe)]
         }
      }
      error.f = 0
      for (m in (n-n.p):(n-1)){
         X.mv = matrix(0, d.ffpe, 1)
         for (i in 1:p.ffpe){
            V = as.matrix(beta.hat.ffpe[[i]]) %*% as.matrix(fpca$scores[m-i+1,1:d.ffpe])
            X.mv = X.mv + V
         }
         X.ff = 0*v.hat[1]
         for (el in 1:d.ffpe){
            X.ff = X.ff + X.mv[el] * v.hat[el]
         }
         eps = inprod(fdata[m+1]-X.ff, fdata[m+1]-X.ff)
         error.f = error + eps 
      }
      err.ffpe = error.f/n.p}
   list(d.new = di.s, error.new=error.new, error.ffpe=err.ffpe, ffpe.d = d.ffpe, ffpe.p=p.ffpe)
}


# Slow decay n=100
N=100
sigma = (1:D)^-1

ds = matrix(0, 1000, 5)
er.ffpe = er.new = d.hat = p.hat = c()
for (M in 1:1000){
   result = simulation(n=N, Pmax=5, sigma=sigma)
   ds[M,] = result$d.new
   er.ffpe[M] = result$error.ffpe
   er.new[M] = result$error.new
   d.hat[M] = result$ffpe.d
   p.hat[M] = result$ffpe.p
}

save(ds, er.ffpe, er.new, d.hat, p.hat , file="fAR.Slow.100.Compare.rda")


