
library(vars)

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

model.select = function(ffdata, Qmax=5, Dmax=10, Pmax){
   
   n = ncol(ffdata$coefs)
   D = nrow(ffdata$coefs)
   fdata = center.fd(ffdata)
   
   n.p = n/10
   
   #fPCA
   fpca = pca.fd(fdata, nharm=D, centerfns=TRUE)
   v.hat = fpca$harmonics
   
   alpha = seq(0, Qmax, 0.5)
   ell = seq(0 , 1, length=Pmax)
   values = matrix(0, Dmax, length(alpha))
   
   for (d in 1:Dmax){
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
         values[d, l] = error
      }
   }
   res = which(values == min(values), arr.ind = TRUE)[1,]
   d.hat = res[1]
   alpha.hat = res[2]
   di.s = round(g(alpha[alpha.hat], ell)*d.hat)
   return(di.s)
}
   

   






