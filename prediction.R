library(vars)

prediction.new = function(ffdata){
   
   n = ncol(ffdata$coef)
   D = nrow(ffdata$coef)
   mu = mean.fd(ffdata)
   fdata = center.fd(ffdata)
   
   n.p = n/10
   
   #fPCA
   fpca = pca.fd(fdata, nharm=D, centerfns=TRUE)
   v.hat = fpca$harmonics
   
   Qmax = 5
   Dmax = 10
   Pmax=5
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
   
   d = di.s[1]
   
   scores = as.matrix(fpca$scores[,1:d]) # full projected data
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
   X.mv = matrix(0, d, 1)
   for (i in 1:Pmax){
      V = beta.hat[[i]] %*% reduction(scores[n-i+1,],di.s[i])
      X.mv = X.mv + V
   }
   X.f = 0*v.hat[1]
   for (el in 1:d){
      X.f = X.f + X.mv[el] * v.hat[el]
   }
   return(X.f+mu)
}
   
   
prediction.ffpe = function(fdata){
   n = ncol(fdata$coef)
   D = nrow(fdata$coef)
   
   #center the data
   mu = mean(fdata)
   data = center.fd(fdata)
   
   # ffpe Criterion
   ffpe = fFPE(fdata, Pmax=5)
   d.hat = ffpe[1]
   p.hat = ffpe[2]
   
   #fPCA
   fpca = pca.fd(data, nharm=D, centerfns=TRUE)

      scores = fpca$scores[, 1:d.hat]
         if(d.hat==1)
         {
            VAR.pre = predict(arima(scores, order=c(p.hat,0,0)), n.ahead=1)$pred[1]
            yhat = VAR.pre
         }
         else
         {
            if(p.hat==0)
            {
               yhat = colMeans(scores)
            }
            else
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
      far.pre = FAR.pre + mu
      return(far.pre)      
}
   
   