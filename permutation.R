
library(fda)

Cov.h = function(X, h){
       n = ncol(X)
       D = nrow(X)
       A = matrix(0, D, D)
       for (i in 1:(n-h)){
            A = A + X[,i] %*% t(X[,i+h])
          }
      E = eigen(A/n)
       Lam = E$value
       V = E$vector
       list(vector = V, value = Lam)
}

 C_h = function(data, h, x){
          D = length(x)
          CC = Cov.h(data, h)
          lambda = CC$value
          v_ell = CC$vector
          a = sum(sapply(1:D, function(l) lambda[l]^2*(t(x)%*%v_ell[,l])^2))
          a
 }
 g = function(x, y){exp(-x*y)}
 reduction = function(data, d){
    dd = length(data)
    if (d==0){
       out = rep(0, dd)
    } else {
       out = c(data[1:d], rep(0, dd-d))
    }
    out
 }
 
 # ---------------------------------------------------------------------------------------------------------------------
 # ---------------------------------------------------------------------------------------------------------------------
 # ---------------------------------------------------------------------------------------------------------------------
 # ---------------------------------------------------------------------------------------------------------------------
 # ---------------------------------------------------------------------------------------------------------------------

 Predict_new = function(fdata, Pmax=5 , Qmax=5, Dmax=10 ){
    # desciptive st, and center the functional time series
    n = ncol(fdata$coef)
    D = nrow(fdata$coef)
    mu = mean.fd(fdata)
    fdata = center.fd(fdata)
    
    # size for for APE
    n.p = n/10
    
    # fPCA
    fpca = pca.fd(fdata, nharm=D, centerfns=TRUE)
    v.hat = fpca$harmonics
    v_0 = v.hat$coefs
    lambda_0 = fpca$values
    scores = t(fpca$scores)
    
    alpha = seq(0, Qmax, 0.5)
    ell = seq(0 , 1, length=Pmax)
    values = matrix(0, Dmax, length(alpha))
    
    # matrix of computed covariance at lag of eigenfunctions ||C_h(v_l)||
    res = matrix(0, Pmax, D)
    for (h in 1:Pmax){
       for (dn in 1:D){
          res[h,dn]=C_h(fdata$coefs, h, v_0[,dn])
       }
    }
    # rank of res
    A = matrix(0, Pmax,D)
    for (i in 1:Pmax){
       A[i, ] = rank(-res[i, ])
    }
    # which eigenfunction to use, rows=lags, columns = chosen eigenfunctions
    B = A
    for (i in 1:Pmax){
       for (dd in 1:D){
          B[i,dd] = which(A[i,]==dd)
       }
    }
    chosen.scores = list()
    for (d in 1:Dmax){
       Proj = t(B[,1:d])
       for (i in 1:ceiling(n)){
          Proj = cbind(Proj, t(B[,1:d]))
       }
       
       # eigenfunctions to use at each projection
       if (d==1){new.scores = t(as.matrix(scores[1:d,]))}else{
          new.scores = scores[1:d,]
       }
       
       for (j in 1:n){
          for (i in 1:d){
             new.scores[i,j] = scores[Proj[i,j],j]
          }
       }
       #full.scores = as.matrix(fpca$scores[,1:d]) # full projected data
       training.scores = new.scores[,1:(n-n.p)] # training set for estimation
       beta.hat = list()
       if(d==1)
       {
          ar = arima(training.scores, order = c(Pmax, 0, 0), method="ML")
          betas = ar$coef
          for (j in 1:Pmax){
             beta.hat[[j]] = betas[((j-1)*d+1):(j*d)]
          }
       } else {
          # give colnames to the scores to avoid warnings 
          rownames(training.scores) = as.character(seq(1:d))
          var.p = VAR(t(training.scores), p=Pmax, type="const")
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
             X.mv = 0*fdata[1]
             for (i in 1:Pmax){
                V = beta.hat[[i]] %*% reduction(new.scores[,m-i+1],d.prime[i]) 
                X.f = 0*v.hat[1]
                for (el in 1:d){
                   X.f = X.f + V[el] * v.hat[Proj[el,i]]
                }
                X.mv = X.mv + X.f
             }
             eps = inprod(fdata[m+1]-X.mv, fdata[m+1]-X.mv)
             error = error + eps 
          }
          values[d, l] = error
       }
       chosen.scores[[d]] = new.scores
    }
    res = which(values == min(values), arr.ind = TRUE)[1,]
    d.hat = res[1]
    alpha.hat = res[2]
    di.s = round(g(alpha[alpha.hat], ell)*d.hat)
    directions = data.frame(B)
    colnames(directions) = 1:D
    
    
    ######################################################################
    # Prediction:
    d = di.s[1]
    
    Scores = chosen.scores[[d]] # full projected data
    beta.hat = list()
    if(d==1)
    {
       ar = arima(t(Scores), order = c(Pmax, 0, 0), method="ML")
       betas = ar$coef
       for (j in 1:Pmax){
          beta.hat[[j]] = betas[((j-1)*d+1):(j*d)]
       }
    } else {
       # give colnames to the scores to avoid warnings 
       rownames(Scores) = as.character(seq(1:d))
       var.p = VAR(t(Scores), p=Pmax, type="const")
       betas = t(sapply(1:d, function(k) coef(var.p)[[k]][,1]))
       for (j in 1:Pmax){
          beta.hat[[j]] = betas[, ((j-1)*d+1):(j*d)]
       }
    }
    
    for (i in 1:Pmax){
       V = beta.hat[[i]] %*% reduction(Scores[,n-i+1],di.s[i]) 
       X.f = 0*v.hat[1]
       for (el in 1:d){
          X.f = X.f + V[el] * v.hat[Proj[el,i]]
       }
       X.mv = X.mv + X.f
    }
    X.hat = X.mv+mu
    out = list(di = di.s,
               Proj.Matrix = directions,
               pred = X.hat)
    return(out)
 }
 

