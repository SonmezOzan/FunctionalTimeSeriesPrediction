
library(fda)
library(vars)


D=21
BASIS = create.fourier.basis(rangeval = c(0,1), nbasis=D)



fun.prediction = function(scores, d, p, eigenfun){
	n = length(scores)

	if(d==1)
         {
            pred = predict(arima(scores, order=c(p,0,0)), n.ahead=1)$pred[1]
         }
         else
         {
            if(p==0)
            {
               pred = t(matrix(rep(colMeans(scores),n),d))
            }
            else
            {
               n = nrow(scores)
               # we have to give colnames to the scores, 
               # otherwise we get warnings from vars resid function
               colnames(scores) = as.character(seq(1:d))
               VAR.pre= predict(VAR(scores, p), n.ahead=1, type="const")$fcst
               pred = sapply(1:d, function(i) VAR.pre[[i]][1])
            }
        }

        f.pred = eigenfun[1] * 0
         for (i in 1:d){
         	f.pred = f.pred + pred[i] * eigenfun[i]
         }
         f.pred
}




permute.pred = function(fdata, Pmax=4){

	 n = ncol(fdata$coefs)
     D = nrow(fdata$coefs)
    
     m = n/10 

	 P = matrix(0, D, n)
     P[, 1] = c(1:D)
     for (i in 2:n){
	     P[, i] = sample(c(1:D), D)
        }

       M = 100

       Res = matrix(0, 4, M)

      FPCA = pca.fd(fdata, nharm=D, centerfns=TRUE)
      scores_0 = FPCA$scores
      v.hat = FPCA$harmonics
      Lambda = FPCA$values

      

     for (per in 1:M){
     	  values = matrix(0,D,(Pmax+1))
	      Score = scores_0[,P[,per]]
	      v_hat = v.hat[P[,per]]

	        for (d in 1:D){

		      score.d = as.matrix(Score[, 1:d])

		         for (p in 0:Pmax){
            
                    Error = c()
                    for (ell in (n-m):(n-1)){
                       pred.m = fun.prediction(score.d[1:ell, ], d, p, eigenfun=v_hat)
                       err = inprod(fdata[ell+1]-pred.m, fdata[ell+1]-pred.m)
                       Error[ell-(n-m)+1] = err
                           }
                    values[d, p+1] = mean(Error)
		        }
	        }
         hat.p = (which.min(values)-1)%/%D
         hat.d = which.min(values)%%D
         if(hat.d==0) hat.d = hat.d + D
         Res[1, per] = per
         Res[2, per] = min(values)
         Res[3, per] = hat.p
         Res[4, per] = hat.d
        }
        Res
}
         




per = 78
d = 6
p = 4
Score = scores_0[,P[,per]]
v_hat = v.hat[P[,per]]
scores = Score[, 1:d]
colnames(scores) = as.character(seq(1:d))
VAR.pre= predict(VAR(scores, p), n.ahead=1, type="const")$fcst
pred = sapply(1:d, function(i) VAR.pre[[i]][1])
f.pred = v_hat[1] * 0
         for (i in 1:d){
         	f.pred = f.pred + pred[i] * v_hat[i]
         }
ffpe.pred = prediction.ffpe(fdata)

plot(fdata, col="grey")
lines(f.pred, col="red")
lines(ffpe.pred, col="blue")
lines(fdata[100], col="black")

