
project = function(fdata, d){
   
   mu = mean(fdata)
   cdata =center.fd(fdata)
   n = ncol(fdata$coef)
   D = nrow(fdata$coef)
   
   fpca = pca.fd(fdata, nharm=D, centerfns=TRUE)
   v = fpca$harmonics
   par(mfrow=c(1,2))
   plot(fdata, col="black",ylim=c(-4,4))
   plot(10,ann=FALSE, ylim=c(-4,4), xlim=c(0,1))
   for (i in 1:n){
      y = 0*v[1]
      for (j in 1:d){
         a = inprod(fdata[i], v[j])*v[j]
         y = y +a + mu
      }
      lines(y)
   }
}


library(fds)
library(fda)
dat = Australiafertility$y
D=21
BASIS = create.fourier.basis(rangeval = c(0,1), nbasis=D)
fdata = Data2fd(argvals=seq(0, to=1, length.out = 35) , y=dat, basisobj = BASIS)


dat = ECBYieldcurve$y
fdata = Data2fd(argvals=seq(0, to=1, length.out = nrow(dat)) , y=dat, basisobj = BASIS)


dat = Fatspectrum$y
fdata = Data2fd(argvals=seq(0, to=1, length.out = nrow(dat)) , y=dat, basisobj = BASIS)
