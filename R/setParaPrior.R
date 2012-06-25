# Default prior values for PING
#priottype=1 means use the new prior regulate both forward/reverse peaks, otherwise old prior is used to regulate sum of precision of peaks
#paraPriorH<-list(xi=200,rho=1,alpha=20,beta=40000,lambda=-0.00001,dMu=200,priortype=1)
#paraPriorTF<-list(xi=200,rho=1,alpha=20,beta=40000,lambda=0,dMu=0,priortype=1)
#I changed the value hyperparameter, 
# If both forward/reverse precision become centered at previous center of sum precision, 
# then rho should change from 1 to 0.5 to keep variance of delta not change.
# The previous variance of precision was too small, I change alpha, beta to make it value of sigma_f and sigma_r center @ 50, and CI= c(1,100)


##I changed the value of beta, since now the gamma distribution is for precision of F or R reads, but not their sum.
# following is 95% bound of sigma based on Gamma(20,20000) and its median
##MNase
#> (qgamma(c(0.025,0.5,0.975),20,20000))^(-.5)
#[1] 40.46143 31.88882 25.96271
##Sonication
#> (qgamma(c(0.025,0.5,0.975),10,20000))^(-.5)
#[1] 64.58075 45.48107 34.21448

# If two sigma on each side, the F/R peak is about 100 to 160 bps wide, and median is 120 bps wide.
## I use rho=0.8 so that 95% bound and median of varaince of delta is given as follows
##MNase
#> sqrt((3)^(-1)/(qgamma(c(0.025,0.5,0.975),20,20000)*2))
#[1] 16.51831 13.01856 10.59923
##Sonication
#>  sqrt((1.2)^(-1)/(qgamma(c(0.025,0.5,0.975),10,20000)*2))
#[1] 41.68670 29.35790 22.08535
##According to ovservation of estimated delta, we find they are centered at 150 about the length of nucleosome, 
# so we change xi=160 for MNase data, we might use larger value of xi for sonication data
# using this xi and variance
# wide peak can have delta 95% range as 98-222, and narrow peak can has "delta" 95% range 60-260


#paraPriorH<-list(xi=150,rho=0.8,alpha=20,beta=20000,lambda=-0.000064,dMu=200) # for MNase data
paraPriorH<-list(xi=150,rho=1.2,alpha=10,beta=20000,lambda=-0.000064,dMu=200) # for sonication data
paraPriorTF<-list(xi=200,rho=0.5,alpha=20,beta=20000,lambda=0,dMu=0)


setParaPrior<-function(xi=150,rho=1.2,alpha=10,beta=20000,lambda=-0.000064,dMu=200)
{
  dataType="H"
  if(dataType!="TF" & dataType!="H")
  {
    stop("Object 'dataType' must be either 'TF' or 'H'")
  }
  if(!is.finite(xi))
  {
    stop("'xi' must be a numeric value")
  }
  if(!is.finite(rho) & rho<=0)
  {
    stop("'rho' must be a positive number")
  }
  if(!is.finite(alpha) & alpha<=0)
  {
    stop("'alpha' must be a positive number")
  }
  if(!is.finite(beta) & beta<=0)
  {
    stop("'beta' must be a positive number")
  }
  if(!is.finite(lambda) & lambda<=0)
  {
    stop("'lambda' must be a positive number")
  }
  if(!is.finite(dMu) & dMu<=0)
  {
    stop("'dMu' must be a positive number")
  }
  
  if(dataType=="TF")
  {
  unlockBinding("paraPriorTF", environment(PING))
  assign("paraPriorTF",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=environment(PING))
  assign("paraPriorTF",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=.GlobalEnv)
  lockBinding("paraPriorTF", environment(PING))
  }
  else if(dataType=="H")
  {
    if(!(lambda<0))
    {
      warning("You have selected the Histone option but lambda=0")
    }
    unlockBinding("paraPriorH", environment(PING))
    assign("paraPriorH",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=environment(PING))
    assign("paraPriorH",list(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu), envir=.GlobalEnv)
    lockBinding("paraPriorH", environment(PING))
  }
}
