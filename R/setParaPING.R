#setParaPrior with PING default values
setParaPriorPING<-function(xi=150,rho=0.8,alpha=20,beta=20000,lambda=-0.000064,dMu=200)
{
  #default for PING MNase
  #(xi=150,rho=0.8,alpha=20,beta=20000,lambda=-0.000064,dMu=200)
  #default for PING sonication
  #(xi=150,rho=1.2,alpha=10,beta=20000,lambda=-0.000064,dMu=200)

  #Call to function imported from PICS
  paraList<-setParaPrior(xi=xi,rho=rho,alpha=alpha,beta=beta,lambda=lambda,dMu=dMu)
  return(paraList)
}

#setParaEM with PING default values
setParaEMPING<-function(minK=0,maxK=0,tol=1e-4,B=100,mSelect="AIC3",mergePeaks=TRUE,mapCorrect=TRUE)
{
  #Call to function imported from PICS
  paraList<-setParaEM(minK=minK, maxK=maxK, tol=tol, B=B, mSelect=mSelect, mergePeaks=mergePeaks, mapCorrect=mapCorrect)
  return(paraList)
}
