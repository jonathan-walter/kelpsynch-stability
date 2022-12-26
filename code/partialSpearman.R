## Compute partial spearman correlation sensu Ghosh et al. 2020 Ecosphere

#u and v are numeric vectors
#range is a length-2 numeric vector giving the range of quantiles over which to
#compute the partial spearman correlation

partialSpearman<-function(u, v, bounds){
  
  if(!is.numeric(u) | !is.numeric(v)){
    stop("u and v must be numeric")
  }
  if(length(u)!=length(v)){
    stop("u and v have different length")
  }
  if(min(bounds)<0 | max(bounds)>1){
    stop("bounds must be between 0 and 1")
  }
  
  n<-length(u)
  x<-rank(u)/(n+1)
  y<-rank(v)/(n+1)
  
  x_b<-x[x+y > (2*min(bounds)) & x+y < (2*max(bounds))]
  y_b<-y[x+y > (2*min(bounds)) & x+y < (2*max(bounds))]
  
  r<-sum((x_b-mean(x))*(y_b-mean(y)))/((n-1)*sqrt(var(x)*var(y)))
  
  return(r)
}

#tests

# library(mvtnorm)
# rho=0.1
# sigmat<-matrix(c(1,rho,rho,1),2,2)
# xx<-rmvnorm(1000, sigma=sigmat)
# u=xx[,1]
# v=xx[,2]
# 
# cor(u,v,method="spearman")
# partialSpearman(u,v,c(0,1))
# partialSpearman(u,v,c(0,0.5))
# partialSpearman(u,v,c(0.5,1))