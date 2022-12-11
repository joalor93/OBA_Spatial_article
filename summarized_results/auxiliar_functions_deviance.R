#### auxiliar function
deviance = function(theta,coords,X,y,cov.model = "matern",kappa = 0.5){
  p = ncol(X)
  beta = theta[1:p]
  cov.pars = theta[(p+1):(length(theta)-1)]
  nu = theta[length(theta)]
  covmatrix = varcov.spatial(coords = coords,cov.model = cov.model,
                             kappa = kappa,nugget =0,cov.pars=cov.pars)
  covmatrix = covmatrix$varcov
  predlin = X%*%beta
  
  y = as.matrix(y)
  dens = dmvt(t(y),mu = t(predlin), S = covmatrix, df =nu,log = TRUE)
  return(-2*dens)
}


DICinf= function(obj){
  X = obj$X
  theta = obj$theta
  y = obj$y
  kappa =obj$kappa
  coords = obj$coords
  
  vec = 0
  for( i in 1: dim(obj$dist)[1]){
    vec[i] = deviance(coords = coords,X = X,theta =obj$dist[i,],y = y,cov.model = "matern",kappa = kappa)
    print(i)
  }
  
  #result = apply(obj$dist,MARG=1,FUN=deviance,coords = coords,X = X,y = y,cov.model = "matern",kappa = kappa)
  
  dev = deviance(coords = coords,X = X,theta = theta,y = y,cov.model = "matern",kappa = kappa)
  expecteddev = mean(vec)
  return(list(dev = dev,expdev=expecteddev,DIC=2*expecteddev-dev))
}
