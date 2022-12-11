##########################################################################################################
###################### First and second derivative of spatial correlation function #######################
##########################################################################################################

derivcormatrix1<- function(H, phi, kappa=0, type="exponential"){

  if (type=="exponential"){
    H1 <- (abs(H)/phi^2)*exp(-(abs(H)/phi))
    H2 <- abs(H)*(abs(H)-2*phi)*exp(-(abs(H)/phi))/(phi^4)
  }

  if (type=="gaussian"){
    H1 <- (2*abs(H)^2/phi^3)*exp(-(abs(H)/phi)^2)
    H2 <- (4*abs(H)^4 - 6*abs(H)^2*phi^2)*exp(-(abs(H)/phi)^2)/(phi^6)
  }

  if (type=="matern"){
    H[H==0]<-1
    Ak <- besselK(abs(H)/phi,(kappa-1)) + besselK(abs(H)/phi,(kappa+1))
    Bk <- besselK(abs(H)/phi,(kappa-2)) + 2*besselK(abs(H)/phi,kappa) + besselK(abs(H)/phi,(kappa+2))
    H1 <- -1/((2^kappa)*(phi^2)*gamma(kappa))*(abs(H)/phi)^kappa*(2*kappa*phi*besselK(abs(H)/phi,kappa) - abs(H)*Ak)
    H2 <- (abs(H)^kappa)/(2^(kappa+1)*gamma(kappa)*phi^(kappa+4))*(4*kappa*(kappa+1)*phi^2*besselK(abs(H)/phi,kappa) - 4*(kappa+1)*phi*abs(H)*Ak + abs(H)^2*Bk)
  }

  if (type=="pow.exp"){
    H1 <- (kappa/phi)*(abs(H)/phi)^(kappa)*exp(-(abs(H)/phi)^(kappa))
    H2 <- H1*(kappa*abs(H)^kappa/(phi^(kappa+1)) - (kappa+1)/phi)
  }

  if (type=="spherical"){
    H1 <- 1.5*(abs(H)/phi^2) - 1.5*(abs(H)^3/phi^4)
    Haux <- (abs(H)>phi) + 0
    H1[Haux==1] <- 0
    H2 <- 6*(abs(H)^3)/(phi^5) - 3*abs(H)/(phi^3)
    H2[Haux==1] <- 0
  }
  diag(H1) <- 0 # First derivative correlation matrix
  diag(H2) <- 0 # Second derivative correlation matrix

  return(list(dev1=H1, dev2=H2))
}


########################################################
####### Ccovariance Matrix##############################
########################################################



varcov.spatial=function(coordsi=F,coords,H,cov.model="cauchy",cov.pars,nugget,kappa){

phi=cov.pars[2]
sigma=cov.pars[1]

if(coordsi==T){
H=as.matrix(dist(coords,upper=T,diag=T))
}

#matern
if(cov.model=="matern"){
if(kappa==.5){
cov.model="exponential"
}
else{
c1=(0.5^(kappa-1))/gamma(kappa)
c2=(abs(H)/phi)^kappa
c3=besselK(abs(H)/phi,kappa)
R=c1*c2*c3
diag(R)=1
}
}
#power exponential
if(cov.model=="power.exponential"||cov.model=="pow.exp"){
if(kappa==1){
cov.model="exponential"
}

else{
R=exp(-(abs(H)/phi)^kappa)
}
}

#exponential
if(cov.model=="exponential"||cov.model=="exp"){
R=exp(-abs(H)/phi)
}

#spherical
if(cov.model=="spherical"){

A=abs(H)/phi
R=1-(1.5*A)+(0.5*(A^3))
 Raux <- (abs(H)>phi) + 0
    R[Raux==1] <- 0
}


#Cauchy
if(cov.model=="cauchy"){
A=abs(H)/phi
R=(1+ A^2)^(-kappa)
}
m=dim(R)[1]
varcova= sigma*R + nugget*diag(m)
return(varcova)
}


###############################################################
######## Encontrando o valor de phi############################
###############################################################


find.phi <- function(d, kappa=0.5,range=c(1e-04,100),cut=0.05,cov.model){
if(cov.model=="matern"){
if(kappa==0.5) cov.model=="exponential"

f <- function(x,d,kappa,cut){
     out <- 1/(2^(kappa-1)*gamma(kappa)) * (d/x)^kappa * besselK(d/x,kappa)
     out <- out - cut
   }
}

#power exponential
if(cov.model=="power.exponential"||cov.model=="pow.exp"){
if(kappa==1) cov.model=="exponential"
f <- function(x,d,kappa,cut){
     out <- exp(-(abs(d)/x)^kappa)
     out <- out - cut
   }

}

if(cov.model=="cauchy"){
f <- function(x,d,kappa,cut){
 A=abs(d)/x
 out <-(1+ A^2)^(-kappa)
   out <- out - cut
}
}

if(cov.model=="exponential"){
f <- function(x,d,kappa,cut){
out= exp(-abs(d)/x)
 out <- out - cut
}
}

  out <- uniroot(f=f,interval=range,d=d,kappa=kappa,cut=cut)
   out$root
}



#######################################################
#########densidade priori cjta (phi,nu)=theta1#########
#######################################################

prioriref2=function(sigma,x,H,kappa,cov.model,tau2,xmat){
covpars=c(sigma,x[1])
p=ncol(xmat)
ndata=nrow(xmat)
nu1=ndata-p
nucand=x[2]
cand=covpars
sigma=covpars[1]
phi=covpars[2]
phicand=phi
covacand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
covacand=covacand+(10e-10*diag(1,ndata))

derivcovacand=sigma*derivcormatrix1(H,phi=phicand,kappa=kappa,type=cov.model)$dev1
covacandinv=solve(covacand,tol=10e-20)
Pcand=diag(ndata)-(xmat%*%solve(t(xmat)%*%covacandinv%*%xmat)%*%t(xmat)%*%covacandinv)
Wcand=derivcovacand%*%covacandinv%*%Pcand

###comp sigma2
acand=(nu1*nucand)/(nu1+nucand+2)

b11cand= (nucand*(nu1))/(nucand-2)
b12cand=sum(diag(Wcand))
b13cand=((nucand+nu1)/(nucand+nu1+2))-1


#####comp phi sigma
b1cand=b11cand*b12cand*b13cand


aest11cand=(nucand^2)/((nucand-2)*(nucand-4))
aest12cand= 2*sum(diag(Wcand%*%Wcand))
aest13cand=sum(diag(Wcand))^2
aestcand= aest11cand*(aest12cand+aest13cand)

cons1cand=(nu1)/(nucand*(nucand+nu1+2))
cons2cand=(nucand+2)/(nucand-2)

####comp phi2

aestcand2= ((cons1cand+1)*aestcand) - (cons2cand*aest13cand)



#c1cand=nu1*(nucand+ndata)/(2*nucand*(nucand+nu1))
#c2cand= log((nucand+nu1)/nucand)/2
c3cand= trigamma((nucand+nu1)/2)-trigamma(nucand/2)
c4cand= ((2*nu1)/(nucand))*((nucand+nu1+4)/((nucand+nu1+2)*(nucand+nu1)))



###esperanca nu2######

ccand=-c3cand-c4cand

### comp nu sigma

ba12cand= -(nu1)/((nucand+nu1+2)*(nucand+nu1))


####

#d11cand=(1/(nucand-2))*(1+ log((nucand+nu1)/nucand))
d12cand=(nu1)/((nu1+nucand)*(nu1+nucand+2)*(nucand-2))

##### comp nu phi

dcand= d12cand*b12cand

############

#res1cand= acand*aestcand2*ccand + b1cand*dcand*ba12cand
#res2cand =(ba12cand^2)*aestcand2 + acand*(dcand^2)+ (b1cand^2)*ccand
#rescand=res1cand- 0.5*res2cand

res1cand= acand*aestcand2*ccand + (16*b1cand*dcand*ba12cand)
res2cand =(8*(ba12cand^2)*aestcand2) + (16*acand*(dcand^2))+ (0.5*(b1cand^2)*ccand)
rescand=res1cand-res2cand


priorcand=sqrt(rescand)


return(priorcand)
}





#######################################################
#########densidade priori cjta (phi,nu)=theta1#########
#######################################################

prioriref2mat=function(covacand,covacandinv,vbetacand,phi,nu,H,kappa,cov.model,tau2,xmat){
  #covpars=c(sigma,x[1])
  p=ncol(xmat)
  ndata=nrow(xmat)
  nu1=ndata-p
  nucand=nu
  sigma=1
  #cand=covpars
  #sigma=covpars[1]
  #phi=covpars[2]
  phicand=phi
  #covacand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
  #covacand=covacand+(10e-10*diag(1,ndata))

  derivcovacand=sigma*derivcormatrix1(H,phi=phicand,kappa=kappa,type=cov.model)$dev1
  #covacandinv=solve(covacand,tol=10e-20)
  #Pcand=diag(ndata)-(xmat%*%solve(t(xmat)%*%covacandinv%*%xmat)%*%t(xmat)%*%covacandinv)
  Pcand=diag(ndata)-(xmat%*%vbetacand%*%t(xmat)%*%covacandinv)
  Wcand=derivcovacand%*%covacandinv%*%Pcand

  ###comp sigma2
  acand=(nu1*nucand)/(nu1+nucand+2)

  b11cand= (nucand*(nu1))/(nucand-2)
  b12cand=sum(diag(Wcand))
  b13cand=((nucand+nu1)/(nucand+nu1+2))-1


  #####comp phi sigma
  b1cand=b11cand*b12cand*b13cand


  aest11cand=(nucand^2)/((nucand-2)*(nucand-4))
  aest12cand= 2*sum(diag(Wcand%*%Wcand))
  aest13cand=sum(diag(Wcand))^2
  aestcand= aest11cand*(aest12cand+aest13cand)

  cons1cand=(nu1)/(nucand*(nucand+nu1+2))
  cons2cand=(nucand+2)/(nucand-2)

  ####comp phi2

  aestcand2= ((cons1cand+1)*aestcand) - (cons2cand*aest13cand)



  #c1cand=nu1*(nucand+ndata)/(2*nucand*(nucand+nu1))
  #c2cand= log((nucand+nu1)/nucand)/2
  c3cand= trigamma((nucand+nu1)/2)-trigamma(nucand/2)
  c4cand= ((2*nu1)/(nucand))*((nucand+nu1+4)/((nucand+nu1+2)*(nucand+nu1)))



  ###esperanca nu2######

  ccand=-c3cand-c4cand

  ### comp nu sigma

  ba12cand= -(nu1)/((nucand+nu1+2)*(nucand+nu1))


  ####

  #d11cand=(1/(nucand-2))*(1+ log((nucand+nu1)/nucand))
  d12cand=(nu1)/((nu1+nucand)*(nu1+nucand+2)*(nucand-2))

  ##### comp nu phi

  dcand= d12cand*b12cand

  ############

  #res1cand= acand*aestcand2*ccand + b1cand*dcand*ba12cand
  #res2cand =(ba12cand^2)*aestcand2 + acand*(dcand^2)+ (b1cand^2)*ccand
  #rescand=res1cand- 0.5*res2cand

  res1cand= acand*aestcand2*ccand + (16*b1cand*dcand*ba12cand)
  res2cand =(8*(ba12cand^2)*aestcand2) + (16*acand*(dcand^2))+ (0.5*(b1cand^2)*ccand)
  rescand=res1cand-res2cand


  priorcand=sqrt(rescand)


  return(priorcand)
}













#####################################################################33
###########A priori Jeffreys 2#######################################
#####################################################################


jefpriort2=function(x,sigma,H,kappa,cov.model,tau2,xmat){
covpars=c(sigma,x[1])
p=ncol(xmat)
n=nrow(xmat)
ndata=n
nu1=ndata-p
nucand=x[2]
nu0=nucand
cand=covpars
sigma=covpars[1]
phi=covpars[2]
phicand=phi

cova=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covpars,nugget=tau2,kappa=kappa)
#cova=varcov.spatial(coords,cov.model=cov.model,cov.pars=covpars,nugget=tau2,kappa=kappa)$varcov

cova=cova+(10e-4*diag(1,ndata))

#R=(cova/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
#R=(cova/sigma)-((tau2/sigma)*diag(n)) +(10e-7*diag(1,n))
#Rinv=solve(R,tol=10e-20)
covainv=solve(cova,tol=10e-20)



derivcova=sigma*derivcormatrix1(H,phi=phi,kappa=kappa,type=cov.model)$dev1


matu=(covainv%*%derivcova)

compbeta=((nu0+n)/(nu0+n+2))*(t(xmat)%*%covainv%*%xmat)


trmatu=sum(diag(matu))
trmatu2=sum(diag(matu%*%matu))
atr1=((nu0+n)/(nu0+n+2))*0.5
atr2=(0.5/(nu0+n+2))

compphi2=(atr1*trmatu2)-(atr2*((trmatu)^2))





#compsigma2= ( (((nu0+n)/(nu0+n+2))*(0.5*n))- ((0.5/(nu0+n+2))*(n^2)))
compsigma2=0.5*((n*nu0)/(nu0+n+2))
compsigma2=compsigma2


c1=trigamma((n+nu0)/2)- trigamma(nu0/2)
#c2= n/(nu0*(nu0+n))
#c3=1/(nu0+n)
#c4= (nu0+2)/(nu0*(nu0+n+2))
c2=(n*(nu0+n+4))/(nu0*(nu0+n)*(nu0+n+2))

#compnu2= -0.5*((0.5*c1)+c2-c3+c4)

compnu2= -0.5*((0.5*c1)+c2)




# (0.5)*(((nu0+n)/(nu0+n+2))-(n/(nu0+n+2)))*(trmatu)

csigmaa=(nu0)/(2*(nu0+n+2))

compphisigma=csigmaa*(trmatu)


compphinu= -(1/((nu0+n+2)*(nu0+n)))*trmatu

compnusigma= -(n/((nu0+n+2)*(nu0+n)))


#mat1=c(compsigma2,compphisigma,compnusigma,compphisigma,compphi2,compphinu,compnusigma,compphinu,compnu2)

#mat2=matrix(t(mat1),3,3)

expresion= (compsigma2*compphi2*compnu2)+(compphisigma*compphinu*compnusigma)+ (compnusigma*compphisigma*compphinu)-
((compnusigma^2)*compphi2)- ((compphinu^2)*compsigma2)-((compphisigma^2)*(compnu2))

detcompbeta=det(compbeta)
#detmat2=det(mat2)



if(sqrt(expresion*detcompbeta)<1e-323){
priornuphijef= 1e-323
}else{
priornuphijef= sqrt(detcompbeta*expresion)
}

return(priornuphijef)
}

########priorivaga####

vagprior=function(x,sigma,anu,bnu,aphi,bphi){

covpars=c(sigma,x[1])
nu0=x[2]
taxnucand=x[3]
dens=dunif(covpars[2],aphi,bphi)*dtrunc(nu0,spec="exp",rate=taxnucand,a=anu,b=bnu)*dunif(taxnucand,0.02,0.5)
return(dens)
}

#####################################################################################
############# A priori ind de jeffreys 2####################################
##################################################################


jefpriortind2=function(x,sigma,H,kappa,cov.model,tau2){
ndata=dim(H)[1]
n=ndata
nu0=x[2]
sigma=sigma
phi=x[1]
covpars=c(sigma,x[1])
cova=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covpars,nugget=tau2,kappa=kappa)
#cova=varcov.spatial(coords,cov.model=cov.model,cov.pars=covpars,nugget=tau2,kappa=kappa)$varcov

cova=cova+(10e-4*diag(1,ndata))

covainv=solve(cova,tol=10e-20)

derivcova=sigma*derivcormatrix1(H,phi=phi,kappa=kappa,type=cov.model)$dev1


matu=(covainv%*%derivcova)

trmatu=sum(diag(matu))
trmatu2=sum(diag(matu%*%matu))

compphi2=(((nu0+n)/(nu0+n+2))*0.5*trmatu2)-((0.5/(nu0+n+2))*((trmatu)^2))

compsigma2= ( (((nu0+n)/(nu0+n+2))*(0.5*n))- ((0.5/(nu0+n+2))*(n^2)))
compsigma2=compsigma2/(sigma^2)

c1=trigamma((n+nu0)/2)- trigamma(nu0/2)
#c2= n/(nu0*(nu0+n))
#c3=1/(nu0+n)
#c4= (nu0+2)/(nu0*(nu0+n+2))
c2=(n*(nu0+n+4))/(nu0*(nu0+n)*(nu0+n+2))

#compnu2= -0.5*((0.5*c1)+c2-c3+c4)

compnu2= -0.5*((0.5*c1)+c2)


# (0.5)*(((nu0+n)/(nu0+n+2))-(n/(nu0+n+2)))*(trmatu)

compphisigma=(0.5/sigma)*(nu0/(nu0+n+2))*(trmatu)


compphinu= -(1/((nu0+n+2)*(nu0+2)))*trmatu

compnusigma= -(n/(sigma*(nu0+n+2)*(nu0+2)))


mat1=c(compsigma2,compphisigma,compnusigma,compphisigma,compphi2,compphinu,compnusigma,compphinu,compnu2)

mat2=matrix(t(mat1),3,3)
detmat2=abs(det(mat2))


if(detmat2<1e-323){
priornuphijef= 0
}else{
priornuphijef= sqrt(detmat2)
}

return(priornuphijef)
}



###########################################################################
######## obtendo um valor da posteriori de referencia cjta ################
###########################################################################


posteriorref=function(y,cons,sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat){
ycomp=y
n=nrow(xmat)
p=ncol(xmat)
v=n-p
last=c(sigma,x[1])
nulast=x[2]
Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilastaux=varcov.spatial(coords,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)$varcov
#Psilast=Psilastaux$varcov
#Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



reflast=prioriref2(sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)

detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}

#priorphilast=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast*dlnorm(nucand,meanlog=log(nulast),sdlog=sdlog)
posterior=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast
return(posterior/cons)
}


#posteriorref(y=yobs,cons=consm1$integral,coords=coordsobs,sigma=sigma,x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmatobs1)

#y=yobs;cons=consm1$integral;coords=coordsobs;sigma=100;theta1=theta1;H=H;kappa=kappa;cov.model=cov.model;tau2=tau2;x=xobs1


#####obtendo um valor da aposteriori jeffreys########

posteriorjef=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat){
ycomp=y
n=nrow(xmat)
p=ncol(xmat)
v=n-p
last=c(sigma,x[1])
nulast=x[2]
Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilastaux=varcov.spatial(coords=coords,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux$varcov
#Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



jeflast=jefpriort2(x=x,H=H,sigma=sigma,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)

detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}


a= (p/2) + 1
cons1=((nulast)^(a-1))*gamma((nulast/2)-(a-1))/gamma(nulast/2)
if(nulast>40){
cons1= 2^(1-a)
}

posterior=cons1*(((S2last)^(-v/2-a+1))*sqrt(det(vbetalast))/sqrt(detrlast))*jeflast
return(posterior/cons)

}


#####################################################
#####obtendo um valor da aposteriori Jeffreys########
#####################################################

posteriorjef=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat){
ycomp=y
n=nrow(xmat)
p=ncol(xmat)
v=n-p
last=c(sigma,x[1])
nulast=x[2]
Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilastaux=varcov.spatial(coords=coords,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux$varcov
#Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



jeflast=jefpriort2(x=x,H=H,sigma=sigma,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)

detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}


a= (p/2) + 1
cons1=((nulast)^(a-1))*gamma((nulast/2)-(a-1))/gamma(nulast/2)
if(nulast>40){
cons1= 2^(1-a)
}

posterior=cons1*(((S2last)^(-v/2-a+1))*sqrt(det(vbetalast))/sqrt(detrlast))*jeflast
return(posterior/cons)

}




################################################333333
######### A posteriori jef ind#######################
###################################################3

#x=x;y=y;cons=cons;sigma=sigma;H=H;kappa=kappa;cov.model=cov.model;tau2=tau2;xmat=xmat

posteriorjefind=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat){
ycomp=y
ndata=nrow(H)
n=ndata
p=ncol(xmat)
v=n-p
last=c(sigma,x[1])
nulast=x[2]
Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilastaux=varcov.spatial(coords=coords,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux$varcov
#Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)



jeflast=jefpriortind2(x=x,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,sigma=sigma)

detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}


posterior=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*jeflast
return(posterior/cons)

}


######################################################
######obtendo um valor da aposteriori vaga############
######################################################

posteriorvag=function(y,cons,sigma,x,H,kappa,cov.model,tau2,xmat,aphi,bphi,anu,bnu,asigma){
ycomp=y
n=nrow(H)
p=ncol(xmat)
v=n-p
last=c(sigma,x[1])
nulast=x[2]
phi=x[1]
taxnucand=x[3]
Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilastaux=varcov.spatial(coords=coords,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
#Psilast=Psilastaux$varcov
#Psilast=Psilastaux

Rlast=(Psilast/sigma)-((tau2/sigma)*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)


#taxnucand=runif(1,0.02,0.5)
vaglast=dtrunc(nulast,spec="exp",rate=taxnucand,a=anu,b=bnu)*dunif(phi,aphi,bphi)*dunif(taxnucand,0.02,0.5)

detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}


a= asigma
cons1=((nulast)^(a-1))*gamma((nulast/2)-(a-1))/gamma(nulast/2)
if(nulast>40){
cons1= 2^(1-a)
}

posterior=cons1*(((S2last)^(-v/2 -a+1))*sqrt(det(vbetalast))/sqrt(detrlast))*vaglast
return(posterior/cons)

}


#######################################
###Integral M caso T###################
#######################################

intm1tref=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval){
consm1mat2=hcubature(prioriref2, lowerLimit=c(0,4.1), upperLimit=c(Inf,Inf),H=H,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
consm1posmat2=hcubature(posteriorref, lowerLimit=c(0,4+1e-10), upperLimit=c(Inf,Inf),H=H,y=y,cons=consm1mat2$integral,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)

ms=c(consm1posmat2$integral)
return(ms)
}

intm1tjef=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval){
consm1mat2=hcubature(jefpriort2, lowerLimit=c(0,4.1), upperLimit=c(Inf,Inf),H=H,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
consm1posmat2=hcubature(posteriorjef, lowerLimit=c(0,4.1), upperLimit=c(Inf,Inf),H=H,y=y,cons=consm1mat2$integral,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)

ms=c(consm1posmat2$integral)
return(ms)
}


intm1tjefind=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval){
consm1mat2=hcubature(jefpriortind2, lowerLimit=c(0,4.1), upperLimit=c(Inf,Inf),H=H,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,maxEval=maxEval)
consm1posmat2=hcubature(posteriorjefind, lowerLimit=c(0,4.1), upperLimit=c(Inf,Inf),H=H,y=y,cons=consm1mat2$integral,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval)
ms=c(consm1posmat2$integral)
return(ms)
}

intm1tvag=function(kappa,H,sigma,cov.model,tau2,xmat,y,maxEval,asigma,aphi,bphi,anu,bnu){
consm1posmat2=hcubature(posteriorvag, lowerLimit=c(aphi,anu,0.02), upperLimit=c(bphi,bnu,0.5),H=H,y=y,cons=1,kappa=kappa,sigma=sigma,cov.model=cov.model,tau2=tau2,xmat=xmat,maxEval=maxEval,aphi=aphi,bphi=bphi,anu=anu,bnu=bnu,asigma=asigma)
ms=c(consm1posmat2$integral)
return(ms)
}





#############################################################################
##############amostrando da posteriori T################################
#############################################################################

#############################################
####Reference and jeffreys independent priors
#############################################

bayspaestT1=function(sdnu=1,method,xmat,y,coords,covini,nuini,tau2,kappa,cov.model,aphi,bphi,anu,bnu,burn,iter,thin,prior="refprior",cprop){
p=ncol(xmat)
n=nrow(xmat)
H=as.matrix(dist(coords,upper=T,diag=T))
distancias=unique(as.vector(H))
Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covini,nugget=tau2,kappa=kappa)
#Psi1=varcov.spatial(coords,cov.model=cov.model,cov.pars=covini,nugget=tau2,kappa=kappa)$varcov
Psi1=(Psi1+t(Psi1))/2
beta1=solve(t(xmat)%*%solve(Psi1)%*%xmat)%*%t(xmat)%*%solve(Psi1)%*%y

betaF=matrix(0,iter+1,p)
betaF[1,]=beta1
sigmaF=c(covini[1],rep(0,iter))
phiF=c(covini[2],rep(0,iter))
nuF=c(nuini,rep(0,iter))

count=0
countiter=0
if(n<=501){
for(j in 2:iter){

ycomp=y

n=nrow(xmat)
p=ncol(xmat)
v=n-p
philast=phiF[j-1]
#phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
aphiproplast = min(aphi,max(philast-cprop,aphi))
bphiproplast = max(bphi,min(philast + cprop,bphi))

phicand=runif(1,aphiproplast,bphiproplast)
aphipropcand = min(aphi,max(aphi,phicand-cprop))
bphipropcand = max(bphi,min(phicand + cprop,bphi))


nulast=nuF[j-1]
#nucand=rtrunc(1,spec="exp",rate=1/nulast,a=anu,b=bnu)
nucand=rtrunc(1,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)

cand=c(sigmaF[j-1],phicand)
Psicand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)

Rcand=(Psicand/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))

Rcandinv=solve(Rcand,tol=10e-20)
vbetacand<- solve(t(xmat)%*%Rcandinv%*%xmat)
     vbetacand<-(vbetacand+t(vbetacand))/2

mubetacand<-   vbetacand%*%t(xmat)%*%Rcandinv%*%ycomp
S2cand=t(ycomp-xmat%*%mubetacand)%*%Rcandinv%*%(ycomp-xmat%*%mubetacand)/(n-p)

detrcand=det(Rcand)
if(detrcand<1e-323){
detrcand=1e-323
}

philast=phiF[j-1]
last=c(sigmaF[j-1],phiF[j-1])
Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)

Rlast=(Psilast/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast,tol=10e-20)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2

mubetalast<-   vbetalast%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)


detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}

#if(prior=="reference"){
#refcand=prioriref2(x=c(phicand,nucand),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
#reflast=prioriref2(x=c(philast,nulast),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
  refcand=prioriref2mat(covacand=Rcand,phi=phicand,vbetacand=vbetacand,covacandinv=Rcandinv,nu=nucand,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
  reflast=prioriref2mat(covacand=Rlast,phi=philast,vbetacand=vbetalast,covacandinv=Rlastinv,nu=nulast,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
  a=1
#}



#if(prior=="jef.ind"){
#refcand=jefpriortind2(x=c(phicand,nucand),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
#reflast=jefpriortind2(x=c(philast,nulast),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
#a=1
#}

priorphicand=(((S2cand)^(-v/2))*sqrt(det(vbetacand))/sqrt(detrcand))*refcand*dtrunc(nulast,spec="lnorm",meanlog=log(nucand),sdlog=sdnu,a=anu,b=bnu)*dunif(philast,aphipropcand,bphipropcand)
#priorphicand=(((S2cand)^(-v/2))*sqrt(det(vbetacand))/sqrt(detrcand))*refcand*dtrunc(nulast,spec="exp",rate=1/nulast,a=anu,b=bnu)



priorphilast=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dunif(phicand,aphiproplast,bphiproplast)
#priorphilast=(((S2last)^(-v/2))*sqrt(det(vbetalast))/sqrt(detrlast))*reflast*dtrunc(nucand,spec="exp",rate=1/nulast,a=anu,b=bnu)




Num=priorphicand
Den=priorphilast

unif=runif(1)
  if(unif<Num/Den)
  {
next.pi2F <- cand
mubeta=mubetacand
vbeta=vbetacand
S2=S2cand
next.nuF <-nucand
count=count+1
}else
  {
    next.pi2F <- last
next.nuF <- nulast
  mubeta=mubetalast
vbeta=vbetalast
S2=S2last
  }



betaF[j,]=rmvt(1,mu=t(mubeta),S=as.numeric(S2)*vbeta,df=v)
sigmaaux2=rf(1,next.nuF,v)
sigmaF[j]=sigmaaux2*S2
phiF[j]=next.pi2F[2]
nuF[j]=next.nuF
covinii=c(sigmaF[j],phiF[j])
#Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covinii,nugget=tau2,kappa=kappa)
#Psi1=(Psi1+t(Psi1))/2
countiter=countiter+1
cat("Iteration ",countiter," of ",iter,"\r")
#print(c(betaF[j,],sigmaF[j],phiF[j],nuF[j],j))
}
}

else{
  for(j in 2:iter){
    ycomp=y

    n=nrow(xmat)
    p=ncol(xmat)
    v=n-p

    philast=phiF[j-1]
    #phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
    aphiproplast = min(aphi,max(philast-cprop,aphi))
    bphiproplast = max(bphi,min(philast + cprop,bphi))

    phicand=runif(1,aphiproplast,bphiproplast)
    aphipropcand = min(aphi,max(aphi,phicand-cprop))
    bphipropcand = max(bphi,min(phicand + cprop,bphi))


    #phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=0,b=Inf)
    nulast=nuF[j-1]
    #nucand=rtrunc(1,spec="exp",rate=1/nulast,a=anu,b=bnu)
    nucand=rtrunc(1,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)
    last=c(sigmaF[j-1],phiF[j-1])

    Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
    Rlast=(Psilast/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))
    Rlastinv=solve(Rlast,tol=10e-20)
    vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
    vbetalast<-(vbetalast+t(vbetalast))/2
    mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
    S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)

    sigmacand=S2last*rf(1,nulast,v)
    sigmalast=sigmaF[j-1]
    betalast=betaF[j-1,]
    cand=c(sigmacand,phicand)

    Psicand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
    Rcand=(Psicand/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))

    Rcandinv=solve(Rcand,tol=10e-20)
    vbetacand<- solve(t(xmat)%*%Rcandinv%*%xmat)
    vbetacand<-(vbetacand+t(vbetacand))/2

    #mubetacand<-   solve(t(xmat)%*%Rcandinv%*%xmat)%*%t(xmat)%*%Rcandinv%*%ycomp
    #S2cand=t(ycomp-xmat%*%mubetacand)%*%Rcandinv%*%(ycomp-xmat%*%mubetacand)/(n-p)
    betacand=t(rmvt(1,mu=t(mubetalast),S=as.numeric(S2last)*vbetalast,df=v))

    covcand=c(sigmacand,phicand)


    #detrcand=det(Rcand)
    # if(detrcand<1e-300){
    #  detrcand=1e-300
    #}

    detrlast=det(Rlast)
    if(detrlast<1e-300){
      detrlast=1e-300
    }


    sigmalast=last[1]
    philast=last[2]
    covlast=c(sigmalast,philast)


    refcand=prioriref2mat(covacand=Rcand,phi=phicand,vbetacand=vbetacand,covacandinv=Rcandinv,nu=nucand,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
    reflast=prioriref2mat(covacand=Rlast,phi=philast,vbetacand=vbetalast,covacandinv=Rlastinv,nu=nulast,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
    a=1

    #print(c(detrcand,detrlast,refcand,reflast))

    xbetacand=xmat%*%betacand
    priorphicand=(dmvt(t(ycomp),mu=t(xbetacand),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nucand),sdlog=sdnu,a=anu,b=bnu)*dunif(philast,aphipropcand,bphipropcand)*S2last*df(sigmalast,nulast,v)*dmvt(betalast,mu=t(mubetalast),S=as.numeric(S2last)*vbetalast,df=v)
    #priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)*dtrunc(philast, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
    #priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmalast,log(sigmalast),sdsigma)

    #priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)
    # priorphicand=priorphicand1*refcand*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmalast,nulast,n)
    #priorphicand=priorphicand1*refcand*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmalast,log(sigmalast),sdsigma)
    ##################################
    #consint= (t1cand^v1cand)
    #print(priorphicand)

    if(priorphicand<1e-323){
      priorphicand=1e-323
    }

    xbetalast=xmat%*%betalast
    priorphilast=(dmvt(t(ycomp),mu=t(xbetalast),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dunif(phicand,aphiproplast,bphiproplast)*S2last*df(sigmacand,nulast,v)*dmvt(betalast,mu=t(mubetalast),S=as.numeric(S2last)*vbetalast,df=v)
    #priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2cand*df(sigmacand,nucand,v)*dtrunc(philast, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
    #priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmacand,log(sigmalast),sdsigma)
    #priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmacand,nulast,v)

    #print(c(priorphicand,priorphilast))
    Num=priorphicand
    Den=priorphilast

    unif=runif(1)
    if(unif<Num/Den)
    {
      next.pi2F <- cand
      next.nuF <-nucand
      next.sigma=sigmacand
      betanew=betacand
      count=count+1
    }else
    {
      next.pi2F <- last
      next.nuF <- nulast
      next.sigma=sigmalast
      betanew=betalast
    }

    betaF[j,]=betanew
    sigmaF[j]=next.sigma
    phiF[j]=next.pi2F[2]
    nuF[j]=next.nuF
    covinii=c(sigmaF[j],phiF[j])

    #    Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covinii,nugget=tau2,kappa=kappa)
    #   Psi1=(Psi1+t(Psi1))/2

    #print(c(betaF[j,],sigmaF[j],phiF[j],nuF[j],j))
    countiter=countiter+1
    cat("Iteration ",countiter," of ",iter,"\r")
  }
}





betaburn=as.matrix(betaF[(burn+1):iter,])
betaval=as.matrix(betaburn[seq((burn+1),iter-burn,thin),])
phiburn=phiF[(burn+1):iter]
phival=phiburn[seq((burn+1),iter-burn,thin)]
sigmaburn=sigmaF[(burn+1):iter]
sigmaval=sigmaburn[seq((burn+1),iter-burn,thin)]
nuburn=nuF[(burn+1):iter]
nuval=nuburn[seq((burn+1),iter-burn,thin)]
dist=cbind(betaval,sigmaval,phival,nuval)

###estimadores########

if(method=="mode"){
modebeta=apply(betaval,2,mlv)
modephi=mlv(phival)
modesigma=mlv(sigmaval)
mediannu=mlv(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}

if(method=="mean"){
modebeta=apply(betaval,2,mean)
modephi=mean(phival)
modesigma=mean(sigmaval)
mediannu=mean(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}


if(method=="median"){
modebeta=apply(betaval,2,median)
modephi=median(phival)
modesigma=median(sigmaval)
mediannu=median(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}


dist=cbind(betaval,sigmaval,phival,nuval)
return(list(prob=count/iter,dist=dist,betaF=betaval,sigmaF=sigmaval,phiF=phival,nuF=nuval,coords=coords,nugget=tau2,kappa=kappa,X=xmat,type=cov.model,theta=theta,y=ycomp))
}


##########################################################
#####Amostrando da A posteriori de Jeffreys##############
#########################################################


bayspaestTjef=function(candpar, prior,sdnu, method,xmat,y,coords,covini,nuini,tau2,kappa,cov.model,aphi,bphi,anu,bnu,burn,iter,thin,cprop){

p=ncol(xmat)
n=nrow(xmat)

H=as.matrix(dist(coords,upper=T,diag=T))
Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covini,nugget=tau2,kappa=kappa)

Psi1=(Psi1+t(Psi1))/2
beta1=solve(t(xmat)%*%solve(Psi1)%*%xmat)%*%t(xmat)%*%solve(Psi1)%*%y

betaF=matrix(0,iter+1,p)
betaF[1,]=beta1
sigmaF=c(covini[1],rep(0,iter))
phiF=c(covini[2],rep(0,iter))
nuF=c(nuini,rep(0,iter))
count=0
countiter=0

if(n<=200){
for(j in 2:iter){
ycomp=y

n=nrow(xmat)
p=ncol(xmat)
v=n-p
philast=phiF[j-1]
aphiproplast = min(aphi,max(philast-cprop,aphi))
bphiproplast = max(bphi,min(philast + cprop,bphi))

phicand=runif(1,aphiproplast,bphiproplast)
aphipropcand = min(aphi,max(aphi,phicand-cprop))
bphipropcand = max(bphi,min(phicand + cprop,bphi))

#phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=0,b=Inf)
nulast=nuF[j-1]
#nucand=rtrunc(1,spec="exp",rate=1/nulast,a=anu,b=bnu)
nucand=rtrunc(1,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)
last=c(sigmaF[j-1],phiF[j-1])

Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
Rlast=(Psilast/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast,tol=10e-20)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2
mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)

sigmacand=S2last*rf(1,nulast,v)
sigmalast=sigmaF[j-1]
#sigmacand=rlnorm(log(sigmalast),sdsigma)
cand=c(sigmacand,phicand)

Psicand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
Rcand=(Psicand/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))

Rcandinv=solve(Rcand,tol=10e-20)
vbetacand<- solve(t(xmat)%*%Rcandinv%*%xmat)
     vbetacand<-(vbetacand+t(vbetacand))/2

mubetacand<-   solve(t(xmat)%*%Rcandinv%*%xmat)%*%t(xmat)%*%Rcandinv%*%ycomp
S2cand=t(ycomp-xmat%*%mubetacand)%*%Rcandinv%*%(ycomp-xmat%*%mubetacand)/(n-p)

covcand=c(sigmacand,phicand)


detrcand=det(Rcand)
if(detrcand<1e-300){
detrcand=1e-300
}

detrlast=det(Rlast)
if(detrlast<1e-300){
  detrlast=1e-300
}


sigmalast=last[1]
philast=last[2]
covlast=c(sigmalast,philast)


if(prior=="jef.ind"){
  refcand=jefpriortind2(x=c(phicand,nucand),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
  reflast=jefpriortind2(x=c(philast,nulast),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
  a=1
}




if(prior=="jef.rul"){
  refcand=jefpriort2(x=c(phicand,nucand),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
  reflast=jefpriort2(x=c(philast,nulast),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
  a=1+p/2
}








#print(c(detrcand,detrlast,refcand,reflast))

#xbeta=xmat%*%betaF[j-1,]
#priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)
#priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)*dtrunc(philast, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
#priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmalast,log(sigmalast),sdsigma)

#priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)

#############################
t1cand= nucand + (v*S2cand/sigmacand)
v1cand=(nucand+v)/2

if(v1cand>=160){
  v1cand=160
  gamma1cand=sqrt(2*pi*(v1cand-1))*(((v1cand-1)/exp(1))^(v1cand-1))
}
else{
  gamma1cand=gamma(v1cand)
}

t2cand=gamma1cand/gamma(nucand/2)
t3cand=sigmacand^(-((v/2)+a))
priorphicand1= t2cand* (t3cand/(t1cand^v1cand))*sqrt(det(vbetacand)/detrcand)*(nucand^(nucand/2))

if(priorphicand1<1e-323){
  priorphicand1=1e-323
}

priorphicand=priorphicand1*refcand*dtrunc(nulast,spec="lnorm",meanlog=log(nucand),sdlog=sdnu,a=anu,b=bnu)*dunif(philast,aphipropcand,bphipropcand)*S2last*df(sigmalast,nulast,n)
#priorphicand=priorphicand1*refcand*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmalast,log(sigmalast),sdsigma)
##################################
#consint= (t1cand^v1cand)
#print(priorphicand)

if(priorphicand<1e-323){
priorphicand=1e-323
}


#priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2cand*df(sigmacand,nulast,v)
#priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2cand*df(sigmacand,nucand,v)*dtrunc(philast, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
#priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmacand,log(sigmalast),sdsigma)
#priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmacand,nulast,v)


#############################
t1last= nulast + (v*S2last/sigmalast)
v1last=(nulast+v)/2

if(v1last>=160){
  v1last=160
  gamma1last=sqrt(2*pi*(v1last-1))*(((v1last-1)/exp(1))^(v1last-1))
}
else{
  gamma1last=gamma(v1last)
}
t2last=gamma(v1last)/gamma(nulast/2)
t3last=sigmalast^(-((v/2)+a))
priorphilast1= t2last* (t3last/(t1last^v1last))*sqrt(det(vbetalast)/detrlast)*(nulast^(nulast/2))

if(priorphilast1<1e-323){
  priorphilast1=1e-323
}


priorphilast=priorphilast1*reflast*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dunif(phicand,aphiproplast,bphiproplast)*S2last*df(sigmacand,nulast,n)
#priorphilast=priorphilast1*reflast*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmalast,log(sigmalast),sdsigma)
##################################

if(priorphilast<1e-323){
priorphilast=1e-323
}

#print(c(priorphicand,priorphilast))




Num=priorphicand
Den=priorphilast

unif=runif(1)
  if(unif<Num/Den)
  {
next.pi2F <- cand
mubeta=mubetacand
vbeta=vbetacand
S2=S2cand
next.nuF <-nucand
next.sigma=sigmacand
#betanew=betacand
count=count+1
  }else
  {
    next.pi2F <- last
next.nuF <- nulast
next.sigma=sigmalast
  mubeta=mubetalast
vbeta=vbetalast
S2=S2last
 }

betaF[j,]=rmvt(1,mu=t(mubeta),S=as.numeric(S2)*vbeta,df=v)
sigmaF[j]=next.sigma
phiF[j]=next.pi2F[2]
nuF[j]=next.nuF
covinii=c(sigmaF[j],phiF[j])

Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covinii,nugget=tau2,kappa=kappa)
Psi1=(Psi1+t(Psi1))/2

#print(c(betaF[j,],sigmaF[j],phiF[j],nuF[j],j,"J"))
countiter=countiter+1
cat("Iteration ",countiter," of ",iter,"\r")
}
}
else{
  for(j in 2:iter){
    ycomp=y

    n=nrow(xmat)
    p=ncol(xmat)
    v=n-p
    philast=phiF[j-1]
    aphiproplast = min(aphi,max(philast-cprop,aphi))
    bphiproplast = max(bphi,min(philast + cprop,bphi))

    phicand=runif(1,aphiproplast,bphiproplast)
    aphipropcand = min(aphi,max(aphi,phicand-cprop))
    bphipropcand = max(bphi,min(phicand + cprop,bphi))

    #phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=0,b=Inf)
    nulast=nuF[j-1]
    #nucand=rtrunc(1,spec="exp",rate=1/nulast,a=anu,b=bnu)
    nucand=rtrunc(1,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)
    last=c(sigmaF[j-1],phiF[j-1])

    Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
    Rlast=(Psilast/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))
    Rlastinv=solve(Rlast,tol=10e-20)
    vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
    vbetalast<-(vbetalast+t(vbetalast))/2
    mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
    S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)

    sigmacand=S2last*rf(1,nulast,v)
    sigmalast=sigmaF[j-1]
    betalast=betaF[j-1,]
    cand=c(sigmacand,phicand)

    Psicand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
    #Rcand=(Psicand/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))

    #Rcandinv=solve(Rcand,tol=10e-20)
    #vbetacand<- solve(t(xmat)%*%Rcandinv%*%xmat)
    #vbetacand<-(vbetacand+t(vbetacand))/2

    #mubetacand<-   solve(t(xmat)%*%Rcandinv%*%xmat)%*%t(xmat)%*%Rcandinv%*%ycomp
    #S2cand=t(ycomp-xmat%*%mubetacand)%*%Rcandinv%*%(ycomp-xmat%*%mubetacand)/(n-p)
    betacand=t(rmvt(1,mu=t(mubetalast),S=as.numeric(S2last)*vbetalast,df=v))

    covcand=c(sigmacand,phicand)


    #detrcand=det(Rcand)
     # if(detrcand<1e-300){
    #  detrcand=1e-300
    #}

    detrlast=det(Rlast)
    if(detrlast<1e-300){
      detrlast=1e-300
    }


    sigmalast=last[1]
    philast=last[2]
    covlast=c(sigmalast,philast)


    if(prior=="jef.ind"){
      refcand=jefpriortind2(x=c(phicand,nucand),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
      reflast=jefpriortind2(x=c(philast,nulast),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2)
      a=1
    }




    if(prior=="jef.rul"){
      refcand=jefpriort2(x=c(phicand,nucand),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
      reflast=jefpriort2(x=c(philast,nulast),sigma=1,H=H,kappa=kappa,cov.model=cov.model,tau2=tau2,xmat=xmat)
      a=1+p/2
    }








    #print(c(detrcand,detrlast,refcand,reflast))

    xbetacand=xmat%*%betacand
    priorphicand=(dmvt(t(ycomp),mu=t(xbetacand),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nucand),sdlog=sdnu,a=anu,b=bnu)*dunif(philast,aphipropcand,bphipropcand)*S2last*df(sigmalast,nulast,v)*dmvt(betalast,mu=t(mubetalast),S=as.numeric(S2last)*vbetalast,df=v)
    #priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)*dtrunc(philast, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
    #priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmalast,log(sigmalast),sdsigma)

    #priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)
   # priorphicand=priorphicand1*refcand*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmalast,nulast,n)
    #priorphicand=priorphicand1*refcand*dtrunc(nulast,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmalast,log(sigmalast),sdsigma)
    ##################################
    #consint= (t1cand^v1cand)
    #print(priorphicand)

    if(priorphicand<1e-323){
      priorphicand=1e-323
    }

    xbetalast=xmat%*%betalast
    priorphilast=(dmvt(t(ycomp),mu=t(xbetalast),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dunif(phicand,aphiproplast,bphiproplast)*S2last*df(sigmacand,nulast,v)*dmvt(betalast,mu=t(mubetalast),S=as.numeric(S2last)*vbetalast,df=v)
    #priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2cand*df(sigmacand,nucand,v)*dtrunc(philast, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
    #priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dlnorm(sigmacand,log(sigmalast),sdsigma)
    #priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*S2last*df(sigmacand,nulast,v)

#print(c(priorphicand,priorphilast))
    Num=priorphicand
    Den=priorphilast

    unif=runif(1)
    if(unif<Num/Den)
    {
      next.pi2F <- cand
      next.nuF <-nucand
      next.sigma=sigmacand
      betanew=betacand
      count=count+1
    }else
    {
      next.pi2F <- last
      next.nuF <- nulast
      next.sigma=sigmalast
      betanew=betalast
    }

    betaF[j,]=betanew
    sigmaF[j]=next.sigma
    phiF[j]=next.pi2F[2]
    nuF[j]=next.nuF
    covinii=c(sigmaF[j],phiF[j])

#    Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covinii,nugget=tau2,kappa=kappa)
 #   Psi1=(Psi1+t(Psi1))/2

    #print(c(betaF[j,],sigmaF[j],phiF[j],nuF[j],j,"J"))
    countiter=countiter+1
    cat("Iteration ",countiter," of ",iter,"\r")
  }
}

betaburn=as.matrix(betaF[(burn+1):iter,])
betaval=as.matrix(betaburn[seq((burn+1),iter-burn,thin),])
phiburn=phiF[(burn+1):iter]
phival=phiburn[seq((burn+1),iter-burn,thin)]
sigmaburn=sigmaF[(burn+1):iter]
sigmaval=sigmaburn[seq((burn+1),iter-burn,thin)]
nuburn=nuF[(burn+1):iter]
nuval=nuburn[seq((burn+1),iter-burn,thin)]

dist=cbind(betaval,sigmaval,phival,nuval)

if(method=="mode"){
modebeta=apply(betaval,2,mlv)
modephi=mlv(phival)
modesigma=mlv(sigmaval)
mediannu=median(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}


if(method=="mean"){
modebeta=apply(betaval,2,mean)
modephi=mean(phival)
modesigma=mean(sigmaval)
mediannu=median(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}


if(method=="median"){
modebeta=apply(betaval,2,median)
modephi=median(phival)
modesigma=median(sigmaval)
mediannu=median(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}




return(list(prob=count/iter,dist=dist,betaF=betaval,sigmaF=sigmaval,phiF=phival,nuF=nuval,coords=coords,X=xmat,nugget=tau2,kappa=kappa,type=cov.model,theta=theta,y=ycomp))
}



##########################################################
#####Amostrando da A posteriori de vaga##############
#########################################################



bayspaestTvag=function(candpar,method,xmat,y,coords,covini,nuini,tau2,kappa,cov.model,aphi,bphi,asigma,anu,bnu,sdnu,burn,iter,thin,cprop){

p=ncol(xmat)
n=nrow(xmat)

H=as.matrix(dist(coords,upper=T,diag=T))
Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covini,nugget=tau2,kappa=kappa)

Psi1=(Psi1+t(Psi1))/2
beta1=solve(t(xmat)%*%solve(Psi1)%*%xmat)%*%t(xmat)%*%solve(Psi1)%*%y

betaF=matrix(0,iter+1,p)
betaF[1,]=beta1
sigmaF=c(covini[1],rep(0,iter))
phiF=c(covini[2],rep(0,iter))
nuF=c(nuini,rep(0,iter))
count=0
countiter=0
for(j in 2:iter){
ycomp=y

n=nrow(xmat)
p=ncol(xmat)
v=n-p

philast=phiF[j-1]
aphiproplast = min(aphi,max(philast-cprop,aphi))
bphiproplast = max(bphi,min(philast + cprop,bphi))

phicand=runif(1,aphiproplast,bphiproplast)
aphipropcand = min(aphi,max(aphi,phicand-cprop))
bphipropcand = max(bphi,min(phicand + cprop,bphi))

#phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=0,b=Inf)
nulast=nuF[j-1]
#nucand=rtrunc(1,spec="exp",rate=1/nulast,a=anu,b=bnu)
nucand=rtrunc(1,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)
#phicand=rtrunc(1, spec="gamma",shape=phiF[j-1],scale=candpar,a=aphi,b=Inf)
last=c(sigmaF[j-1],phiF[j-1])

Psilast=varcov.spatial(H=H,cov.model=cov.model,cov.pars=last,nugget=tau2,kappa=kappa)
Rlast=(Psilast/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))
Rlastinv=solve(Rlast)
vbetalast<- solve(t(xmat)%*%Rlastinv%*%xmat)
     vbetalast<-(vbetalast+t(vbetalast))/2
mubetalast<-   solve(t(xmat)%*%Rlastinv%*%xmat)%*%t(xmat)%*%Rlastinv%*%ycomp
S2last=t(ycomp-xmat%*%mubetalast)%*%Rlastinv%*%(ycomp-xmat%*%mubetalast)/(n-p)

sigmacand=0.5*S2last*rf(1,nulast,n)
sigmalast=sigmaF[j-1]
betacand=rmvt(1,mu=t(mubetalast),S=as.numeric(S2last)*vbetalast,df=v)
betalast=betaF[j-1,]
cand=c(sigmacand,phicand)

Psicand=varcov.spatial(H=H,cov.model=cov.model,cov.pars=cand,nugget=tau2,kappa=kappa)
Rcand=(Psicand/sigmaF[j-1])-((tau2/sigmaF[j-1])*diag(n)) +(10e-4*diag(1,n))

Rcandinv=solve(Rcand)
vbetacand<- solve(t(xmat)%*%Rcandinv%*%xmat)
     vbetacand<-(vbetacand+t(vbetacand))/2

mubetacand<-   solve(t(xmat)%*%Rcandinv%*%xmat)%*%t(xmat)%*%Rcandinv%*%ycomp
S2cand=t(ycomp-xmat%*%mubetacand)%*%Rcandinv%*%(ycomp-xmat%*%mubetacand)/(n-p)

covcand=c(sigmacand,phicand)



detrcand=det(Rcand)
if(detrcand<1e-323){
detrcand=1e-323
}



sigmalast=last[1]
philast=last[2]
covlast=c(sigmalast,philast)

detrlast=det(Rlast)
if(detrlast<1e-323){
detrlast=1e-323
}



xbeta=xmat%*%betaF[j-1,]
a=asigma
#priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)*refcand/(sigmacand^a))*dtrunc(nulast,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2last*df(sigmalast,nulast,v)
priorphicand=(dmvt(t(ycomp),mu=t(xbeta),S=Psicand,df=nucand,log=F)/(sigmacand^a))*S2last*df(sigmalast,nucand,n)*dtrunc(nulast,spec="lnorm",meanlog=log(nucand),sdlog=sdnu,a=anu,b=bnu)*dunif(philast,aphipropcand,bphipropcand)



if(priorphicand<1e-323){
priorphicand=1e-323
}


#priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)*reflast/(sigmalast^a))*dtrunc(nucand,spec="exp",rate=1/nulast,a=anu,b=bnu)*S2cand*df(sigmacand,nucand,v)
priorphilast=(dmvt(t(ycomp),mu=t(xbeta),S=Psilast,df=nulast,log=F)/(sigmalast^a))*S2last*df(sigmacand,nulast,n)*dtrunc(nucand,spec="lnorm",meanlog=log(nulast),sdlog=sdnu,a=anu,b=bnu)*dunif(phicand,aphiproplast,bphiproplast)

if(priorphilast<1e-323){
priorphilast=1e-323
}





Num=priorphicand
Den=priorphilast

unif=runif(1)
  if(unif<Num/Den)
  {
next.pi2F <- cand
mubeta=mubetacand
vbeta=vbetacand
S2=S2cand
next.nuF <-nucand
next.sigma=sigmacand
#betanew=betacand
count=count+1
  }else
  {
    next.pi2F <- last
next.nuF <- nulast
next.sigma=sigmalast
  mubeta=mubetalast
vbeta=vbetalast
S2=S2last
 }

betaF[j,]=rmvt(1,mu=t(mubeta),S=as.numeric(S2)*vbeta,df=v)
sigmaF[j]=next.sigma
phiF[j]=next.pi2F[2]
nuF[j]=next.nuF
covinii=c(sigmaF[j],phiF[j])

Psi1=varcov.spatial(H=H,cov.model=cov.model,cov.pars=covinii,nugget=tau2,kappa=kappa)
Psi1=(Psi1+t(Psi1))/2

#print(c(betaF[j,],sigmaF[j],phiF[j],nuF[j],j,"J"))
countiter=countiter+1
cat("Iteration ",countiter," of ",iter,"\r")
}

betaburn=as.matrix(betaF[(burn+1):iter,])
betaval=as.matrix(betaburn[seq((burn+1),iter-burn,thin),])
phiburn=phiF[(burn+1):iter]
phival=phiburn[seq((burn+1),iter-burn,thin)]
sigmaburn=sigmaF[(burn+1):iter]
sigmaval=sigmaburn[seq((burn+1),iter-burn,thin)]
nuburn=nuF[(burn+1):iter]
nuval=nuburn[seq((burn+1),iter-burn,thin)]

dist=cbind(betaval,sigmaval,phival,nuval)

if(method=="mode"){
modebeta=apply(betaval,2,mlv)
modephi=mlv(phival)
modesigma=mlv(sigmaval)
mediannu=median(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}


if(method=="mean"){
modebeta=apply(betaval,2,mean)
modephi=mean(phival)
modesigma=mean(sigmaval)
mediannu=median(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}


if(method=="median"){
modebeta=apply(betaval,2,median)
modephi=median(phival)
modesigma=median(sigmaval)
mediannu=median(nuval)

theta=c(modebeta,modesigma,modephi,mediannu)
}



return(list(prob=count/iter,dist=dist,betaF=betaval,sigmaF=sigmaval,phiF=phival,nuF=nuval,coords=coords,X=xmat,nugget=tau2,kappa=kappa,type=cov.model,theta=theta,y=ycomp))
}

###obj(tsroba)

prediction1=function(z,xpred,coordspred,xobs,coordsobs,theta,cov.model,tau2,kappa){
  #####values of obj (class(tsroba))
  #include coordspred.cols,covar.col
  coords=rbind(coordsobs,coordspred)
  H=as.matrix(dist(coords,upper=T,diag=T))
  xpred=as.matrix(xpred)
  xobs=as.matrix(xobs)
  type=cov.model
  n=dim(xobs)[1]+ dim(xpred)[1]
  npred=dim(xpred)[1]
  s=dim(xobs)[2]
  beta=as.matrix(theta[1:s])
  phi=theta[(s+2)]
  sigma2=theta[(s+1)]
  nu=theta[(s+3)]
  nobs=dim(xobs)[1]
  pred.row= (nobs+1):(nobs+npred)
  pred=rep(0,n)
  pred[pred.row]=1

  Sigma=varcov.spatial(H=H,cov.model=type,cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)

   Sigmaa2=Sigma[pred==1,pred==0]
  Sigmaa3=solve(Sigma[pred==0,pred==0])
  Sigmaa4=Sigma[pred==1,pred==1]
  Sigmaa5=Sigma[pred==0,pred==1]
  dista=z-(xobs%*%beta)
  prediction=(xpred%*%beta) + (Sigmaa2%*%Sigmaa3%*%dista)

  delta=t(dista)%*%Sigmaa3%*%(dista)

  cons= (nu+delta)/(nu+nobs)
  cons=as.numeric(cons)

  S=cons*(Sigmaa4-(Sigmaa2%*%Sigmaa3%*%Sigmaa5))
  S1=t(S)
  S[lower.tri(S)]=S1[lower.tri(S1)]
  cons1=(nu+nobs)/(nu+nobs-2)

  #sdpred=sqrt(diag(cons1*S))
  coordspred=coords[pred==1,]
  coordsobs=coords[pred==0,]
  #predict=cbind(prediction,sdpred)
  amos=rmvt(n=1, S= S, df = nu+nobs, mu =t(prediction))
  #amos=rmvt(n=1, S = S, df = nu+nobs, mu =t(prediction))
  #return(amos)
  return(amos)
}


###############7.Prediction for the SAEM algorithm###################

predictionnorm=function(z,xpred,coordspred,xobs,coordsobs,theta,cov.model,tau2,kappa){
  xobs=as.matrix(xobs)
  xpred=as.matrix(xpred)
  colnames(coordsobs)=c("x","y")
  colnames(coordspred)=c("x","y")
  coords=rbind(coordsobs,coordspred)
  H=as.matrix(dist(coords,upper=T,diag=T))
  n=dim(xobs)[1]+ dim(xpred)[1]
  npred=dim(xpred)[1]
  s=dim(xobs)[2]
  nobs=dim(xobs)[1]
  pred.row= (nobs+1):(nobs+npred)
  beta=theta[1:s]
  phi=theta[(s+2)]
  sigma2=theta[(s+1)]
  Sigma=varcov.spatial(H=H,cov.model=cov.model,cov.pars=c(sigma2,phi),nugget=tau2,kappa=kappa)
  pred=rep(0,n)
  pred[pred.row]=1
  prediction=(xpred%*%beta)+ ( (Sigma[pred==1,pred==0]%*%solve(Sigma[pred==0,pred==0]))%*%(z-(xobs%*%beta)))
  S=Sigma[pred==1,pred==1]-(Sigma[pred==1,pred==0]%*%solve(Sigma[pred==0,pred==0])%*%Sigma[pred==0,pred==1])
  sdpred=sqrt(diag(S))
  #predict=
  amos=rmvnorm(n=1, sigma= S, mean =t(prediction))
  #amos=rmvn(n=1, Sigma= S, mu =t(prediction))
  coordspred=coords[pred==1,]
  coordsobs=coords[pred==0,]
  return(amos)
}


