library(OBASpatial)
library(geoR)
library(LaplacesDemon)
library(tidyverse)

source("auxiliar_functions_deviance.R")


set.seed(25)
data(dataca20)
d1=dataca20[1:158,]

xpred=model.matrix(calcont~altitude+area,data=dataca20[159:178,])
xobs=model.matrix(calcont~altitude+area,data=dataca20[1:158,])
coordspred=dataca20[159:178,1:2]

####################################################################
######covariance matern: kappa=0.3 prior:reference
res03=readRDS('newcadeia11kappa0.3.rds')
summary(res03)
######covariance matern: kappa=0.3 prior:jef.rul
res103=readRDS('newcadeiakappa0.3jefrul.rds')
######covariance matern: kappa=0.3 prior:jef.ind
res203=readRDS('newcadeiakappa0.3jefind.rds')
######covariance matern: kappa=0.3 prior:vague
res303=readRDS('newcadeiakappa0.3vague.rds')
######################################################################

#############################
##### DIC kappa = 0.3 ######
############################
auxdic03 = DICinf(res03)
auxdic103 = DICinf(res103)
auxdic203 = DICinf(res203)
auxdic303 = DICinf(res303)

### prior:reference
dic03 =auxdic03$DIC
### prior: jef.rul
dic103 =auxdic103$DIC
### prior:jef.ind
dic203 =auxdic203$DIC
### prior:vague
dic303 =auxdic303$DIC



##############################
##### predictions: kappa = 0.3
##############################

### prior:reference
pme03=tsrobapred(res03,xpred=xpred,coordspred=coordspred)$pred

### prior:jef.rul
pme103=tsrobapred(res103,xpred=xpred,coordspred=coordspred)$pred

### prior:jef.ind
pme203=tsrobapred(res203,xpred=xpred,coordspred=coordspred)$pred

### prior:vague
pme303=tsrobapred(res303,xpred=xpred,coordspred=coordspred)$pred

#### MEAN SQUARE OF ERROR Kappa = 0.3

### prior:reference
mse03=mean((pme03-dataca20$calcont[159:178])^2)
### prior:jef.rul
mse103=mean((pme103-dataca20$calcont[159:178])^2)
### prior:jef.ind
mse203=mean((pme203-dataca20$calcont[159:178])^2)
### prior:vague
mse303=mean((pme303-dataca20$calcont[159:178])^2)



##### Parentthesis P(M_i|Y) for kappa = 0.3###############
###(just comparing priors with kappa=0.3).
###the real aplication (see Ordonez et.al) consider kappa=0.3,0.5,0.7.
######### Using reference prior ###########

#m103=intmT(prior="reference",formula=calcont~altitude+area,
 ##        kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000)
######### Using Jeffreys' rule prior ###########
#m1j03=intmT(prior="jef.rul",formula=calcont~altitude+area,
 #         kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000)
######### Using Jeffreys' independent prior ###########
#m1ji03=intmT(prior="jef.ind",formula=calcont~altitude+area
 #          ,kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000)
#m1v03=intmT(prior="vague",formula=calcont~altitude+area
 #         ,kappa=0.3,cov.model="matern",data=dataca20,maxEval=1000,intphi="default")
#tot03=m103+m1j03+m1ji03+m1v03

####posterior probabilities#####
#p103=m103/tot03
#pj03=m1j03/tot03
#pji03=m1ji03/tot03
#pv03=m1v03/tot03






######covariance matern: kappa=0.5 prior:reference
res05=readRDS('newcadeia11kappa0.5.rds')
######covariance matern: kappa=0.5 prior:jef.rul
res105=readRDS('newcadeiakappa0.5jefrul.rds')
######covariance matern: kappa=0.5 prior:jef.ind
res205=readRDS('newcadeiakappa0.5jefind.rds')


######covariance matern: kappa=0.5 prior:vague
res305=readRDS('newcadeiakappa0.5vague.rds')


#############################
##### DIC kappa = 0.5 ######
############################
auxdic05 = DICinf(res05)
auxdic105 = DICinf(res105)
auxdic205 = DICinf(res205)
auxdic305 = DICinf(res305)

### prior:reference
dic05 =auxdic05$DIC
### prior: jef.rul
dic105 =auxdic105$DIC
### prior:jef.ind
dic205 =auxdic205$DIC
### prior:vague
dic305 =auxdic305$DIC


##############################
##### predictions: kappa = 0.5
##############################


### prior:reference
pme05=tsrobapred(res05,xpred=xpred,coordspred=coordspred)$pred
### prior:jef.rul
pme105=tsrobapred(res105,xpred=xpred,coordspred=coordspred)$pred
### prior:jef.ind
pme205=tsrobapred(res205,xpred=xpred,coordspred=coordspred)$pred
### prior:vague
pme305=tsrobapred(res305,xpred=xpred,coordspred=coordspred)$pred

#### MEAN SQUARE OF ERROR Kappa = 0.5

### prior:reference
mse05=mean((pme05-dataca20$calcont[159:178])^2)
### prior:jef.rul
mse105=mean((pme105-dataca20$calcont[159:178])^2)
### prior:jef.ind
mse205=mean((pme205-dataca20$calcont[159:178])^2)
### prior:vague
mse305=mean((pme305-dataca20$calcont[159:178])^2)


##### Parentthesis P(M_i|Y) for kappa = 0.5###############
###(just comparing priors with kappa=0.5).
###the real aplication (see Ordonez et.al) consider kappa=0.3,0.5,0.7.
######### Using reference prior ###########

#m105=intmT(prior="reference",formula=calcont~altitude+area,
##        kappa=0.5,cov.model="matern",data=dataca20,maxEval=1000)
######### Using Jeffreys' rule prior ###########
#m1j05=intmT(prior="jef.rul",formula=calcont~altitude+area,
#         kappa=0.5,cov.model="matern",data=dataca20,maxEval=1000)
######### Using Jeffreys' independent prior ###########
#m1ji05=intmT(prior="jef.ind",formula=calcont~altitude+area
#          ,kappa=0.5,cov.model="matern",data=dataca20,maxEval=1000)
#m1v05=intmT(prior="vague",formula=calcont~altitude+area
#         ,kappa=0.5,cov.model="matern",data=dataca20,maxEval=1000,intphi="default")
#tot05=m105+m1j05+m1ji05+m1v05

####posterior probabilities#####
#p105=m105/tot05
#pj05=m1j05/tot05
#pji05=m1ji05/tot05
#pv05=m1v05/tot05



######covariance matern: kappa=0.7 prior:reference
res07=readRDS('newcadeia11kappa0.7.rds')
#####covariance matern: kappa=0.7 prior:jef.rul
res107=readRDS('newcadeiakappa0.7jefrul.rds')
######covariance matern: kappa=0.7 prior:jef.ind
res207=readRDS('newcadeiakappa0.7jefind.rds')
######covariance matern: kappa=0.7 prior:vague
res307=readRDS('newcadeiakappa0.7vague.rds')

#############################
##### DIC kappa = 0.7 ######
############################
auxdic07 = DICinf(res07)
auxdic107 = DICinf(res107)
auxdic207 = DICinf(res207)
auxdic307 = DICinf(res307)

### prior:reference
dic07 =auxdic07$DIC
### prior: jef.rul
dic107 =auxdic107$DIC
### prior:jef.ind
dic207 =auxdic207$DIC
### prior:vague
dic307 =auxdic307$DIC



##############################
##### predictions: kappa = 0.7
##############################

### prior:reference
pme07=tsrobapred(res07,xpred=xpred,coordspred=coordspred)$pred
### prior:jef.rul
pme107=tsrobapred(res107,xpred=xpred,coordspred=coordspred)$pred
### prior:jef.ind
pme207=tsrobapred(res207,xpred=xpred,coordspred=coordspred)$pred
### prior:vague
pme307=tsrobapred(res307,xpred=xpred,coordspred=coordspred)$pred


#### MEAN SQUARE OF ERROR Kappa = 0.7

### prior:reference
mse07=mean((pme07-dataca20$calcont[159:178])^2)
### prior:jef.rul
mse107=mean((pme107-dataca20$calcont[159:178])^2)
### prior:jef.ind
mse207=mean((pme207-dataca20$calcont[159:178])^2)
### prior:vague
mse307=mean((pme307-dataca20$calcont[159:178])^2)


##### Parentthesis P(M_i|Y) for kappa = 0.7###############
###(just comparing priors with kappa=0.7).
###the real aplication (see Ordonez et.al) consider kappa=0.3,0.7,0.7.
######### Using reference prior ###########

#m107=intmT(prior="reference",formula=calcont~altitude+area,
##        kappa=0.7,cov.model="matern",data=dataca20,maxEval=1000)
######### Using Jeffreys' rule prior ###########
#m1j07=intmT(prior="jef.rul",formula=calcont~altitude+area,
#         kappa=0.7,cov.model="matern",data=dataca20,maxEval=1000)
######### Using Jeffreys' independent prior ###########
#m1ji07=intmT(prior="jef.ind",formula=calcont~altitude+area
#          ,kappa=0.7,cov.model="matern",data=dataca20,maxEval=1000)
#m1v07=intmT(prior="vague",formula=calcont~altitude+area
#         ,kappa=0.7,cov.model="matern",data=dataca20,maxEval=1000,intphi="default")
#tot07=m107+m1j07+m1ji07+m1v07

####posterior probabilities#####
#p107=m107/tot07
#pj07=m1j07/tot07
#pji07=m1ji07/tot07
#pv07=m1v07/tot07





##### joining results kappa = 0.3
refer03=c(res03$theta,dic03,mse03)
vague03=c(res303$theta,dic303,mse303)
jef03=c(res103$theta,dic103,mse103)
ji03=c(res203$theta,dic203,mse203)

result03=rbind(refer03,vague03,jef03,ji03)
rownames(result03)=rep("",4)

##### joining results kappa = 0.5
refer05=c(res05$theta,dic05,mse05)
vague05=c(res305$theta,dic305,mse305)
jef05=c(res105$theta,dic105,mse105)
ji05=c(res205$theta,dic205,mse205)

result05=rbind(refer05,vague05,jef05,ji05)
rownames(result05)=rep("",4)


##### joining results kappa = 0.7
refer07=c(res07$theta,dic07,mse07)
vague07=c(res307$theta,dic307,mse307)
jef07=c(res107$theta,dic107,mse107)
ji07=c(res207$theta,dic207,mse207)

result07=rbind(refer07,vague07,jef07,ji07)
rownames(result07)=rep("",4)

####

rownames(result03)=c("Reference","Vague","Jef.rul","Jef.ind")
rownames(result05)=c("Reference","Vague","Jef.rul","Jef.ind")
rownames(result07)=c("Reference","Vague","Jef.rul","Jef.ind")


#### summarized table ##########
table=rbind(result03,result05,result07)




########################################
### CONVERGENCE SELECTED MODEL #########
########################################


library(coda)
library(ggplot2)
library(tidyverse)
library(patchwork)

cad1 = readRDS('newcadeia11kappa0.3.rds')
cad2 = readRDS('newcadeia21kappa0.3.rds')

#### special requirements for using the function gelman.diag (convergence test for two chains)

##### transforming each parameter chain into an mcmc object
chain1beta0 = mcmc(cad1$betaF[,1],thin=10)
chain2beta0 = mcmc(cad2$betaF[,1],thin=10)

chain1beta1 = mcmc(cad1$betaF[,2],thin=10)
chain2beta1 = mcmc(cad2$betaF[,2],thin=10)

chain1beta2 = mcmc(cad1$betaF[,3],thin=10)
chain2beta2 = mcmc(cad2$betaF[,3],thin=10)

chain1beta3 = mcmc(cad1$betaF[,4],thin=10)
chain2beta3 = mcmc(cad2$betaF[,4],thin=10)


chain1sigma = mcmc(cad1$sigmaF,thin=10)
chain2sigma = mcmc(cad2$sigmaF,thin=10)

chain1phi = mcmc(cad1$phiF,thin=10)
chain2phi = mcmc(cad2$phiF,thin=10)

chain1nu = mcmc(cad1$nuF,thin=10)
chain2nu = mcmc(cad2$nuF,thin=10)


####### joining the two chains ######
listmcmcbeta0 = list(chain1beta0,chain2beta0)
listmcmcbeta0 = mcmc.list(listmcmcbeta0)

listmcmcbeta1 = list(chain1beta1,chain2beta1)
listmcmcbeta1 = mcmc.list(listmcmcbeta1)

listmcmcbeta2 = list(chain1beta2,chain2beta2)
listmcmcbeta2 = mcmc.list(listmcmcbeta2)

listmcmcbeta3 = list(chain1beta3,chain2beta3)
listmcmcbeta3 = mcmc.list(listmcmcbeta3)

listmcmcsigma = list(chain1sigma,chain2sigma)
listmcmcsigma = mcmc.list(listmcmcsigma)

listmcmcphi = list(chain1phi,chain2phi)
listmcmcphi = mcmc.list(listmcmcphi)

listmcmcnu = list(chain1nu,chain2nu)
listmcmcnu = mcmc.list(listmcmcnu)
###########################################

###### convergence criteria ####################
gelman.diag(listmcmcbeta0)
gelman.diag(listmcmcbeta1)
gelman.diag(listmcmcbeta2)
gelman.diag(listmcmcbeta3)
gelman.diag(listmcmcsigma)
gelman.diag(listmcmcphi)
gelman.diag(listmcmcnu)


###### getting the traceplot for each chain #######
data1 = data.frame(cad1$dist,'Chain 1',1:dim(cad1$dist)[1])
colnames(data1) = c(paste0('beta_',0:3),'sigma2','phi','nu','chain','Index')

data2 = data.frame(cad2$dist,'Chain 2',1:dim(cad2$dist)[1])
colnames(data2) = c(paste0('beta_',0:3),'sigma2','phi','nu','chain','Index')

datatot = rbind(data1,data2)



labelsb = c("Chain 1", "Chain 2")

g1 = datatot%>%ggplot(aes(x=Index,y=beta_0,color=factor(chain)))+ geom_line() +
  theme_bw()+
  scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsb)+
  ylab(expression(beta[1])) + xlab("Index")

g2 = datatot%>%ggplot(aes(x=Index,y=beta_1,color=factor(chain)))+ geom_line() +
  theme_bw()+
  scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsb)+
  ylab(expression(beta[2])) + xlab("Index")

g3 = datatot%>%ggplot(aes(x=Index,y=beta_2,color=factor(chain)))+ geom_line() +
  theme_bw()+
  scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsb)+
  ylab(expression(beta[3])) + xlab("Index")

g4 = datatot%>%ggplot(aes(x=Index,y=beta_3,color=factor(chain)))+ geom_line() +
  theme_bw()+
  scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsb)+
  ylab(expression(beta[4])) + xlab("Index")

g5 = datatot%>%ggplot(aes(x=Index,y=sigma2,color=factor(chain)))+ geom_line() +
  theme_bw()+
  scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsb)+
  ylab(expression(sigma^2)) + xlab("Index")

g6 = datatot%>%ggplot(aes(x=Index,y=phi,color=factor(chain)))+ geom_line() +
  theme_bw()+
  scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsb)+
  ylab(expression(phi)) + xlab("Index")


g7 = datatot%>%ggplot(aes(x=Index,y=nu,color=factor(chain)))+ geom_line() +
  theme_bw()+
  scale_color_manual(name = "Chains", values  = rep(c("black","gray")), labels = labelsb)+
  ylab(expression(nu)) + xlab("Index")


##### traceplots of both chains for all the parameters
g1+g2 +g3 +g4 + g5 + g6 + g7 + plot_layout(guides = 'collect')


##### acf plots for all the parameters
library(ggfortify)

p1 <- autoplot(acf(cad1$betaF[,1], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p1 = p1+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[0]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()



p1 <- autoplot(acf(cad1$betaF[,1], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p1 = p1+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[0]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()


p2 <- autoplot(acf(cad1$betaF[,2], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p2 = p2+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[1]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()


p3 <- autoplot(acf(cad1$betaF[,3], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p3 = p3+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[2]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()


p4 <- autoplot(acf(cad1$betaF[,4], plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p4 = p4+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(beta[3]))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()


p5 <- autoplot(acf(cad1$sigmaF, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p5 = p5+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(sigma^2))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()


p6 <- autoplot(acf(cad1$phiF, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p6 = p6+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(phi))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()


p7 <- autoplot(acf(cad1$nuF, plot = FALSE,lag.min=0,lag.max=30), conf.int.value = 0.95
               ,conf.int.fill=1,conf.int.col=1)

p7 = p6+geom_hline(yintercept = 0)+ coord_cartesian(x=c(0,28.5))+
  ggtitle(expression(nu))+
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()

p1 + p2 + p3 + p4 + p5 + p6 + p7







