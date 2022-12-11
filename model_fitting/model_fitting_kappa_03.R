library(OBASpatial)
set.seed(25)
data(dataca20)


###training data###
d1=dataca20[1:158,]

########################
##### priors############
########################
#
#### kappa = 0.3#######

##### reference
res=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
           ini.pars=c(10,390,10),iter=11000,burn=1000,thin=10)


##### jeffreys 

resjef=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
              ini.pars=c(10,390,10),iter=11000,burn=1000,thin=10,
              prior = 'jef.rul')


### independent jeffreys 
resjefind=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
                 ini.pars=c(10,390,10),iter=11000,burn=1000,thin=10,
                 prior = 'jef.ind')

#### vague 
resvague=tsroba(calcont~altitude+area, kappa = 0.3, data=d1,
                ini.pars=c(10,390,10),iter=11000,burn=1000,thin=10,
                prior = 'vague')

