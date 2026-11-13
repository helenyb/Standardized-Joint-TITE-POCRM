###simulations setup
source("JOINT_TITE_POCRM_schedule.R")
library(pocrm)

#define the labelling of doses
dose.ind.mat<-which(matrix(c(1:6),nrow=3)>0,arr.ind=T)



tox_scen1<-c(0.1,
             0.15,
             0.2,
             0.15,
             0.2,
             0.3)
tox_scen2<-c(0.1,
             0.2,
             0.4,
             0.15,
             0.3,
             0.5
)

tox_scen3<-c(0.1,
             0.2,
             0.3,
             0.4,
             0.45,
             0.5
)
tox_scen4<-c(0.1,
             0.3,
             0.4,
             0.2,
             0.45,
             0.5
)

tox_scen5<-c(0.3,
             0.4,
             0.45,
             0.45,
             0.5,
             0.6
)


eff_scen1<-c(0.2,
             0.25,
             0.4,
             0.3,
             0.35,
             0.45
)
eff_scen2<-c(0.05,
             0.2,
             0.3,
             0.1,
             0.4,
             0.5
)
eff_scen3<-c(0.05,
             0.1,
             0.3,
             0.15,
             0.4,
             0.5
)
eff_scen4<-c(0.05,
             0.1,
             0.15,
             0.2,
             0.3,
             0.4
)

ndoses<-6



#how to split the activity across cycles for data generation
eff.pattern1<-rep(1,3)/3
# eff.pattern2<-c(1:3)/6
# eff.pattern3<-c(3:1)/6
eff.pattern<-1
ncycles<-3

tox_scens<-c(1:5)
eff_scens<-c(1:4)
#calculate parameters for data generation
for(pattern in 1){
  for(eff_scen in eff_scens){
    parMat<-matrix(nrow=2,ncol=ndoses)
    cycleMat<-     cyc_func_eff_v2(cyc_all_vec = get(paste(c("eff_scen",eff_scen),collapse="")), 
                                   split_vec = get(paste(c("eff.pattern",pattern),collapse="")))
    for(j in 1:ncol(cycleMat)){
      parMat[,j]<- find_lognormal_parms3(p1=cycleMat[1,j],p3=sum(cycleMat[,j]),int2=seq(0.01,10,0.01))
    }
    assign(paste(c("eff_scen",eff_scen,".pars"),collapse=""),
           parMat)
  }
}


for(tox in tox_scens){
  parMat<-matrix(nrow=2,ncol=ndoses)
  cycleMat<-cyc_func_tox_v2(cyc_all_vec=get(paste(c("tox_scen",tox),collapse="")),cyc1_prop=0.75)
  for(j in 1:ncol(cycleMat)){
    parMat[,j]<-  find_lognormal_parms3(p1=cycleMat[1,j],p3=cycleMat[2,j],int2=seq(0.01,10,0.01))
  }
  assign(paste(c("tox_scen",tox,".pars"),collapse=""),
         parMat)
}



prior_hyp2<-matrix(nrow=10,ncol=3)

prior_hyp2[1]<- log(1/16) #mu_alphaT
prior_hyp2[2]<-1 #tau_alphaT
prior_hyp2[3]<- log(1/4) #mu_betaT
prior_hyp2[4]<- 1/2 #tau_betaT


prior_hyp2[5]<- -3 #mu_alphaE
prior_hyp2[6]<-1 #tau_alphaE
prior_hyp2[7]<- -.2 #mu_betaE
prior_hyp2[8]<- 1 #tau_betaE


#relationship
prior_hyp2[9,]<-0 #mu_phi
prior_hyp2[10,]<-0.01 #tau_phi


generate.psdose.matrix<-function(min.psdose,max.psdose,realdose.1,realdose.2){
  psdose.matrix<-matrix(nrow=7,ncol=6)
  
  A.val<-(realdose.2+max.psdose)/2
  B.val<-realdose.1+(realdose.2-realdose.1)/3
  E.val<-(realdose.1+min.psdose)/2
  F.val<-realdose.1+2*(realdose.2-realdose.1)/3
  # real doses
  psdose.matrix[,2]<-realdose.1
  psdose.matrix[,5]<-realdose.2
  
  psdose.matrix[,1]<-min.psdose
  psdose.matrix[,6]<-max.psdose
  
  #order 1
  psdose.matrix[1,c(3,4)]<-c(B.val,F.val) #cols
  psdose.matrix[2,c(3,4)]<-c(A.val,E.val) #rows
  psdose.matrix[3,c(3,4)]<-c(F.val,B.val) #up
  psdose.matrix[4,c(3,4)]<-c(B.val,E.val) #down-up1
  psdose.matrix[5,c(3,4)]<-c(F.val,E.val) #down-up2
  psdose.matrix[6,c(3,4)]<-c(A.val,B.val) #up-down1
  psdose.matrix[7,c(3,4)]<-c(A.val,F.val) #up-downp2
  
  return(psdose.matrix)
}

reals<-c(4.5,6)

psdose.matrix.prop.prA<- generate.psdose.matrix(min.psdose=reals[1]-1.5,max.psdose=reals[2]+1.5,realdose.1=reals[1],realdose.2=reals[2])
psdose.matrix.prop.prB<- generate.psdose.matrix(min.psdose=reals[1]-1.9,max.psdose=reals[2]+1.9,realdose.1=reals[1],realdose.2=reals[2])
psdose.matrix.prop.prC<- generate.psdose.matrix(min.psdose=reals[1]-2.4,max.psdose=reals[2]+2.4,realdose.1=reals[1],realdose.2=reals[2])

psdose.matrix.cal.prA<- generate.psdose.matrix(min.psdose=reals[1]-1.5,max.psdose=reals[2]+1.5,realdose.1=reals[1],realdose.2=reals[2])
psdose.matrix.cal.prB<- generate.psdose.matrix(min.psdose=reals[1]-1.5,max.psdose=reals[2]+1.5,realdose.1=reals[1],realdose.2=reals[2])
psdose.matrix.cal.prC<- generate.psdose.matrix(min.psdose=reals[1]-1,max.psdose=reals[2]+1,realdose.1=reals[1],realdose.2=reals[2])


#number of simulations
nsims<-1000 

#escalation only allowed in 1 direction at a time
ds_rule<-"ON"

#default ordering
def.order<-c(1:6)

#prior probs assigned to working models for part 1
prior.probsT<-prior.probsE<-rep(1/nrow(psdose.matrix.prop.prA),nrow(psdose.matrix.prop.prA))
load("Previous_Data_v2.RData")



save.image("standardzed_sims_setup.RData")
