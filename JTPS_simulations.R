
load("standardzed_sims_setup.RData")

library(doParallel)
library(mvtnorm)
library(rjags)
library(cubature)
#library(profvis)
registerDoParallel(cores=30)
print(timestamp())
counter<-0

for(prev.data.i in c("A","B","C")){   
  for(eff.scen in eff_scens){   
    for(tox.scen in tox_scens){  
      assign(paste(c("JOINTTITEPOCRM_stan_cal.pr",prev.data.i,"eff",eff.scen,"_tox",tox.scen,""),collapse=""),
             
            foreach(i=1:nsims, combine = list) %dopar% {
            #     foreach(i=1:nsims, combine = list) %do% {
               ##function
               JOINT.TITE.POCRM_schedule(seed=i,tru.E.pars = get(paste(c("eff_scen",eff.scen,".pars"),collapse="")),
                                             tru.T.pars=get(paste(c("tox_scen",tox.scen,".pars"),collapse="")),tru.corET=-0.5,
                                             co_size=3,ncohorts=20 ,target=0.3,
                                             ncycles=ncycles,dose.skipping.rule=ds_rule, 
                                             prior_vec2=prior_hyp2,
                                             sufficient.information=T,sufficient.information.lim=9,hard.safety.rule=95,safety.stopping.low.unsafe=T,
                                             safety.stopping.high.toosafe=F,initial.one.cycle=T,
                                             w1=0.33, w2=1.09, upper.tox=0.301,C_eff=0.2,C_tox=0.2,effbound=0.2,toxbound=0.3,
                                             TITE=T,dose.indices=dose.ind.mat,
                                             wmE=get(paste(c("psdose.matrix.cal.pr",prev.data.i),collapse="")),wmT=get(paste(c("psdose.matrix.cal.pr",prev.data.i),collapse="")),prior.probs = "EQUAL",
                                             prior.o.T=prior.probsT,prior.o.E=prior.probsE,
                                             default.order=def.order,start.dose=5,pause=T,a.stop.bound = 0,
                                             not.fixed = c(3,4), prev.data = get(paste(c("prev.data.",prev.data.i),collapse="")),
                                             order.E.index = c(1:3, 4, 4, 5, 5), order.T.index = c(1:3,4, 4, 5, 5))
               
             }
             
             
      )
      
      
      save.image("JOINT_TITEPOCRM_stan_sims.RData")
      print(timestamp())
      counter<-counter+1
      print(paste(c(
        "Standardized, Calibrated Spacing. Scenario T",tox.scen,".A",eff.scen, ", Previous Dataset ",prev.data.i,"(",counter,"/",3*length(tox_scens)*length(eff_scens),")"),collapse=""))
    }
  }
}
