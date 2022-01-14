rm(list=ls())
library(Hmisc) # cut2
library(survey) # svycoxph
library(doParallel) # mclapply
library(nnet)# multinom
library(kyotil);           stopifnot(packageVersion("kyotil")>="2021.2-2")
library(marginalizedRisk); stopifnot(packageVersion("marginalizedRisk")>="2021.2-4")
library(DengueTrialsYF)
#
B=1e3 # number of available cores
numCores=30 # bootstrap replicates
time.start=Sys.time()


####################################################################################################
# vaccine arm

for(setting in c("cont","cat")) {    
    res=sapply(c("cyd14","cyd15"), simplify="array", function (trial) {    
    
# setting="cat"; trial="cyd15"
    
        dat=make.m13.dat(trial, stype=0)
        dat=subset(dat, trt==1)
        dat$wt=1/dat$sampling.p
        
        # adjust for protocol-specified age categories, sex, and country
        if (trial=="cyd14") {
            f.0=Surv(X,d) ~ old + little + gender + MYS + PHL + THA + VNM
            #q.a=c(-Inf,log10(58),log10(266),Inf)# cut points from Moodie et al (2018)
        } else {
            f.0=Surv(X,d) ~ old + gender + COL + HND + MEX + PRI 
            #q.a=c(-Inf,log10(135),log10(631),Inf)# cut points from Moodie et al (2018)
        }
        
        if (setting=="cat") {
            q.a=c(-Inf,quantile(dat$titer[dat$fasi=="Y"], c(1/3,2/3), na.rm=T),Inf)                
            dat$s = factor(cut2(dat$titer,cuts=q.a)) 
            print(table(dat$s))
            ss=NULL
        } else {
            dat$s = dat$titer
            ss=quantile(dat$titer[dat$fasi=="Y"], seq(.05,.95,by=0.01), na.rm=TRUE) 
            myprint(ss)        
        }
        
        #hist(dat$titer[dat$fasi=="Y"], main=trial)
        #with(data, table(fasi, indicators==1, d, useNA="ifany"))
            
        t0=365; myprint(max(dat$X[dat$d==1]))#363 
        
        get.marginalized.risk=function(dat){
            # risk regression
            dat.design=twophase(id=list(~1,~1),strata=list(NULL,~d),subset=~indicators, data=dat)
            fit.risk = svycoxph(update(f.0, ~.+s), design=dat.design)
#            # marker regression
#            fit.s=if(setting=="cat") nnet::multinom(update(f.0, s~.), dat[dat$fasi=="Y",], trace=FALSE) else lm(update(f.0, s~.), dat[dat$fasi=="Y",]) 
            # marginal risk estimation
            marginalized.risk(fit.risk, "s", data=subset(dat, indicators==1), categorical.s=setting=="cat", weights=subset(dat, indicators==1, wt, drop=T), t=t0, ss=ss)    
        }
        
        prob=get.marginalized.risk(dat)
        
        # preparing for bootstrap
        # ptids to bootstrap
        ptids.cases=   subset(dat, d==1, ptid, drop=TRUE)
        ptids.controls=subset(dat, d==0 & indicators==1, ptid, drop=TRUE)
        ptids.rest=    subset(dat, d==0 & indicators==0, ptid, drop=TRUE)    
            
        # bootstrapping
        # store the current rng state 
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
        #   
        out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
            set.seed(seed) 
            ## three-step bootstrap process
            tmp=list()
            # 1. bootstrap the cases 
            tmp[[1]]=sample(ptids.cases, replace=TRUE)        
            # 2. bootstrap the controls
            tmp[[2]]=sample(ptids.controls, replace=TRUE)        
            # 3. add rest, only to get the counts right
            tmp[[3]]=ptids.rest
            idxes=do.call(c, tmp)
            dat.b=dat[match(idxes, dat$ptid),]
            
            get.marginalized.risk(dat.b)
        })
        boot=do.call(cbind, out)
        # restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)    
        
        if (setting=="cat") {
            c(prob[3]/prob[1], quantile(boot[3,]/boot[1,], c(.025,.975)))            
        } else {
            cbind(marker=ss, prob, boot)            
        }    
    })
    if(setting=="cat") rownames(res)[1]="est"
    if (!dir.exists("input")) dir.create("input")
    save(res, file=paste0("input/res_", setting, ".Rdata"))
}



####################################################################################################
# placebo arm

# only implement continuous markers and not categorical markers
res.placebo.cont=sapply(c("cyd14","cyd15"), simplify="array", function (trial) {    
# trial="cyd15"
    
    dat=make.m13.dat(trial, stype=0)
    dat=subset(dat, trt==0)
    n=nrow(dat)
    
    # adjust for protocol-specified age categories, sex, and country
    if (trial=="cyd14") {
        f.0=Surv(X,d) ~ old + little + gender + MYS + PHL + THA + VNM
        #q.a=c(-Inf,log10(58),log10(266),Inf)# cut points from Moodie et al (2018)
    } else {
        f.0=Surv(X,d) ~ old + gender + COL + HND + MEX + PRI 
        #q.a=c(-Inf,log10(135),log10(631),Inf)# cut points from Moodie et al (2018)
    }
    
    t0=365; myprint(max(dat$X[dat$d==1]))#363 
    
    get.marginalized.risk=function(dat){
        fit.risk = coxph(f.0, dat, model=T) # model=T is required because the type of prediction requires it, see Note on ?predict.coxph
        dat$X=t0
        risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
        mean(risks)
    }
    
    prob=get.marginalized.risk(dat)
    
    # bootstrapping
    # store the current rng state 
    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
    if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }   
    out=mclapply(1:B, mc.cores = numCores, FUN=function(seed) {   
        set.seed(seed) 
        dat.b=dat[sample.int(n, replace=T),]            
        get.marginalized.risk(dat.b)
    })
    boot=do.call(cbind, out)
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    c(est=prob, boot)
})
save(res.placebo.cont, file=paste0("input/res_placebo_cont.Rdata"))

print(Sys.time()-time.start)
