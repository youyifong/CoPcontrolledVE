# nonparametric estimation of placebo risk for revision
library(kyotil)
library(survtmle)
library(DengueTrialsYF)


# process arguments
Args <- commandArgs(trailingOnly=TRUE) 
if (length(Args)==0) {
    Args=c(batch.size="1",batch.number="0") 
}
myprint(Args)
i=0;
i=i+1; batch.size=as.numeric(Args[i])
i=i+1; batch=as.numeric(Args[i])
seeds=1:batch.size+batch.size*(batch-1); names(seeds)=seeds
myprint(batch, batch.size)

trials=c("cyd14", "cyd15"); names(trials)=trials
t0=365

begin=Sys.time()
res=
sapply(seeds, simplify="array", function (seed) {
    
    t.0=Sys.time()   
    myprint(seed)
    set.seed(seed) 
    
    out=sapply(trials, function(trial) {
    #trial="cyd14"
        myprint(trial)
        dat=make.m13.dat(trial, stype=0)
        if (seed>0) {
            dat.0=subset(dat, trt==0)
            dat.1=subset(dat, trt==1)
            dat=rbind(dat.0[sample.int(nrow(dat.0), replace=T),], dat.1[sample.int(nrow(dat.1), replace=T),])
        }
        #myprint(max(dat$X[dat$d==1]))
        
        ftime=dat$X
        ftype=dat$d
        trt=dat$trt
        if (trial=="cyd14") {
            adjustVars=subset(dat, select=c(old , little , gender , MYS , PHL , THA , VNM))
            fstr="trt + old + little + gender + MYS + PHL + THA + VNM"
        } else {
            adjustVars=subset(dat, select=c(old , gender , COL , HND , MEX , PRI ))
            fstr="trt + old + gender + COL + HND + MEX + PRI"
        }
            
#        fit <- try(survtmle(ftime = ftime, ftype = ftype,
#                         trt = trt, adjustVars = adjustVars,
#                         glm.ftime = fstr,
#                         glm.ctime = fstr,
#                         method = "mean", t0 = t0, verbose=T))
        
        fit <- survtmle(ftime = ftime, ftype = ftype,
                         trt = trt, adjustVars = adjustVars,
                         SL.ftime = c("SL.glm","SL.mean","SL.step"),
                         glm.ctime = fstr,
                         method = "mean", t0 = t0, verbose=F)
        
        if(!inherits(fit,"try-error")) fit$est[,1] else c(NA,NA)
    })    
    
    gc()# there are some memory problem, seem to quit automatically
    print("time used under this seed: "%.%format(Sys.time()-t.0))      
    out
                
})

# save results
save (res, file="input/tmlest2/tmlest_batch"%.%formatInt(batch, 3)%.%".Rdata")
# note time passed
done = Sys.time()
body1=format(done-begin)
print(date())
print("time used: "%.%body1)
