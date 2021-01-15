rm(list=ls())
library(kyotil);       stopifnot(packageVersion("kyotil")>="2020.11.20")
library(marginalRisk); stopifnot(packageVersion("marginalRisk")>="2021.1.7") # bias.factor, E.value, controlled.risk.bias.factor
library(DengueTrialsYF) # need data to estimate overall attack rate to get s.ref
load(file=paste0("res_placebo_cont.Rdata")) # placebo arm results res.placebo.cont
RRud=RReu=4
bias.factor=bias.factor(RRud, RReu)
    

for(setting in c("cont","cat")) {
# setting="cont"; trial="cyd15"
    
    load(file=paste0("res_", setting, ".Rdata"))    
    
    if (setting=="cat") {
    
        # Table 1
        for (trial in c("cyd14","cyd15")) {
            write(paste0(
                # marginal RR and ci
                formatDouble(res[1,trial],2,remove.leading0=F), "&", formatDouble(res[2,trial],2,remove.leading0=F), "--", formatDouble(res[3,trial],2,remove.leading0=F)
                , "&" ,
                # causal RR and ci
                formatDouble(res[1,trial]*bias.factor,2,remove.leading0=F), "&", formatDouble(res[2,trial]*bias.factor,2,remove.leading0=F), "--", formatDouble(res[3,trial]*bias.factor,2,remove.leading0=F)
                , "&" ,
                # E-value and ub
                formatDouble(E.value(res[1,trial]),1), "&", formatDouble(E.value(res[3,trial]),1)
            ), file="input/CoPVeryHighVE_"%.%trial)    
        }
        
    } else {
    # continous S
        
        # these two reference quantiles are used in the next two blocks of code
        s2="85%"; s1="15%"        
                
        # results in the text is from comparison of s1 and s1
        for (trial in c("cyd14","cyd15")) {
            val=c(res[s2,"prob",trial]/res[s1,"prob",trial], quantile(res[s2,3:dim(res)[2],trial]/res[s1,3:dim(res)[2],trial], c(.025,.975)))
            # write to file, % is appended to avoid extra space
            # marginal RR and ci
            write(formatDouble(val[1],2,remove.leading0=F)%.%"%", file="input/"%.%trial%.%"_mrr")
            write(paste0(formatDouble(val[2],2,remove.leading0=F), "--", formatDouble(val[3],2,remove.leading0=F))%.%"%", file="input/"%.%trial%.%"_mrrci")
            # causal RR and ci
            write(formatDouble(val[1]*bias.factor,2,remove.leading0=F)%.%"%", file="input/"%.%trial%.%"_crr")
            write(paste0(formatDouble(val[2]*bias.factor,2,remove.leading0=F), "--", formatDouble(val[3]*bias.factor,2,remove.leading0=F))%.%"%", file="input/"%.%trial%.%"_crrci")
            # E-value and ub
            write(formatDouble(E.value(val[1]),1)%.%"%", file="input/"%.%trial%.%"_e")
            write(formatDouble(E.value(val[3]),1)%.%"%", file="input/"%.%trial%.%"_eul")
        }            
        
        # Fig 2
        ylim=c(0,0.055)
        mypdf(onefile=F, file=paste0("input/CoPveryhighVE_Fig2"), mfrow=c(1,2), oma=c(0,0,1,0))
            for (trial in c("cyd14","cyd15")) {        
                lwd=2
                
                # choose a reference marker value
                # use median 
                #s.ref=res["50%","marker",trial]
                # use the alternative 
                dat=make.m13.dat(trial, stype=0)
                dat=subset(dat, trt==1)
                which=which.min(abs(res[,"prob",trial]-mean(dat$d))); print(which)
                s.ref=res[which,"marker",trial]
    
                Bias=controlled.risk.bias.factor(ss=res[,"marker",trial], s.cent=s.ref, s1=res[s1,"marker",trial], s2=res[s2,"marker",trial], RRud) 
                
                ci.band=apply(res[,,trial], 1, function (x) quantile(x[3:length(x)], c(.025,.975)))
                mymatplot(res[,"marker",trial], t(rbind(res[,"prob",trial], ci.band)),      type="l", lty=c(1,2,2), col=4, lwd=lwd, make.legend=F, xlab="Month 13 Log10 Average Neutralizing Antibody Titer", ylab="Probability of VCD", main=toupper(trial), ylim=ylim)
                mymatplot(res[,"marker",trial], t(rbind(res[,"prob",trial], ci.band))*Bias, type="l", lty=c(1,2,2), col=3, lwd=lwd, make.legend=F, add=T)
                title(main="Marginalized and Controlled Risk of Virologically Confirmed Dengue by Antibody Titer", outer=T)
                mylegend(x=3,legend=c("Marginalized risk", "Controlled risk (conservative)"), lty=1, col=c(4,3), lwd=2, cex=.8)
            }
        dev.off()    
                

        # Fig 3
        ylim=NULL
        mypdf(onefile=F, file=paste0("input/CoPveryhighVE_Fig3"), mfrow=c(1,2), oma=c(0,0,1,0))
            for (trial in c("cyd14","cyd15")) {        
                lwd=2
                
                # choose a reference marker value
                # use median 
                #s.ref=res["50%","marker",trial]
                # use the alternative 
                dat=make.m13.dat(trial, stype=0)
                dat=subset(dat, trt==1)
                mean(dat$d)
                which=which.min(abs(res[,"prob",trial]-mean(dat$d))); print(which)
                s.ref=res[which,"marker",trial]
                # equation (6): Bias * marginalized risk
                Bias=controlled.risk.bias.factor(ss=res[,"marker",trial], s.cent=s.ref, s1=res[s1,"marker",trial], s2=res[s2,"marker",trial], RRud) 
                
                # CVE
                est = 1 - res[,"prob",trial]*Bias/res.placebo.cont["est",trial]
                boot = 1 - t(t(res[,3:dim(res)[2],trial]*Bias)/res.placebo.cont[2:dim(res.placebo.cont)[1],trial])                         
                ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
                mymatplot(res[,"marker",trial], t(rbind(est, ci.band)), type="l", lty=c(1,2,2), col=4, lwd=lwd, make.legend=F, xlab="Month 13 Log10 Average Neutralizing Antibody Titer", ylab="VE", main=toupper(trial), ylim=ylim)
                title(main="Controlled VE of Virologically Confirmed Dengue by Antibody Titer", outer=T)
                #mylegend(x=3,legend=c("Marginalized risk", "Controlled risk (conservative)"), lty=1, col=c(4,3), lwd=2, cex=.8)
            }
        dev.off()    

    }
    
}


# fix the column names
