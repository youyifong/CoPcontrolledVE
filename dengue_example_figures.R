rm(list=ls())
library(kyotil);       stopifnot(packageVersion("kyotil")>="2020.11.20")
library(DengueTrialsYF) # need data to estimate overall attack rate to get s.ref


for(setting in c("cat","cont")) {
# setting="cont"; trial="cyd14"
    load(file=paste0("res_", setting, ".Rdata"))
    ###############################
    # compute sensitivity measures
    
    RRud=RReu=4
    bias.factor=RRud*RReu/(RRud+RReu-1); myprint(bias.factor)
    E.value=function(rr) { rr1=1/rr; rr1 + sqrt(rr1 * (rr1-1)) }
    
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
            ), file="texinputs/CoPVeryHighVE_"%.%trial)    
        }
        
    } else {
    # if setting=="cont"
        
        # these two reference quantiles are used in the next two blocks of code
        s2="85%"; s1="15%"        
                
        # results in the text is from comparison of s1 and s1
        for (trial in c("cyd14","cyd15")) {
            val=c(res[s2,"prob",trial]/res[s1,"prob",trial], quantile(res[s2,3:dim(res)[2],trial]/res[s1,3:dim(res)[2],trial], c(.025,.975)))
            # write to file, % is appended to avoid extra space
            # marginal RR and ci
            write(formatDouble(val[1],2,remove.leading0=F)%.%"%", file="texinputs/"%.%trial%.%"_mrr")
            write(paste0(formatDouble(val[2],2,remove.leading0=F), "--", formatDouble(val[3],2,remove.leading0=F))%.%"%", file="texinputs/"%.%trial%.%"_mrrci")
            # causal RR and ci
            write(formatDouble(val[1]*bias.factor,2,remove.leading0=F)%.%"%", file="texinputs/"%.%trial%.%"_crr")
            write(paste0(formatDouble(val[2]*bias.factor,2,remove.leading0=F), "--", formatDouble(val[3]*bias.factor,2,remove.leading0=F))%.%"%", file="texinputs/"%.%trial%.%"_crrci")
            # E-value and ub
            write(formatDouble(E.value(val[1]),1)%.%"%", file="texinputs/"%.%trial%.%"_e")
            write(formatDouble(E.value(val[3]),1)%.%"%", file="texinputs/"%.%trial%.%"_eul")
        }            
        
        # Fig 2
        mypdf(onefile=F, file=paste0("CoPveryhighVE_Fig2"), mfrow=c(1,2), oma=c(0,0,1,0))
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
    
                # compute RR(s.ref, s) according to the display after equation (5) in the manuscript
                RR.ref.s <- exp( (res[,"marker",trial]-s.ref) / (res[s2,"marker",trial]-res[s1,"marker",trial]) * log(RRud) )
                RR.s.ref <- exp( (s.ref-res[,"marker",trial]) / (res[s2,"marker",trial]-res[s1,"marker",trial]) * log(RRud) )
                # compute B accoring to the display before equation (5) in the manuscript
                B.s.ref=RR.s.ref*RR.s.ref/(RR.s.ref+RR.s.ref-1)
                B.ref.s=RR.ref.s*RR.ref.s/(RR.ref.s+RR.ref.s-1)
                # compute bias factor according to equation (5)
                Bias=ifelse(res[,"marker",trial]>=s.ref, B.ref.s, 1/B.s.ref)
                
                ci.band=apply(res[,,trial], 1, function (x) quantile(x[3:length(x)], c(.025,.975)))
                mymatplot(res[,"marker",trial], t(rbind(res[,"prob",trial], ci.band)),      type="l", lty=c(1,2,2), col=4, lwd=lwd, make.legend=F, xlab="Month 13 Log10 Average Neutralizing Antibody Titer", ylab="Probability of VCD", main=toupper(trial))
                mymatplot(res[,"marker",trial], t(rbind(res[,"prob",trial], ci.band))*Bias, type="l", lty=c(1,2,2), col=3, lwd=lwd, make.legend=F, add=T)
                title(main="Marginalized and Controlled Effect Risk of Virologically Confirmed Dengue by Antibody Titer", outer=T)
                mylegend(x=3,legend=c("Marginalized risk", "Controlled effect risk (conservative)"), lty=1, col=c(4,3), lwd=2, cex=.8)
            }
        dev.off()    
    }
    
}
