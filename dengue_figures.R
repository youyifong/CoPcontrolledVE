rm(list=ls())
library(kyotil);           stopifnot(packageVersion("kyotil")>="2021.2-2")
library(marginalizedRisk); stopifnot(packageVersion("marginalizedRisk")>="2021.2-4") # bias.factor, E.value, controlled.risk.bias.factor
library(DengueTrialsYF) # need data to estimate overall attack rate to get s.ref
load(file=paste0("input/res_placebo_cont.Rdata")) # placebo arm results res.placebo.cont
RRud=RReu=4
bias.factor=bias.factor(RRud, RReu)
    

####################################################################################################
# categorical markers
    
setting="cat"
load(file=paste0("input/res_", setting, ".Rdata"))    
    
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
    


####################################################################################################
# continuous markers

setting="cont"
load(file=paste0("input/res_", setting, ".Rdata"))    
s2="85%"; s1="15%" # these two reference quantiles are used in the next two blocks of code

        
# generate text results comparing s1 and s1
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
mypdf(onefile=F, file=paste0("input/CoPveryhighVE_Fig2"), mfrow=c(1,2), oma=c(0,0,1,0), width=11, height=5)
    par(mar=c(4,5,3,2), las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    for (trial in c("cyd14","cyd15")) {        
    #trial="cyd14"
        lwd=2
        
        # choose a reference marker value
        # use median 
        #s.ref=res["50%","marker",trial]
        # use the alternative 
        dat=make.m13.dat(trial, stype=0)
        dat=subset(dat, trt==1)
        which=which.min(abs(res[,"prob",trial]-mean(dat$d)))
        s.ref=res[which,"marker",trial]
            
        Bias=controlled.risk.bias.factor(ss=res[,"marker",trial], s.cent=s.ref, s1=res[s1,"marker",trial], s2=res[s2,"marker",trial], RRud) 
        
        xlim=log10(c(20,3000))
        
        ci.band=apply(res[,,trial], 1, function (x) quantile(x[3:length(x)], c(.025,.975)))
        mymatplot(res[,"marker",trial], t(rbind(res[,"prob",trial], ci.band)),      type="l", lty=c(1,2,2), col=4, lwd=lwd, make.legend=F, xlab="Month 13 Average nAb Titer of Vaccinees", ylab="Probability of VCD", main=toupper(trial), ylim=ylim, xaxt="n", draw.x.axis=F, xlim=xlim)
        mymatplot(res[,"marker",trial], t(rbind(res[,"prob",trial], ci.band))*Bias, type="l", lty=c(1,2,2), col=3, lwd=lwd, make.legend=F, add=T)
        tmp=c(30,100,300,1000,3000)
        axis(side=1,at=log10(tmp),labels=tmp)
        title(main="Marginalized and Controlled Risk of Virologically Confirmed Dengue by Antibody Titer", outer=T)
        mylegend(x=3,legend=c("Marginalized risk", "Controlled risk (conservative)"), lty=1, col=c(4,3), lwd=2, cex=.8)
    
        # add histogram
        par(new=TRUE) 
        col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
        col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
        tmp=hist(dat$titer[dat$d==0],plot=F,breaks=15)    
        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,ylim=c(0,1.25*max(tmp$density)), xlim=xlim)    
        #axis(side=4, at=axTicks(side=4)[1:5])
        #mtext("Density", side=4, las=0, line=2, cex=1, at=.3)  
        #mylegend(x=6, fill=col, border=col, legend="Vaccine Group", bty="n", cex=0.7)  
    }
dev.off()    
        


# Fig 3 controlled VE curves
mypdf(onefile=F, file=paste0("input/CoPveryhighVE_Fig3"), mfrow=c(1,2), oma=c(0,0,1,0), width=11, height=5)
    par(mar=c(4,5,3,2), las=1, cex.axis=0.9, cex.lab=1)# axis label orientation
    for (trial in c("cyd14","cyd15")) {        
    # trial="cyd14"
        lwd=2.5
        
        # compute Bias as a vector, which is a function of s. equation (6) = Bias * marginalized risk
        # choose a reference marker value
        # use median 
        #s.ref=res["50%","marker",trial]
        # use the alternative 
        dat=make.m13.dat(trial, stype=0)
        dat=subset(dat, trt==1)
        mean(dat$d)
        which=which.min(abs(res[,"prob",trial]-mean(dat$d)))
        s.ref=res[which,"marker",trial]
        Bias=controlled.risk.bias.factor(ss=res[,"marker",trial], s.cent=s.ref, s1=res[s1,"marker",trial], s2=res[s2,"marker",trial], RRud) 
    
        ylim=c(0, max(hist(dat$titer[dat$d==0],breaks=10,plot=FALSE)$density, 1))
        xlim=log10(c(20,3000))
    
        # CVE
        est = 1 - res[,"prob",trial]*Bias/res.placebo.cont["est",trial]
        boot = 1 - t(  t( res[,3:ncol(res),trial]*Bias )/res.placebo.cont[2:nrow(res.placebo.cont),trial])                         
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
        mymatplot(res[,"marker",trial], t(rbind(est, ci.band)), type="l", lty=c(1,2,2), col="red", lwd=lwd, make.legend=F, xlab="Month 13 Average nAb Titer of Vaccinees", ylab="Vaccine Efficacy", 
            main=toupper(trial), ylim=ylim, xlim=xlim, yaxt="n", xaxt="n", draw.x.axis=F)
        tmp=c(30,100,300,1000,3000)
        axis(side=1,at=log10(tmp),labels=tmp)
        axis(side=2,at=seq(0,1,by=.25),labels=c('0','25','50','75','100')%.%"%")
        title(main="Controlled Vaccine Efficacy against Virologically Confirmed Dengue by Antibody Titer", outer=T)
        cve = cbind(10**res[,"marker",trial], t(rbind(est, ci.band))); cve
        # after comparing both trials, we choose 36 and 1200 (2 significant digits) as two titers and choose percentile based on that
        if (trial=="cyd15") {
            titer.1="5%"
            titer.2="79%"
        } else {
            titer.1="8%"
            titer.2="90%"
        }
        tmp=formatDouble(cve[titer.1, -1]*100, 1)
        write(paste0(tmp[1], "\\% (95\\% CI ", tmp[2], "--", tmp[3], ")"), file="input/fig3_low_"%.%trial)    
        tmp=formatDouble(cve[titer.2, -1]*100, 1)
        write(paste0(tmp[1], "\\% (95\\% CI ", tmp[2], "--", tmp[3], ")"), file="input/fig3_high_"%.%trial)    
        
        
        # VE
        est = 1 - res[,"prob",trial]/res.placebo.cont["est",trial]
        boot = 1 - t(  t( res[,3:ncol(res),trial] )/res.placebo.cont[2:nrow(res.placebo.cont),trial])                         
        ci.band=apply(boot, 1, function (x) quantile(x, c(.025,.975)))
        mymatplot(res[,"marker",trial], t(rbind(est, ci.band)), type="l", lty=c(1,2,2), col="pink", lwd=lwd, make.legend=F, add=T)
        mylegend(x=1,legend=c("Controlled VE Sens. Analysis","Controlled VE"), lty=1, col=c("red","pink"), lwd=2, cex=.8)
            
        # add histogram
        par(new=TRUE) 
        col <- c(col2rgb("olivedrab3")) # orange, darkgoldenrod2
        col <- rgb(col[1], col[2], col[3], alpha=255*0.4, maxColorValue=255)
        tmp=hist(dat$titer[dat$d==0],breaks=15,plot=F)    
        plot(tmp,col=col,axes=F,labels=F,main="",xlab="",ylab="",border=0,freq=F,ylim=ylim,xlim=xlim)    
        #axis(side=4, at=axTicks(side=4)[1:5]) # note that this is not correct if ylim is set
        #mtext("Density", side=4, las=0, line=2, cex=1, at=.3)  
        #mylegend(x=6, fill=col, border=col, legend="Vaccine Group", bty="n", cex=0.7)  
    
    }
dev.off()    
