# nonparametric estimation of placebo risk for revision
library(kyotil)

label="tmlestsl"# tmlestglm, tmlestsl
res=get.sim.res("input/"%.%label, verbose=1)

res.placebo.cont = t(res[1,,])
rownames(res.placebo.cont)[1]="est"
save(res.placebo.cont, file="input/res_placebo_cont_"%.%label%.%".Rdata")

res.vaccine.cont = t(res[2,,])
rownames(res.vaccine.cont)[1]="est"
save(res.vaccine.cont, file="input/res_vaccine_cont_"%.%label%.%".Rdata")

est = 1 - res.vaccine.cont/res.placebo.cont
est[1,]
ci=apply(est[-1,], 2, function (x) quantile(x, c(.025,.975)))*100
ci
s1=paste0(formatDouble(est[1,"cyd14"]*100,1,remove.leading0=F), "\\% (95\\% CI ", formatDouble(ci[1,"cyd14"],1,remove.leading0=F), "--", formatDouble(ci[2,"cyd14"],1,remove.leading0=F), ")"); s1
write(s1, file="input/cyd14_ve_"%.%label)
s2=paste0(formatDouble(est[1,"cyd15"]*100,1,remove.leading0=F), "\\% (95\\% CI ", formatDouble(ci[1,"cyd15"],1,remove.leading0=F), "--", formatDouble(ci[2,"cyd15"],1,remove.leading0=F), ")"); s2
write(s2, file="input/cyd15_ve_"%.%label)
