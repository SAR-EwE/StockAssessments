#===============================================================================
# LEGACY BASELINE SCRIPT -- do not develop further.
#
# This is "BAM2EwE age0 corr.R" copied from the OneDrive BAM2EwE folder
# (last modified 2025-03-25), with only the minimal edits needed to run it
# end-to-end from the repo root and capture baseline outputs for regression
# testing the refactored pipeline. Edits are marked "BASELINE EDIT".
# Known bugs (t0 assignment, linear indexing in the total BA rows, etc.) are
# intentionally left in place so the baseline reflects the original code.
#
# Run from the repo root:  Rscript legacy/BAM2EwE_age0_corr.R
# Outputs are written to tests/baseline/
#===============================================================================
# BASELINE EDIT: removed interactive setup (rm(list=ls()), windows(), etc.)
library('xlsx')
library('plyr')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# setup----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# BASELINE EDIT: repo-relative input/output folders instead of getwd()
data.dir = "data/bam_rdats"
out.dir = "tests/baseline"
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)

#get list of rdat files in the working directory
bam.spp = data.frame(filename=list.files(path=data.dir,pattern="\\.rdat$",full.names = T))
bam.spp$spp = gsub(".rdat","",basename(bam.spp$filename))
bam.spp$EwEages = NA
bam.spp$enyr = bam.spp$styr = NA
for(i in 1:nrow(bam.spp)){
  #i=1
  dat.i = dget(bam.spp$filename[i])
  bam.spp$styr[i] = dat.i$parms$styr
  bam.spp$enyr[i] = dat.i$parms$endyr
  rm(dat.i);gc()
}
print(bam.spp)
print(basename(bam.spp[,1]))

#assign EwE stanzas to each species
# BASELINE EDIT: assign by file name instead of row position, because the
# mutton snapper (ASAP) file is not in data/bam_rdats and the positions shift.
# Same species/stanza assignments as the original.
bam.spp$EwEages[bam.spp$spp %in% c("GAG_SEDAR71","SNWY_SEDAR36_update")] = paste(0,1,5,sep=',')  #gag, snowy grouper
bam.spp$EwEages[bam.spp$spp %in% c("RSN_SEDAR73_update","VSN_SEDAR55_FESupdate","RGRP_SEDAR53_FESupdate","BSB_SEDAR76","GTLF_SEDAR66")] = paste(0,1,4,sep=',')  #red snapper, vermilion, red grouper, black seabass, golden tilefish
bam.spp$EwEages[bam.spp$spp %in% c("GAJ_SEDAR59","RPGY_SEDAR60")] = paste(0,1,3,sep=",") #gaj, red porgy
bam.spp$EwEages[bam.spp$spp %in% c("GTF_SEDAR82")] = 0
print(bam.spp)
bam.spp <- bam.spp[!is.na(bam.spp$EwEages),]
rownames(bam.spp) = 1:nrow(bam.spp)
print(bam.spp)

units.df = data.frame()
for(i in 1:nrow(bam.spp)){
  #i=7
  dat = dget(bam.spp$filename[i])
  tmp = (as.data.frame(dat$info[which(substr(names(dat$info),1,5)=='units')]))
  tmp$species = dat$info$species
  tmp = tmp[,c(ncol(tmp),2:(ncol(tmp)-1))]
  units.df = rbind.fill(units.df,tmp)
}
rm(tmp)
print(units.df)
#units.df$units.landings[1] = "1000 lb whole"
bam.spp <- cbind(bam.spp,units.df)
print(bam.spp)
write.csv(bam.spp,file.path(out.dir,'BAM list.csv'),row.names=F)  # BASELINE EDIT: out.dir

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# species loop----
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
spp.list = list()
indices.out <- data.frame(year=1995:2022)
indices.head <- data.frame()

for(s in 1:nrow(bam.spp)){
  # s=2  # BASELINE EDIT: commented out; it forced every iteration to species 2
  do.spp = s

##read rdat report file-----------------------------------------------------------------------------
dat = dget(bam.spp$filename[do.spp])
names(dat$size.age.fishery) = tolower(names(dat$size.age.fishery))
names(dat$t.series) = tolower(names(dat$t.series))
names(dat$sel.age) = tolower(names(dat$sel.age))

sppname = dat$info$species
units = (as.data.frame(dat$info[which(substr(names(dat$info),1,5)=='units')]))
yrs = dat$parms$styr:dat$parms$endyr
Nages = ncol(dat$N.age)
Nyrs = length(yrs)

eweages = as.numeric(unlist(strsplit(bam.spp$EwEages[do.spp],split=',')))
idx.stanzas = rep(eweages,diff(c(eweages,Nages)))
dat$a.series
#this is because BAM doesn't include age-0
#if(!sppname %in% c('SA gray triggerfish','Black sea bass')){
if(min(dat$a.series$age)==1 & length(eweages)>1){
  eweages = eweages[-1]
  idx.stanzas = rep(eweages,diff(c(eweages,Nages+1)))
}



##get biomass-at-age and sum over EwE ages------------------------------------------
#BAM B-at-age is Jan 1 numbers times spawning weight at age
#I calculate mid year B-at-age using 
names(dat)
naaj = dat$N.age #jan1 numbers
naa = dat$N.age.mdyr  #mid year numbers
waa = dat$a.series$wgt.mt #weight at age in mt
baa = sweep(naa,MARGIN=2,waa,'*') #mid year biomass in mt
baaj = sweep(naaj,MARGIN=2,waa,"*")
baa.stanzas = t(rowsum(t(baa),group=idx.stanzas)) 
baaj.stanzas = t(rowsum(t(baaj),group=idx.stanzas))


# naa.vec = colMeans(naa[which(yrs%in%c(1994,1995,1996)),])
# naa.vec2 = colMeans(naa)
# baa.vec = colMeans(baa[which(yrs%in%c(1994,1995,1996)),])
# baa.vec2 = colMeans(baa)
# plot(12*(1:Nages),naa.vec/max(naa.vec),type='l',col='red',ylab='Normalize Value',main=sppname,xlab='Age in months')
# lines(12*(1:Nages),naa.vec2/max(naa.vec2),col='red',lty=2)
# lines(12*(1:Nages),baa.vec/max(baa.vec),col='black')
# lines(12*(1:Nages),baa.vec2/max(baa.vec2),col='black',lty=2)
# lines(12*(1:Nages),dat$a.series$wgt.mt/max(dat$a.series$wgt.mt),col='blue')
# legend('topright',legend=c('Naa 3yr avg','Naa all avg','Baa 3yr avg','Baa all avg','Waa'),col=c('red','red','black','black','blue'),lty=c(1,2,1,2,1))


##get landings at age by fleet--------------------------------------------------------------
units.df
names(dat)
names(dat$CLD.est.mats)
flts.land = unique(gsub("Lw.","",gsub("Ln.","",names(dat$CLD.est.mats)[which(substr(names(dat$CLD.est.mats),1,1)=="L")])))
flts.land = flts.land[-which(flts.land=='total')]
flts.disc = unique(gsub("Dw.","",gsub("Dn.","",names(dat$CLD.est.mats)[which(substr(names(dat$CLD.est.mats),1,1)=="D")])))
flts.disc = flts.disc[-which(flts.disc=='total')]

if(units.df$units.landings[do.spp]=='1000 lb gutted'){
  gw2w.ratio = dat$a.series$weight/dat$a.series$weight.gut
  if(sppname == "SA Gag Grouper"){
    gw2w.ratio = rep(mean(gw2w.ratio[gw2w.ratio>1]),Nages)
  }
} else {gw2w.ratio = rep(1,Nages)}

land.age.klb = matrix(0, nrow=Nyrs, ncol=Nages)
land.stanzas.flt = data.frame(year=yrs)
land.stanzas = data.frame(matrix(0,nrow=Nyrs,ncol=length(eweages)))
dat$info$units.landings
for(f in 1:length(flts.land)){
  #f=1
  flt = flts.land[f]
  flt.m = as.data.frame(dat$CLD.est.mats[which(names(dat$CLD.est.mats)==paste0('Lw.',flt))])
  land.age.klb = land.age.klb+flt.m
  flt.m = sweep(flt.m,MARGIN=2,gw2w.ratio,"*")  #convert gutted weight to whole weight
  tmp = t(rowsum(t(flt.m),group=idx.stanzas))
  tmp = tmp*0.453592  #convert klb to mt
  colnames(tmp) = paste0(flt,'.',eweages)
  land.stanzas = land.stanzas+tmp
  land.stanzas.flt = cbind(land.stanzas.flt,tmp)
}
land.stanzas.flt = land.stanzas.flt[,-1]

disc.age.klb = matrix(0, nrow=Nyrs, ncol=Nages)
disc.stanzas.flt = data.frame(year=yrs)
disc.stanzas = data.frame(matrix(0,nrow=Nyrs,ncol=length(eweages)))
if(!sppname %in% c("Tilefish","SA Snowy Grouper")){
for(f in 1:length(flts.disc)){
  #f=1
  flt = flts.disc[f]
  flt.m = as.data.frame(dat$CLD.est.mats[which(names(dat$CLD.est.mats)==paste0('Dw.',flt))])
  disc.age.klb = disc.age.klb+flt.m
  flt.m = sweep(flt.m,MARGIN=2,gw2w.ratio,"*")  #convert gutted weight to whole weight
  tmp = t(rowsum(t(flt.m),group=idx.stanzas))
  tmp = tmp*0.453592  #convert klb to mt
  colnames(tmp)= paste0(flt,'.',eweages)
  disc.stanzas = disc.stanzas+tmp
  disc.stanzas.flt = cbind(disc.stanzas.flt,tmp)
}
}
disc.stanzas.flt = disc.stanzas.flt[,-1]

if(sppname %in% c("Tilefish","SA Snowy Grouper")){
  disc.stanzas.flt = land.stanzas.flt
  disc.stanzas.flt[,] = 0
}

##F-at-stanza as (Landings + dead discards)/mid year Biomass at stanza--------------------------------------------
faa.stanzas = (land.stanzas+disc.stanzas)/baa.stanzas
rownames(faa.stanzas) = yrs

##M-at-stanza--------------------------------------------------------
dat$a.series
Maa = dat$a.series$M
m.stanzas = matrix(NA,nrow(baa),ncol=length(eweages))
for(i in 1:nrow(m.stanzas)){
    #i=10
    df = data.frame(m = as.numeric(Maa),w = as.numeric(naa[i,]),grp=idx.stanzas)
    m.stanzas[i,] = sapply(split(df, df$grp), function(x){weighted.mean(x$m, w=x$w)}) 
}

##Z-at-age---------------------------------------------------------------
z.stanzas = faa.stanzas + m.stanzas
row.names(z.stanzas) = row.names(m.stanzas) = yrs

#name columns
colnames(baa.stanzas) = colnames(faa.stanzas) = colnames(z.stanzas) = colnames(m.stanzas) = colnames(land.stanzas) = colnames(disc.stanzas) = paste0('age',eweages)

##get other life history stuff----------------------------------------------------------
vbk = dat$parms$K
if(sppname %in% c('Red Porgy')){
  Wmat = min(dat$a.series$weight[dat$a.series$mat.fem.endyr>=0.5])
  Winf = max(dat$a.series$weight[dat$a.series$mat.fem.endyr>=0.5])
  Pmat = round(aggregate(mat.fem.endyr~idx.stanzas,dat=dat$a.series,mean),2)
}else{
  Wmat = min(dat$a.series$weight[dat$a.series$mat.female>=0.5])
  Winf = max(dat$a.series$weight[dat$a.series$mat.female>=0.5])
  Pmat = round(aggregate(mat.female~idx.stanzas,dat=dat$a.series,mean),2)
}
Wmat.Winf = Wmat/Winf

multistanza = data.frame(par=c('vbK','Wmat','Winf','Wmat.Winf',paste0('Pmat',eweages)),val=c(vbk,Wmat,Winf,Wmat.Winf,Pmat[,2]))

##check calcs -----------------------------------------------
names(dat$t.series)
names(dat)
t(units)
chk = data.frame(year = yrs)

#biomass
chk$ts.b.mt = dat$t.series$b[1:Nyrs]
chk$rowsum.B.age.mt = rowSums(dat$B.age[1:Nyrs,])
if(length(eweages)==1){
  chk$rowsum.Bcalc.syr.mt = baaj.stanzas[1:Nyrs,1]
  chk$rowsum.Bcalc.myr.mt = baa.stanzas[1:Nyrs,1]
} else{
  chk$rowsum.Bcalc.syr.mt = rowSums(baaj.stanzas[1:Nyrs,])
  chk$rowsum.Bcalc.myr.mt = rowSums(baa.stanzas[1:Nyrs,])
}
   
#landings
chk$ts.tot.l.klb = dat$t.series$total.l.klb[1:Nyrs]
if(!sppname %in% c("SA gray triggerfish")){
chk$rowsum.L.age.pred.klb = rowSums(dat$L.age.pred.klb)
}
chk$rowsum.L.CDL.total.klb    = rowSums(dat$CLD.est.mats$Lw.total)
chk$rowsum.calc.flt = rowSums(land.stanzas.flt)/.453592
chk$rowsum.calc.flt.mt = rowSums(land.stanzas.flt)

#discards
if(!sppname %in% c("Tilefish","SA Snowy Grouper")){
  chk$ts.tot.d.klb = dat$t.series$total.d.klb[1:Nyrs]
  chk$rowsum.D.CDL.total = rowSums(dat$CLD.est.mats$Dw.total)
  chk$rowsum.calc.disc.flt.klb =  rowSums(disc.stanzas.flt)/0.453592
  chk$rowsum.calc.disc.flt.mt = rowSums(disc.stanzas.flt)
}

##get indices--------------------------------------------------------------------
sppabb <- substr(basename(bam.spp$spp[s]),1,unlist(gregexpr("_",basename(bam.spp$spp[s])))[1]-1)
ts <- dat$t.series
names(ts) = tolower(names(ts))
ts.U <- ts[,c(1,which(substr(names(ts),1,1)=='u'),which(substr(names(ts),1,4)=='cv.u'))]
ts.U[ts.U<0] <- NA
ts.U.ob <- ts.U[,c(1,which(substr(names(ts.U),nchar(names(ts.U))-1,nchar(names(ts.U)))=='ob'))]
ts.U.cv <- ts.U[,c(1,which(substr(names(ts.U),1,4)=='cv.u'))]
#matplot(ts.U.ob[,1], ts.U.ob[,-1],type='l',col=1:4, lty=1, lwd=2)  
#legend('bottomleft',legend=names(ts.U.ob)[-1],col=1:4, lty=1, lwd=2)
names(ts.U.ob); names(ts.U.cv)
ageselex <- dat$sel.age
names(ageselex) = tolower(names(ageselex))
names(ageselex)[which(names(ageselex)=='sel.v.mcvt')] = 'sel.v.cvid'
ageselex <- ageselex[which(gsub("sel.m.","",gsub('sel.v.',"",names(ageselex))) %in% gsub("u.","",gsub(".ob","",names(ts.U.ob)[-1])))]

tmp.head <- data.frame(name=paste0(sppabb,"_",names(ts.U.ob)[-1]), cv=colMeans(ts.U.cv[,-1],na.rm=T), Asel50=NA, Asel95=NA, Aselmax=NA, sel50=NA, sel95=NA, selmax=NA)
for(i in 1:length(ageselex)){
  #i=2
  if(is.matrix(ageselex[[i]])) sel = colMeans(ageselex[[i]],na.rm=T)
  if(is.vector(ageselex[[i]])) sel = ageselex[[i]]
  sel
  selmax = sel[which.min(abs(1-sel) )]
  sel2 = sel[1:min(which(sel==selmax))]
  sel50 = sel2[which.min(abs(.5-sel2))]
  sel95 = sel2[which.min(abs(.95-sel2))]
  Aselmax = as.numeric(names(sel)[min(which(sel==selmax))])
  Asel50 = as.numeric(names(sel)[min(which(sel==sel50))])
  Asel95 = as.numeric(names(sel)[min(which(sel==sel95))])
  tmp.head$Asel50[i] = Asel50
  tmp.head$Aselmax[i] = Aselmax
  tmp.head$sel50[i] = sel50
  tmp.head$selmax[i] = selmax
  tmp.head$sel95[i] = sel95
  tmp.head$Asel95[i] = Asel95
  
}

names(ts.U.ob)[-1] <- paste0(sppabb,"_",names(ts.U.ob)[-1])
indices.out <- merge(indices.out,ts.U.ob,all=T)
indices.head <- rbind(indices.head,tmp.head)

##save in species list----
spp.out.tmp = list(B_mt=baa.stanzas, F=faa.stanzas, M=m.stanzas, L_mt=land.stanzas, L_mt_flt=land.stanzas.flt, 
               D_mt=disc.stanzas, D_mt_flt=disc.stanzas.flt,
               indices=ts.U.ob, indices.cv=tmp.head, multistanza=multistanza)
spp.list[[s]] = spp.out.tmp
names(spp.list)[s] = sppabb

#export----------------------------------------------------------------------------
if(min(dat$a.series$age)==1) outfile = file.path(out.dir,paste0(sppname,"_no age0.xlsx"))  # BASELINE EDIT: out.dir
if(min(dat$a.series$age)==0) outfile = file.path(out.dir,paste0(sppname,".xlsx"))  # BASELINE EDIT: out.dir
write.xlsx(baa.stanzas,outfile,sheetName='Biomass mt',row.names=T)
write.xlsx(faa.stanzas,outfile,sheetName='F rates',row.names=T,append=T)
write.xlsx(z.stanzas,outfile,sheetName='Z rates',row.names=T,append=T)
write.xlsx(m.stanzas,outfile,sheetName='M rates',row.names=T,append=T)
write.xlsx(land.stanzas,outfile,sheetName='landings mt',row.names=T,append=T)
write.xlsx(land.stanzas.flt,outfile,sheetName='landings mt fleet',row.names=T,append=T)
#if(!sppname %in% c("Tilefish")){
write.xlsx(disc.stanzas.flt,outfile,sheetName='discards mt fleet',row.names=T,append=T)
#}
write.xlsx(multistanza,outfile,sheetName='multistanza',row.names=F,append=T)
write.xlsx(chk,outfile,sheetName='check',row.names=F,append=T)
rm(land.stanzas.flt,disc.stanzas.flt,disc.stanzas,land.stanzas)
if(s==nrow(bam.spp)){
  write.csv(indices.out,file.path(out.dir,'indices.csv'),row.names=F)  # BASELINE EDIT: out.dir
  write.csv(indices.head,file.path(out.dir,'indices_header.csv'),row.names=F)  # BASELINE EDIT: out.dir
}
}

#get Ecopath inputs---------------------------------------------------------------------------------
eco.in <- data.frame()
baseyear<-1995
modelarea_km2 = 336192
for(i in 1:length(spp.list)){
  #i = 1
  spp = spp.list[[i]]
  spp.abb = names(spp.list)[i]
  baseyear.idx = which(rownames(spp$B_mt)==baseyear)
  nstz = ncol(spp$B_mt)
  B.stz.bar = spp$B_mt[,1:nstz]
  F.stz.bar = spp$F[,1:nstz]
  L.stz = spp$L_mt[,1:nstz]
  M.stz = spp$M[,1:nstz]
  D.stz = spp$D_mt[,1:nstz]
  
  if(nstz==1){
    spp.in <- data.frame(species = spp.abb,
                         stanza = 'total',
                         biomass_mt = B.stz.bar[baseyear.idx],
                         biomass_mtkm2 = B.stz.bar[baseyear.idx]/modelarea_km2,
                         F = F.stz.bar[baseyear.idx],
                         M = M.stz[baseyear.idx],
                         Z = F.stz.bar[baseyear.idx]+M.stz[baseyear.idx],
                         landings_mt = L.stz[baseyear.idx],
                         landings_mtkm2 = L.stz[baseyear.idx]/modelarea_km2,
                         discards_mt = D.stz[baseyear.idx],
                         discards_mtkm2 = D.stz[baseyear.idx]/modelarea_km2,
                         BA.1985 = (B.stz.bar[baseyear.idx+1]/B.stz.bar[baseyear.idx])-1,
                         BA.3yr = sum(B.stz.bar[baseyear.idx:(baseyear.idx+2)])/sum(B.stz.bar[(baseyear.idx-min(3,baseyear.idx-1)):(baseyear.idx-min(1,baseyear.idx-1))])-1)
    
  }else{
    spp.in <- data.frame(species = rep(spp.abb,nstz+1),
                         stanza=c(1:nstz,'total'),
                         biomass_mt=c(B.stz.bar[baseyear.idx,],sum(B.stz.bar[baseyear.idx,])),
                         biomass_mtkm2=c(B.stz.bar[baseyear.idx,],sum(B.stz.bar[baseyear.idx,]))/modelarea_km2,
                         F=as.numeric(c(F.stz.bar[baseyear.idx,],sum(L.stz[baseyear.idx,],D.stz[baseyear.idx,])/sum(B.stz.bar[baseyear.idx,]))),
                         M=as.numeric(c(M.stz[baseyear.idx,],NA)),
                         Z=as.numeric(c(F.stz.bar[baseyear.idx,]+M.stz[baseyear.idx,],NA)),
                         landings_mt=as.numeric(c(L.stz[baseyear.idx,],sum(L.stz[baseyear.idx,]))),
                         landings_mtkm2=as.numeric(c(L.stz[baseyear.idx,],sum(L.stz[baseyear.idx,])))/modelarea_km2,
                         discards_mt=as.numeric(c(D.stz[baseyear.idx,],sum(D.stz[baseyear.idx,]))),
                         discards_mtkm2=as.numeric(c(D.stz[baseyear.idx,],sum(D.stz[baseyear.idx,])))/modelarea_km2,
                         BA.1985=c(apply(B.stz.bar,2,FUN=function(x){x[baseyear.idx+1]/x[max(1,baseyear.idx)]-1}),
                                   sum(B.stz.bar[baseyear.idx+1])/sum(B.stz.bar[baseyear.idx])-1),
                         BA.3yr=c(apply(B.stz.bar,2,FUN=function(x){sum(x[(baseyear.idx):(baseyear.idx+2)])/sum(x[(baseyear.idx-min(3,baseyear.idx-1)):(baseyear.idx-min(1,baseyear.idx-1))])-1}),
                                  sum(B.stz.bar[(baseyear.idx):(baseyear.idx+2)])/sum(B.stz.bar[(baseyear.idx-min(3,baseyear.idx-1)):(baseyear.idx-min(1,baseyear.idx-1))])-1))
  }
  eco.in <- rbind(eco.in,spp.in)
  
}
# BASELINE EDIT: dir.base was undefined; write to out.dir without a date stamp
write.csv(eco.in,file.path(out.dir,"stock assessment Ecopath inputs.csv"),row.names=F)



###################################################################################################
# loop through species and do the alternative M calculations
###################################################################################################
print(units.df)
for(s in 1:nrow(bam.spp)){
  #s=7
  do.spp = s
  eweages = as.numeric(unlist(strsplit(bam.spp$EwEages[do.spp],split=',')))
  #read rdat report file-----------------------------------------------------------------------------
  dat = dget(bam.spp$filename[do.spp])
  names(dat$size.age.fishery) = tolower(names(dat$size.age.fishery))
  names(dat$t.series) = tolower(names(dat$t.series))
  names(dat$sel.age) = tolower(names(dat$sel.age))
  
  sppname = dat$info$species
  units = (as.data.frame(dat$info[which(substr(names(dat$info),1,5)=='units')]))
  Nages = ncol(dat$N.age)

  #new M calc based on monthly ages
  names(dat)
  age.mo = 0:(Nages*12)
  age.yr = age.mo/12
  Lmax = dat$parms$Linf/10
  vbk = dat$parms$K
  t0 = ifelse(dat$parms$t0<-1,-0.5,dat$parms$t0)
  lwa = dat$parms$wgt.a
  lwb = dat$parms$wgt.b
  len.cm = Lmax*(1-exp(-vbk*(age.yr-t0)))
  
  if(sppname %in% c('Black sea bass')){
    Mtarget = round(c(tail(dat$a.series$M,1),tail(dat$a.series$M,1)/12),3)
    wt.kg = lwa*len.cm^lwb
  } else if(sppname %in% c('Tilefish')){
    Mtarget = round(c(dat$parms$M.constant,dat$parms$M.constant/12),3)
    wt.kg = 1e6*lwa*len.cm^lwb
  } else if(sppname %in% c('SA Greater Amberjack')){
    #Mtarget = round(c(dat$parms$M.constant,dat$parms$M.constant/12),3)
    Mtarget = round(c(tail(dat$a.series$M,1),tail(dat$a.series$M,1)/12),3)
    wt.kg = lwa*len.cm^lwb
  } else if(sppname %in% c('SA Vermilion Snapper')){
    Mtarget = round(c(dat$parms$M.constant,dat$parms$M.constant/12),3)
    wt.kg = lwa*len.cm^lwb
  } else{
    Mtarget = round(c(dat$parms$M.constant,dat$parms$M.constant/12),3)
    wt.kg = 1000*lwa*len.cm^lwb
  }
  tc = head(which(colMeans(dat$F.age)/tail(colMeans(dat$F.age),1)>=0.5),1)
  tc = c(tc,tc*12)
  tmax = tail(which(colMeans(dat$F.age)/tail(colMeans(dat$F.age),1)>=0.5),1)
  tmax = c(tmax,tmax*12)
  nc = tmax-tc

  M.yr = 3*(wt.kg*1000)^-0.288
  M.mo = M.yr/12
  M.mo.scaled = M.yr.scaled = rep(NA,length(age.mo))
  for(i in 1:length(M.mo.scaled)){
    M.mo.scaled[i] = Mtarget[2]*((nc[2]*M.mo[i])/sum(M.mo[age.mo>=tc[2]&age.mo<=tmax[2]]))
    #M.yr.scaled[i] = Mtarget[1]*((nc[1]*M.mo[i])/sum(M.mo[age.yr>=tc[1]&age.yr<=tmax[1]]))
  }
  M.yr.scaled = M.mo.scaled*12
  
  survship = survship.scaled = rep(NA,length(age.mo))
  survship[1] = survship.scaled[1] = 1
  for(i in 2:length(survship)){
    survship[i] = survship[i-1]*exp(-M.mo[i-1])
    survship.scaled[i] = survship.scaled[i-1]*exp(-M.mo.scaled[i-1])
  }
  
  stanza.start = c(eweages*12,max(age.yr)*12)
  eweages.yr = stanza.start/12
  midage = M.wtavg.yr = M.styr = M.midyr = M.wtavg.yr.scaled = M.styr.scaled = M.midyr.scaled = rep(NA,length(eweages))
  for(i in 1:length(eweages)){
    #i=3
    midage[i] = eweages.yr[i]+(eweages.yr[i+1]-eweages.yr[i])/2
    idx.ages = which(age.yr>=eweages.yr[i]&age.yr<eweages.yr[i+1])
    #unscaled
    M.wtavg.yr[i] = weighted.mean(M.yr[idx.ages],w=survship[idx.ages])
    M.styr[i] = M.yr[which(age.yr==eweages.yr[i])]
    M.midyr[i] = M.yr[which(age.yr==midage[i])]
    #scaled
    M.wtavg.yr.scaled[i] = weighted.mean(M.yr.scaled[idx.ages],w=survship.scaled[idx.ages])
    M.styr.scaled[i] = M.yr.scaled[which(age.yr==eweages.yr[i])]
    M.midyr.scaled[i] = M.yr.scaled[which(age.yr==midage[i])]
  }
  Mcalc.out = data.frame(startage.mo=stanza.start[1:length(eweages)],startage=eweages.yr[1:length(eweages)],
                         midage,M.wtavg.yr,M.styr,M.midyr,M.wtavg.yr.scaled,M.styr.scaled,M.midyr.scaled)
  
  write.xlsx(Mcalc.out,file.path(out.dir,'alternative M calcs.xlsx'),sheetName=sppname,row.names=F,append=ifelse(s==1,F,T))  # BASELINE EDIT: out.dir
}
  





      
