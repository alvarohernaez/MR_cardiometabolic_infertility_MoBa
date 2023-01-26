rm(list=ls())

library(devtools)
library(googleAuthR)
library(MendelianRandomization)
library(mr.raps)
library(meta)
library(MRPRESSO)
library(MRInstruments)
library(MRMix)
library(RadialMR)
library(ieugwasr)
library(TwoSampleMR)


### SELF-MADE FUNCTIONS TO OBTAIN CLEAN ESTIMATES ###
#####################################################

guapa<-function(x)
{
  redondeo<-ifelse(abs(x)<0.00001,signif(x,1),
                   ifelse(abs(x)<0.0001,signif(x,1),
                          ifelse(abs(x)<0.001,signif(x,1),
                                 ifelse(abs(x)<0.1,sprintf("%.3f",round(x,3)),
                                        ifelse(abs(x)<1,sprintf("%.2f",round(x,2)),
                                               ifelse(abs(x)<10,sprintf("%.2f",round(x,2)),
                                                      ifelse(abs(x)<100,sprintf("%.1f",round(x,1)),
                                                             ifelse(abs(x)>=100,round(x,0),round(x,0)))))))))
  return(redondeo)
}

ic_guapa<-function(x,y,z)
{
  ic<-paste(x," [",y,"; ",z,"]",sep="")
  return(ic)
}

pval_guapa<-function(x)
{
  pval<-ifelse(x<0.00001,"<0.00001",
               ifelse(x<0.001,"<0.001",round(x,3)))
  return(pval)
}

mean_ic_guapa <- function(x, na.rm=FALSE) 
{
  if (na.rm) x <- na.omit(x)
  se<-sqrt(var(x)/length(x))
  z<-qnorm(1-0.05/2)
  media<-mean(x)
  ic95a<-guapa(media-(z*se))
  ic95b<-guapa(media+(z*se))
  media<-guapa(media)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

mean_sd_guapa <- function(x) 
{
  media<-guapa(mean(x, na.rm=TRUE))
  sd<-guapa(sd(x, na.rm=TRUE))
  end<-paste(media," (",sd,")",sep="")
  return(end)
}

beta_se_ic_guapa <- function(x, y) 
{
  z<-qnorm(1-0.05/2)
  ic95a<-guapa(x-(z*y))
  ic95b<-guapa(x+(z*y))
  media<-guapa(x)
  ic_ok<-ic_guapa(media,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-guapa(exp(x))
  ic95a<-guapa(exp(x-(z*y)))
  ic95b<-guapa(exp(x+(z*y)))
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

risk_se_ic_guapa2 <- function(x,y) 
{
  z<-qnorm(1-0.05/2)
  hr<-round(exp(x),5)
  ic95a<-round(exp(x-(z*y)),5)
  ic95b<-round(exp(x+(z*y)),5)
  ic_ok<-ic_guapa(hr,ic95a,ic95b)
  return(ic_ok)
}

header.true <- function(df)
{
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}

closest<-function(xv,sv){
  xv[which(abs(xv-sv)==min(abs(xv-sv)))] }

setwd("N:/data/durable/Projects/Hernaez_MR_CVRF")


################################
### DATABASE FOR DESCRIPTION ###
################################

### MERGING GRSs ###

load("N:/data/durable/Projects/Hernaez_MR_BMI/R/MoBa_raw.RData")
dat$sentrixid_mom<-NULL
dat$sentrixid_dad<-NULL


### BMI AND HEIGHT MEASURED VALUES ###

dat$bmi_mom<-dat$aa85/((dat$aa87/100)^2)
dat$bmi_dad<-dat$aa89/((dat$aa88/100)^2)
dat$bmi_mom<-with(dat,ifelse(bmi_mom<13 | bmi_mom>60,NA,bmi_mom))
dat$bmi_dad<-with(dat,ifelse(bmi_dad<13 | bmi_dad>60,NA,bmi_dad))


### HIGHEST EDUCATIONAL LEVEL, COMPLETED OR ONGOING ###

# Mothers #

dat$aa1125<-with(dat,ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125>aa1124),aa1125,
                            ifelse(!is.na(aa1124) & !is.na(aa1125) & (aa1125<aa1124),NA,aa1125)))
dat$aa1128x<-with(dat,ifelse(aa1128==1 | aa1129==1,1,0))
dat$aa1128x<-with(dat,ifelse(is.na(aa1128x),0,aa1128x))

dat$edu_mom<-with(dat,ifelse(aa1124==1,1,
                             ifelse(aa1124==2,2,
                                    ifelse(aa1124==3,3,
                                           ifelse(aa1124==4,4,
                                                  ifelse(aa1124==5,5,
                                                         ifelse(aa1124==6,6,NA)))))))
dat$edu_mom<-with(dat,ifelse(is.na(aa1125),edu_mom,
                             ifelse(aa1125==1,1,
                                    ifelse(aa1125==2,2,
                                           ifelse(aa1125==3,3,
                                                  ifelse(aa1125==4,4,
                                                         ifelse(aa1125==5,5,
                                                                ifelse(aa1125==6,6,NA))))))))
dat$edu_mom<-with(dat,ifelse(is.na(edu_mom),0,edu_mom))

dat$edu_mom<-with(dat,ifelse(edu_mom==0 & aa1128x==1,3,
                             ifelse((edu_mom>0 & edu_mom<3) & aa1128x==1,3,
                                    ifelse(edu_mom>=3 & aa1128x==1,edu_mom,
                                           ifelse(edu_mom==0 & aa1128x==0,NA,
                                                  ifelse((edu_mom>0 & edu_mom<3) & aa1128x==0,edu_mom,
                                                         ifelse(edu_mom>=3 & aa1128x==0,edu_mom,NA)))))))

dat$eduyears_mom<-with(dat,ifelse(edu_mom==1,10,
                                  ifelse(edu_mom==2,10,
                                         ifelse(edu_mom==3,13,
                                                ifelse(edu_mom==4,13,
                                                       ifelse(edu_mom==5,19,
                                                              ifelse(edu_mom==6,20,
                                                                     ifelse(is.na(edu_mom),NA,NA))))))))

dat$edu_mom<-with(dat,ifelse(edu_mom==1,1,
                             ifelse(edu_mom==2,1,
                                    ifelse(edu_mom==3,2,
                                           ifelse(edu_mom==4,2,
                                                  ifelse(edu_mom==5,3,
                                                         ifelse(edu_mom==6,4,
                                                                ifelse(is.na(edu_mom),NA,NA))))))))

# Fathers #

dat$aa1127<-with(dat,ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127>aa1126),aa1127,
                            ifelse(!is.na(aa1126) & !is.na(aa1127) & (aa1127<aa1126),NA,aa1127)))
dat$aa1130x<-with(dat,ifelse(aa1130==1 | aa1131==1,1,0))
dat$aa1130x<-with(dat,ifelse(is.na(aa1130x),0,aa1130x))

dat$edu_dad<-with(dat,ifelse(aa1126==1,1,
                             ifelse(aa1126==2,2,
                                    ifelse(aa1126==3,3,
                                           ifelse(aa1126==4,4,
                                                  ifelse(aa1126==5,5,
                                                         ifelse(aa1126==6,6,NA)))))))
dat$edu_dad<-with(dat,ifelse(is.na(aa1127),edu_dad,
                             ifelse(aa1127==1,1,
                                    ifelse(aa1127==2,2,
                                           ifelse(aa1127==3,3,
                                                  ifelse(aa1127==4,4,
                                                         ifelse(aa1127==5,5,
                                                                ifelse(aa1127==6,6,NA))))))))
dat$edu_dad<-with(dat,ifelse(is.na(edu_dad),0,edu_dad))

dat$edu_dad<-with(dat,ifelse(edu_dad==0 & aa1130x==1,3,
                             ifelse((edu_dad>0 & edu_dad<3) & aa1130x==1,3,
                                    ifelse(edu_dad>=3 & aa1130x==1,edu_dad,
                                           ifelse(edu_dad==0 & aa1130x==0,NA,
                                                  ifelse((edu_dad>0 & edu_dad<3) & aa1130x==0,edu_dad,
                                                         ifelse(edu_dad>=3 & aa1130x==0,edu_dad,NA)))))))

dat$eduyears_dad<-with(dat,ifelse(edu_dad==1,10,
                                  ifelse(edu_dad==2,10,
                                         ifelse(edu_dad==3,13,
                                                ifelse(edu_dad==4,13,
                                                       ifelse(edu_dad==5,19,
                                                              ifelse(edu_dad==6,20,
                                                                     ifelse(is.na(edu_dad),NA,NA))))))))
dat$edu_dad<-with(dat,ifelse(edu_dad==1,1,
                             ifelse(edu_dad==2,1,
                                    ifelse(edu_dad==3,2,
                                           ifelse(edu_dad==4,2,
                                                  ifelse(edu_dad==5,3,
                                                         ifelse(edu_dad==6,4,
                                                                ifelse(is.na(edu_dad),NA,NA))))))))


### DEFINITION OF PARITY (number of previous deliveries) ###

dat$parity<-with(dat,ifelse(paritet_5==0,0,
                            ifelse(paritet_5==1,1,
                                   ifelse(paritet_5==2,2,
                                          ifelse(paritet_5==3,3,
                                                 ifelse(paritet_5==4,3,NA))))))

length(which(is.na(dat$parity)))


### DEFINITION OF EVER SMOKERS (smkinit) ###

# Mothers #
dat$smkinit_mom<-with(dat,ifelse(aa1355==2 | aa1356>1 | aa1357!=0 | aa1358!=0 | aa1359>1 | aa1360!=0 | aa1361!=0 | 
                                   !is.na(aa1362) | aa1363==2 | !is.na(aa1364) | !is.na(aa1365),1,0))
dat$smkinit_mom<-with(dat,ifelse(is.na(smkinit_mom),0,smkinit_mom))
dat$smkinit_mom<-with(dat,ifelse(is.na(aa1355) & is.na(aa1356) & is.na(aa1357) & is.na(aa1358) & is.na(aa1359) & is.na(aa1360) & is.na(aa1361) & 
                                   is.na(aa1362) & is.na(aa1363) & is.na(aa1364) & is.na(aa1365),NA,smkinit_mom))

# Fathers #
dat$smkinit_dad<-with(dat,ifelse(aa1353==2 | aa1354==2 | ff214==2 | ff215>1 | ff216!=0 | ff217!=0 | ff218>1 | ff219!=0 | ff220!=0,1,0))
dat$smkinit_dad<-with(dat,ifelse(is.na(smkinit_dad),0,smkinit_dad))
dat$smkinit_dad<-with(dat,ifelse(is.na(aa1353) & is.na(aa1354) & is.na(ff214) & is.na(ff215) & is.na(ff216) & is.na(ff217) & is.na(ff218) & is.na(ff219) & is.na(ff220),NA,smkinit_dad))


### DEFINITION OF SUBFERTILITY ###

dat$exclude<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))
dat$art<-with(dat,ifelse(!is.na(art) & art>0,1,0))
dat$preg_plan<-with(dat,ifelse(is.na(aa46),9,
                               ifelse(aa46==0,9,
                                      ifelse(aa46==1,0,
                                             ifelse(aa46==2,1,9)))))
attributes(dat$preg_plan)$value.label<-c("0=No","1=Yes","9=NA")
attributes(dat$preg_plan)$label<-c("Planned pregnancy")

dat$subf<-with(dat,ifelse(aa48<12,0,
                          ifelse(aa48>=12,1,NA)))
dat$subf<-with(dat,ifelse(art==0,subf,
                          ifelse(art==1,1,NA)))
dat$subf<-with(dat,ifelse(exclude==0 & is.na(subf),0,
                          ifelse(exclude==0 & subf==0,0,
                                 ifelse(exclude==0 & subf==1,1,
                                        ifelse(exclude==1 & is.na(subf),NA,
                                               ifelse(exclude==1 & subf==0,NA,
                                                      ifelse(exclude==1 & subf==1,NA,NA)))))))


### FINAL FORMAT OF DATA ###

dat<-rename.vars(dat,
                 from=c("mors_alder","fars_alder"),
                 to=c("agedelivery_mom","agedelivery_dad"))
attributes(dat$agedelivery_mom)$label<-c("agedelivery_mom")
attributes(dat$agedelivery_dad)$label<-c("agedelivery_dad")

dat$exclude<-with(dat,ifelse(flerfodsel==1,1,0))
dat$exclude<-with(dat,ifelse(is.na(exclude),0,exclude))
table(dat$exclude)
dat<-subset2(dat,"dat$exclude==0")

dat$exclude<-with(dat,ifelse(versjon_skjema1_tbl1=="",1,0))
dat$exclude<-with(dat,ifelse(is.na(exclude),0,exclude))
table(dat$exclude)
dat<-subset2(dat,"dat$exclude==0")

nam_ok<-c("sentrixid_mom","sentrixid_dad","m_id_2374","f_id_2374","preg_id_2374",
          "agedelivery_mom","agedelivery_dad","bmi_mom","bmi_dad","eduyears_mom","eduyears_dad","smkinit_mom","smkinit_dad","parity",
          "subf","preg_plan","art")
dat<-dat[,nam_ok]

attributes(dat$bmi_mom)$label<-c("bmi_mom")
attributes(dat$bmi_dad)$label<-c("bmi_dad")


excl<-spss.get("N:/data/durable/RAW/1.MoBa_questionnaires_MBRN/INFO_files_start_your_analysis_here/PDB2374_SV_INFO_V12_20220824.sav",
               use.value.labels=FALSE,to.data.frame=TRUE,allow="_")
names(excl)<-tolower(names(excl))
excl$consent<-1
excl<-excl[,c("preg_id_2374","consent")]
dat<-merge2(dat,excl,by.id=c("preg_id_2374"),all.x=TRUE,sort=FALSE)
dat$consent<-with(dat,ifelse(is.na(consent),0,consent))
table(dat$consent)
dat<-dat[dat$consent==1,]

save(dat,file="N:/data/durable/Projects/Hernaez_MR_CVRF/MoBa_cvrf.RData")


### DEFINITION OF STUDY FLOW CHART ###

# Total singleton pregnancies
length(which(!is.na(dat$preg_id_2374)))

# Unique IDs of mothers/fathers
length(unique(dat[which(!is.na(dat$subf==0) & dat$preg_plan==1),c("m_id_2374")]))
length(unique(dat[which(!is.na(dat$subf==0) & dat$preg_plan==1),c("f_id_2374")]))

# IDs of mothers/fathers with genotype data
length(dat[which(!is.na(dat$tc_grs_mom) & dat$preg_plan==1),c("sentrixid_mom")])
length(dat[which(!is.na(dat$tc_grs_dad) & dat$preg_plan==1),c("sentrixid_dad")])

# IDs of unique mothers/fathers with genotype data
length(unique(dat[which(!is.na(dat$tc_grs_mom) & dat$preg_plan==1),c("sentrixid_mom")]))
length(unique(dat[which(!is.na(dat$tc_grs_dad) & dat$preg_plan==1),c("sentrixid_dad")]))




#####################
### MAIN ANALYSES ###
#####################

setwd("N:/data/durable/Projects/Hernaez_MR_CVRF")

dir.create("./Outputs")
dir.create("./Outputs/descriptive")
dir.create("./Outputs/results")

setwd("N:/data/durable/Projects/Hernaez_MR_CVRF/Outputs")


### POPULATION DESCRIPTION ###
##############################

load("N:/data/durable/Projects/Hernaez_MR_CVRF/MoBa_cvrf.RData")
dat$parity<-with(dat,ifelse(parity==0,1,2))

datx<-subset2(dat,"!is.na(dat$tc_grs_mom) & preg_plan==1")
xxx<-datx[,c("agedelivery_mom","eduyears_mom","bmi_mom","smkinit_mom","parity","subf")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.
                               -subf,
                               xxx, method=c("bmi_mom"=2,"smkinit_mom"=3,"parity"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=0)

comp<-NULL
comp<-createTable(compareGroups(subf~.
                                -sel,
                                xxx, method=c("bmi_mom"=2,"smkinit_mom"=3,"parity"=3)),
                  show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=0)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1],comp$descr))
colnames(tab1)<-c("Mothers-All","Mothers-Non-subfertile","Mothers-Subfertile","Mothers-P-value","Mothers-N")
write.table(tab1,file="./descriptive/descriptive_mothers.csv",sep=";",col.names=NA)


datx<-subset2(dat,"!is.na(dat$tc_grs_dad) & preg_plan==1")
xxx<-datx[,c("agedelivery_dad","eduyears_dad","bmi_dad","smkinit_dad","parity","subf")]
xxx$sel<-1

all<-NULL
all<-createTable(compareGroups(sel~.
                               -subf,
                               xxx, method=c("bmi_dad"=2,"smkinit_dad"=3,"parity"=3)),
                 show.n=TRUE, show.p.overall=FALSE, show.p.trend=FALSE, hide.no=NA)

comp<-NULL
comp<-createTable(compareGroups(subf~.
                                -sel,
                                xxx, method=c("bmi_dad"=2,"smkinit_dad"=3,"parity"=3)),
                  show.n=TRUE, show.p.overall=TRUE, show.p.trend=FALSE, hide.no=NA)

tab1<-NULL
tab1<-as.data.frame(cbind(all$descr[,1],comp$descr))
colnames(tab1)<-c("Fathers-All","Fathers-Non-subfertile","Fathers-Subfertile","Fathers-P-value","Fathers-N")
write.table(tab1,file="./descriptive/descriptive_fathers.csv",sep=";",col.names=NA)



###############################################
### SENSITIVITY ANALYSES: STEIGER FILTERING ###
###############################################

setwd("N:/data/durable/Projects/Hernaez_MR_CVRF")

### GETTING READY THE TOP HIT DATABASES FOR THE STEIGER FILTERING ###
#####################################################################

### LDL CHOLESTEROL ###

ldlc<-read.csv2("./gwas/ldlc_graham2021/ldlc_graham2021.csv",header=TRUE,sep=";",dec=".")
ldlc<-rename.vars(ldlc,
                from=c("rsid","reference_allele","af","n"),
                to=c("SNP","other_allele","eaf","samplesize"))
ldlc$Phenotype<-c("LDL cholesterol")
ldlc$units<-c("units")
ldlc$id<-c("Graham 2021")
ldlc$se<-ldlc$sd/sqrt(ldlc$samplesize)
ldlc<-ldlc[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
              "Phenotype","units","samplesize","id")]
write.table(ldlc,"./gwas/ldlc_graham2021/ldlc_graham2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### HDL CHOLESTEROL ###

hdlc<-read.csv2("./gwas/hdlc_graham2021/hdlc_graham2021.csv",header=TRUE,sep=";",dec=".")
hdlc<-rename.vars(hdlc,
                from=c("ï..rsid","ref","alt","af","n"),
                to=c("SNP","other_allele","effect_allele","eaf","samplesize"))
hdlc$Phenotype<-c("HDL cholesterol")
hdlc$units<-c("units")
hdlc$id<-c("Graham 2021")
hdlc$se<-hdlc$sd/sqrt(hdlc$samplesize)
hdlc<-hdlc[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
              "Phenotype","units","samplesize","id")]
write.table(hdlc,"./gwas/hdlc_graham2021/hdlc_graham2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### TRIGLYCERIDES ###

tg<-read.csv2("./gwas/tg_graham2021/tg_graham2021.csv",header=TRUE,sep=";",dec=".")
tg<-rename.vars(tg,
                  from=c("rsid","reference_allele","af","n"),
                  to=c("SNP","other_allele","eaf","samplesize"))
tg$Phenotype<-c("Triglycerides")
tg$units<-c("units")
tg$id<-c("Graham 2021")
tg$se<-tg$sd/sqrt(tg$samplesize)
tg<-tg[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
          "Phenotype","units","samplesize","id")]
write.table(tg,"./gwas/tg_graham2021/tg_graham2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### SYSTOLIC BLOOD PRESSURE ###

sbp<-read.csv2("./gwas/sbp_evangelou2018/sbp_evangelou2018.csv",header=TRUE,sep=";",dec=".")
sbp<-rename.vars(sbp,
                 from=c("rsid"),
                 to=c("SNP"))
sbp$Phenotype<-c("Systolic blood pressure")
sbp$units<-c("units")
sbp$id<-c("Evangelou 2018")
sbp$samplesize<-757601
sbp<-sbp[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
            "Phenotype","units","samplesize","id")]
write.table(sbp,"./gwas/sbp_evangelou2018/sbp_evangelou2018.txt",sep="\t",row.names=FALSE, quote=FALSE)

sbp<-read.csv2("./gwas/sbp_evangelou2018/sbp_georgiopoulos2021.csv",header=TRUE,sep=";",dec=".")
sbp<-rename.vars(sbp,
                 from=c("rsid"),
                 to=c("SNP"))
sbp$Phenotype<-c("Systolic blood pressure")
sbp$units<-c("units")
sbp$id<-c("Georgiopoulos 2021")
sbp$samplesize<-757601
sbp<-sbp[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
            "Phenotype","units","samplesize","id")]
write.table(sbp,"./gwas/sbp_evangelou2018/sbp_georgiopoulos2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### DIASTOLIC BLOOD PRESSURE ###

dbp<-read.csv2("./gwas/dbp_evangelou2018/dbp_evangelou2018.csv",header=TRUE,sep=";",dec=".")
dbp<-rename.vars(dbp,
                 from=c("rsid"),
                 to=c("SNP"))
dbp$Phenotype<-c("Diastolic blood pressure")
dbp$units<-c("units")
dbp$id<-c("Evangelou 2018")
dbp$samplesize<-757601
dbp<-dbp[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
            "Phenotype","units","samplesize","id")]
write.table(dbp,"./gwas/dbp_evangelou2018/dbp_evangelou2018.txt",sep="\t",row.names=FALSE, quote=FALSE)

dbp<-read.csv2("./gwas/dbp_evangelou2018/dbp_georgiopoulos2021.csv",header=TRUE,sep=";",dec=".")
dbp<-rename.vars(dbp,
                 from=c("rsid"),
                 to=c("SNP"))
dbp$Phenotype<-c("Diastolic blood pressure")
dbp$units<-c("units")
dbp$id<-c("Georgiopoulos 2021")
dbp$samplesize<-757601
dbp<-dbp[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
            "Phenotype","units","samplesize","id")]
write.table(dbp,"./gwas/dbp_evangelou2018/dbp_georgiopoulos2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### FASTING GLUCOSE ###

fg<-read.csv2("./gwas/fg_chen2021/fg_chen2021.csv",header=TRUE,sep=";",dec=".")
fg<-rename.vars(fg,
                from=c("rsid","pval_all"),
                to=c("SNP","pval"))
fg$Phenotype<-c("Fasting glucose")
fg$units<-c("units")
fg$id<-c("Chen 2021")
fg$samplesize<-200602
fg<-fg[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
          "Phenotype","units","samplesize","id")]
write.table(fg,"./gwas/fg_chen2021/fg_chen2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### FASTING GLUCOSE ###

hb1ac<-read.csv2("./gwas/hb1ac_chen2021/hb1ac_chen2021.csv",header=TRUE,sep=";",dec=".")
hb1ac<-rename.vars(hb1ac,
                from=c("rsid","pval_all"),
                to=c("SNP","pval"))
hb1ac$Phenotype<-c("Glycated hemoglobin")
hb1ac$units<-c("units")
hb1ac$id<-c("Chen 2021")
hb1ac$samplesize<-146806
hb1ac<-hb1ac[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
          "Phenotype","units","samplesize","id")]
write.table(hb1ac,"./gwas/hb1ac_chen2021/hb1ac_chen2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### FASTING INSULIN ###

fins<-read.csv2("./gwas/fins_chen2021/fins_chen2021.csv",header=TRUE,sep=";",dec=".")
fins<-rename.vars(fins,
                  from=c("rsid","pval_all"),
                  to=c("SNP","pval"))
fins$Phenotype<-c("Fasting insulin")
fins$units<-c("units")
fins$id<-c("Chen 2021")
fins$samplesize<-151013
fins<-fins[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
              "Phenotype","units","samplesize","id")]
write.table(fins,"./gwas/fins_chen2021/fins_chen2021.txt",sep="\t",row.names=FALSE, quote=FALSE)


### BODY MASS INDEX ###

bmi<-read.csv2("./gwas/bmi_yengo2018/bmi_yengo2018.csv",header=TRUE,sep=";",dec=".")
bmi<-rename.vars(bmi,
                 from=c("rsid","effect_allele_freq"),
                 to=c("SNP","eaf"))
bmi$Phenotype<-c("Body mass index")
bmi$units<-c("units")
bmi$id<-c("Yengo 2018")
bmi$samplesize<-681275
bmi<-bmi[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
            "Phenotype","units","samplesize","id")]
write.table(bmi,"./gwas/bmi_yengo2018/bmi_yengo2018.txt",sep="\t",row.names=FALSE, quote=FALSE)


### SMKINIT ###

smkinit<-read.csv2("./gwas/smkinit_liu2019/smkinit_liu2019.csv",header=TRUE,sep=";",dec=".")
smkinit<-rename.vars(smkinit,
                 from=c("rsid","effect_allele_freq"),
                 to=c("SNP","eaf"))
smkinit$Phenotype<-c("Smoking initiation")
smkinit$units<-c("units")
smkinit$id<-c("Liu 2019")
smkinit$samplesize<-1232091
smkinit<-smkinit[,c("SNP","chr","pos","other_allele","effect_allele","eaf","beta","se","pval",
            "Phenotype","units","samplesize","id")]
write.table(smkinit,"./gwas/smkinit_liu2019/smkinit_liu2019.txt",sep="\t",row.names=FALSE, quote=FALSE)



### GETTING READY THE SUMMARY DATABASES FOR THE STEIGER FILTERING ###

### LDL CHOLESTEROL ###

ldlc<-fread("./gwas/ldlc_graham2021/summary_ldlc_graham2021.results",header=TRUE,sep="\t",sep2="\t")
names(ldlc)<-tolower(names(ldlc))
ldlc<-rename.vars(ldlc,
                from=c("rsid","chrom","pos_b37","ref","alt","pooled_alt_af","n","effect_size","pvalue"),
                to=c("SNP","chr","pos","other_allele","effect_allele","eaf","samplesize","beta","pval"))
ldlc$Units<-c("Units")
ldlc$Phenotype<-c("Total cholesterol")
ldlc<-ldlc[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(ldlc,"./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
ldlc<-NULL


### HDL CHOLESTEROL ###

hdlc<-fread("./gwas/hdlc_graham2021/summary_hdlc_graham2021.results",header=TRUE,sep="\t",sep2="\t")
names(hdlc)<-tolower(names(hdlc))
hdlc<-rename.vars(hdlc,
                from=c("rsid","chrom","pos_b37","ref","alt","pooled_alt_af","n","effect_size","pvalue"),
                to=c("SNP","chr","pos","other_allele","effect_allele","eaf","samplesize","beta","pval"))
hdlc$Units<-c("Units")
hdlc$Phenotype<-c("Total cholesterol")
hdlc<-hdlc[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(hdlc,"./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
hdlc<-NULL


### TRIGLYCERIDES ###

tg<-fread("./gwas/tg_graham2021/summary_tg_graham2021.results",header=TRUE,sep="\t",sep2="\t")
names(tg)<-tolower(names(tg))
tg<-rename.vars(tg,
                from=c("rsid","chrom","pos_b37","ref","alt","pooled_alt_af","n","effect_size","pvalue"),
                to=c("SNP","chr","pos","other_allele","effect_allele","eaf","samplesize","beta","pval"))
tg$Units<-c("Units")
tg$Phenotype<-c("Total cholesterol")
tg<-tg[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(tg,"./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
tg<-NULL


### SYSTOLIC BLOOD PRESSURE ###

moba_all<-fread("N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_translate.txt", header=TRUE, sep="\t", sep2="\t")
names(moba_all)<-tolower(names(moba_all))
hunt_all<-fread("N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/hunt_translate.txt", header=TRUE, sep="\t", sep2="\t")
names(hunt_all)<-tolower(names(hunt_all))
alspac_all<-fread("N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/alspac_translate.txt", header=TRUE, sep="\t", sep2="\t")
names(alspac_all)<-tolower(names(alspac_all))
graham_all<-fread("N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/graham2021_translate.txt", header=TRUE, sep="\t", sep2="\t")
names(graham_all)<-tolower(names(graham_all))

sbp<-fread("./gwas/sbp_evangelou2018/summary_sbp_evangelou2018.txt",header=TRUE,sep=" ",sep2=" ")
names(sbp)<-tolower(names(sbp))
sbp$chr<-as.numeric(sub("\\:.*", "", sbp$markername))
sbp$pos<-as.numeric(gsub(".*:(.*)\\:.*", "\\1", sbp$markername))
sbp$markername<-paste(paste(sbp$chr,sbp$pos,sep=":"),paste(toupper(sbp$allele1),toupper(sbp$allele2),sep="/"),sep="_")

sbp<-merge(sbp,moba_all,by="markername",all.x=TRUE)
sbp$markername<-with(sbp,ifelse(!is.na(rsid),rsid,markername))
sbp$rsid<-NULL
sbp<-merge(sbp,hunt_all,by="markername",all.x=TRUE)
sbp$markername<-with(sbp,ifelse(!is.na(rsid),rsid,markername))
sbp$rsid<-NULL
sbp<-merge(sbp,alspac_all,by="markername",all.x=TRUE)
sbp$markername<-with(sbp,ifelse(!is.na(rsid),rsid,markername))
sbp$rsid<-NULL
sbp<-merge(sbp,graham_all,by="markername",all.x=TRUE)
sbp$markername<-with(sbp,ifelse(!is.na(rsid),rsid,markername))
sbp$rsid<-NULL
sbp<-sbp[startsWith(sbp$markername,"rs"),]

sbp<-rename.vars(sbp,
                 from=c("markername","allele2","allele1","freq1","n_effective","effect","stderr","p"),
                 to=c("SNP","other_allele","effect_allele","eaf","samplesize","beta","se","pval"))
sbp$Units<-c("Units")
sbp$Phenotype<-c("Systolic blood pressure")
sbp<-sbp[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(sbp,"./gwas/sbp_evangelou2018/summaryclean_sbp_evangelou2018.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
sbp<-NULL


### DIASTOLIC BLOOD PRESSURE ###

dbp<-fread("./gwas/dbp_evangelou2018/summary_dbp_evangelou2018.txt",header=TRUE,sep=" ",sep2=" ")
names(dbp)<-tolower(names(dbp))
dbp$chr<-as.numeric(sub("\\:.*", "", dbp$markername))
dbp$pos<-as.numeric(gsub(".*:(.*)\\:.*", "\\1", dbp$markername))
dbp$markername<-paste(paste(dbp$chr,dbp$pos,sep=":"),paste(toupper(dbp$allele1),toupper(dbp$allele2),sep="/"),sep="_")

dbp<-merge(dbp,moba_all,by="markername",all.x=TRUE)
dbp$markername<-with(dbp,ifelse(!is.na(rsid),rsid,markername))
dbp$rsid<-NULL
dbp<-merge(dbp,hunt_all,by="markername",all.x=TRUE)
dbp$markername<-with(dbp,ifelse(!is.na(rsid),rsid,markername))
dbp$rsid<-NULL
dbp<-merge(dbp,alspac_all,by="markername",all.x=TRUE)
dbp$markername<-with(dbp,ifelse(!is.na(rsid),rsid,markername))
dbp$rsid<-NULL
dbp<-merge(dbp,graham_all,by="markername",all.x=TRUE)
dbp$markername<-with(dbp,ifelse(!is.na(rsid),rsid,markername))
dbp$rsid<-NULL
dbp<-dbp[startsWith(dbp$markername,"rs"),]

dbp<-rename.vars(dbp,
                 from=c("markername","allele2","allele1","freq1","n_effective","effect","stderr","p"),
                 to=c("SNP","other_allele","effect_allele","eaf","samplesize","beta","se","pval"))
dbp$Units<-c("Units")
dbp$Phenotype<-c("Systolic blood pressure")
dbp<-dbp[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(dbp,"./gwas/dbp_evangelou2018/summaryclean_dbp_evangelou2018.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
dbp<-NULL


### FASTING GLUCOSE ###

fg<-fread("./gwas/fg_chen2021/summary_fg_chen2021",header=TRUE,sep="\t",sep2="\t")
names(fg)<-tolower(names(fg))
fg$markername<-paste(paste(fg$chromosome,fg$base_pair_location,sep=":"),paste(toupper(fg$other_allele),toupper(fg$effect_allele),sep="/"),sep="_")

fg<-merge(fg,moba_all,by="markername",all.x=TRUE)
fg$markername<-with(fg,ifelse(!is.na(rsid),rsid,markername))
fg$rsid<-NULL
fg<-merge(fg,hunt_all,by="markername",all.x=TRUE)
fg$markername<-with(fg,ifelse(!is.na(rsid),rsid,markername))
fg$rsid<-NULL
fg<-merge(fg,alspac_all,by="markername",all.x=TRUE)
fg$markername<-with(fg,ifelse(!is.na(rsid),rsid,markername))
fg$rsid<-NULL
fg<-merge(fg,graham_all,by="markername",all.x=TRUE)
fg$markername<-with(fg,ifelse(!is.na(rsid),rsid,markername))
fg$rsid<-NULL
fg<-fg[startsWith(fg$markername,"rs"),]

fg<-rename.vars(fg,
                from=c("markername","chromosome","base_pair_location","effect_allele_frequency","standard_error","p_value","sample_size"),
                to=c("SNP","chr","pos","eaf","se","pval","samplesize"))
fg$Units<-c("Units")
fg$Phenotype<-c("Fasting glucose")
fg<-fg[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(fg,"./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
fg<-NULL


### FASTING INSULIN ###

fins<-fread("./gwas/fins_chen2021/summary_fins_chen2021",header=TRUE,sep="\t",sep2="\t")
names(fins)<-tolower(names(fins))
fins$markername<-paste(paste(fins$chromosome,fins$base_pair_location,sep=":"),paste(toupper(fins$other_allele),toupper(fins$effect_allele),sep="/"),sep="_")

fins<-merge(fins,moba_all,by="markername",all.x=TRUE)
fins$markername<-with(fins,ifelse(!is.na(rsid),rsid,markername))
fins$rsid<-NULL
fins<-merge(fins,hunt_all,by="markername",all.x=TRUE)
fins$markername<-with(fins,ifelse(!is.na(rsid),rsid,markername))
fins$rsid<-NULL
fins<-merge(fins,alspac_all,by="markername",all.x=TRUE)
fins$markername<-with(fins,ifelse(!is.na(rsid),rsid,markername))
fins$rsid<-NULL
fins<-merge(fins,graham_all,by="markername",all.x=TRUE)
fins$markername<-with(fins,ifelse(!is.na(rsid),rsid,markername))
fins$rsid<-NULL
fins<-fins[startsWith(fins$markername,"rs"),]

fins<-rename.vars(fins,
                from=c("markername","chromosome","base_pair_location","effect_allele_frequency","standard_error","p_value","sample_size"),
                to=c("SNP","chr","pos","eaf","se","pval","samplesize"))
fins$Units<-c("Units")
fins$Phenotype<-c("Fasting insulin")
fins<-fins[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(fins,"./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
fins<-NULL


### GLYCATED HEMOGLOBIN ###

hb1ac<-fread("./gwas/hb1ac_chen2021/summary_hb1ac_chen2021",header=TRUE,sep="\t",sep2="\t")
names(hb1ac)<-tolower(names(hb1ac))
hb1ac$markername<-paste(paste(hb1ac$chromosome,hb1ac$base_pair_location,sep=":"),paste(toupper(hb1ac$other_allele),toupper(hb1ac$effect_allele),sep="/"),sep="_")

hb1ac<-merge(hb1ac,moba_all,by="markername",all.x=TRUE)
hb1ac$markername<-with(hb1ac,ifelse(!is.na(rsid),rsid,markername))
hb1ac$rsid<-NULL
hb1ac<-merge(hb1ac,hunt_all,by="markername",all.x=TRUE)
hb1ac$markername<-with(hb1ac,ifelse(!is.na(rsid),rsid,markername))
hb1ac$rsid<-NULL
hb1ac<-merge(hb1ac,alspac_all,by="markername",all.x=TRUE)
hb1ac$markername<-with(hb1ac,ifelse(!is.na(rsid),rsid,markername))
hb1ac$rsid<-NULL
hb1ac<-merge(hb1ac,graham_all,by="markername",all.x=TRUE)
hb1ac$markername<-with(hb1ac,ifelse(!is.na(rsid),rsid,markername))
hb1ac$rsid<-NULL
hb1ac<-hb1ac[startsWith(hb1ac$markername,"rs"),]

hb1ac<-rename.vars(hb1ac,
                from=c("markername","chromosome","base_pair_location","effect_allele_frequency","standard_error","p_value","sample_size"),
                to=c("SNP","chr","pos","eaf","se","pval","samplesize"))
hb1ac$Units<-c("Units")
hb1ac$Phenotype<-c("Glycated hemoglobin")
hb1ac<-hb1ac[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(hb1ac,"./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
hb1ac<-NULL
moba_all<-NULL
hunt_all<-NULL
alspac_all<-NULL
graham_all<-NULL


### BMI ###

bmi<-fread("N:/data/durable/Projects/Hernaez_MR_CVRF/2smr/bmi.txt",header=TRUE,sep="\t",sep2="\t")
names(bmi)<-tolower(names(bmi))
bmi<-rename.vars(bmi,
                 from=c("snp","tested_allele","p","n","freq_tested_allele_in_hrs"),
                 to=c("SNP","effect_allele","pval","samplesize","eaf"))
bmi$Units<-c("Units")
bmi$Phenotype<-c("Body mass index")
bmi<-bmi[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(bmi,"./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
bmi<-NULL


### SMKINIT ###

smkinit<-fread("./gwas/smkinit_liu2019/summary_smkinit_liu2019.txt",header=TRUE,sep="\t",sep2="\t")
names(smkinit)<-tolower(names(smkinit))
smkinit<-rename.vars(smkinit,
                     from=c("rsid","chrom","ref","alt","af","pvalue","effective_n"),
                     to=c("SNP","chr","other_allele","effect_allele","eaf","pval","samplesize"))
smkinit$Units<-c("Units")
smkinit$Phenotype<-c("Smoking initiation")
smkinit<-smkinit[,c("SNP","chr","pos","effect_allele","other_allele","eaf","beta","se","pval","samplesize","Units","Phenotype")]
fwrite(smkinit,"./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)
smkinit<-NULL



### STEIGER FILTERING - LDL CHOLESTEROL ###
###########################################

ldlc<-read_exposure_data(
  filename = "./gwas/ldlc_graham2021/ldlc_graham2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

ldlc_orig<-ldlc

hdlc <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, hdlc, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

tg <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, tg, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

sbp <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/sbp_evangelou2018/summaryclean_sbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, sbp, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

dbp <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/dbp_evangelou2018/summaryclean_dbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, dbp, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

fg <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, fg, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

fins <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, fins, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

hb1ac <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, hb1ac, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

bmi <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, bmi, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

smkinit <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

ldlc<-harmonise_data(ldlc, smkinit, action = 1)
ldlc<-steiger_filtering(ldlc)
ldlc<-subset2(ldlc,"ldlc$steiger_dir==TRUE")
ldlc<-ldlc_orig[ldlc_orig$SNP%in%ldlc$SNP,]
dim(ldlc)[1]

fwrite(as.data.frame(ldlc$SNP),"./gwas/ldlc_graham2021/rsid_steiger_ldlc_graham2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)


### STEIGER FILTERING - HDL CHOLESTEROL ###
###########################################

hdlc<-read_exposure_data(
  filename = "./gwas/hdlc_graham2021/hdlc_graham2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

hdlc_orig<-hdlc

ldlc <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, ldlc, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

tg <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, tg, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

sbp <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/sbp_evangelou2018/summaryclean_sbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, sbp, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

dbp <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/dbp_evangelou2018/summaryclean_dbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, dbp, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

fg <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, fg, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

fins <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, fins, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

hb1ac <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, hb1ac, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

bmi <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, bmi, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

smkinit <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hdlc<-harmonise_data(hdlc, smkinit, action = 1)
hdlc<-steiger_filtering(hdlc)
hdlc<-subset2(hdlc,"hdlc$steiger_dir==TRUE")
hdlc<-hdlc_orig[hdlc_orig$SNP%in%hdlc$SNP,]
dim(hdlc)[1]

fwrite(as.data.frame(hdlc$SNP),"./gwas/hdlc_graham2021/rsid_steiger_hdlc_graham2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=FALSE)


### STEIGER FILTERING - TRIGLYCERIDES ###
#########################################

tg<-read_exposure_data(
  filename = "./gwas/tg_graham2021/tg_graham2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

tg_orig<-tg

ldlc <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, ldlc, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

hdlc <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, hdlc, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

sbp <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/sbp_evangelou2018/summaryclean_sbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, sbp, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

dbp <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/dbp_evangelou2018/summaryclean_dbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, dbp, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

fg <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, fg, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

fins <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, fins, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

hb1ac <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, hb1ac, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

bmi <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, bmi, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

smkinit <- read_outcome_data(
  snps = tg$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

tg<-harmonise_data(tg, smkinit, action = 1)
tg<-steiger_filtering(tg)
tg<-subset2(tg,"tg$steiger_dir==TRUE")
tg<-tg_orig[tg_orig$SNP%in%tg$SNP,]
dim(tg)[1]

fwrite(as.data.frame(tg$SNP),"./gwas/tg_graham2021/rsid_steiger_tg_graham2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=FALSE)


### STEIGER FILTERING - SYSTOLIC BLOOD PRESSURE ###
###################################################

sbp<-read_exposure_data(
  filename = "./gwas/sbp_evangelou2018/sbp_georgiopoulos2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

sbp_orig<-sbp

ldlc <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, ldlc, action = 1)
dim(sbp)[1]
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

hdlc <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, hdlc, action = 1)
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

tg <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, tg, action = 1)
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

fg <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, fg, action = 1)
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

fins <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, fins, action = 1)
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

hb1ac <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, hb1ac, action = 1)
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

bmi <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, bmi, action = 1)
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

smkinit <- read_outcome_data(
  snps = sbp$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

sbp<-harmonise_data(sbp, smkinit, action = 1)
sbp<-steiger_filtering(sbp)
sbp<-subset2(sbp,"sbp$steiger_dir==TRUE")
sbp<-sbp_orig[sbp_orig$SNP%in%sbp$SNP,]
dim(sbp)[1]

fwrite(as.data.frame(sbp$SNP),"./gwas/sbp_evangelou2018/rsid_steiger_sbp_georgiopoulos2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)


### STEIGER FILTERING - DIASTOLIC BLOOD PRESSURE ###
####################################################

dbp<-read_exposure_data(
  filename = "./gwas/dbp_evangelou2018/dbp_georgiopoulos2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

dbp_orig<-dbp

ldlc <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, ldlc, action = 1)
dim(dbp)[1]
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

hdlc <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, hdlc, action = 1)
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

tg <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, tg, action = 1)
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

fg <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, fg, action = 1)
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

fins <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, fins, action = 1)
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

hb1ac <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, hb1ac, action = 1)
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

bmi <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, bmi, action = 1)
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

smkinit <- read_outcome_data(
  snps = dbp$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dbp<-harmonise_data(dbp, smkinit, action = 1)
dbp<-steiger_filtering(dbp)
dbp<-subset2(dbp,"dbp$steiger_dir==TRUE")
dbp<-dbp_orig[dbp_orig$SNP%in%dbp$SNP,]
dim(dbp)[1]

fwrite(as.data.frame(dbp$SNP),"./gwas/dbp_evangelou2018/rsid_steiger_dbp_georgiopoulos2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=FALSE)


### STEIGER FILTERING - FASTING GLUCOSE ###
###########################################

fg<-read_exposure_data(
  filename = "./gwas/fg_chen2021/fg_chen2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

fg_orig<-fg

ldlc <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, ldlc, action = 1)
dim(fg)[1]
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

hdlc <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, hdlc, action = 1)
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

tg <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, tg, action = 1)
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

sbp <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/sbp_evangelou2018/summaryclean_sbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, sbp, action = 1)
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

dbp <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/dbp_evangelou2018/summaryclean_dbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, fins, action = 1)
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

hb1ac <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, hb1ac, action = 1)
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

bmi <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, bmi, action = 1)
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

smkinit <- read_outcome_data(
  snps = fg$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fg<-harmonise_data(fg, smkinit, action = 1)
fg<-steiger_filtering(fg)
fg<-subset2(fg,"fg$steiger_dir==TRUE")
fg<-fg_orig[fg_orig$SNP%in%fg$SNP,]
dim(fg)[1]

fwrite(as.data.frame(fg$SNP),"./gwas/fg_chen2021/rsid_steiger_fg_chen2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=FALSE)


### STEIGER FILTERING - FASTING INSULIN ###
###########################################

fins<-read_exposure_data(
  filename = "./gwas/fins_chen2021/fins_chen2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

fins_orig<-fins

ldlc <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, ldlc, action = 1)
dim(fins)[1]
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

hdlc <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, hdlc, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

tg <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, tg, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

sbp <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/sbp_evangelou2018/summaryclean_sbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, sbp, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

dbp <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/dbp_evangelou2018/summaryclean_dbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, dbp, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

fg <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, fg, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

hb1ac <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/hb1ac_chen2021/summaryclean_hb1ac_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, hb1ac, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

bmi <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, bmi, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

smkinit <- read_outcome_data(
  snps = fins$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

fins<-harmonise_data(fins, smkinit, action = 1)
fins<-steiger_filtering(fins)
fins<-subset2(fins,"fins$steiger_dir==TRUE")
fins<-fins_orig[fins_orig$SNP%in%fins$SNP,]
dim(fins)[1]

fwrite(as.data.frame(fins$SNP),"./gwas/fins_chen2021/rsid_steiger_fins_chen2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=TRUE)


### STEIGER FILTERING - GLYCATED HEMOGLOBIN ###
###############################################

hb1ac<-read_exposure_data(
  filename = "./gwas/hb1ac_chen2021/hb1ac_chen2021.txt",
  clump = FALSE,
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  ncase_col = "ncase",
  ncontrol_col = "ncontrol",
  samplesize_col = "samplesize",
  gene_col = "gene",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE
)

hb1ac_orig<-hb1ac

ldlc <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/ldlc_graham2021/summaryclean_ldlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, ldlc, action = 1)
dim(hb1ac)[1]
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

hdlc <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/hdlc_graham2021/summaryclean_hdlc_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, hdlc, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

tg <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/tg_graham2021/summaryclean_tg_graham2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, tg, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

sbp <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/sbp_evangelou2018/summaryclean_sbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, sbp, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

dbp <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/dbp_evangelou2018/summaryclean_dbp_evangelou2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, dbp, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

fg <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/fg_chen2021/summaryclean_fg_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, fg, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

fins <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/fins_chen2021/summaryclean_fins_chen2021.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, fins, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

bmi <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/bmi_yengo2018/summaryclean_bmi_yengo2018.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, bmi, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

smkinit <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "./gwas/smkinit_liu2019/summaryclean_smkinit_liu2019.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

hb1ac<-harmonise_data(hb1ac, smkinit, action = 1)
hb1ac<-steiger_filtering(hb1ac)
hb1ac<-subset2(hb1ac,"hb1ac$steiger_dir==TRUE")
hb1ac<-hb1ac_orig[hb1ac_orig$SNP%in%hb1ac$SNP,]
dim(hb1ac)[1]

fwrite(as.data.frame(hb1ac$SNP),"./gwas/hb1ac_chen2021/rsid_steiger_hb1ac_chen2021.txt",
       append=FALSE, sep="\t", sep2=c("\t","|","\t"), row.names=FALSE, col.names=FALSE)




##########################################
### TWO-SAMPLE MENDELIAN RANDOMIZATION ###
##########################################

setwd("N:/data/durable/Projects/Hernaez_MR_CVRF/")

### WOMEN ###
#############

### EXPOSURE: LDL CHOLESTEROL ###

ldlc<-read_exposure_data(
  filename = "./gwas/ldlc_graham2021/ldlc_graham2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(ldlc, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_ldlc_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/ldlc_subf_mom.RData")

ldlc_steiger<-fread("./gwas/ldlc_graham2021/rsid_steiger_ldlc_graham2021.txt",header=FALSE,sep="\t",sep2="\t")
names(ldlc_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%ldlc_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_ldlc_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/ldlc_steiger_subf_mom.RData")


### EXPOSURE: HDL CHOLESTEROL ###

hdlc<-read_exposure_data(
  filename = "./gwas/hdlc_graham2021/hdlc_graham2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(hdlc, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hdlc_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hdlc_subf_mom.RData")

hdlc_steiger<-fread("./gwas/hdlc_graham2021/rsid_steiger_hdlc_graham2021.txt",header=FALSE,sep="\t",sep2="\t")
names(hdlc_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%hdlc_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hdlc_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hdlc_steiger_subf_mom.RData")


### EXPOSURE: TRIGLYCERIDES ###

tg<-read_exposure_data(
  filename = "./gwas/tg_graham2021/tg_graham2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = tg$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(tg, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_tg_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/tg_subf_mom.RData")

tg_steiger<-fread("./gwas/tg_graham2021/rsid_steiger_tg_graham2021.txt",header=FALSE,sep="\t",sep2="\t")
names(tg_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%tg_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_tg_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/tg_steiger_subf_mom.RData")


### EXPOSURE: FASTING GLUCOSE ###

fg<-read_exposure_data(
  filename = "./gwas/fg_chen2021/fg_chen2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = fg$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(fg, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fg_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fg_subf_mom.RData")

fg_steiger<-fread("./gwas/fg_chen2021/rsid_steiger_fg_chen2021.txt",header=FALSE,sep="\t",sep2="\t")
names(fg_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%fg_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fg_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fg_steiger_subf_mom.RData")


### EXPOSURE: FASTING GLUCOSE ###

fins<-read_exposure_data(
  filename = "./gwas/fins_chen2021/fins_chen2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = fins$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(fins, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fins_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fins_subf_mom.RData")

fins_steiger<-fread("./gwas/fins_chen2021/rsid_steiger_fins_chen2021.txt",header=FALSE,sep="\t",sep2="\t")
names(fins_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%fins_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fins_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fins_steiger_subf_mom.RData")


### EXPOSURE: GLYCATED HEMOGLOBIN ###

hb1ac<-read_exposure_data(
  filename = "./gwas/hb1ac_chen2021/hb1ac_chen2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(hb1ac, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hb1ac_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hb1ac_subf_mom.RData")

hb1ac_steiger<-fread("./gwas/hb1ac_chen2021/rsid_steiger_hb1ac_chen2021.txt",header=FALSE,sep="\t",sep2="\t")
names(hb1ac_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%hb1ac_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hb1ac_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hb1ac_steiger_subf_mom.RData")


### EXPOSURE: SYSTOLIC BLOOD PRESSURE ###

sbp<-read_exposure_data(
  filename = "./gwas/sbp_evangelou2018/sbp_evangelou2018.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = sbp$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(sbp, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_sbp_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/sbp_subf_mom.RData")

sbp_steiger<-fread("./gwas/sbp_evangelou2018/rsid_steiger_sbp_georgiopoulos2021.txt",header=FALSE,sep="\t",sep2="\t")
names(sbp_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%sbp_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_sbp_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/sbp_steiger_subf_mom.RData")


### EXPOSURE: DIASTOLIC BLOOD PRESSURE ###

dbp<-read_exposure_data(
  filename = "./gwas/dbp_evangelou2018/dbp_evangelou2018.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = dbp$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_mom_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(dbp, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_dbp_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/dbp_subf_mom.RData")

dbp_steiger<-fread("./gwas/dbp_evangelou2018/rsid_steiger_dbp_georgiopoulos2021.txt",header=FALSE,sep="\t",sep2="\t")
names(dbp_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%dbp_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_dbp_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/dbp_steiger_subf_mom.RData")


### MEN ###
###########

### EXPOSURE: LDL CHOLESTEROL ###

ldlc<-read_exposure_data(
  filename = "./gwas/ldlc_graham2021/ldlc_graham2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = ldlc$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(ldlc, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_ldlc_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/ldlc_subf_dad.RData")

ldlc_steiger<-fread("./gwas/ldlc_graham2021/rsid_steiger_ldlc_graham2021.txt",header=FALSE,sep="\t",sep2="\t")
names(ldlc_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%ldlc_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_ldlc_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/ldlc_steiger_subf_dad.RData")


### EXPOSURE: HDL CHOLESTEROL ###

hdlc<-read_exposure_data(
  filename = "./gwas/hdlc_graham2021/hdlc_graham2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = hdlc$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(hdlc, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hdlc_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hdlc_subf_dad.RData")

hdlc_steiger<-fread("./gwas/hdlc_graham2021/rsid_steiger_hdlc_graham2021.txt",header=FALSE,sep="\t",sep2="\t")
names(hdlc_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%hdlc_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hdlc_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hdlc_steiger_subf_dad.RData")


### EXPOSURE: TRIGLYCERIDES ###

tg<-read_exposure_data(
  filename = "./gwas/tg_graham2021/tg_graham2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = tg$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(tg, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_tg_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/tg_subf_dad.RData")

tg_steiger<-fread("./gwas/tg_graham2021/rsid_steiger_tg_graham2021.txt",header=FALSE,sep="\t",sep2="\t")
names(tg_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%tg_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_tg_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/tg_steiger_subf_dad.RData")


### EXPOSURE: FASTING GLUCOSE ###

fg<-read_exposure_data(
  filename = "./gwas/fg_chen2021/fg_chen2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = fg$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(fg, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fg_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fg_subf_dad.RData")

fg_steiger<-fread("./gwas/fg_chen2021/rsid_steiger_fg_chen2021.txt",header=FALSE,sep="\t",sep2="\t")
names(fg_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%fg_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fg_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fg_steiger_subf_dad.RData")


### EXPOSURE: FASTING GLUCOSE ###

fins<-read_exposure_data(
  filename = "./gwas/fins_chen2021/fins_chen2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = fins$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(fins, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fins_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fins_subf_dad.RData")

fins_steiger<-fread("./gwas/fins_chen2021/rsid_steiger_fins_chen2021.txt",header=FALSE,sep="\t",sep2="\t")
names(fins_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%fins_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_fins_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/fins_steiger_subf_dad.RData")


### EXPOSURE: GLYCATED HEMOGLOBIN ###

hb1ac<-read_exposure_data(
  filename = "./gwas/hb1ac_chen2021/hb1ac_chen2021.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = hb1ac$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(hb1ac, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hb1ac_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hb1ac_subf_dad.RData")

hb1ac_steiger<-fread("./gwas/hb1ac_chen2021/rsid_steiger_hb1ac_chen2021.txt",header=FALSE,sep="\t",sep2="\t")
names(hb1ac_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%hb1ac_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_hb1ac_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/hb1ac_steiger_subf_dad.RData")


### EXPOSURE: SYSTOLIC BLOOD PRESSURE ###

sbp<-read_exposure_data(
  filename = "./gwas/sbp_evangelou2018/sbp_evangelou2018.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = sbp$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(sbp, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_sbp_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/sbp_subf_dad.RData")

sbp_steiger<-fread("./gwas/sbp_evangelou2018/rsid_steiger_sbp_georgiopoulos2021.txt",header=FALSE,sep="\t",sep2="\t")
names(sbp_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%sbp_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_sbp_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/sbp_steiger_subf_dad.RData")


### EXPOSURE: DIASTOLIC BLOOD PRESSURE ###

dbp<-read_exposure_data(
  filename = "./gwas/dbp_evangelou2018/dbp_evangelou2018.txt",
  clump = FALSE,
  sep="\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  units_col = "units",
  samplesize_col = "samplesize",
  id_col = "id",
  min_pval = 1e-200,
  log_pval = FALSE)

subf <- read_outcome_data(
  snps = dbp$SNP,
  filename = "N:/data/durable/Projects/Magnus_GWAS_Infertility/Metaanalysis/moba_dad_2smr.txt",
  sep = "\t",
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  pval_col = "pval",
  units_col = "Units",
  samplesize_col = "samplesize"
)

dat<-harmonise_data(dbp, subf, action = 2)
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_dbp_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/dbp_subf_dad.RData")

dbp_steiger<-fread("./gwas/dbp_evangelou2018/rsid_steiger_dbp_georgiopoulos2021.txt",header=FALSE,sep="\t",sep2="\t")
names(dbp_steiger)<-c("SNP")
dat<-dat[dat$SNP%in%dbp_steiger$SNP,]
xxx<-dat[,c("SNP","mr_keep")]
xxx$mr_keep<-with(xxx,ifelse(mr_keep=="TRUE",1,0))
xxx<-subset2(xxx,"xxx$mr_keep==1")
write.table(xxx$SNP,file="./Outputs/descriptive/snps_dbp_steiger_2smr.csv",sep=";",col.names=FALSE,row.names=FALSE)
save(dat,file="./2smr/dbp_steiger_subf_dad.RData")


##############################
### TWO-SAMPLE MR ANALYSES ###
##############################

setwd("N:/data/durable/Projects/Hernaez_MR_CVRF/Outputs/results/")

z<-qnorm(1-0.05/2)
Isq = function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

# F>10: good IVW performance
# Isq close to 1: good MR-Egger performance

vars01<-c("ldlc_","ldlc_steiger_","ldlc_","ldlc_steiger_","hdlc_","hdlc_steiger_","hdlc_","hdlc_steiger_",
          "tg_","tg_steiger_","tg_","tg_steiger_","sbp_","sbp_steiger_","sbp_","sbp_steiger_","dbp_","dbp_steiger_","dbp_","dbp_steiger_",
          "fg_","fg_steiger_","fg_","fg_steiger_","fins_","fins_steiger_","fins_","fins_steiger_","hb1ac_","hb1ac_steiger_","hb1ac_","hb1ac_steiger_")
vars02<-c("mom","mom","dad","dad","mom","mom","dad","dad","mom","mom","dad","dad",
          "mom","mom","dad","dad","mom","mom","dad","dad","mom","mom","dad","dad","mom","mom","dad","dad","mom","mom","dad","dad")

tab<-NULL
for(i in 1:length(vars01))
  
{
  namedat<-paste("N:/data/durable/Projects/Hernaez_MR_CVRF/2smr/",vars01[i],"subf_",vars02[i],".RData",sep="")
  load(namedat)
  
  mr_results<-mr(dat,method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))
  mr_results$pval<-pval_guapa(mr_results$pval)
  mr_results$beta<-paste(guapa(mr_results$b)," (",guapa(mr_results$se),")",sep="")
  mr_results$or<-risk_se_ic_guapa(mr_results$b,mr_results$se)
  ivw_or<-round(exp(mr_results$b[1]),6)
  ivw_ic95lo<-round(exp(mr_results$b[1]-(z*mr_results$se[1])),6)
  ivw_ic95hi<-round(exp(mr_results$b[1]+(z*mr_results$se[1])),6)
  mr_raps<-mr.raps(b_exp=dat$beta.exposure,b_out=dat$beta.outcome,se_exp=dat$se.exposure,se_out=dat$se.outcome,diagnosis=FALSE)
  mr_raps_beta<-paste(guapa(mr_raps$beta.hat)," (",guapa(mr_raps$beta.se),")",sep="")
  mr_raps_or<-risk_se_ic_guapa(mr_raps$beta.hat,mr_raps$beta.se)
  mr_raps_pval<-pval_guapa(mr_raps$beta.p.value)
  q1<-paste(guapa(mr_heterogeneity(dat)$Q[2])," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[2]),")",sep="")
  q2<-paste(guapa(mr_heterogeneity(dat)$Q[1])," (P=",pval_guapa(mr_heterogeneity(dat)$Q_pval[1]),")",sep="")
  plei<-pval_guapa(mr_pleiotropy_test(dat)$pval)
  fstat<-guapa(mean((abs(dat$beta.exposure))^2/dat$se.exposure^2,na.rm=TRUE))
  unIsq<-guapa(Isq(abs(dat$beta.exposure),dat$se.exposure))
  
  tab<-rbind(tab,cbind(mr_results$nsnp[1],
                       mr_results$beta[1],mr_results$or[1],mr_results$pval[1],
                       mr_results$beta[2],mr_results$or[2],mr_results$pval[2],
                       mr_results$beta[3],mr_results$or[3],mr_results$pval[3],
                       mr_results$beta[4],mr_results$or[4],mr_results$pval[4],
                       mr_raps_beta,mr_raps_or,mr_raps_pval,plei,q1,q2,ivw_or,ivw_ic95lo,ivw_ic95hi,fstat,unIsq))
}

colnames(tab)<-c("SNPs","IVW_b","IVW_or","IVW_p","Egger_b","Egger_or","Egger_p",
                 "WMe_b","WMe_or","Wme_p","WMo_b","WMo_or","WMo_p",
                 "RAPS_b","RAPS_or","RAPS_p","Egger_pleio","Cochran_Q","Rucker_Q",
                 "ivw_or","ivw_ic95lo","ivw_ic95hi","F-stat","I2")
rownames(tab)<-paste(vars01,"subf_",vars02,sep="")
write.table(tab,file="./2smr.csv",sep=";",col.names=NA)

