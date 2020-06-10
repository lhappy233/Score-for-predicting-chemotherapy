rm(list = ls())
getwd()
setwd("/Users/linhaiping/Documents/Breast Cancer_HER2-HR+")

library(readstata13)
library(cmprsk)
Data <- read.dta13("./From Fu/age.progress.hrpostive.dta")

#####Data prepare#####
#年份
table(Data$year)
Data <- subset(Data,Data$year<=2014)

#年龄
table(Data$age)
table(Data$age4)
Data$age4 <- as.factor(Data$age4)

#种族
table(Data$race)
Data$race <- as.factor(Data$race)

#婚姻
table(Data$marry)
Data$marry <- as.factor(Data$marry)

#组织学
table(Data$histology)
Data$histology <- as.factor(Data$histology)

#部位
table(Data$location)
Data$location <- as.factor(Data$location)

#分级
table(Data$Grade)
Data$Grade <- factor(Data$Grade,levels = c("Well differentiated; Grade I",
                                           "Moderately differentiated; Grade II","Poorly differentiated; Grade III"))

#HR all HR+
table(Data$HR)
Data$HR <- as.factor(Data$HR)
Data$ER <- as.factor(Data$ER)
Data$PR <- as.factor(Data$PR)

#HER2 all HER2-
table(Data$her2)
Data$her2 <- as.factor(Data$her2)

#T
table(Data$stageT)
Data$stageT <- as.factor(Data$stageT)

#N
table(Data$N)
Data$N <- as.factor(Data$N)

#化疗
table(Data$chemotherapy)
Data$chemotherapy <- as.factor(Data$chemotherapy)

#放疗
table(Data$radiotherapy)
Data$radiotherapy <- as.factor(Data$radiotherapy)

#生存时间
table(Data$Survivalmonths)
Data <- subset(Data,Data$Survivalmonths>1)

#死因
table(Data$death)
table(Data$death3)
#Data$death <- as.factor(Data$death)
#Data$death3 <- as.factor(Data$death3)

#####PSM#####
library(MatchIt)
set.seed(2000)
attach(Data)
data_m <- matchit(chemotherapy~age4+race+marry+histology+Grade+stageT+N,data=Data,ratio=1,caliper=0.05)
detach(Data)
summary(data_m)
#save(data_m,file='./data_m.Rdata')
load("~/Documents/Breast Cancer_HER2-HR+/data_m.Rdata")

data_match <- match.data(data_m)
#plot(data_m,type="jitter")
plot(data_m,type="hist")
library(cobalt)
love.plot(data_m, stats = c("mean.diffs", "variance.ratios"), threshold = c(m = 0.1, v = 2), abs = TRUE, binary = "std", var.order = "unadjusted")

#####Data analysis#####
data_match$death <- as.numeric(data_match$death)
data_match$death3 <- as.numeric(data_match$death3)

attach(data_match)
mod.o_chemom <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_match)
summary(mod.o_chemom)#
detach(data_match)

#after matching
attach(data_match)
mod.m_chemo <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_match)
summary(mod.m_chemo)#
detach(data_match)

#multivariable
attach(data_match)
mod.m <- coxph(Surv(Survivalmonths,death)~age4+race+marry+histology+Grade+stageT+N+chemotherapy,data_match)
summary(mod.m)
detach(data_match)

#subgroup
data_matchage1 <- subset(data_match,data_match$age4=="1")
attach(data_matchage1)
mod.m_age1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchage1)
summary(mod.m_age1)
detach(data_matchage1)

data_matchage2 <- subset(data_match,data_match$age4=="2")
attach(data_matchage2)
mod.m_age2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchage2)
summary(mod.m_age2)
detach(data_matchage2)

data_matchage3 <- subset(data_match,data_match$age4=="3")
attach(data_matchage3)
mod.m_age3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchage3)
summary(mod.m_age3)
detach(data_matchage3)

data_matchage4 <- subset(data_match,data_match$age4=="4")
attach(data_matchage4)
mod.m_age4 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchage4)
summary(mod.m_age4)
detach(data_matchage4)

data_matchrace1 <- subset(data_match,data_match$race=="1")
attach(data_matchrace1)
mod.m_race1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchrace1)
summary(mod.m_race1)
detach(data_matchrace1)

data_matchrace2 <- subset(data_match,data_match$race=="2")
attach(data_matchrace2)
mod.m_race2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchrace2)
summary(mod.m_race2)
detach(data_matchrace2)

data_matchrace3 <- subset(data_match,data_match$race=="3")
attach(data_matchrace3)
mod.m_race3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchrace3)
summary(mod.m_race3)
detach(data_matchrace3)

data_matchmarri1 <- subset(data_match,data_match$marry=="1")
attach(data_matchmarri1)
mod.m_marri1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchmarri1)
summary(mod.m_marri1)
detach(data_matchmarri1)

data_matchmarri2 <- subset(data_match,data_match$marry=="2")
attach(data_matchmarri2)
mod.m_marri2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchmarri2)
summary(mod.m_marri2)
detach(data_matchmarri2)

data_matchmarri3 <- subset(data_match,data_match$marry=="3")
attach(data_matchmarri3)
mod.m_marri3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchmarri3)
summary(mod.m_marri3)
detach(data_matchmarri3)

data_matchhisto1 <- subset(data_match,data_match$histology=="1")
attach(data_matchhisto1)
mod.m_histo1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchhisto1)
summary(mod.m_histo1)
detach(data_matchhisto1)

data_matchhisto2 <- subset(data_match,data_match$histology=="2")
attach(data_matchhisto2)
mod.m_histo2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchhisto2)
summary(mod.m_histo2)
detach(data_matchhisto2)

data_matchhisto3 <- subset(data_match,data_match$histology=="3")
attach(data_matchhisto3)
mod.m_histo3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchhisto3)
summary(mod.m_histo3)
detach(data_matchhisto3)

data_matchgrade1 <- subset(data_match,data_match$Grade=="Well differentiated; Grade I")
attach(data_matchgrade1)
mod.m_grade1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchgrade1)
summary(mod.m_grade1)
detach(data_matchgrade1)

data_matchgrade2 <- subset(data_match,data_match$Grade=="Moderately differentiated; Grade II")
attach(data_matchgrade2)
mod.m_grade2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchgrade2)
summary(mod.m_grade2)
detach(data_matchgrade2)

data_matchgrade3 <- subset(data_match,data_match$Grade=="Poorly differentiated; Grade III")
attach(data_matchgrade3)
mod.m_grade3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchgrade3)
summary(mod.m_grade3)
detach(data_matchgrade3)

data_matchT1 <- subset(data_match,data_match$stageT=="1")
attach(data_matchT1)
mod.m_T1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchT1)
summary(mod.m_T1)
detach(data_matchT1)

data_matchT2 <- subset(data_match,data_match$stageT=="2")
attach(data_matchT2)
mod.m_T2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchT2)
summary(mod.m_T2)
detach(data_matchT2)

data_matchT3 <- subset(data_match,data_match$stageT=="3")
attach(data_matchT3)
mod.m_T3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchT3)
summary(mod.m_T3)
detach(data_matchT3)

data_matchT4 <- subset(data_match,data_match$stageT=="4")
attach(data_matchT4)
mod.m_T4 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchT4)
summary(mod.m_T4)
detach(data_matchT4)

data_matchN0 <- subset(data_match,data_match$N=="0")
attach(data_matchN0)
mod.m_N0 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchN0)
summary(mod.m_N0)
detach(data_matchN0)

data_matchN1 <- subset(data_match,data_match$N=="1")
attach(data_matchN1)
mod.m_N1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchN1)
summary(mod.m_N1)
detach(data_matchN1)

data_matchN2 <- subset(data_match,data_match$N=="2")
attach(data_matchN2)
mod.m_N2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchN2)
summary(mod.m_N2)
detach(data_matchN2)

data_matchN3 <- subset(data_match,data_match$N=="3")
attach(data_matchN3)
mod.m_N3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,data_matchN3)
summary(mod.m_N3)
detach(data_matchN3)

#####Forestplot#####
#forest_match <- read.csv("/Users/linhaiping/Documents/Breast Cancer_HER2-HR+/Forestplot_after_cox.csv",header = TRUE,na.strings = "NA")
library(forestplot)
head(forest_match,50)
labeltext_match <- as.matrix(forest_match[,c(1:3,7)])
attach(forest_match)
forestplot(labeltext_match,
           HR,
           Hrdown,
           Hrup,
           align = "l",
           graph.pos = 3,
           title="Hazard Ratio Plot",
           xlab="<---Favor chemotherapy---                      ---Favor non-chemotherapy--->",
           is.summary= c(T,T,F,T,
                         T,F,F,F,F,F,
                         T,F,F,F,F,
                         T,F,F,F,F,
                         T,F,F,F,F,
                         T,F,F,F,F,
                         T,F,F,F,F,F,
                         T,F,F,F,F,F),
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex=0.8),
                          title=gpar(cex=1.0)),
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           lwd.ci=1, boxsize=0.2,
           ci.vertices=TRUE, ci.vertices.height = 0.2,
           xlog=FALSE
)

#####P for trend#####
### z = coef / se（coef）
mod.m_age1 <- summary(mod.m_age1)
mod.m_age2 <- summary(mod.m_age2)
mod.m_age3 <- summary(mod.m_age3)
mod.m_age4 <- summary(mod.m_age4)

coef_age <- c(mod.m_age1$coefficients[1],mod.m_age2$coefficients[1],mod.m_age3$coefficients[1],mod.m_age4$coefficients[1])
se_age <- c(mod.m_age1$coefficients[3],mod.m_age2$coefficients[3],mod.m_age3$coefficients[3],mod.m_age4$coefficients[3])
n_age <- c(mod.m_age1$n,mod.m_age2$n,mod.m_age3$n,mod.m_age4$n)
age.gen <- metagen(coef_age, se_age, n.e = n_age, sm = "HR")
age.gen #< 0.0001
forest(age.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))


mod.m_race1 <- summary(mod.m_race1)
mod.m_race2 <- summary(mod.m_race2)
mod.m_race3 <- summary(mod.m_race3)

coef_race <- c(mod.m_race1$coefficients[1],mod.m_race2$coefficients[1],mod.m_race3$coefficients[1])
se_race <- c(mod.m_race1$coefficients[3],mod.m_race2$coefficients[3],mod.m_race3$coefficients[3])
n_race <- c(mod.m_race1$n,mod.m_race2$n,mod.m_race3$n)
race.gen <- metagen(coef_race, se_race, n.e = n_race, sm = "HR")
race.gen #0.0684
forest(race.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))


mod.m_marri1 <- summary(mod.m_marri1)
mod.m_marri2 <- summary(mod.m_marri2)
mod.m_marri3 <- summary(mod.m_marri3)

coef_marri <- c(mod.m_marri1$coefficients[1],mod.m_marri2$coefficients[1],mod.m_marri3$coefficients[1])
se_marri <- c(mod.m_marri1$coefficients[3],mod.m_marri2$coefficients[3],mod.m_marri3$coefficients[3])
n_marri <- c(mod.m_marri1$n,mod.m_marri2$n,mod.m_marri3$n)
marri.gen <- metagen(coef_marri, se_marri, n.e = n_marri, sm = "HR")
marri.gen #0.0013
forest(marri.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))


mod.m_histo1 <- summary(mod.m_histo1)
mod.m_histo2 <- summary(mod.m_histo2)
mod.m_histo3 <- summary(mod.m_histo3)

coef_histo <- c(mod.m_histo1$coefficients[1],mod.m_histo2$coefficients[1],mod.m_histo3$coefficients[1])
se_histo <- c(mod.m_histo1$coefficients[3],mod.m_histo2$coefficients[3],mod.m_histo3$coefficients[3])
n_histo <- c(mod.m_histo1$n,mod.m_histo2$n,mod.m_histo3$n)
histo.gen <- metagen(coef_histo, se_histo, n.e = n_histo, sm = "HR")
histo.gen #0.1539
forest(histo.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))


mod.m_grade1 <- summary(mod.m_grade1)
mod.m_grade2 <- summary(mod.m_grade2)
mod.m_grade3 <- summary(mod.m_grade3)

coef_grade <- c(mod.m_grade1$coefficients[1],mod.m_grade2$coefficients[1],mod.m_grade3$coefficients[1])
se_grade <- c(mod.m_grade1$coefficients[3],mod.m_grade2$coefficients[3],mod.m_grade3$coefficients[3])
n_grade <- c(mod.m_grade1$n,mod.m_grade2$n,mod.m_grade3$n)
grade.gen <- metagen(coef_grade, se_grade, n.e = n_grade, sm = "HR")
grade.gen #0.7187
forest(grade.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))


mod.m_T1 <- summary(mod.m_T1)
mod.m_T2 <- summary(mod.m_T2)
mod.m_T3 <- summary(mod.m_T3)
mod.m_T4 <- summary(mod.m_T4)

coef_T <- c(mod.m_T1$coefficients[1],mod.m_T2$coefficients[1],mod.m_T3$coefficients[1],mod.m_T4$coefficients[1])
se_T <- c(mod.m_T1$coefficients[3],mod.m_T2$coefficients[3],mod.m_T3$coefficients[3],mod.m_T4$coefficients[3])
n_T <- c(mod.m_T1$n,mod.m_T2$n,mod.m_T3$n,mod.m_T4$n)
T.gen <- metagen(coef_T, se_T, n.e = n_T, sm = "HR")
T.gen #< 0.0001
forest(T.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))


mod.m_N0 <- summary(mod.m_N0)
mod.m_N1 <- summary(mod.m_N1)
mod.m_N2 <- summary(mod.m_N2)
mod.m_N3 <- summary(mod.m_N3)

coef_N <- c(mod.m_N0$coefficients[1],mod.m_N1$coefficients[1],mod.m_N2$coefficients[1],mod.m_N3$coefficients[1])
se_N <- c(mod.m_N0$coefficients[3],mod.m_N1$coefficients[3],mod.m_N2$coefficients[3],mod.m_N3$coefficients[3])
n_N <- c(mod.m_N0$n,mod.m_N1$n,mod.m_N2$n,mod.m_N3$n)
N.gen <- metagen(coef_N, se_N, n.e = n_N, sm = "HR")
N.gen #< 0.0001
forest(N.gen, leftcols = c("studlab", "n.e", "TE", "seTE"))


#####Score#####
coef_age
coef_marri
coef_T
coef_N

age.beta<-coef_age
age.beta.min<-age.beta[which.min(age.beta)]
age.beta1<-age.beta-age.beta.min

marri.beta<-coef_marri
marri.beta.min<-marri.beta[which.min(marri.beta)]
marri.beta1<-marri.beta-marri.beta.min

Tstage.beta<-coef_T
Tstage.beta.min<-Tstage.beta[which.min(Tstage.beta)]
Tstage.beta1<-Tstage.beta-Tstage.beta.min

Nstage.beta<-coef_N
Nstage.beta.min<-Nstage.beta[which.min(Nstage.beta)]
Nstage.beta1<-Nstage.beta-Nstage.beta.min

coef <- c(age.beta1,marri.beta1,Tstage.beta1,Nstage.beta1)

lmaxbeta<-which.max(abs(coef))
maxbeta<-abs(coef[lmaxbeta])


data_match$agers<-data_match$age4
age.scale<-(age.beta1/maxbeta*100)
data_match$agers <- as.numeric(data_match$agers)
data_match$agers[data_match$agers==1] <- age.scale[1]
data_match$agers[data_match$agers==2] <- age.scale[2]
data_match$agers[data_match$agers==3] <- age.scale[3]
data_match$agers[data_match$agers==4] <- age.scale[4]

data_match$marrirs<-data_match$marry
marri.scale<-(marri.beta1/maxbeta*100)
data_match$marrirs <- as.numeric(data_match$marrirs)
data_match$marrirs[data_match$marrirs==1] <- marri.scale[1]
data_match$marrirs[data_match$marrirs==2] <- marri.scale[2]
data_match$marrirs[data_match$marrirs==3] <- marri.scale[3]

data_match$Tstagers<-data_match$stageT
Tstage.scale<-(Tstage.beta1/maxbeta*100)
data_match$Tstagers <- as.numeric(data_match$Tstagers)
data_match$Tstagers[data_match$Tstagers==1] <- Tstage.scale[1]
data_match$Tstagers[data_match$Tstagers==2] <- Tstage.scale[2]
data_match$Tstagers[data_match$Tstagers==3] <- Tstage.scale[3]
data_match$Tstagers[data_match$Tstagers==4] <- Tstage.scale[4]

data_match$Nstagers<-data_match$N
Nstage.scale<-(Nstage.beta1/maxbeta*100)
data_match$Nstagers <- as.numeric(data_match$Nstagers)
data_match$Nstagers[data_match$Nstagers==0] <- Nstage.scale[0]
data_match$Nstagers[data_match$Nstagers==1] <- Nstage.scale[1]
data_match$Nstagers[data_match$Nstagers==2] <- Nstage.scale[2]
data_match$Nstagers[data_match$Nstagers==3] <- Nstage.scale[3]

data_match$rs <- data_match$agers+data_match$marrirs+data_match$Tstagers+data_match$Nstagers
table(data_match$rs) #The lower point, the more favorable chemotherapy.

breakpoint <- fivenum(data_match$rs)
breakpoint <- c(0,171.4203,202.9886,245.7139,999)
labels <- c("low", "medium", "high", "max")
data_match$scoregroup <-cut(data_match$rs,breakpoint,labels,ordered_result = T)

#####Validate score#####
scoregroup1 <- subset(data_match,data_match$scoregroup=="low")
attach(scoregroup1)
mod.scoregroup1 <- coxph(Surv(Survivalmonths,death)~chemotherapy,scoregroup1)
summary(mod.scoregroup1)
detach(scoregroup1)

scoregroup2 <- subset(data_match,data_match$scoregroup=="medium")
attach(scoregroup2)
mod.scoregroup2 <- coxph(Surv(Survivalmonths,death)~chemotherapy,scoregroup2)
summary(mod.scoregroup2)
detach(scoregroup2)

scoregroup3 <- subset(data_match,data_match$scoregroup=="high")
attach(scoregroup3)
mod.scoregroup3 <- coxph(Surv(Survivalmonths,death)~chemotherapy,scoregroup3)
summary(mod.scoregroup3)
detach(scoregroup3)

scoregroup4 <- subset(data_match,data_match$scoregroup=="max")
attach(scoregroup4)
mod.scoregroup4 <- coxph(Surv(Survivalmonths,death)~chemotherapy,scoregroup4)
summary(mod.scoregroup4)
detach(scoregroup4)

#forest_score <- read.csv("/Users/linhaiping/Documents/Breast Cancer_HER2-HR+/Forest_score.csv",header = TRUE,na.strings = "NA")
head(forest_score,50)
labeltext_match <- as.matrix(forest_score[,c(1:3,7)])
attach(forest_score)
forestplot(labeltext_match,
           HR,
           Hrdown,
           Hrup,
           align = "l",
           graph.pos = 3,
           title="Hazard Ratio Plot",
           xlab="<---Favor chemotherapy---                      ---Favor non-chemotherapy--->",
           is.summary= c(T,T,F,T,
                         T,F,F,F,F),
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8),
                          xlab=gpar(cex=0.8),
                          title=gpar(cex=1.0)),
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           zero=1,
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           lwd.ci=1, boxsize=0.2,
           ci.vertices=TRUE, ci.vertices.height = 0.2,
           xlog=FALSE
)
