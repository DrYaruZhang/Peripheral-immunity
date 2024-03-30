library(dplyr)
library(survival)
library(ggfortify)
library(survminer)
library(forcats)

####data process####
Blood_UKB_bl$WBC_30000_0.0[Blood_UKB_bl$WBC_30000_0.0<=0]<-NA
Blood_UKB_bl$Basophill_30160_0.0[Blood_UKB_bl$Basophill_30160_0.0<=0]<-NA
Blood_UKB_bl$Eosinophill_30150_0.0[Blood_UKB_bl$Eosinophill_30150_0.0<=0]<-NA
Blood_UKB_bl$Neutrophill_30140_0.0[Blood_UKB_bl$Neutrophill_30140_0.0<=0]<-NA
Blood_UKB_bl$Monocyte_30130_0.0[Blood_UKB_bl$Monocyte_30130_0.0<=0]<-NA
Blood_UKB_bl$Lymphocyte_30120_0.0[Blood_UKB_bl$Lymphocyte_30120_0.0<=0]<-NA
Blood_UKB_bl$Platelet_30080_0.0[Blood_UKB_bl$Platelet_30080_0.0<=0]<-NA

#De extremum
Blood_UKB_bl$WBC_30000_0.0[Blood_UKB_bl$WBC_30000_0.0<(mean(Blood_UKB_bl$WBC_30000_0.0,na.rm=TRUE)-3*sd(Blood_UKB_bl$WBC_30000_0.0,na.rm=TRUE))
                           | Blood_UKB_bl$WBC_30000_0.0>(mean(Blood_UKB_bl$WBC_30000_0.0,na.rm=TRUE)+3*sd(Blood_UKB_bl$WBC_30000_0.0,na.rm=TRUE))]<-NA
Blood_UKB_bl$Basophill_30160_0.0[Blood_UKB_bl$Basophill_30160_0.0<(mean(Blood_UKB_bl$Basophill_30160_0.0,na.rm=TRUE)-3*sd(Blood_UKB_bl$Basophill_30160_0.0,na.rm=TRUE))
                                 | Blood_UKB_bl$Basophill_30160_0.0>(mean(Blood_UKB_bl$Basophill_30160_0.0,na.rm=TRUE)+3*sd(Blood_UKB_bl$Basophill_30160_0.0,na.rm=TRUE))]<-NA
Blood_UKB_bl$Eosinophill_30150_0.0[Blood_UKB_bl$Eosinophill_30150_0.0<(mean(Blood_UKB_bl$Eosinophill_30150_0.0,na.rm=TRUE)-3*sd(Blood_UKB_bl$Eosinophill_30150_0.0,na.rm=TRUE))
                                   | Blood_UKB_bl$Eosinophill_30150_0.0>(mean(Blood_UKB_bl$Eosinophill_30150_0.0,na.rm=TRUE)+3*sd(Blood_UKB_bl$Eosinophill_30150_0.0,na.rm=TRUE))]<-NA
Blood_UKB_bl$Neutrophill_30140_0.0[Blood_UKB_bl$Neutrophill_30140_0.0<(mean(Blood_UKB_bl$Neutrophill_30140_0.0,na.rm=TRUE)-3*sd(Blood_UKB_bl$Neutrophill_30140_0.0,na.rm=TRUE))
                                   | Blood_UKB_bl$Neutrophill_30140_0.0>(mean(Blood_UKB_bl$Neutrophill_30140_0.0,na.rm=TRUE)+3*sd(Blood_UKB_bl$Neutrophill_30140_0.0,na.rm=TRUE))]<-NA
Blood_UKB_bl$Monocyte_30130_0.0[Blood_UKB_bl$Monocyte_30130_0.0<(mean(Blood_UKB_bl$Monocyte_30130_0.0,na.rm=TRUE)-3*sd(Blood_UKB_bl$Monocyte_30130_0.0,na.rm=TRUE))
                                | Blood_UKB_bl$Monocyte_30130_0.0>(mean(Blood_UKB_bl$Monocyte_30130_0.0,na.rm=TRUE)+3*sd(Blood_UKB_bl$Monocyte_30130_0.0,na.rm=TRUE))]<-NA
Blood_UKB_bl$Lymphocyte_30120_0.0[Blood_UKB_bl$Lymphocyte_30120_0.0<(mean(Blood_UKB_bl$Lymphocyte_30120_0.0,na.rm=TRUE)-3*sd(Blood_UKB_bl$Lymphocyte_30120_0.0,na.rm=TRUE))
                                  | Blood_UKB_bl$Lymphocyte_30120_0.0>(mean(Blood_UKB_bl$Lymphocyte_30120_0.0,na.rm=TRUE)+3*sd(Blood_UKB_bl$Lymphocyte_30120_0.0,na.rm=TRUE))]<-NA
Blood_UKB_bl$Platelet_30080_0.0[Blood_UKB_bl$Platelet_30080_0.0<(mean(Blood_UKB_bl$Platelet_30080_0.0,na.rm=TRUE)-3*sd(Blood_UKB_bl$Platelet_30080_0.0,na.rm=TRUE))
                                | Blood_UKB_bl$Platelet_30080_0.0>(mean(Blood_UKB_bl$Platelet_30080_0.0,na.rm=TRUE)+3*sd(Blood_UKB_bl$Platelet_30080_0.0,na.rm=TRUE))]<-NA

#Z = (log10(value)-mean(log10(value))/SD(log10(value))
Blood_UKB_bl$log_WBC<-log10(Blood_UKB_bl$WBC_30000_0.0)
Blood_UKB_bl$Z_WBC<-((Blood_UKB_bl$log_WBC-mean(Blood_UKB_bl$log_WBC,na.rm=TRUE)))/sd(Blood_UKB_bl$log_WBC,na.rm=TRUE)

Blood_UKB_bl$log_Basophill<-log10(Blood_UKB_bl$Basophill_30160_0.0)
Blood_UKB_bl$Z_Basophill<-((Blood_UKB_bl$log_Basophill-mean(Blood_UKB_bl$log_Basophill,na.rm=TRUE)))/sd(Blood_UKB_bl$log_Basophill,na.rm=TRUE)

Blood_UKB_bl$log_Eosinophill<-log10(Blood_UKB_bl$Eosinophill_30150_0.0)
Blood_UKB_bl$Z_Eosinophill<-((Blood_UKB_bl$log_Eosinophill-mean(Blood_UKB_bl$log_Eosinophill,na.rm=TRUE)))/sd(Blood_UKB_bl$log_Eosinophill,na.rm=TRUE)

Blood_UKB_bl$log_Neutrophill<-log10(Blood_UKB_bl$Neutrophill_30140_0.0)
Blood_UKB_bl$Z_Neutrophill<-((Blood_UKB_bl$log_Neutrophill-mean(Blood_UKB_bl$log_Neutrophill,na.rm=TRUE)))/sd(Blood_UKB_bl$log_Neutrophill,na.rm=TRUE)

Blood_UKB_bl$log_Monocyte<-log10(Blood_UKB_bl$Monocyte_30130_0.0)
Blood_UKB_bl$Z_Monocyte<-((Blood_UKB_bl$log_Monocyte-mean(Blood_UKB_bl$log_Monocyte,na.rm=TRUE)))/sd(Blood_UKB_bl$log_Monocyte,na.rm=TRUE)

Blood_UKB_bl$log_Lymphocyte<-log10(Blood_UKB_bl$Lymphocyte_30120_0.0)
Blood_UKB_bl$Z_Lymphocyte<-((Blood_UKB_bl$log_Lymphocyte-mean(Blood_UKB_bl$log_Lymphocyte,na.rm=TRUE)))/sd(Blood_UKB_bl$log_Lymphocyte,na.rm=TRUE)

Blood_UKB_bl$log_Platelet<-log10(Blood_UKB_bl$Platelet_30080_0.0)
Blood_UKB_bl$Z_Platelet<-((Blood_UKB_bl$log_Platelet-mean(Blood_UKB_bl$log_Platelet,na.rm=TRUE)))/sd(Blood_UKB_bl$log_Platelet,na.rm=TRUE)

Blood_UKB_bl$NLR<-Blood_UKB_bl$Neutrophill_30140_0.0/Blood_UKB_bl$Lymphocyte_30120_0.0
Blood_UKB_bl$log_NLR<-log10(Blood_UKB_bl$NLR)
Blood_UKB_bl$Z_NLR<-((Blood_UKB_bl$log_NLR-mean(Blood_UKB_bl$log_NLR,na.rm=TRUE)))/sd(Blood_UKB_bl$log_NLR,na.rm=TRUE)

Blood_UKB_bl$PLR<-Blood_UKB_bl$Platelet_30080_0.0/Blood_UKB_bl$Lymphocyte_30120_0.0
Blood_UKB_bl$log_PLR<-log10(Blood_UKB_bl$PLR)
Blood_UKB_bl$Z_PLR<-((Blood_UKB_bl$log_PLR-mean(Blood_UKB_bl$log_PLR,na.rm=TRUE)))/sd(Blood_UKB_bl$log_PLR,na.rm=TRUE)

Blood_UKB_bl$SII<-(Blood_UKB_bl$Neutrophill_30140_0.0*Blood_UKB_bl$Platelet_30080_0.0)/Blood_UKB_bl$Lymphocyte_30120_0.0
Blood_UKB_bl$log_SII<-log10(Blood_UKB_bl$SII)
Blood_UKB_bl$Z_SII<-((Blood_UKB_bl$log_SII-mean(Blood_UKB_bl$log_SII,na.rm=TRUE)))/sd(Blood_UKB_bl$log_SII,na.rm=TRUE)

Blood_UKB_bl$LMR<-Blood_UKB_bl$Lymphocyte_30120_0.0/Blood_UKB_bl$Monocyte_30130_0.0
Blood_UKB_bl$log_LMR<-log10(Blood_UKB_bl$LMR)
Blood_UKB_bl$Z_LMR<-((Blood_UKB_bl$log_LMR-mean(Blood_UKB_bl$log_LMR,na.rm=TRUE)))/sd(Blood_UKB_bl$log_LMR,na.rm=TRUE)

u1<-subset(Blood_UKB_bl, is.na(Blood_UKB_bl$WBC_30000_0.0) & is.na(Blood_UKB_bl$Basophill_30160_0.0) & is.na(Blood_UKB_bl$Eosinophill_30150_0.0)
           & is.na(Blood_UKB_bl$Neutrophill_30140_0.0) & is.na(Blood_UKB_bl$Monocyte_30130_0.0) & is.na(Blood_UKB_bl$Lymphocyte_30120_0.0)
           & is.na(Blood_UKB_bl$Platelet_30080_0.0))
u2<-subset(Blood_UKB_bl,!Blood_UKB_bl$eid %in% u1$eid)
#write.csv(u2,"Blood_UKB_bl_Zscore.csv")

#exclusion of dementia patients at baseline
ACD_losefl<-subset(ACD_AD_VD_covariants,is.na(ACD_AD_VD_covariants$dementia_days)|(ACD_AD_VD_covariants$dementia_days<=0 & ACD_AD_VD_covariants$dementia_status==0))
ACD_bl<-subset(ACD_AD_VD_covariants,ACD_AD_VD_covariants$dementia_days<=0 & ACD_AD_VD_covariants$dementia_status==1)
sum(ACD_bl$AD_days<=0 & ACD_bl$AD_status==1)
sum(ACD_bl$VD_days<=0 & ACD_bl$VD_status==1)
ACD_AD_VD_covariants$dementia_status[ACD_AD_VD_covariants$dementia_days<=0]<-NA
x1<-subset(ACD_AD_VD_covariants,!is.na(ACD_AD_VD_covariants$dementia_status))

AD_losefl<-subset(x1,is.na(x1$AD_days)|(x1$AD_days<=0 & x1$AD_status==0))
AD_bl<-subset(x1,x1$AD_days<=0 & x1$AD_status==1)
x1$AD_status[x1$AD_days<=0]<-NA

VD_losefl<-subset(x1,is.na(x1$VD_days)|(x1$VD_days<=0 & x1$VD_status==0))
VD_bl<-subset(x1,x1$VD_days<=0 & x1$VD_status==1)
x1$VD_status[x1$VD_days<=0]<-NA

dementia_bl<-merge(ACD_bl,AD_bl,by="eid",all=TRUE)
dementia_losefl<-ACD_losefl
u2_dementia_bl<-subset(u2,u2$eid %in% dementia_bl$eid)
u2_dementia_losefl<-subset(u2,u2$eid %in% dementia_losefl$eid)

x2<-subset(x1,!is.na(x1$AD_status))
sum(is.na(x2$AD_status))
sum(is.na(x2$VD_status))

x2$AD_status[x2$AD_status==0 & x2$dementia_status==1]<-NA
x2$VD_status[x2$VD_status==0 & x2$dementia_status==1]<-NA
sum(is.na(x2$AD_status))
sum(is.na(x2$VD_status))

a<-merge(u2,x2,by="eid")
x<-subset(a,a$dementia_date=="2021-04-07")
a$dementia_date[a$dementia_date=="2037-07-07"]<-"2021-04-07"
a$AD_date[a$AD_date=="2037-07-07"]<-"2021-04-07"
a$VD_date[a$VD_date=="2037-07-07"]<-"2021-04-07"
a$dementia_days[a$dementia_date=="2021-04-07"]<-x$dementia_days
a$AD_days[a$AD_date=="2021-04-07"]<-x$dementia_days
a$VD_days[a$VD_date=="2021-04-07"]<-x$dementia_days
a$dementia_months[a$dementia_date=="2021-04-07"]<-x$dementia_months
a$AD_months[a$AD_date=="2021-04-07"]<-x$dementia_months
a$VD_months[a$VD_date=="2021-04-07"]<-x$dementia_months
a$dementia_years[a$dementia_date=="2021-04-07"]<-x$dementia_years
a$AD_years[a$AD_date=="2021-04-07"]<-x$dementia_years
a$VD_years[a$VD_date=="2021-04-07"]<-x$dementia_years
name_a<-colnames(a[,c(25,11,17,15,13,19,23,21,28,31,34,37)])
name_a
name_a<-colnames(a[,c(17,19,23,21,28,31,34,37)])
name_a


#exclusion of participants with morbidities that could influence leukocyte differential counts at baseline
Exc_a<-subset(a,a$eid %in% Excluded_population$eid)
a1<-subset(a,!a$eid %in% Excluded_population$eid)
a<-a1

####demographic####
library(dplyr)
mean(a$dementia_years)
median(a$dementia_years)
mean(a$Age)
sd(a$Age)
sum(a$Sex==1)
sum(a$Sex==1)/sum(!is.na(a$Sex))
sum(a$dementia_status==1)
sum(a$dementia_status==0)
X<-subset(a,a$AD_status==1)
X<-subset(a,a$VD_status==1)

#ACD VS CN
a %>%
  group_by(dementia_status)%>%
  summarise(followup = mean(dementia_years, na.rm=T),
            Age = mean(Age, na.rm = T),BMI = mean(BMI, na.rm = T),Townsend = mean(Townsend, na.rm = T),
            neutrophils = mean(Neutrophill_30140_0.0, na.rm = T),
            monocytes = mean(Monocyte_30130_0.0, na.rm = T),platelets = mean(Platelet_30080_0.0, na.rm = T),
            lymphocytes = mean(Lymphocyte_30120_0.0, na.rm = T),NLR = mean(NLR, na.rm = T),PLR = mean(PLR, na.rm = T),
            SII = mean(SII, na.rm = T),LMR = mean(LMR, na.rm = T))
a %>%
  group_by(dementia_status)%>%
  summarise(followup = sd(dementia_years, na.rm=T),
            Age = sd(Age, na.rm = T),BMI = sd(BMI, na.rm = T),Townsend = sd(Townsend, na.rm = T),
            neutrophils = sd(Neutrophill_30140_0.0, na.rm = T),
            monocytes = sd(Monocyte_30130_0.0, na.rm = T),platelets = sd(Platelet_30080_0.0, na.rm = T),
            lymphocytes = sd(Lymphocyte_30120_0.0, na.rm = T),NLR = sd(NLR, na.rm = T),PLR = sd(PLR, na.rm = T),
            SII = sd(SII, na.rm = T),LMR = sd(LMR, na.rm = T))

b<-subset(a,a$dementia_status==0 & !is.na(a$Sex))
c<-subset(a,a$dementia_status==1 & !is.na(a$Sex))
table(b$Sex)
table(c$Sex)

sum(!is.na(a$APOE4))
X_CN<-subset(a,!is.na(a$APOE4) & a$dementia_status==0)
X_ACD<-subset(a,!is.na(a$APOE4) & a$dementia_status==1)
X_AD<-subset(a,!is.na(a$APOE4) & a$AD_status==1)
X_VD<-subset(a,!is.na(a$APOE4) & a$VD_status==1)
sum(!is.na(a$APOE4))
d<-subset(a,a$dementia_status==0 & !is.na(a$APOE4))
e<-subset(a,a$dementia_status==1 & !is.na(a$APOE4))
table(d$APOE4)
table(e$APOE4)

sum(!is.na(a$Qualification))
X_CN<-subset(a,!is.na(a$Qualification) & a$dementia_status==0)
X_ACD<-subset(a,!is.na(a$Qualification) & a$dementia_status==1)
X_AD<-subset(a,!is.na(a$Qualification) & a$AD_status==1)
X_VD<-subset(a,!is.na(a$Qualification) & a$VD_status==1)
a$Qualification_group[a$Qualification==1 | a$Qualification==6]<-1
a$Qualification_group[a$Qualification==2 | a$Qualification==3 | a$Qualification==4 | a$Qualification==5]<-0
f<-subset(a,a$dementia_status==0 & !is.na(a$Qualification_group))
g<-subset(a,a$dementia_status==1 & !is.na(a$Qualification_group))
table(f$Qualification_group)
table(g$Qualification_group)

sum(!is.na(a$race_group))
X_CN<-subset(a,!is.na(a$race_group) & a$dementia_status==0)
X_ACD<-subset(a,!is.na(a$race_group) & a$dementia_status==1)
X_AD<-subset(a,!is.na(a$race_group) & a$AD_status==1)
X_VD<-subset(a,!is.na(a$race_group) & a$VD_status==1)
h<-subset(a,a$dementia_status==0 & !is.na(a$race_group))
i<-subset(a,a$dementia_status==1 & !is.na(a$race_group))
table(h$race_group)
table(i$race_group)

sum(!is.na(a$Smoking))
X_CN<-subset(a,!is.na(a$Smoking) & a$dementia_status==0)
X_ACD<-subset(a,!is.na(a$Smoking) & a$dementia_status==1)
X_AD<-subset(a,!is.na(a$Smoking) & a$AD_status==1)
X_VD<-subset(a,!is.na(a$Smoking) & a$VD_status==1)
j<-subset(a,a$dementia_status==0 & !is.na(a$Smoking))
k<-subset(a,a$dementia_status==1 & !is.na(a$Smoking))
table(j$Smoking)
table(k$Smoking)

sum(!is.na(a$Alcohol))
X_CN<-subset(a,!is.na(a$Alcohol) & a$dementia_status==0)
X_ACD<-subset(a,!is.na(a$Alcohol) & a$dementia_status==1)
X_AD<-subset(a,!is.na(a$Alcohol) & a$AD_status==1)
X_VD<-subset(a,!is.na(a$Alcohol) & a$VD_status==1)
l<-subset(a,a$dementia_status==0 & !is.na(a$Alcohol))
m<-subset(a,a$dementia_status==1 & !is.na(a$Alcohol))
table(l$Alcohol)
table(m$Alcohol)

#AD VS CN
a %>%
  group_by(AD_status)%>%
  summarise(followup = mean(AD_years, na.rm=T),
            Age = mean(Age, na.rm = T),BMI = mean(BMI, na.rm = T),Townsend = mean(Townsend, na.rm = T),
            neutrophils = mean(Neutrophill_30140_0.0, na.rm = T),
            monocytes = mean(Monocyte_30130_0.0, na.rm = T),platelets = mean(Platelet_30080_0.0, na.rm = T),
            lymphocytes = mean(Lymphocyte_30120_0.0, na.rm = T),NLR = mean(NLR, na.rm = T),PLR = mean(PLR, na.rm = T),
            SII = mean(SII, na.rm = T),LMR = mean(LMR, na.rm = T))
a %>%
  group_by(AD_status)%>%
  summarise(followup = sd(AD_years, na.rm=T),
            Age = sd(Age, na.rm = T),BMI = sd(BMI, na.rm = T),Townsend = sd(Townsend, na.rm = T),
            neutrophils = sd(Neutrophill_30140_0.0, na.rm = T),
            monocytes = sd(Monocyte_30130_0.0, na.rm = T),platelets = sd(Platelet_30080_0.0, na.rm = T),
            lymphocytes = sd(Lymphocyte_30120_0.0, na.rm = T),NLR = sd(NLR, na.rm = T),PLR = sd(PLR, na.rm = T),
            SII = sd(SII, na.rm = T),LMR = sd(LMR, na.rm = T))

c<-subset(a,a$AD_status==1 & !is.na(a$Sex))
table(c$Sex)
e<-subset(a,a$AD_status==1 & !is.na(a$APOE4))
table(e$APOE4)
g<-subset(a,a$AD_status==1 & !is.na(a$Qualification_group))
table(g$Qualification_group)
i<-subset(a,a$AD_status==1 & !is.na(a$race_group))
table(i$race_group)
k<-subset(a,a$AD_status==1 & !is.na(a$Smoking))
table(k$Smoking)
m<-subset(a,a$AD_status==1 & !is.na(a$Alcohol))
table(m$Alcohol)

#VD VS CN
a %>%
  group_by(VD_status)%>%
  summarise(followup = mean(VD_years, na.rm=T),
            Age = mean(Age, na.rm = T),BMI = mean(BMI, na.rm = T),Townsend = mean(Townsend, na.rm = T),
            neutrophils = mean(Neutrophill_30140_0.0, na.rm = T),
            monocytes = mean(Monocyte_30130_0.0, na.rm = T),platelets = mean(Platelet_30080_0.0, na.rm = T),
            lymphocytes = mean(Lymphocyte_30120_0.0, na.rm = T),NLR = mean(NLR, na.rm = T),PLR = mean(PLR, na.rm = T),
            SII = mean(SII, na.rm = T),LMR = mean(LMR, na.rm = T))
a %>%
  group_by(VD_status)%>%
  summarise(followup = sd(VD_years, na.rm=T),
            Age = sd(Age, na.rm = T),BMI = sd(BMI, na.rm = T),Townsend = sd(Townsend, na.rm = T),
            neutrophils = sd(Neutrophill_30140_0.0, na.rm = T),
            monocytes = sd(Monocyte_30130_0.0, na.rm = T),platelets = sd(Platelet_30080_0.0, na.rm = T),
            lymphocytes = sd(Lymphocyte_30120_0.0, na.rm = T),NLR = sd(NLR, na.rm = T),PLR = sd(PLR, na.rm = T),
            SII = sd(SII, na.rm = T),LMR = sd(LMR, na.rm = T))


c<-subset(a,a$VD_status==1 & !is.na(a$Sex))
table(c$Sex)
e<-subset(a,a$VD_status==1 & !is.na(a$APOE4))
table(e$APOE4)
g<-subset(a,a$VD_status==1 & !is.na(a$Qualification_group))
table(g$Qualification_group)
i<-subset(a,a$VD_status==1 & !is.na(a$race_group))
table(i$race_group)
k<-subset(a,a$VD_status==1 & !is.na(a$Smoking))
table(k$Smoking)
m<-subset(a,a$VD_status==1 & !is.na(a$Alcohol))
table(m$Alcohol)

####COX####
####_____unadjusted####
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a,!is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i)))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a,!is.na(AD_status),!is.na(AD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i)))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a,!is.na(VD_status),!is.na(VD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i)))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

####_____adjusted for Age+Sex+ApoE4+education####
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

####_____adjusted for Age+Sex+ApoE4+education+race+smoking+alcohol+BMI+Townsend####
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

####COX SUBGROUP####
####_____subgroup_age####
#55_65 years
a_55<-subset(a,a$Age<55)
a_55_65<-subset(a,a$Age>=55 & a$Age<65)
a_65<-subset(a,a$Age>=65)

#<55
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_55, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_55[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_55, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_55[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_55, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_55[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#55_65
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_55_65, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_55_65[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_55_65, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_55_65[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_55_65, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_55_65[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#65
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_65, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_65[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_65, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_65[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_65, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_65[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

####_____subgroup_sex####
a_Male<-subset(a,a$Sex==1)
a_Female<-subset(a,a$Sex==0)

#Male
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_Male, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_Male[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_Male, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_Male[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_Male, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_Male[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#Female
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_Female, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_Female[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_Female, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_Female[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_Female, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_Female[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

apoe4_20210906 <- read.csv("D:/Dr.Yu/UKB/Sleep-dementia/data/leinao/apoe4_20210906.csv", sep="")
colnames(apoe4_20210906)[1] <- "eid"
aa1<-subset(apoe4_20210906,apoe4_20210906$eid>0)
aaa1<-merge(a,aa1,by="eid",all = TRUE)
aaaa1<-subset(aaa1, aaa1$eid %in% a$eid)
a<-aaaa1
a_noe2e4<-subset(a,!a$apoe=="e2_e4")

####_____subgroup_APOE####
a_none_APOE4<-subset(a,a$APOE4==0)
a_one_APOE4<-subset(a,a$APOE4==1 & a$apoe == "e3_e4")
a_two_APOE4<-subset(a,a$APOE4==2)

#none_APOE4
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_none_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_none_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_none_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_none_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_none_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_none_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#one_APOE4
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_one_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_one_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_one_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_one_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_one_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_one_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#two_APOE4
#adjusted for Age+Sex+ApoE4+education
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_two_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a_two_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_two_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(AD_status),!is.na(AD_years),!is.na(a_two_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a_two_APOE4, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(VD_status),!is.na(VD_years),!is.na(a_two_APOE4[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

####_____subgroup_Followupyears####
#Followup>5 years
a_follow_over5years<-subset(a,a$dementia_years>5)
a<-a_follow_over5years
#adjusted for Age+Sex+ApoE4+education+race+smoking+alcohol+BMI+Townsend
#ACD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(dementia_years,dementia_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#AD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(AD_years,AD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2

#VD
coeff.old <- data.frame()
for (i in name_a) {
  a_noNA <- filter(a, !is.na(Age),!is.na(Sex),!is.na(APOE4),!is.na(Qualification),!is.na(race_group),!is.na(Smoking),!is.na(Alcohol),!is.na(BMI),!is.na(Townsend),
                   !is.na(dementia_status),!is.na(dementia_years),!is.na(a[,which(names(a)==i)])) 
  FML <- as.formula(paste0('Surv(VD_years,VD_status)~',paste(i,'Age','Sex','APOE4','Qualification','race_group','Smoking','Alcohol','BMI','Townsend',sep = "+")))
  pheno_cox_AD <- coxph(FML,data=a_noNA)
  
  varid <- i
  coeff <- summary(pheno_cox_AD)$conf.int[1,1]
  lower95 <- summary(pheno_cox_AD)$conf.int[1,3]
  upper95 <- summary(pheno_cox_AD)$conf.int[1,4]
  pvalue <- summary(pheno_cox_AD)$coefficients[1,5]
  
  samplesize <- nrow(a_noNA)
  case_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==1)
  control_r <- filter(a_noNA,a_noNA[,which(names(a_noNA)==i)]==0)
  coeff.new <- data.frame(varid=varid,
                          samplesize=paste0(nrow(a_noNA),'(',
                                            nrow(case_r),'/',
                                            nrow(control_r),')'),
                          coeff=coeff,
                          lower95=lower95,
                          upper95=upper95,
                          pvalue=pvalue
  )
  
  coeff.old <- rbind.data.frame(coeff.old,coeff.new)
}

coeff.2 <- coeff.old
coeff.2