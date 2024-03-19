##survival analysis here is needed
suvrdata <- data.frame(y = d$Neuropsych.MMSE.total,t = d$Age,projid = d$AIBL.Id,
                       fu_year = d$Age_diff)

##pre-process for datafit
DataFit1 <- DataFit[DataFit$Statu != "AD",]
DataFit2 <-  DataFit %>%
group_by(projid ) %>% 
subset(DataFit$Statu == "AD") %>%
distinct( projid, .keep_all = TRUE)
DataFit3 <- rbind(DataFit1, DataFit2) %>% 
group_by(projid)
DataFit <- DataFit3[order(DataFit3$projid),]

suvrdata$cog <- DataFit$cogAge
suvrdata$sex <- d$Demographic.Sex
suvrdata$cogage <- DataFit$cogAge

#For cognitive state
suvrdata$Classification <- d$Neuropsych.Simple.Classification
suvrdata$censored[suvrdata$Classification == "AD"] <- 1
suvrdata$censored[suvrdata$Classification == "MCI"] <- 0
suvrdata$censored[suvrdata$Classification == "HC"] <- 0

#For APOE gene type
suvrdata$Gene <- d$Demographic.ApoE.genotype
subsuvrdata <- subset(suvrdata, suvrdata$Gene=="E2/E2"|suvrdata$Gene=="E3/E2"|suvrdata$Gene=="E3/E3"
                      |suvrdata$Gene=="E4/E2"|suvrdata$Gene=="E4/E3"|suvrdata$Gene=="E4/E4")
subsuvrdata$Gene[subsuvrdata$Gene == "E4/E4"|
                   subsuvrdata$Gene == "E4/E3"|subsuvrdata$Gene == "E4/E2"]<- "E4 Carrier"
subsuvrdata$Gene[subsuvrdata$Gene == "E2/E2"|
                   subsuvrdata$Gene == "E3/E2"|subsuvrdata$Gene == "E3/E3"]<- "E4 non-Carrier"

#For Education
suvrdata$Edu <- d$Demographic.Years.of.Education.Exact
suvrdataedu <- subset(suvrdata, is.na(suvrdata$Edu)!= TRUE)
suvrdataedu$Edu[suvrdataedu$Edu < 15]<- 0
suvrdataedu$Edu[suvrdataedu$Edu > 15|suvrdataedu$Edu == 15]<- 1
suvrdataedu$Edu[suvrdataedu$Edu == 0] <- "low edu"
suvrdataedu$Edu[suvrdataedu$Edu != "low edu"] <- "higher edu"

#For site
suvrdata$site <- d$Demographic.Site
subsuvrdatasite <- subset(suvrdata, is.na(suvrdata$site)!=TRUE)

library(survival)
#Total survival analysis
km_fit <- survfit(Surv(t,censored) ~ 1,data=suvrdata, type="kaplan-meier")
#Sex difference
km_fit2 <- survfit(Surv(cog,censored) ~ sex,data=suvrdata, type="kaplan-meier")
#APOE4 difference
km_fit3 <- survfit(Surv(cog,censored) ~ Gene,
                   data=subsuvrdata, 
                   type="kaplan-meier")

#Site difference
km_fit5 <- survfit(Surv(cog,censored) ~ site,data=subsuvrdatasite, type="kaplan-meier")
#Edu difference
km_fit6 <- survfit(Surv(cog,censored) ~ Edu,data=suvrdataedu, type="kaplan-meier")
library(survminer)
ggsurvplot(km_fit, data=suvrdata,xlim=c(70,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve",ggtheme=theme_minimal())
ggsurvplot(km_fit2, data=suvrdata,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for sex", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
ggsurvplot(km_fit3, data=subsuvrdata,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for Gene", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
ggsurvplot(km_fit5, data=subsuvrdatasite,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for site",pval = TRUE, pval.method = TRUE,
           ggtheme=theme_minimal())
ggsurvplot(km_fit6, data=suvrdataedu,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for Edu", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())

##survival for medical history (cancer)
suvrdata$cancer <- d$Medical.History.Cancer
suvrdatacancer <- subset(suvrdata, suvrdata$cancer=="Yes"|suvrdata$cancer=="No")
suvrdatacancer$cancer[suvrdatacancer$cancer == "Yes"]<- 1
suvrdatacancer$cancer[suvrdatacancer$cancer == "No"]<- 0
km_fitcancer <-  survfit(Surv(cog,censored) ~ cancer,data=suvrdatacancer, type="kaplan-meier")
ggsurvplot(km_fitcancer, data=suvrdatacancer,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for cancer", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
km_fitcancer


##survival for Abeta
suvrdata$pet <- d$Image.PET.Amyloid.Status
suvrdatapet <- subset(suvrdata, suvrdata$pet=="Negative"|suvrdata$pet=="Moderate"
                       |suvrdata$pet=="High"|suvrdata$pet=="VeryHigh")
suvrdatapet$pet[suvrdatapet$pet=="Negative"] <- 0
suvrdatapet$pet[suvrdatapet$pet=="Moderate"
                      |suvrdatapet$pet=="High"|suvrdatapet$pet=="VeryHigh"]<- 1
km_fitpet <-  survfit(Surv(cog,censored) ~ pet,data=suvrdatapet, type="kaplan-meier")
ggsurvplot(km_fitpet, data=suvrdatapet,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for cancer", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
km_fitpet
surv_pvalue(km_fitpet, suvrdatapet)

##survival for Medical History.Hypertension
suvrdata$Hypertension <- d$Medical.History.Hypertension
suvrdataHypertension <- subset(suvrdata, suvrdata$Hypertension=="Yes"|suvrdata$Hypertension=="No")
suvrdataHypertension$Hypertension[suvrdataHypertension$Hypertension == "Yes"]<- 1
suvrdataHypertension$Hypertension[suvrdataHypertension$Hypertension == "No"]<- 0
km_fitHypertension <-  survfit(Surv(cog,censored) ~ Hypertension,data=suvrdataHypertension, type="kaplan-meier")
ggsurvplot(km_fitHypertension, data=suvrdataHypertension,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for Hypertension", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
km_fitHypertension

##survival for medical history (stroke)
suvrdata$stroke <- d$Medical.History.Stroke
suvrdatastroke <- subset(suvrdata, suvrdata$stroke=="Yes"|suvrdata$stroke=="No")
suvrdatastroke$stroke[suvrdatastroke$stroke == "Yes"]<- 1
suvrdatastroke$stroke[suvrdatastroke$stroke == "No"]<- 0
km_fitstroke <-  survfit(Surv(cog,censored) ~ stroke,data=suvrdatastroke, type="kaplan-meier")
ggsurvplot(km_fitstroke, data=suvrdatastroke,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for stroke", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
km_fitstroke
surv_pvalue(km_fitstroke, suvrdatastroke)

##survival for medical history (Gastric)
suvrdata$Gas <- d$Medical.History.Gastric.Complaints
suvrdataGas <- subset(suvrdata, suvrdata$Gas=="Yes"|suvrdata$Gas=="No")
suvrdataGas$Gas[suvrdataGas$Gas == "Yes"]<- 1
suvrdataGas$Gas[suvrdataGas$Gas == "No"]<- 0
km_fitGas <-  survfit(Surv(cog,censored) ~ Gas,data=suvrdataGas, type="kaplan-meier")
ggsurvplot(km_fitGas, data=suvrdataGas,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for Gas", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
km_fitGas

##survival for medical history (Neuropsychological Disorder other than AD)
suvrdata$Disorder <- d$Medical.History.Neurological.Disorders
suvrdataDisorder <- subset(suvrdata, suvrdata$Disorder=="Yes"|suvrdata$Disorder=="No")
suvrdataDisorder$Disorder[suvrdataDisorder$Disorder == "Yes"]<- 1
suvrdataDisorder$Disorder[suvrdataDisorder$Disorder == "No"]<- 0
km_fitDisorder <-  survfit(Surv(cog,censored) ~ Disorder,data=suvrdataDisorder, type="kaplan-meier")
ggsurvplot(km_fitDisorder, data=suvrdataDisorder,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for Disorder", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
km_fitDisorder

##survival for medical history (Psychiatric)
suvrdata$Psychiatric <- d$Medical.History.Psychiatric.Disorders
suvrdataPsychiatric <- subset(suvrdata, suvrdata$Psychiatric=="Yes"|suvrdata$Psychiatric=="No")
suvrdataPsychiatric$Psychiatric[suvrdataPsychiatric$Psychiatric == "Yes"]<- 1
suvrdataPsychiatric$Psychiatric[suvrdataPsychiatric$Psychiatric == "No"]<- 0
km_fitPsychiatric <-  survfit(Surv(cog,censored) ~ Psychiatric,data=suvrdataPsychiatric, type="kaplan-meier")
ggsurvplot(km_fitPsychiatric, data=suvrdataPsychiatric,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
           legend.title = "Survival Curve for Psychiatric Disorder", pval=TRUE,pval.method = TRUE,
           ggtheme=theme_minimal())
km_fitPsychiatric








library(timeROC)
library(survival)
##For AD predict, chro v.s. cog
ROC.bili.marginal<-timeROC(T=suvrdata$cogage,
                           delta=suvrdata$censored,marker=suvrdata$cog,
                           cause=1,weighting="marginal",
                           times=c(8,9,10,11,12),ROC=TRUE)
ROC.bili.marginal

plot(ROC.bili.marginal,time=8)        
plot(ROC.bili.marginal,time=10,add=TRUE,col="blue") 
plot(ROC.bili.marginal,time=12,add=TRUE,col="grey50") 
legend("bottomright",c("Y-8","Y-10","Y-12"),col=c("red","blue","grey50"),lty=1,lwd=2)

##For mci predict, chro v.s. cog
mcidata <- subset(suvrdata,suvrdata$Classification=="HC"|suvrdata$Classification=="MCI"
                  |suvrdata$Classification=="AD")
mcidata$censored[mcidata$Classification == "MCI"] <- 1
mcidata$censored[mcidata$Classification == "HC"] <- 0
mcidata$censored[mcidata$Classification == "AD"] <- 0
ROC.bili.marginal1<-timeROC(T=mcidata$fu_year,
                           delta=mcidata$censored,marker=mcidata$cog,
                           cause=1,weighting="marginal",
                           times=c(1,3,5,7,9),ROC=TRUE,iid=TRUE)
ROC.bili.marginal1
ROC.bili.marginal2<-timeROC(T=mcidata$fu_year,
                            delta=mcidata$censored,marker=mcidata$t,
                            cause=1,weighting="marginal",
                            times=c(1,3,5,7,9),ROC=TRUE,iid=TRUE)
ROC.bili.marginal2
compare(ROC.bili.marginal1,ROC.bili.marginal2,adjusted=TRUE) #compute p-values of comparison tests
plotAUCcurve(ROC.bili.marginal1,conf.int=TRUE,col="red")
plotAUCcurve(ROC.bili.marginal2,conf.int=TRUE,col="blue",add=TRUE)
legend("topleft",c("cognitive age","chronological age"),col=c("red","blue"),lty=1,lwd=2)

##For AD predict, chro v.s. cog
ADdata <- subset(suvrdata,suvrdata$Classification=="HC"|suvrdata$Classification=="AD"|
                    suvrdata$Classification=="MCI")
ADdata$censored[ADdata$Classification == "AD"] <- 1
ADdata$censored[ADdata$Classification == "HC"] <- 0
ADdata$censored[ADdata$Classification == "MCI"] <- 0
ROC.bili.marginal11<-timeROC(T=ADdata$fu_year,
                            delta=ADdata$censored,marker=ADdata$cog,
                            cause=1,weighting="marginal",
                            times=c(1,3,5,7,9),ROC=TRUE,iid=TRUE)
ROC.bili.marginal11
ROC.bili.marginal22<-timeROC(T=ADdata$fu_year,
                            delta=ADdata$censored,marker=ADdata$t,
                            cause=1,weighting="marginal",
                            times=c(1,3,5,7,9),ROC=TRUE,iid=TRUE)
ROC.bili.marginal22
compare(ROC.bili.marginal11,ROC.bili.marginal22,adjusted=TRUE) #compute p-values of comparison tests
plotAUCcurve(ROC.bili.marginal11,conf.int=TRUE,col="red")
plotAUCcurve(ROC.bili.marginal22,conf.int=TRUE,col="blue",add=TRUE)
legend("left",c("cognitive age","chronological age"),col=c("red","blue"),lty=1,lwd=2)

