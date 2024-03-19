#######################################################################
### Cognitive Clock Estimation
#######################################################################
library(nlme)
library(lme4)
library(splines)
library(splines2)
library(gmodels)
library(pracma)
library(dplyr,tidyr)
source('C:/Users/11373/Desktop/code/Functions2.1.R')


#### Data available upon formal request to RADC hub
# DataFit = read.csv("PR934r1updated.csv",header = T, na.strings = list("NA",'.'))
# DataFit$y = DataFit$cts_mmse30
# DataFit$t = DataFit$age_at_visit
# DataFit = DataFit[!is.na(DataFit$y + DataFit$t),]

#### Data details
# DataFit is a long form data frame with four entries
# DataFit$y is the longitudinal cognition outcome
# DataFit$t is the measurement time
# DataFit$projid is the id number of each subject
# DataFit$fu_year is the index of follow-up year (The baseline is coded as fu_year=0)

#######################
setwd('C:/Users/11373/Desktop/code')
d0 <- read.csv("C:/Users/11373/Desktop/code/ADNI_merge.csv") 


# d1 <- d0 %>%
#   group_by(AIBL.Id) %>%
#   mutate(Age_diff = Age - first(Age))

# used for Adni data set 
d0$AGE <- round(d0$AGE)
d0 <- distinct(d0, RID, AGE, .keep_all = T)

# Reorder the records for each participant from small to big and then big to small
d0 <- d0 %>%
  arrange(d0$RID, d0$AGE, desc(d0$AGE))

d1 <- d0 %>%
  group_by(RID) %>%
  mutate(Age_diff = AGE - first(AGE))

DataFit <- data.frame(y = d1$MMSE,t = round(d1$AGE),projid = d1$RID,
                     fu_year = d1$Age_diff)
# removing duplicated rows, 1 year multiple measurements
DataFit <- distinct(DataFit,t, projid, .keep_all = T )







####################


#### Estimate age-specific MMSE mean trajectory by B-spline with equal-spaced interior knots
## Number of interior knots
K = 4
## Define the B-spline basis functions
knot_BD = c(55,100)
knot_All = seq(knot_BD[1],knot_BD[2],by=(knot_BD[2]-knot_BD[1])/(K+1))
knot_Int = knot_All[c(-1,-length(knot_All))]
BS_Fit = bSpline(DataFit$t,knots=knot_Int,degree = 4L,intercept=TRUE,Boundary.knots=knot_BD)
DataFit$BS_Fit = BS_Fit

## Estimation by mixed-effects model
lme.fit = lme(y ~ BS_Fit-1, 
              random = ~t|projid, 
              data = DataFit,
              method = 'REML',
              na.action = na.omit,
              control = lmeControl(
                MaxIter=200, msMaxIter=200, msMaxEval=1000, niterEM=100,
                msVerbose=T, opt = 'optim'))
summary(lme.fit)
pBS = ncol(BS_Fit)
lme.fit.BS = lme.fit$coefficients$fixed[1:pBS]
ymean = BS_Fit%*%lme.fit.BS

## Ensure monotonicity
lme.fit.BS.d = diff(lme.fit.BS)
lme.fit.BS.d = lme.fit.BS.d*I(lme.fit.BS.d<=0) - 0.01
pBS.hf = floor(pBS/2+1)
b.mmse = lme.fit.BS
b.mmse[pBS.hf:pBS] = cumsum(c(lme.fit.BS[pBS.hf],lme.fit.BS.d[pBS.hf:(pBS-1)]))
b.mmse[pBS.hf:1] = cumsum(c(lme.fit.BS[pBS.hf],-lme.fit.BS.d[(pBS.hf-1):1]))

# Plot the longitudinal MMSE trajectories and the estimated mean trajectory
t_plot = 60:100
BS_plot = bSpline(t_plot,knots=knot_Int,degree = 4L,intercept=T,Boundary.knots=knot_BD)
y_plot = BS_plot%*%b.mmse


#### Estimate the cognitive clock
## Define the common shape function
u0 = b.mmse
AGE0 = 60
align.fun.BS <- function(t,a,b,c,u=u0)
{
  t = (t)*exp(-b) + a + AGE0
  n = length(t)
  
  af = rep(0,n)
  I.in = as.logical((t < max(knot_BD))*(t > min(knot_BD)))
  if (sum(I.in)>=1){
    BS.in = bSpline(t[I.in],knots=knot_Int,degree = 4L,intercept=T,Boundary.knots=knot_BD)
    f.in = BS.in%*%u
    af[I.in] = f.in
  }
  
  delta.BD = max(knot_BD) - min(knot_BD)
  I.left = (t <= min(knot_BD))
  if (sum(I.left)>=1){
    tmin0 = min(knot_BD)
    tmin1 = min(knot_BD)+0.2*delta.BD
    BS.left = bSpline(c(tmin0,tmin1),knots=knot_Int,degree = 4L,intercept=T,Boundary.knots=knot_BD)
    fBS.left = BS.left%*%u
    f.left = 1*(31-t[I.left])^4*(31>t[I.left]) + 
      fBS.left[1] + (t[I.left]-tmin0)/(tmin1-tmin0)*(fBS.left[2]-fBS.left[1])
    af[I.left] = f.left
  }
  
  I.right = (t >= max(knot_BD))
  if (sum(I.right)>=1){
    tmax0 = max(knot_BD)
    tmax1 = max(knot_BD)-0.2*delta.BD
    BS.right = bSpline(c(tmax0,tmax1),knots=knot_Int,degree = 4L,intercept=T,Boundary.knots=knot_BD)
    fBS.right = BS.right%*%u
    f.right = 0*(max(knot_BD)-t[I.right])^4 + 
      fBS.right[1] + (t[I.right]-tmax0)/(tmax1-tmax0)*(fBS.right[2]-fBS.right[1])
    af[I.right] = f.right
  }
  return(af+c)
}
DataFit$t_AGE0 = DataFit$t - AGE0

## Estimate the shape invariance model
lcr.fit<- nlme(y~align.fun.BS(t_AGE0,a,b,c),
               data = DataFit, 
               fixed = c ~ 1,
               random = a+b ~ 1|projid,
               start = c(c=0),
               na.action = na.omit,
               control = nlmeControl(maxIter=100,niterEM=100,
                                     msMaxIter=400,msMaxEval=2000,
                                     tolerance=1e-4,
                                     returnObject=TRUE,
                                     eval.max=1000, opt = 'nlm'))
summary(lcr.fit)
lcr.fit.fixed = lcr.fit$coefficients$fixed
lcr.fit.random = as.data.frame(lcr.fit$coefficients$random[[1]])
summary(lcr.fit.random)



####### BLUP estimation of the cognitive age
DataFit_BL = DataFit[DataFit$fu_year==0,]
abBLUP = LCR_FixC_BLUP(data=list(Data_Long = DataFit[!is.na(DataFit$y),], Data_Uni = DataFit_BL),
                         lcr.fit=lcr.fit)

DataFit_BL$time2ad = 0 #
DataFit_BL$time2mci_i = 0 #
DataFit_BL$age_bl = DataFit_BL$t

CogAge.Fit = CogAge_Est(data=list(Data_Long = DataFit, Data_Uni = DataFit_BL),
                        ab.est=abBLUP,AGE0=60)

DataFit$cogAge = CogAge.Fit$cogAge



convery_toAge <- function(ini, lcrfit, med = 90 ){
  # a function 
  # testing----
  # ini <- DataFit
  # lcrfit <- lcr.fit.random
  
  ini <- ini %>%
    group_by(projid) %>%
    mutate(Age0 =  first(t))
  
  ab = lcrfit
  ab$projid <- as.integer(rownames(lcrfit))
  info <- full_join(ini,ab, by = "projid")
  info$aprx_real_age <- (med - 60 - info$a)*exp(info$b) + 60
  info$state = d0$DX
  
  return(info)
}


# without considering the illness-------

  info <- convery_toAge(DataFit, lcr.fit.random, 90)
  info <- na.omit(info)
  
  ad <- info[info$state == "Dementia" & info$aprx_real_age >= 60& info$cogAge >= 87,] # for other, change into "AD"
  ad <- ad %>% 
    group_by(projid)%>%
    mutate(ra = first(aprx_real_age) , ta= first(t))
  
  ad <- distinct(data.frame(ra = ad$ra,ta = ad$ta, id = ad$projid))
  
  library(Metrics)
  # root mean square error
  rmse(ad$ra,ad$ta)
  # absolute mean error
  mean(abs(ad$ra - ad$ta))
  
  # predicted age's mean, sd, range
  summary(ad$ra)
  sd(ad$ra)
  summary(ad$ta)
  sd(ad$ta)
  
  plot(ad$ra , ad$ta, xlab = "Predicted age", ylab = "Observed age",xlim = c(60,100),ylim = c(60,100))
  abline(a = 0,b=1, col = "red")
  
  
  
  # considering for the illness --------
  
  ill <- data.frame(projid = d0$RID, disease = d0$MHPSYCH,t = d0$AGE)
  info <- left_join(info, ill, by = c("projid","t") )
  
  
  
  dis <- function(info, disage, ndage){
    # disage : people with disease's age mean 
    # ndage : without disease's age mean 
    
    dt <- matrix(nrow = 2,ncol = 2)
    info <- convery_toAge(DataFit, lcr.fit.random, disage)
    info <- left_join(info, ill, by = c("projid","t") )
    info <- na.omit(info)
    d <- info[info$state == "Dementia" & info$aprx_real_age >= 60& info$cogAge >= 87 & info$disease  != 0,] # for other, change into "AD"
    d <- d%>% 
      group_by(projid)%>%
      mutate(ra = first(aprx_real_age) , ta= first(t))
    
    d <- distinct(data.frame(ra = d$ra,ta = d$ta, id = d$projid))
    n1 <- dim(d)[1]
    dt[1,] <- c( rmse(d$ra,d$ta),mean(abs(d$ra - d$ta)))

     info <- convery_toAge(DataFit, lcr.fit.random, ndage)
     info <- left_join(info, ill, by = c("projid","t") )
     info <- na.omit(info)
     nd <- info[info$state == "Dementia" & info$aprx_real_age >= 60& info$cogAge >= 87 & info$disease  == 0,]
     nd <- nd%>% 
       group_by(projid)%>%
       mutate(ra = first(aprx_real_age) , ta= first(t))
     
     nd <- distinct(data.frame(ra = nd$ra,ta = nd$ta, id = nd$projid))
     n2 <- dim(nd)[1]
     dt[2,] <- c( rmse(nd$ra,nd$ta),mean(abs(nd$ra - nd$ta)))
     rownames(dt) <- c("disease","nondisease")   
     colnames(dt) <- c("rmse", "abse")
     
     newmae <- (dt[2,2]*n2+dt[1,2]*n1)/(n1+n2)
     newrmse <- ((n2*(dt[2,1])^2+n1*(dt[1,1])^2)/(n1+n2))^(0.5)
     
     print(c(newmae,newrmse))
                              
  }
  
  dis(info, 95.7, 89.9 )
  

  

##number count
  last_records <- DataFit %>%
    group_by(projid) %>%
    slice_head(n = 1)
  lastrisk <- subset(last_records,last_records$cogAge > 80|last_records$cogAge == 80)





##time-depemndent AUC
  DataFit$state <- d0$DX
  DataFit$censored[DataFit$state == "Dementia"] <- 1
  DataFit$censored[DataFit$state == "MCI"] <- 0
  DataFit$censored[DataFit$state == "CN"] <- 0
  library(timeROC)
  library(survival)
  ##For AD predict, chro v.s. cog
  ROC.bili.marginal<-timeROC(T=DataFit$fu_year,
                             delta=DataFit$censored,marker=DataFit$cogAge,
                             cause=1,weighting="marginal",
                             times=c(5,6,7,8,9),ROC=TRUE)
  ROC.bili.marginal
  
  plot(ROC.bili.marginal,time=8)        
  plot(ROC.bili.marginal,time=7,add=TRUE,col="blue") 
  plot(ROC.bili.marginal,time=6,add=TRUE,col="grey50") 
  legend("bottomright",c("Y-8","Y-10","Y-12"),col=c("red","blue","grey50"),lty=1,lwd=2)
  
  
  ##For AD predict, chro v.s. cog
  ROC.bili.marginal11<-timeROC(T=DataFit$fu_year,
                               delta=DataFit$censored,marker=DataFit$cogAge,
                               cause=1,weighting="marginal",
                               times=c(1,3,5,7,9),ROC=TRUE,iid=TRUE)

 
  plotAUCcurve(ROC.bili.marginal11,conf.int=TRUE,col="red")
  legend("topleft",c("AUC of cognitive age"),col=c("red"),lty=1,lwd=2)










