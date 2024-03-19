#######################################################################
### Cognitive Clock Estimation
#######################################################################
library(nlme)
library(lme4)
library(splines)
library(splines2)
library(gmodels)
library(pracma)
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

# 
DataFit <- data.frame(y = BLdataset$Neuropsych.MMSE.total,t = BLdataset$Age,projid = BLdataset$AIBL.Id,
                      fu_year = BLdataset$Age_diff)

####################


#### Estimate age-specific MMSE mean trajectory by B-spline with equal-spaced interior knots
## Number of interior knots
K = 4
## Define the B-spline basis functions
knot_BD = c(55,105)
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
print(ymean)
## Ensure monotonicity
lme.fit.BS.d = diff(lme.fit.BS)
lme.fit.BS.d = lme.fit.BS.d*I(lme.fit.BS.d<=0) - 0.0001
pBS.hf = floor(pBS/2+1)
b.mmse = lme.fit.BS
b.mmse[pBS.hf:pBS] = cumsum(c(lme.fit.BS[pBS.hf],lme.fit.BS.d[pBS.hf:(pBS-1)]))
b.mmse[pBS.hf:1] = cumsum(c(lme.fit.BS[pBS.hf],-lme.fit.BS.d[(pBS.hf-1):1]))

# Plot the longitudinal MMSE trajectories and the estimated mean trajectory
t_plot = 60:100
BS_plot = bSpline(t_plot,knots=knot_Int,degree = 4L,intercept=T,Boundary.knots=knot_BD)
y_plot = BS_plot%*%b.mmse

dev.new()
mplot(x = t, y = y, id = projid, data = DataFit, col = projid, las = 1,
      ylab = 'MMSE',xlab = 'Chronological age',
      main = 'MMSE trajectries',
      xlim=c(60,100),ylim = c(0,30))
lines(t_plot,y_plot,lwd=4)


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
DataFit$Statu = BLdataset$Neuropsych.Simple.Classification



mplot(x = cogAge, y = y, id = projid, data = DataFit, col = projid, las = 1,
      ylab = 'MMSE',xlab = 'Cognitive Age',ylim = c(0,30),xlim=c(65,105),
      main = 'Estimated Cognitive Clock')
lines(65:105,
      align.fun.BS((65:105)-AGE0,0,0,lcr.fit.fixed['c']),lwd=3)




  

# plotting each state ----- 

dt <- DataFit[,c(1:3,7:8)]
  
hc <- dt %>% group_by(projid) %>% filter(first(Statu) == "HC" &last(Statu) == "HC")
mci <- dt %>% group_by(projid) %>% filter(first(Statu) == "HC" &last(Statu) == "MCI")
ad <- dt %>% group_by(projid) %>% filter(first(Statu) == "HC" &last(Statu) == "AD")

par(mfrow = c(1,1))

num_colors <- length(unique(hc$projid))
colors <- rainbow(num_colors)



## Plot with rainbow colors based on projid, with grid turned off and thinner lines ----

### usder for consstructing the abline---
x_values <- c(65:105)
y_values <- align.fun.BS((65:105) - AGE0, 0, 0, lcr.fit.fixed['c'])

# ggplot ------
library(ggplot2)

ggplot() +
  geom_line(data = hc, aes(x = cogAge, y = y, group = projid, color = projid)) + 
  scale_colour_gradientn(colours = rainbow(100)) +
  geom_line(data = data.frame(x = x_values, y = y_values), aes(x = x, y = y), color = "black", size = 1, show.legend = FALSE) +
  coord_cartesian(xlim = c(60, 100), ylim = c(0, 30)) +
  labs(x = "Cogage", y = "MMSE") 


ggplot() +
  geom_line(data = mci, aes(x = cogAge, y = y, group = projid, color = projid)) + 
  scale_colour_gradientn(colours = rainbow(100)) +
  geom_line(data = data.frame(x = x_values, y = y_values), aes(x = x, y = y), color = "black", size = 1, show.legend = FALSE) +
  coord_cartesian(xlim = c(60, 100), ylim = c(0, 30)) +
  labs(x = "Cogage", y = "MMSE")

  
ggplot() +
  geom_line(data = ad, aes(x = cogAge, y = y, group = projid, color = projid)) + 
  scale_colour_gradientn(colours = rainbow(100)) +
  geom_line(data = data.frame(x = x_values, y = y_values), aes(x = x, y = y), color = "black", size = 1, show.legend = FALSE) +
  coord_cartesian(xlim = c(60, 100), ylim = c(0, 30)) +
  labs(x = "Cogage", y = "MMSE")


##chro v.s. cog for individual
x <- c(88,89,91,92)
y <- c(11,9,6,0)
x1 <- c(98.16282,98.75861,99.05650,99.65228)
plot(x,y,col = "red",xlab = "chronological age", ylab = "MMSE",xlim=c(65,100),ylim = c(0,30))
lines(t_plot,y_plot,lwd=2,col = "blue")
lines(x,y,lwd=1, col = "red")
# Add legend
legend("topright", legend = c("observed for AIBL id = 7"), col = c("red"), 
       pch = c(16))


x <- c(89,91,92,94,95,97)
y <- c(29,28,28,27,30,24)
x1 <- c(74.64693,76.30862,77.13947,78.80116,79.63201,81.29370)
plot(x,y,col = "red",xlab = "chronological age", ylab = "MMSE",xlim=c(65,100),ylim = c(0,30))
lines(t_plot,y_plot,lwd=2,col = "blue")
lines(x,y,lwd=1, col = "red")
##lines(65:105,
  ##    align.fun.BS((65:105)-AGE0,0,0,lcr.fit.fixed['c']),lwd=2,col = "blue")
# Add legend
legend("bottomleft", legend = c("observed for AIBL id = 11"), col = c("red"), 
       pch = c(16))









##Cox PH model
##cox_reg1 <- coxph(Surv(cog,censored) ~ Gene, data=subsuvrdata)
##summary(cox_reg1)
##cox_reg2 <- coxph(Surv(cog,censored) ~ sex, data=suvrdata)
##summary(cox_reg2)
##cox_reg3 <- coxph(Surv(cog,censored) ~ sex+Gene, data=subsuvrdata)
##summary(cox_reg3)

#AFT model
#library(flexsurv)

#Fit a AFT survival model
#aft_fit <- flexsurvreg(Surv(cog, censored) ~ 1, data = suvrdata, dist = "gamma")
#ggsurvplot(aft_fit, data=suvrdata,xlim=c(75,100), risk.table = TRUE, conf.int = TRUE,
#legend.title = "Survival Curve",ggtheme=theme_minimal())

# Calculate Harrell's C-index











