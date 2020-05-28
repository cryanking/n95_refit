library(magrittr)
library(tidyverse)
library(survival)
library("interval")
library(cgam)
library(icenReg)
library("km.ci")
library(flexsurv)
library(tableone)

missing_packages <- setdiff(c("caret", "tableone", "rmarkdown", "knitr", "bookdown") , rownames(installed.packages() ) ) 

if(length(missing_packages) > 0) {
stop(paste("install" , paste0(missing_packages, collapse =" , ") ))
}
##############
#### data load and procesing
##############
prelim.df <- read_csv(file="Final Pilot Data.csv")

## this variable controls the minimum amount of time a mask is assumed to function; the analysis assumes ALL masks functioned before this time since all users had passed fit tests and masks undergo factory qc. it is expressed in days
immortal_time <- 0.5

## drop the redundant column
names(prelim.df) %<>% make.names
prelim.df %<>% rename(fit.fail = Qualitative.Fit)
prelim.df %<>% mutate(fit.fail = fit.fail == "Fail")
prelim.df %<>% mutate_if( .predicate=is.character, factor)
prelim.df %<>% mutate_at(vars(one_of("fit.fail","Days.Worn", "Sterilizations", "Donned")), as.integer  )
prelim.df %<>% mutate(left=if_else(fit.fail==1, immortal_time, as.numeric(Days.Worn)  ) , right=if_else(fit.fail==1, as.numeric(Days.Worn), as.numeric(max(Days.Worn )*3 )  )   )

##############
### cross-sectional analysis
##############

## this is my preferred estimate
cross_fit <- cgam(fit.fail ~ s.incr(Days.Worn) , family=binomial(), data=prelim.df )
dummy_data <- data.frame(Days.Worn = 4:20)
dummy_data <- cbind(dummy_data, predict(cross_fit, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )
dummy_data %<>% filter(Days.Worn %in% c(5,10,15) )


## supplement plot comparing methods
## black = non-parametric monotone decreasing (no smoothness constrain), identical to the nonparametric survival estimator
## red = smooth monotone decreasing (spline basis, constrained)
## blue = smooth no monotone requirement (spline bases)
## conclusion: extremely similar, observed fit failure pretty flat across time
cross_fit <- cgam(fit.fail ~ incr(Days.Worn) , family=binomial(), data=prelim.df )
dummy_data <- data.frame(Days.Worn = 4:20)
dummy_data <- cbind(dummy_data, predict(cross_fit, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )
plot(dummy_data$Days.Worn, dummy_data$fit , xlab="Days Used", ylab="Failure probability", ylim=c(0,1), xlim=c(4,20), type="b", pch=19)
lines(dummy_data$Days.Worn, dummy_data$lower, lty=2)
lines(dummy_data$Days.Worn, dummy_data$upper, lty=2)

cross_fit <- cgam(fit.fail ~ s.incr(Days.Worn) , family=binomial(), data=prelim.df )
dummy_data <- data.frame(Days.Worn = 4:20)
dummy_data <- cbind(dummy_data, predict(cross_fit, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )
points(dummy_data$Days.Worn, dummy_data$fit , col='red', type="b", pch=19)
lines(dummy_data$Days.Worn, dummy_data$lower, col='red', lty=2)
lines(dummy_data$Days.Worn, dummy_data$upper, col='red', lty=2)


cross_fit <- glm(fit.fail ~ ns(Days.Worn) , family=binomial(), data=prelim.df )
dummy_data <- data.frame(Days.Worn = 4:20)
dummy_data <- cbind(dummy_data, predict(cross_fit, newdata=dummy_data, se.fit=TRUE)[c("fit", "se.fit")] %>% as.data.frame )

points(dummy_data$Days.Worn, cross_fit$family$linkinv(dummy_data$fit) , col='blue', type="b", pch=19)
dummy_data %<>% mutate( lower=cross_fit$family$linkinv(fit - 1.96*se.fit), upper=cross_fit$family$linkinv(fit+1.96*se.fit))
lines(dummy_data$Days.Worn, dummy_data$lower, col='blue', lty=2)
lines(dummy_data$Days.Worn, dummy_data$upper, col='blue', lty=2)




## supplement plot comparing methods
## black = nonparametric survival estimator with bootstrap (interval package)
## red = parametric survival (generalized gamma)
## blue = monotone smooth spline
## conclusion: survival model forces the early survival upwards
npm_out <- icfit( Surv( left, right, type = "interval2") ~ 1, data = prelim.df, conf.int=TRUE, control = icfitControl(B=2000,seed=1234))
fs_1 <- flexsurvreg(Surv(left, right, type="interval2")~1, data=prelim.df, dist="gengamma")
cross_fit <- cgam(I(1L-fit.fail) ~ s.decr(Days.Worn) , family=binomial(), data=prelim.df )
# fs_2 <- flexsurvspline(Surv(left, right, type="interval2")~1, data=prelim.df, timescale="log", k=1L) ## opaque error



npm_out %>% plot(xlim=c(4,20))
lines(fs_1)
dummy_data <- data.frame(Days.Worn = 4:20)
dummy_data <- cbind(dummy_data, predict(cross_fit, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )
lines(dummy_data$Days.Worn, dummy_data$fit , col='blue', lwd=2)
lines(dummy_data$Days.Worn, dummy_data$lower, col='blue', lty=2)
lines(dummy_data$Days.Worn, dummy_data$upper, col='blue', lty=2)

##ICsurv method - no direct method for se of output (though I think I could get it), is similar
if(FALSE) {
exp_res<-fast.PH.ICsurv.EM(d1=rep(0L, nrow(prelim.df)), d3=rep(0L, nrow(prelim.df)), d2=rep(1L, nrow(prelim.df)), Li=prelim.df$left, Ri=prelim.df$right, Xp=matrix(runif(nrow(prelim.df)), ncol=1) , n.int=3, order=3, g0=rep(1, 3+3) , b0=rep(0,1L), tol=0.001, equal=TRUE, t.seq=4:20)

points(dummy_data$Days.Worn, 1.-exp_res$hz, col="green", type="b", pch=19)
}

## same calculation, different package to check
if(FALSE) {
alt_np <- ic_np(cbind(left, right)~0, data=prelim.df )
getSCurves(alt_np )
}



## supplement plot: comparison of parametric survival families
## conclusion : very similar , generalized gamma probably the best fitting

## to show added uncertainty from the distribution choice
npm_out %>% plot(xlim=c(4,20))
lines(fs_1, type = "survival", col="red")

fs2 <- flexsurvreg(Surv(left, right, type="interval2")~1, data = prelim.df, dist="gompertz")
lines(fs2, type = "survival", col="blue")

fs2 <- flexsurvreg(Surv(left, right, type="interval2")~1, data = prelim.df, dist="weibull")
lines(fs2, type = "survival", col="green")



## supplement plot: comparison of interval censored and "traditional" survival analysis
## conclusion: traditional survival badly biased (as expected)
npm_out %>% plot(xlim=c(4,20))
lines(km.ci(survfit(Surv(time=Days.Worn, event=fit.fail)~1, data=prelim.df) ), col='red')




## analytical claims:
## 1: # steralizations doesn't matter conditional on number of days
## survival semi-para: 0.17 +- 0.27 log scale
## cross-sectional semi-para: .34 +- .34
np_reg <- ic_sp(cbind(left, right)~Sterilizations, data=prelim.df, bs_samples = 1000  )
np_reg %>% summary

cross_fit <- glm(fit.fail ~ ns(Days.Worn) , family=binomial(), data=prelim.df )
cross_fit2 <- glm(fit.fail ~ ns(Days.Worn) +Sterilizations  , family=binomial(), data=prelim.df )
anova(cross_fit2, cross_fit, test="LRT")
cross_fit2 %>% summary
glm(fit.fail ~ Sterilizations , family=binomial(), data=prelim.df ) %>% summary


## 2: number of times donned doesn't matter
##  survival semi-para: -0.005 +- 0.23 log scale
##  cross sectional: -0.0007505  0.0243780
np_reg2a <- ic_sp(cbind(left, right)~Donned, data=prelim.df, bs_samples = 1000  )
np_reg2a %>% summary

cross_fit <- glm(fit.fail ~ ns(Days.Worn) , family=binomial(), data=prelim.df )
cross_fit2 <- glm(fit.fail ~ ns(Days.Worn) +Donned  , family=binomial(), data=prelim.df )
anova(cross_fit2, cross_fit, test="LRT")
cross_fit2 %>% summary
glm(fit.fail ~ Donned , family=binomial(), data=prelim.df ) %>% summary

## 2b per day
## 0.01764 0.1975
## -0.05921    0.22283
np_reg2b <- ic_sp(cbind(left, right)~I(Donned/Days.Worn), data=prelim.df, bs_samples = 1000  )
np_reg2b %>% summary

cross_fit <- glm(fit.fail ~ ns(Days.Worn) , family=binomial(), data=prelim.df )
cross_fit2 <- glm(fit.fail ~ ns(Days.Worn) +I(Donned/Days.Worn)  , family=binomial(), data=prelim.df )
anova(cross_fit2, cross_fit, test="LRT")
cross_fit2 %>% summary
glm(fit.fail ~ I(Donned/Days.Worn) , family=binomial(), data=prelim.df ) %>% summary

## 3: assessment of fit probably does matter
##  -1.071  2.412 decent sized effect!
##  -1.8021     0.8696 p 0.02249 by lrt
##  about the same in marginal model -1.8137     0.8574    0.0344 * (p=.033 by fet)

np_reg3 <- ic_sp(cbind(left, right)~Fits.well, data=prelim.df, bs_samples = 1000  )
np_reg3 %>% summary

cross_fit <- glm(fit.fail ~ ns(Days.Worn) , family=binomial(), data=prelim.df )
cross_fit2 <- glm(fit.fail ~ ns(Days.Worn) +Fits.well  , family=binomial(), data=prelim.df )
anova(cross_fit2, cross_fit, test="LRT")
cross_fit2 %>% summary
glm(fit.fail ~ Fits.well , family=binomial(), data=prelim.df ) %>% summary
prelim.df %>% group_by(Fits.well) %>% summarize(mean(fit.fail ), n())
prelim.df %>% group_by(fit.fail) %>% summarize(mean(as.integer(Fits.well )-1), n())
fisher.test(prelim.df$Fits.well, prelim.df$fit.fail)

## also available as a test: the method matters, the more conservative one does not reject, the more analytical one is marginal
## p=.052, .145
np_test3 <- ictest(Surv( left, right, type = "interval2") ~ Fits.well, data = prelim.df)
np_test3 <- ictest(Surv( left, right, type = "interval2") ~ Fits.well, data = prelim.df, method="wsr.mc")


#         Asymptotic Logrank two-sample test (permutation form), Sun's scores
# 
# data:  Surv(left, right, type = "interval2") by Fits.well 
# Z = 1.9397, p-value = 0.05242
# alternative hypothesis: survival distributions not equal
# 
#                n Score Statistic*
# Fits.well=No  10         2.977042
# Fits.well=Yes 38        -2.977042
# * like Obs-Exp, positive implies earlier failures than expected
# p-value = 0.1452
# p-value estimated from  Monte Carlo replications
# and 999 permutation resamples

## plot the difference
## no obvious change in time, both are pretty flat. it's just an offset
# alt_np3 <- icfit( Surv( left, right, type = "interval2") ~ Fits.well, data = prelim.df, conf.int=TRUE, control = icfitControl(B=2000,seed=1234))
# alt_np3 %>% plot(xlim=c(4,20)) # this plot is hard to control

cross_fit2 <- cgam(I(1L-fit.fail) ~ s.decr(Days.Worn)*Fits.well  , family=binomial(), data=prelim.df )

npm_out %>% plot(xlim=c(4,20))

dummy_data <- data.frame(Days.Worn = 4:20 , Fits.well =prelim.df$Fits.well %>% unique %>% magrittr::extract(1) )
dummy_data <- cbind(dummy_data, predict(cross_fit2, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )
lines(dummy_data$Days.Worn, dummy_data$fit , col='red', lwd=2)
lines(dummy_data$Days.Worn, dummy_data$lower, col='red', lty=2)
lines(dummy_data$Days.Worn, dummy_data$upper, col='red', lty=2)

dummy_data <- data.frame(Days.Worn = 4:20 , Fits.well =prelim.df$Fits.well %>% unique %>% magrittr::extract(2) )
dummy_data <- cbind(dummy_data, predict(cross_fit2, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )
lines(dummy_data$Days.Worn, dummy_data$fit , col='blue', lwd=2)
lines(dummy_data$Days.Worn, dummy_data$lower, col='blue', lty=2)
lines(dummy_data$Days.Worn, dummy_data$upper, col='blue', lty=2)


# dummy_data <- cbind(dummy_data, predict(cross_fit2, newdata=dummy_data, se.fit=TRUE)[c("fit", "se.fit")] %>% as.data.frame )
# lines(dummy_data$Days.Worn, cross_fit2$family$linkinv(dummy_data$fit) , col='blue')
# dummy_data %<>% mutate( lower=cross_fit$family$linkinv(fit - 1.96*se.fit), upper=cross_fit$family$linkinv(fit+1.96*se.fit))
# lines(dummy_data$Days.Worn, dummy_data$lower, col='blue', lty=2)
# lines(dummy_data$Days.Worn, dummy_data$upper, col='blue', lty=2)

# dummy_data <- data.frame(Days.Worn = 4:20 , Fits.well =prelim.df$Fits.well %>% unique %>% magrittr::extract(2) )
# dummy_data <- cbind(dummy_data, predict(cross_fit2, newdata=dummy_data, se.fit=TRUE)[c("fit", "se.fit")] %>% as.data.frame )
# lines(dummy_data$Days.Worn, cross_fit2$family$linkinv(dummy_data$fit) , col='red')
# dummy_data %<>% mutate( lower=cross_fit$family$linkinv(fit - 1.96*se.fit), upper=cross_fit$family$linkinv(fit+1.96*se.fit))
# lines(dummy_data$Days.Worn, dummy_data$lower, col='red', lty=2)
# lines(dummy_data$Days.Worn, dummy_data$upper, col='red', lty=2)
# 


#############
## plots by # donned
## supplemental: comparison of methods
prelim.df2 <-prelim.df %>% mutate(left=if_else(fit.fail==1, 1, as.numeric(Donned)  ) , right=if_else(fit.fail==1, as.numeric(Donned), as.numeric(max(Donned )*3 )  )   )

npm_out2 <- icfit( Surv( left, right, type = "interval2") ~ 1, data = prelim.df2, conf.int=TRUE, control = icfitControl(B=2000,seed=1234))

fs3 <- flexsurvreg(Surv(left, right, type="interval2")~1, dist="gengamma", data=prelim.df2)
npm_out2 %>% plot( xlim=c(6,50) )
lines(fs3, col="blue")

cross_fit3 <- cgam(I(1L-fit.fail) ~ s.incr(Donned) , family=binomial(), data=prelim.df )
dummy_data <- data.frame(Donned = 6:50)
dummy_data <- cbind(dummy_data, predict(cross_fit3, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )
lines(dummy_data$Donned, dummy_data$fit , col='red', lwd=2)
lines(dummy_data$Donned, dummy_data$lower, col='red', lty=2)
lines(dummy_data$Donned, dummy_data$upper, col='red', lty=2)


save(file="n95_intermediates.Rdata", prelim.df, prelim.df2, npm_out2, npm_out, np_test3, np_reg3, np_reg2b, np_reg2a, np_reg)




prelim.df %>% group_by(fit.fail) %>% summarize(sum(as.integer(Fits.well )-1), n())

temp <- prelim.df %>% filter(fit.fail==1) %>% summarize(sum(as.integer(Fits.well )-1), n()) %>% unlist
prop.test(x=temp[1], n=temp[2]) %>% extract2("conf.int") %>% multiply_by(100) %>% round


temp <- prelim.df %>% filter(fit.fail==1) %>% mutate(Mask.quality=Mask.quality %in% c("Like New", "Good") ) %>% summarize(sum(as.integer(Mask.quality )), n()) %>% unlist
prop.test(x=temp[1], n=temp[2]) %>% extract2("conf.int") %>% multiply_by(100) %>% round


# PropCIs::scoreci(x=temp[1], n=temp[2], conf.level=0.95) %>% extract2("conf.int") %>% multiply_by(100) %>% round

with(prelim.df, cor(Donned, Days.Worn) )
with(prelim.df, summary(lm(Donned~ Days.Worn) ))


prelim.df %>% select(fit.fail, Fits.well ) %>% table
with(prelim.df, caret::confusionMatrix(reference=factor(fit.fail, labels=c("Yes", "No")), data=Fits.well))

prelim.df %>% mutate(Mask.quality=Mask.quality %in% c("Poor") ) %>% select(Mask.quality, fit.fail ) %>% table

with(prelim.df  %>% mutate(Mask.quality=factor(!(Mask.quality %in% c("Poor")) , labels=c("Yes", "No" ) )), caret::confusionMatrix(reference=factor(fit.fail, labels=c("Yes", "No")), data=Mask.quality))

# temp_t1 <- table1(~Sex + Title + respirator + Days.Worn + Sterilizations + Donned + I(Donned/Days.Worn ) +Fits.well + Mask.quality | factor(fit.fail, labels=c("Yes", "No"))  , data=prelim.df)

factor_vars<- c('Sex' , 'Title' , 'respirator', 'Fits.well', 'Mask.quality'  )

con_vars<- c( 'Days.Worn', 'Sterilizations', 'Donned',  'Uses.per.day' )
tab3 <- CreateTableOne(vars = c(con_vars,factor_vars) , strata = "fit.fail" , data = prelim.df%>%  mutate( Uses.per.day=Donned/Days.Worn ) , factorVars = factor_vars)

temp <- print(tab3, showAllLevels = TRUE, contDigits=1)

pdf(file="failure_prob.pdf")

cross_fit <- cgam(fit.fail ~ incr(Days.Worn) , family=binomial(), data=prelim.df )
dummy_data <- data.frame(Days.Worn = 4:20)
dummy_data <- cbind(dummy_data, predict(cross_fit, newData=dummy_data, interval="confidence")[c("fit", "lower", "upper")] %>% as.data.frame )


temp_hist <- hist(prelim.df  %>% select(Days.Worn) %>% unlist ,breaks = seq(from=3.5, to=20.5), plot=T, axes=F, xlab=NULL, ylab=NULL, xlim=c(3.5,20) , main=NULL, col= rgb(0.6,1.0,1.0,alpha=0.5))

par(new=T)

plot(dummy_data$Days.Worn, dummy_data$fit , xlab="Days worn", ylab="Failure probability", ylim=c(0,1), xlim=c(3.5,20), type="l", col="red", lwd=2)
lines(dummy_data$Days.Worn, dummy_data$lower, lty=2, col="red", lwd=2)
lines(dummy_data$Days.Worn, dummy_data$upper, lty=2, col="red", lwd=2)


mtext(side=4, text="Frequency of Day", line=1)
par(new=F)


dev.off()

prelim.df %>% mutate(broken_days = cut(Days.Worn, breaks=c(3, 5, 8, 12, 25)) ) %>% select(broken_days) %>% table

prelim.df %>% mutate(broken_days = cut(Days.Worn, breaks=c(3, 5, 8, 12, 25)) ) %>% group_by(broken_days) %>% summarize(n=n() , mean.failure=mean(fit.fail )%>% round(2), lcb=prop.test(x=sum(fit.fail), n=n())$conf.int[1]%>% round(2), ucb=prop.test(x=sum(fit.fail), n=n())$conf.int[2]%>% round(2) ) %>% {knitr::kable(x=.,format="html")}



with(prelim.df , scatter.smooth(Days.Worn, fit.fail ))


