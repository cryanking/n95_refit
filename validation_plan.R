library(PropCIs)
library(magrittr)

baseline_rate <-0.3




sizes_to_use <- seq(from=30, to=200, by=10)
res_mat <- matrix(NA, nrow=length(sizes_to_use), ncol=4 )
num_sim <- 1000
for( local_size in seq_along(sizes_to_use)) {
local_x <-rbinom(n=num_sim, prob=baseline_rate, size=sizes_to_use[local_size]) %>% scoreci(x=., n= sizes_to_use[local_size], conf.level=.95) %>% extract2("conf.int") %>% matrix(ncol=2)
local_x[,1] <- local_x[,1] *-1
local_x <- rowSums(local_x)
res_mat[local_size,] <- c(sizes_to_use[local_size], quantile(local_x, prob=c(.25, .5, .75) ) )
}

pdf("validation_sample_size.pdf")
plot( res_mat[,1], res_mat[,3], xlab="sample size", ylab="expected CI width", lwd=2, type="l")
lines( res_mat[,1], res_mat[,2], lty=2)
lines( res_mat[,1], res_mat[,4], lty=2)

abline(h=c(.15, .2), col='red')


res_mat[which.max(res_mat[,3] <.20 ),1]
res_mat[which.max(res_mat[,3] <.15 ),1]
dev.off()


num_sim <- 200
library(foreach)
library('doParallel')

registerDoParallel(cores=8)
library(flexsurv)
library(dplyr)
res_mat_surv <- matrix(NA, nrow=length(sizes_to_use), ncol=4 )
daily_fail <- c(.1, .1, .1, .1 , .6)

for( local_size in seq_along(sizes_to_use)) {


  local_res <- foreach(local_index = seq_len(num_sim) , .combine='c', .inorder=FALSE) %dopar% {
    local_x <-sample.int(n=5,size=sizes_to_use[local_size], prob=daily_fail, replace=TRUE) 
    local_surv <- data.frame(obs = (local_x)) %>% mutate(left=obs %>% subtract(1) %>% pmax(0.25), right=if_else(obs==5L, 20, as.numeric(obs)) ) %>% mutate(left = if_else(obs==5L, 4, left))
    fs_1 <- flexsurvreg(Surv(left, right, type="interval2")~1, data=local_surv, dist="gamma")
    local_out <- summary(fs_1, type="survival", t=3, se=TRUE, ci=TRUE)[[1]]
    diff(local_out[,3:4, drop=TRUE] %>% unlist)
  }

  res_mat_surv[local_size,] <- c(sizes_to_use[local_size], quantile(local_res, prob=c(.25, .5, .75) ) )

}


pdf("validation_sample_survival_size.pdf")
plot( res_mat_surv[,1], res_mat_surv[,3], xlab="sample size", ylab="expected CI width", lwd=2, type="l", col='red')
lines( res_mat_surv[,1], res_mat_surv[,2], lty=2, col='red')
lines( res_mat_surv[,1], res_mat_surv[,4], lty=2, col='red')


lines( res_mat[,1], res_mat[,3])
lines( res_mat[,1], res_mat[,2], lty=2)
lines( res_mat[,1], res_mat[,4], lty=2)

dev.off()




library("km.ci")
res_mat_km <- matrix(NA, nrow=length(sizes_to_use), ncol=4 )
daily_fail <- c(.1, .1, .1, .1 , .6)

for( local_size in seq_along(sizes_to_use)) {


  local_res <- foreach(local_index = seq_len(num_sim) , .combine='c', .inorder=FALSE) %dopar% {
  
    local_x <-sample.int(n=5,size=sizes_to_use[local_size], prob=daily_fail, replace=TRUE) 
    local_surv <- data.frame(obs = (local_x)) %>% mutate(left=obs %>% subtract(1) %>% pmax(0.25), right=if_else(obs==5L, 20, as.numeric(obs)) ) %>% mutate(left = if_else(obs==5L, 4, left))
    fs_1 <- km.ci(survfit(Surv(time=obs,event=obs<5, type="right")~1, data=local_surv))
    local_out <- summary(fs_1, type="survival", t=3)
  
    local_out$upper - local_out$lower
  }

  res_mat_km[local_size,] <- c(sizes_to_use[local_size], quantile(local_res, prob=c(.25, .5, .75) ) )

}




pdf("validation_sample_survival_size.pdf")
plot( res_mat_surv[,1], res_mat_surv[,3], xlab="sample size", ylab="expected CI width", lwd=2, type="l", col='red')
lines( res_mat_surv[,1], res_mat_surv[,2], lty=2, col='red')
lines( res_mat_surv[,1], res_mat_surv[,4], lty=2, col='red')


lines( res_mat[,1], res_mat[,3])
lines( res_mat[,1], res_mat[,2], lty=2)
lines( res_mat[,1], res_mat[,4], lty=2)


lines( res_mat_km[,1], res_mat_km[,3] , col='green')
lines( res_mat_km[,1], res_mat_km[,2], lty=2, col='green')
lines( res_mat_km[,1], res_mat_km[,4], lty=2, col='green')


dev.off()





