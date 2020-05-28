library(PropCIs)
library(magrittr)

baseline_rate <-0.3
num_sim <- 1000




sizes_to_use <- seq(from=30, to=200, by=10)
res_mat <- matrix(NA, nrow=length(sizes_to_use), ncol=4 )

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
