#
# Some binomial stats to determine
# the copy number of DAXX/ATRX mutations
#

##Patient 3 DAXX
samples <- c("7370","7372","7379","7382","7385","7386")
mut_counts <- c(22,34,38,53,24,33)
depth <- c(50,54,45,74,72,55)
purity <- c(0.52,0.62,0.85,0.71,0.31,0.7)
tot_cn <- c(rep(2,6))
mle <- c()
conf_lower <- c()
conf_upper <- c()

p_after_loh <- c()
p_before_loh <- c()
for (i in 1:6) {
  test_after <- binom.test(mut_counts[i], depth[i], purity[i]/2)
  p_after_loh[i] <- test_after$p.value
  test_before <- binom.test(mut_counts[i], depth[i], purity[i])
  p_before_loh[i] <- test_before$p.value
  conf_lower[i] <- 2*test_before$conf.int[1]/purity[i]
  conf_upper[i] <- 2*test_before$conf.int[2]/purity[i]
  mle[i]=2*(mut_counts[i]/depth[i])/purity[i]
}
results <- data.frame(samples, mle, conf_lower, conf_upper, p_after_loh, p_before_loh)
colnames(results) <- c("Sample", "CN mutation", "95% lower bound", "95% upper bound", "p-value After LOH", "p-value before LOH")
#write.table(results, file="PanNET3_DAXX.txt", sep="\t", row.names=F)

##Patient 3 ATRX
samples <- c("5621","765-3","765-5")
mut_counts <- c(34,39,36)
depth <- c(37,44,44)
purity <- c(0.84, 0.72, 0.81)
tot_cn <- c(rep(2,3))
mle <- c()
conf_lower <- c()
conf_upper <- c()

p_after_loh <- c()
p_before_loh <- c()
for (i in 1:3) {
  test_after <- binom.test(mut_counts[i], depth[i], purity[i]*1/((1-purity[i])*1+tot_cn[i]*purity[i]))
  p_after_loh[i] <- test_after$p.value
  test_before <- binom.test(mut_counts[i], depth[i], purity[i]*2/((1-purity[i])*1+tot_cn[i]*purity[i]))
  p_before_loh[i] <- test_before$p.value
  conf_lower[i] <- (1*(1-purity[i])+tot_cn[i]*purity[i])*test_before$conf.int[1]/purity[i]
  conf_upper[i] <- (1*(1-purity[i])+tot_cn[i]*purity[i])*test_before$conf.int[2]/purity[i]
  mle[i]=(1*(1-purity[i])+tot_cn[i]*purity[i])*(mut_counts[i]/depth[i])/purity[i]
}
results <- data.frame(samples, mle, conf_lower, conf_upper, p_after_loh, p_before_loh)
colnames(results) <- c("Sample", "CN mutation", "95% lower bound", "95% upper bound", "p-value After LOH", "p-value before LOH")
#write.table(results, file="PanNET3_ATRX.txt", sep="\t", row.names=F)


##Patient 4 DAXX
samples <- c("2580","2584","2586","2597","2600","320-10","320-lever","320","6680")
mut_counts <- c(34,42,41,54,36,36,41,37,34)
depth <- c(49,60,60,67,44,42,54,53,63)
purity <- c(0.68,0.76,0.73,0.81,0.86,0.79,0.75,0.74,0.58)
tot_cn <- c(rep(2,9))
mle <- c()
conf_lower <- c()
conf_upper <- c()

p_after_loh <- c()
p_before_loh <- c()
for (i in 1:9) {
  test_after <- binom.test(mut_counts[i], depth[i], purity[i]*1/((1-purity[i])*2+tot_cn[i]*purity[i]))
  p_after_loh[i] <- test_after$p.value
  test_before <- binom.test(mut_counts[i], depth[i], purity[i]*2/((1-purity[i])*2+tot_cn[i]*purity[i]))
  p_before_loh[i] <- test_before$p.value
  conf_lower[i] <- (2*(1-purity[i])+tot_cn[i]*purity[i])*test_before$conf.int[1]/purity[i]
  conf_upper[i] <- (2*(1-purity[i])+tot_cn[i]*purity[i])*test_before$conf.int[2]/purity[i]
  mle[i]=(2*(1-purity[i])+tot_cn[i]*purity[i])*(mut_counts[i]/depth[i])/purity[i]
}
results <- data.frame(samples, mle, conf_lower, conf_upper, p_after_loh, p_before_loh)
colnames(results) <- c("Sample", "CN mutation", "95% lower bound", "95% upper bound", "p-value After LOH", "p-value before LOH")
#write.table(results, file="PanNET4_DAXX.txt", sep="\t", row.names=F)
