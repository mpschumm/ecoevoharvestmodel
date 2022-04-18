# Get mean genotypes at each time point, after mortality, growth, and reproduction have occurred
# all fish could theoretically mature in a year, so just exclude age="0"
mean_parent_genotypes <- vector(mode="numeric", length=runtime)
for (i in 1:runtime) {
  mean_parent_genotypes[i] <- (age_dist_1[1,i]*genotype1l_bar + age_dist_2[1,i]*genotype2l_bar)/(age_dist_1[1,i] + age_dist_2[1,i]) 
}

# Get mean genotypes of offspring
mean_offspring_genotypes <- vector(mode="numeric", length=runtime)
for (i in 1:runtime) {
  mean_offspring_genotypes[i] <- (age_dist_1[1,i]*genotype1l_bar + age_dist_2[1,i]*genotype2l_bar)/(age_dist_1[1,i] + age_dist_2[1,i])
}

# Find where to start vectors so they sync up with time points where reproduction happens
round(add_genotypes/timescale)
age_dist_1[1,368]

mean_parent_genotypes<-mean_parent_genotypes[seq(368,length(mean_parent_genotypes)-4, by=4)]
mean_offspring_genotypes <- mean_offspring_genotypes[seq(372,length(mean_offspring_genotypes), by=4)]

lines(mean_offspring_genotypes-mean_parent_genotypes,col="blue", lty=1, lwd=1.5)
lines(mean_offspring_genotypes-mean_parent_genotypes, col="red", lty=1, lwd=1.5)

# Make a plot for each of the four models (or 8, if herring)

p_selection_r_nh_upper <- vector(mode = "numeric", length=57)
p_selection_r_nh <- vector(mode = "numeric", length=57)
p_selection_r_nh_lower <- vector(mode = "numeric", length=57)
for (i in 1:57) {
  p_selection_r_nh_upper[i] <- std.error(na.omit(selection_r_nh[seq(i, 568+i, by=57)])) + mean(na.omit(selection_r_nh[seq(i, 568+i, by=57)]))
  p_selection_r_nh[i] <- mean(na.omit(selection_r_nh[seq(i, 568+i, by=57)]))
  p_selection_r_nh_lower[i] <- (-1)*std.error(na.omit(selection_r_nh[seq(i, 568+i, by=57)])) + mean(na.omit(selection_r_nh[seq(i, 568+i, by=57)]))
}

plot(p_selection_r_nh, ylab="Selection differential for maturation",xlab="Year", type="l", lty=2, col="blue", ylim=c(-0.002, 0.003))
lines(p_selection_r_nh_upper, type="l", lty=2, lwd=0.3, col="blue")
lines(p_selection_r_nh_lower, type="l", lty=2, lwd=0.3, col="blue")

p_selection_r_h_upper <- vector(mode = "numeric", length=57)
p_selection_r_h <- vector(mode = "numeric", length=57)
p_selection_r_h_lower <- vector(mode = "numeric", length=57)
for (i in 1:57) {
  p_selection_r_h_upper[i] <- std.error(na.omit(selection_r_h[seq(i, 568+i, by=57)])) + mean(na.omit(selection_r_h[seq(i, 568+i, by=57)]))
  p_selection_r_h[i] <- mean(na.omit(selection_r_h[seq(i, 568+i, by=57)]))
  p_selection_r_h_lower[i] <- (-1)*std.error(na.omit(selection_r_h[seq(i, 568+i, by=57)])) + mean(na.omit(selection_r_h[seq(i, 568+i, by=57)]))
}

lines(p_selection_r_h,  lty=2, col="red")
lines(p_selection_r_h_upper, type="l", lty=2, lwd=0.3, col="red")
lines(p_selection_r_h_lower, type="l", lty=2, lwd=0.3, col="red")
legend("topright",legend=c("No harvest, no var", "Harvest, no var","No harvest, var", "Harvest, var"), lty=c(1,1,2,2), col=c("blue","red","blue","red"), cex=0.5)

x = c(1,2,3,4)
sdev = c(0,0,0.001088462, 0.0005739395)
avg = c(0.1077446, 0.08439146, 0.1010704, 0.08395927)
plot(x, avg,
     ylim=range(c(avg-sdev, avg+sdev)),
     pch=19, xlab="", main="Value of maturation trait\nafter 50 years",
     ylab="Mean value of trait"
)
# hack: we draw arrows but with very special "arrowheads"
arrows(x, avg-sdev, x, avg+sdev, length=0.05, angle=90, code=3)
