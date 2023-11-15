library(cowplot)
library(ggsci)
library(ggpubr)

#Cell type 1
X_ct1 <- rgamma(1000, shape = shape_ct1, rate = rate_ct1)
den_ct1 <- seq(min(X_ct1), max(X_ct1), length.out = 1000)
den_plot_ct1 <- dgamma(den_ct1, shape = shape_ct1, rate = rate_ct1)
h1 <- hist(X_ct1, breaks=30, probability=TRUE, plot= F, ylim = c(0,1.1*max(den_ct1)))
lines(density(den_ct1), col="#4DBBD5FF")
lines(density(X_ct1), lty=2)

#Cell type 2
X_ct2 <- rgamma(1000, shape = shape_ct2, rate = rate_ct2)
den_ct2 <- seq(min(X_ct2), max(X_ct2), length.out = 1000)
den_plot_ct2 <- dgamma(den_ct2, shape = shape_ct2, rate = rate_ct2)
h2 <- hist(X_ct2, breaks=30, probability=TRUE, plot = F, ylim = c(0,1.1*max(den_plot_ct2))) + lines(den_ct2, den_plot_ct2, col="#E64B35FF")
lines(density(X_ct2), lty=2)

plot(h1, col = "#4DBBD599", ylim= c(0,200), main = "histogram of gamma distributions",
     xlab = "")
plot(h2, col = "#E64B3599", add = T, ylim = c(0, 200))
