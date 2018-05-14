library(MCMCpack)
set.seed(1984)
n = 300

true_theta_1 = 0
true_theta_2 = 4
true_sigsq_1 = 2
true_sigsq_2 = 1
true_rho = 0.3


x = c(rnorm(n * true_rho, true_theta_1, sqrt(true_sigsq_1)), rnorm(n * (1 - true_rho), true_theta_2, sqrt(true_sigsq_2)))

par(mfrow = c(1, 1))
hist(x, br = 80) 
#plot


#chains
S = 10000
theta1s = array(NA, S)
theta2s = array(NA, S)
sigsq1s = array(NA, S)
sigsq2s = array(NA, S)
rhos = array(NA, S)
Is = matrix(NA, nrow = n, ncol = S)
#start positions
theta1s[1] = mean(x)
theta2s[1] = mean(x)
sigsq1s[1] = var(x)
sigsq2s[1] = var(x)
rhos[1] = 0.5
Is[, 1] = 0.5 

for (t in 2 : S){
  theta1 = theta1s[t - 1]
  theta2 = theta2s[t - 1]
  sigsq1 = sigsq1s[t - 1]
  sigsq2 = sigsq2s[t - 1]
  rho = rhos[t - 1]
  I = Is[, t - 1]
  
  sum_I = sum(I)
  sum_1_min_I = n - sum_I
  theta1s[t] = rnorm(1, sum(I * x) / sum_I, sqrt(sigsq1 / sum_I))
  theta2s[t] = rnorm(1, sum((1 - I) * x) / sum_1_min_I, sqrt(sigsq2 / sum_1_min_I))
  sigsq1s[t] = rinvgamma(1, sum_I / 2, sum(I * (x - theta1s[t])^2) / 2)
  sigsq2s[t] = rinvgamma(1, sum_1_min_I / 2, sum((1 -I) * (x - theta2s[t])^2) / 2)
  
  for (i in 1 : n){#now draw the Is
    a = rho * dnorm(x[i], theta1s[t], sqrt(sigsq1s[t]))
    b = (1 - rho) * dnorm(x[i], theta2s[t], sqrt(sigsq2s[t]))
    Is[i, t] = rbinom(1, 1, a / (a + b))
    # cat("a =", a, "b = ", b, "p =", a / (a + b), "I =", Is[i, t], "\n")
  }
  rhos[t] = rbeta(1, 1 + sum_I, 1 + sum_1_min_I)
  # cat("t =", t, "sum_I = ", sum_I, "sum_1_min_I =", sum_1_min_I, "rhos[t] =", rhos[t], "\n")
}



###assess convergence
#plot
B = 70
###
par(mfrow = c(5, 1))
S0 = 150
plot(1 : S0, theta1s[1 : S0])
abline(h = mean(theta1s[B : S0]), col = "blue")
# abline(h = true_theta_1, col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, theta2s[1 : S0])
abline(h = mean(theta2s[B : S0]), col = "blue")
# abline(h = true_theta_2, col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, sigsq1s[1 : S0])
abline(h = mean(sigsq1s[B : S0]), col = "blue")
# abline(h = sqrt(true_sigsq_1), col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, sigsq2s[1 : S0])
abline(h = mean(sigsq2s[B : S0]), col = "blue")
# abline(h = sqrt(true_sigsq_2), col = "red")
# abline(v = B, col = "grey")

plot(1 : S0, rhos[1 : S0])
abline(h = mean(rhos[B : S0]), col = "blue")
# abline(h = sqrt(true_rho), col = "red")
# abline(v = B, col = "grey")
#plot

##assess autocorrelation

par(mfrow = c(5, 1))
Kmax = 35
acf(theta1s[B : S], xlim = c(0, Kmax), lag.max = Kmax)
acf(theta2s[B : S], xlim = c(0, Kmax), lag.max = Kmax)
acf(sigsq1s[B : S], xlim = c(0, Kmax), lag.max = Kmax)
acf(sigsq2s[B : S], xlim = c(0, Kmax), lag.max = Kmax)
acf(rhos[B : S], xlim = c(0, Kmax), lag.max = Kmax)
T = 25
#plot

#burn and thin
theta1s = theta1s[B : S]
theta1s = theta1s[seq(1, S - B, by = T)]
theta2s = theta2s[B : S]
theta2s = theta2s[seq(1, S - B, by = T)]
sigsq1s = sigsq1s[B : S]
sigsq1s = sigsq1s[seq(1, S - B, by = T)]
sigsq2s = sigsq2s[B : S]
sigsq2s = sigsq2s[seq(1, S - B, by = T)]
rhos = rhos[B : S]
rhos = rhos[seq(1, S - B, by = T)]


#look at posteriors with post-exp at 95% CI
par(mfrow = c(5, 1))
res = 200

hist(theta1s, br = res)
abline(v = mean(theta1s), col = "blue", lwd = 3)
abline(v = quantile(theta1s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(theta1s, 0.975), col = "grey", lwd = 3)
abline(v = true_theta_1, col = "red", lwd = 3)

hist(theta2s, br = res)
abline(v = mean(theta2s), col = "blue", lwd = 3)
abline(v = quantile(theta2s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(theta2s, 0.975), col = "grey", lwd = 3)
abline(v = true_theta_2, col = "red", lwd = 3)

hist(sigsq1s, br = res)
abline(v = mean(sigsq1s), col = "blue", lwd = 3)
abline(v = quantile(sigsq1s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(sigsq1s, 0.975), col = "grey", lwd = 3)
abline(v = true_sigsq_1, col = "red", lwd = 3)

hist(sigsq2s, br = res)
abline(v = mean(sigsq2s), col = "blue", lwd = 3)
abline(v = quantile(sigsq2s, 0.025), col = "grey", lwd = 3)
abline(v = quantile(sigsq2s, 0.975), col = "grey", lwd = 3)
abline(v = true_sigsq_2, col = "red", lwd = 3)

hist(rhos, br = res)
abline(v = mean(rhos), col = "blue", lwd = 3)
abline(v = quantile(rhos, 0.025), col = "grey", lwd = 3)
abline(v = quantile(rhos, 0.975), col = "grey", lwd = 3)
abline(v = true_rho, col = "red", lwd = 3)
#plot
