W_bin <- W
W_bin[W_bin > 0] <- 1

D <- matrix(nrow = 68, ncol = 68)
for (i in 1:68) {
  for (j in 1:68) {
    D[i, j] <- sqrt((ROI[i, 1] - ROI[j, 1])^2 + (ROI[i, 2] - ROI[j, 2])^2 + 
                      (ROI[i, 3] - ROI[j, 3])^2)/100
  }
}

H <- matrix(0, nrow = 68, ncol = 68)
H[1:34, 1:34] = H[35:68, 35:68] = 1

#CCI <- covariate[, 3]
#CCI[which(is.na(CCI))] <- mean(CCI, na.rm = TRUE)

#CAQ <- covariate[, 8]
#CAQ[which(is.na(CAQ))] <- mean(CAQ, na.rm = TRUE)

#attach(covariate)

edge <- matrix(nrow = 76, ncol = 68*67/2)
for (i in 1:76) {
  B = W_bin[, , i]
  edge[i, ] <- as.vector(B[which(lower.tri(B, diag = FALSE))])
}

dist <- as.vector(D[which(lower.tri(D, diag = FALSE))])
hem <- as.vector(H[which(lower.tri(H, diag = FALSE))])

iq_adj_0 <- iq_adj[1:76]

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

scode_1 <- "
data {
  int<lower=0> I;
  int<lower=0> N;
  int<lower=0,upper=1> x[I, N];
  real<lower=0> dist[N];
  real<lower=0, upper=1> hem[N];
  real<lower=0> y[I];
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real delta;
  real<lower=0> sigma[4];
}
transformed parameters {
  real<lower=0, upper=1> p[I, N];
  for (i in 1:I) {
    for (n in 1:N) {
      p[i, n] = normal_cdf((alpha + beta * dist[n] + gamma * hem[n] + 
                              delta * y[i]), 0, 1);
    }
  }
}
model {
  alpha ~ normal(0, sigma[1]);
  beta ~ normal(0, sigma[2]);
  gamma ~ normal(0, sigma[3]);
  delta ~ normal(0, sigma[4]);
  for (i in 1:I) {
    for (n in 1:N) {
      x[i, n] ~ bernoulli(p[i, n]);
    }
  }
}"

scode_2 <- "
data {
  int<lower=0> I;
  int<lower=0> N;
  int<lower=0,upper=1> x[I, N];
  real<lower=0> dist[N];
  real<lower=0, upper=1> hem[N];
  real<lower=0> y[I];
}
parameters {
  real alpha;
  real beta;
  real gamma;
  real delta;
  real<lower=0> mu[I];
  real<lower=0> sigma;
}
transformed parameters {
  real<lower=0, upper=1> p[I, N];
  for (i in 1:I) {
    for (n in 1:N) {
      p[i, n] = normal_cdf((alpha + beta * dist[n] + gamma * hem[n] + 
                              delta * mu[I]), 0, 1);
    }
  }
}
model {
  for (i in 1:I) {
    y[i] ~ normal(mu[i], sigma);
  }
  for (i in 1:I) {
    for (n in 1:N) {
      x[i, n] ~ bernoulli(p[i, n]);
    }
  }
}"

dat_1 <- list(
  I = 76,
  N = 68*67/2,
  x = edge,
  dist = dist,
  hem = hem,
  y = iq_adj_0
)

dat_2 <- list(
  I = 76,
  N = 68*67/2,
  x = edge,
  dist = dist,
  hem = hem,
  y = iq_adj_0
)

fit1 <- stan(model_name = "fit1", model_code = scode_1, data = dat_1, 
             pars = c("alpha", "beta", "gamma", "delta", "sigma"), 
             include = TRUE, chains = 1, iter = 500)
#print(fit1)

fit2 <- stan(model_name = "fit2", model_code = scode_2, data = dat_2, 
             pars = c("alpha", "beta", "gamma", "delta", "mu", "sigma"), 
             include = TRUE, chains = 1, iter = 500)
#print(fit2)

fit1_summary <- as.matrix(summary(fit1)$summary)
stan_hist(fit1, bins = 20)
par_est <- fit1_summary[1:4, 1]
alpha_est <- par_est[1]
beta_est <- par_est[2]
gamma_est <- par_est[3]
delta_est <- par_est[4]

error_1 <- vector(length = 38)
for (i in 77:114) {
  c_est <- vector(length = 68*67/2)
  for (n in 1:(68*67/2)) {
    p <- pnorm(alpha_est + beta_est * dist[n] + gamma_est * hem[n] + 
                 delta_est * iq_adj[i])
    c_est[n] <- rbinom(1, 1, p)
  }
  
  B = W_bin[, , i]
  c <- as.vector(B[which(lower.tri(B, diag = FALSE))])
  
  diff <- abs(c_est - c)
  error_1[i - 76] <- sum(diff)/length(diff)
}

plot(error_1, xlab = "index of validation test", ylab = "error rate", type = "b")

fit2_summary <- as.matrix(summary(fit2)$summary)
stan_plot(fit2, pars = c("alpha", "beta", "gamma","delta"))
par_est <- fit2_summary[1:81, 1]
alpha_est <- par_est[1]
beta_est <- par_est[2]
gamma_est <- par_est[3]
delta_est <- par_est[4]
mu_est <- par_est[5:80]

error_2 <- vector(length = 76)

for (i in 1:76) {
  p_est <- vector(length = 68*67/2)
  c_est <- vector(length = 68*67/2)
  for (n in 1:(68*67/2)) {
    p_est[n] = pnorm(alpha_est + beta_est * dist[n] + gamma_est * hem[n] + 
                       delta_est * mu_est[i])
    c_est[n] <- rbinom(1, 1, p_est[n])
  }
  
  diff <- abs(c_est - edge[i, ])
  error_2[i] <- sum(diff)/length(diff)
}
