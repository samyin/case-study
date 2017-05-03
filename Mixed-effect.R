W_bin <- W
W_bin[W_bin > 0] = 1

D <- matrix(nrow = 68, ncol = 68)
for (i in 1:68) {
  for (j in 1:68) {
    D[i, j] <- sqrt((ROI[i, 1] - ROI[j, 1])^2 + (ROI[i, 2] - ROI[j, 2])^2 + 
                      (ROI[i, 3] - ROI[j, 3])^2)/100
  }
}

H <- matrix(0, nrow = 68, ncol = 68)
H[1:34, 1:34] = H[35:68, 35:68] = 1

edge <- matrix(nrow = 38, ncol = 68*67/2)
for (i in 1:38) {
  B = W_bin[, , i]
  edge[i, ] <- as.vector(B[which(lower.tri(B, diag = FALSE))])
}

dist <- as.vector(D[which(lower.tri(D, diag = FALSE))])
dist_mean <- mean(dist)
hem <- as.vector(H[which(lower.tri(H, diag = FALSE))])
hem_mean <- mean(hem)

iq_adj_0 <- iq_adj[1:38]

iq_adj_mean <- mean(iq_adj_0)

pop <- c(1, dist_mean, hem_mean, iq_adj_mean)

sub <- array(dim = c(4, 68*67/2, 38))
for (i in 1:38) {
  for (n in 1:(68*67/2)) {
    sub[, n, i] <- c(1, dist[n], hem[n], iq_adj_0[i])
  }
}

scode_3 <- "
data {
  int<lower=0> I;
  int<lower=0> N;
  int<lower=0,upper=1> x[I, N];
  real<lower=0> pop[4];
  real<lower=0> sub[4, N, I];
}
parameters {
  real alpha[4];
  real beta[I, 4];
}
transformed parameters {
  real<lower=0, upper=1> p[I, N];
  for (i in 1:I) {
    for (n in 1:N) {
      p[i, n] = normal_cdf((alpha[1] * pop[1] + alpha[2] * pop[2] + 
                              alpha[3] * pop[3] + alpha[4] * pop[4] + 
                              beta[i, 1] * sub[1, n, i] + 
                              beta[i, 2] * sub[2, n, i] + 
                              beta[i, 3] * sub[3, n, i] + 
                              beta[i, 4] * sub[4, n, i]), 0, 1);
    }
  }
}
model {
  for (i in 1:I) {
    for (n in 1:N) {
      x[i, n] ~ bernoulli(p[i, n]);
    }
  }
}"

dat_3 <- list(
  I = 38,
  N = 68*67/2,
  x = edge,
  pop = pop,
  sub = sub
)

fit3 <- stan(model_name = "fit3", model_code = scode_3, data = dat_3, 
             pars = c("alpha", "beta"), include = TRUE, chains = 1, iter = 500)
fit3_summary <- as.matrix(summary(fit3)$summary)
head(fit3_summary)
