data <- read.csv("case02_data.csv", stringsAsFactors = FALSE)
data$Prenatal_Care <- factor(data$Prenatal_Care, 
                             levels = c("inade", "inter", "adequ", "inten"), 
                             labels = c("1", "2", "3", "4"), ordered = TRUE)
for (i in 1:ncol(data)) {
  data[, i] <- as.numeric(data[, i])
}

library("mice")
data <- complete(mice(data))

model_1 <- glm(formula = Gest_age_.37wks ~ Socioeconomic_Index, 
               data = data, family = binomial(link = "logit"))
summary(model_1)
model_2 <- glm(formula = Gest_age_.37wks ~ DDE + Socioeconomic_Index + Sum_of_PCBS, 
               data = data, family = binomial(link = "logit"))
summary(model_2)
model_3 <- glm(formula = Y ~ X, family = binomial(link = "logit"))
summary(model_3)

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

gest_week <- vector(length = nrow(data))
for (i in 1:length(gest_week)) {
  if (data$Gest_Days[i] <= 217) {
    gest_week[i] <- 1
  }
  else if (data$Gest_Days[i] > 329) {
    gest_week[i] <- 18
  }
  else {
    gest_week[i] <- ceiling(data$Gest_Days[i]/7) - 30
  }
}

I <- sum(gest_week)
X <- matrix(0, nrow = I, ncol = 26)

l <- 0
for (i in 1:length(gest_week)) {
  l <- l + gest_week[i]
  for (j in 1:gest_week[i]) {
    X[l - gest_week[i] + j, j] <- 1
    X[l - gest_week[i] + j, 19:25] <- c(data$Sum_of_PCBS[i], log(data$DDE[i]), 
                                        log(data$Triglycerides[i]), 
                                        data$Socioeconomic_Index[i], 
                                        data$Num_Visits[i], data$Race_Black[i], 
                                        data$BMI_pre_Preg[i])
  }
  X[l, 26] <- 1
}

Y <- X[, 26]
X <- X[, -26]

scode_1 <- "
data {
  int<lower=0> I;
  int<lower=0> N;
  int<lower=0> T_0;
  matrix<lower=0>[I, T_0 + N] X;
  int<lower=0, upper=1> Y[I];
}
parameters {
  vector[T_0] alpha;
  real alpha_0;
  real mu_0;
  real<lower=0> sigma2_0;
  real<lower=0> sigma2;
  real<lower=0> a_0;
  real<lower=0> b_0;
  
  vector[N] beta;
  vector[N] mu_beta;
  vector<lower=0>[N] sigma2_beta;
}
model {
  alpha_0 ~ normal(mu_0, sqrt(sigma2_0));
  sigma2 ~ inv_gamma(a_0, b_0);
  alpha[1] ~ normal(alpha_0, sqrt(sigma2));
  for (t in 2:T_0) 
    alpha[t] ~ normal(alpha[t - 1], sqrt(sigma2));
  
  for (n in 1:N) 
    beta[n] ~ normal(mu_beta[n], sqrt(sigma2_beta[n]));
  
  for (i in 1:I) 
    Y[i] ~ bernoulli_logit(X[i] * append_row(alpha, beta));
}"

dat_1 <- list(
  I = 10000, 
  N = 7, 
  T_0 = 18, 
  X = X[1:10000, ], 
  Y = Y[1:10000]
)

fit_1 <- stan(model_name = "fit_1", model_code = scode_1, data = dat_1, 
              pars = c("alpha", "beta", "sigma2"), include = TRUE, 
              chains = 4, iter = 1000, control = list(adapt_delta = 0.9))
print(fit_1)
stan_plot(fit_1, pars = c("alpha", "beta", "sigma2"), include = TRUE)

fit_1_summary <- as.matrix(summary(fit_1)$summary)
par_est_1 <- fit_1_summary[1:25, 1]
alpha_est_1 <- par_est_1[1:18]
beta_est_1 <- par_est_1[19:25]

X_t <- X[10001:23581, ]
Y_t <- Y[10001:23581]
err_1 <- vector(length = 13581)
for (i in 1:13581) {
  logit <- 1/(1 + exp(-X_t[i, ] * par_est_1))
  y_pre <- rbinom(1, 1, logit)
  err_1[i] <- abs(y_pre - Y_t[i])
}
sum(err_1)/length(err_1)

# Final report #

library("rstan")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

gest_week <- vector(length = nrow(data))
for (i in 1:length(gest_week)) {
  if (data$Gest_Days[i] <= 217) {
    gest_week[i] <- 1
  }
  else if (data$Gest_Days[i] > 329) {
    gest_week[i] <- 18
  }
  else {
    gest_week[i] <- ceiling(data$Gest_Days[i]/7) - 30
  }
}

I <- sum(gest_week)
N <- 7
X <- matrix(0, nrow = I, ncol = N)
Ti <- vector("integer", length = I)
Y <- vector("integer", length = I)

l <- 0
for (i in 1:length(gest_week)) {
  for (j in 1:gest_week[i]) {
    X[l + j, 1:N] <- c(data$Sum_of_PCBS[i], log(data$DDE[i]), 
                       log(data$Triglycerides[i]), data$Socioeconomic_Index[i], 
                       data$Num_Visits[i], data$Race_Black[i], 
                       data$BMI_pre_Preg[i])
    Ti[l + j] <- j
  }
  l <- l + gest_week[i]
  Y[l] <- 1
}

scode_2 <- "
data {
  int<lower=0> I;
  int<lower=0> N;
  int<lower=0> Tn;
  matrix<lower=0>[I, N] X;
  int<lower=1, upper=Tn> Ti[I];
  int<lower=0, upper=1> Y[I];
}
parameters {
  vector[Tn] alpha;
  real alpha_0;
  real mu_0;
  real<lower=0> sigma2_0;
  real<lower=0> sigma2;
  real<lower=0> a_0;
  real<lower=0> b_0;

  matrix[N, Tn] beta;
  matrix[N, Tn] mu_beta;
  vector<lower=0>[N] prec_beta;
}
transformed parameters {
  matrix[N, N] omega_beta;
  omega_beta = diag_matrix(prec_beta);
}
model {
  alpha_0 ~ normal(mu_0, sqrt(sigma2_0));
  sigma2 ~ inv_gamma(a_0, b_0);
  alpha[1] ~ normal(alpha_0, sqrt(sigma2));
  for (t in 2:Tn) 
    alpha[t] ~ normal(alpha[t - 1], sqrt(sigma2));

  for (t in 1:Tn) 
    beta[, t] ~ multi_normal_prec(mu_beta[, t], omega_beta);

  for (i in 1:I) 
    Y[i] ~ bernoulli_logit(alpha[Ti[i]] + X[i] * beta[, Ti[i]]);
}"

sample_1 <- sample.int(I, size = 20000)
X_1 <- matrix(nrow = 20000, ncol = N)
Ti_1 <- vector("integer", length = 20000)
Y_1 <- vector("integer", length = 20000)
for (i in 1:20000) {
  X_1[i, ] <- X[sample_1[i], ]
  Ti_1[i] <- Ti[sample_1[i]]
  Y_1[i] <- Y[sample_1[i]]
}

X_t <- X[-sample_1, ]
Ti_t <- Ti[-sample_1]
Y_t <- Y[-sample_1]

dat_2 <- list(
  I = 20000, 
  N = 7, 
  Tn = 18, 
  X = X_1, 
  Ti = Ti_1, 
  Y = Y_1
)

fit_2 <- stan(model_name = "fit_2", model_code = scode_2, data = dat_2, 
              pars = c("alpha", "beta", "sigma2"), include = TRUE, 
              chains = 1, iter = 1000)
print(fit_2)
stan_rhat(fit_2, pars = "beta")
fit_2_summary <- as.matrix(summary(fit_2)$summary)
par_est_2 <- fit_2_summary[, 1]
alpha_est_2 <- par_est_2[1:18]
beta_est_2 <- par_est_2[19:144]

err_2 <- vector(length = 3581)
Y_pre <- vector(length = 3581)
for (i in 1:3581) {
  logit <- 1/(1 + exp(-(alpha_est_2[Ti_t[i]] + X_t[i, 1] * beta_est_2[Ti_t[i]] + 
                          X_t[i, 2] * beta_est_2[Ti_t[i] + 18] + 
                          X_t[i, 3] * beta_est_2[Ti_t[i] + 36] + 
                          X_t[i, 4] * beta_est_2[Ti_t[i] + 54] + 
                          X_t[i, 5] * beta_est_2[Ti_t[i] + 72] + 
                          X_t[i, 6] * beta_est_2[Ti_t[i] + 90] + 
                          X_t[i, 7] * beta_est_2[Ti_t[i] + 108])))
  Y_pre[i] <- rbinom(1, 1, logit)
  err_2[i] <- abs(y_pre - Y_t[i])
}
sum(err_2)/length(err_2)

library(pROC)
plot.roc(Y_t, Y_pre, print.auc = TRUE)