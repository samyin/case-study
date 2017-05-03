lymph <- read.table(gzfile("lymph.fr.gz"), sep = ",")
gene_names <- vector(length = ncol(lymph) - 2)
for (i in 1:length(gene_names)) {
  gene_names[i] <- paste0("gene_", i)
}
colnames(lymph) <- c("survtime", "status", gene_names)
remove(gene_names)

lymph_uncensored <- lymph[lymph$status == 0, ]
lymph_censored <- lymph[lymph$status == 1, ]
lymph_uncensored <- lymph_uncensored[, -2]
lymph_censored <- lymph_censored[, -2]

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## imputation ##
scode_1 <- "
data {
  int<lower=0> N_cen;
  int<lower=0> N_uncen;
  int<lower=0> V;
  matrix[N_cen, V] X_cen;
  vector[N_cen] Y_cen;
  matrix[N_uncen, V] X_uncen;
  vector[N_uncen] Y_uncen;
}
parameters {
  vector[N_cen] Y_impute;
  real alpha;
  vector[V] theta;
  real<lower=0> sigma2;
  real<lower=0> a;
  real<lower=0> b;
}
model {
  sigma2 ~ inv_gamma(a, b);
  for (i in 1:N_uncen) {
    Y_uncen[i] ~ normal(alpha + X_uncen[i, ] * theta, sqrt(sigma2));}
  for (j in 1:N_cen) {
    Y_impute[j] ~ normal(alpha + X_cen[j, ] * theta, sqrt(sigma2)) T[Y_cen[j], Y_cen[j] + 4];}
}"

var_impute <- as.character(read.csv("Ordered_Pvalues.csv", header = TRUE)[1:500, 1])

dat_1 <- list(
  N_cen = nrow(lymph_censored), 
  N_uncen = nrow(lymph_uncensored), 
  V = length(var_impute), 
  X_cen = lymph_censored[, var_impute], 
  Y_cen = log(lymph_censored[, 1]), 
  X_uncen = lymph_uncensored[, var_impute], 
  Y_uncen = log(lymph_uncensored[, 1])
)

dat_1$Y_cen[dat_1$Y_cen == -Inf] = min(dat_1$Y_cen[dat_1$Y_cen != -Inf])

init <- function(x) {
  list(Y_impute = dat_1$Y_cen + 1, 
       alpha = 1, 
       theta = rep(0, length(var_impute)), 
       sigma2 = 8, 
       a = 2,
       b = 8)
}

n_chains <- 1
init_1 <- lapply(1:n_chains, init)

Y_impute <- matrix(nrow = 138, ncol = 20)
rhat <- matrix(nrow = 138, ncol = 20)
Y_impute_est <- vector(length = 138)
for (i in 1:20) {
  fit_1 <- stan(model_name = "impute_1", model_code = scode_1, data = dat_1, 
                pars = c("Y_impute"), include = TRUE, 
                chains = n_chains, iter = 5000, init = init_1)
  fit_1_summary <- as.matrix(summary(fit_1)$summary)
  Y_impute[, i] <- as.vector(fit_1_summary[1:138, 1])
  rhat[, i] <- as.vector(fit_1_summary[1:138, 10])
}
for (j in 1:138) {
  min <- which.min(rhat[j, ])
  Y_impute_est[j] <- Y_impute[j, min]
}

Y_impute_est <- exp(Y_impute_est)
lymph_censored$survtime <- Y_impute_est
lymph_imputed <- rbind(lymph_censored, lymph_uncensored)
write.csv(lymph_imputed, "lymph_impute.csv")

## latent factor model ##
lymph_uncensored_norm_x <- lymph_uncensored[, var_impute]
for (i in 1:ncol(lymph_uncensored_norm_x)) {
  mean = mean(lymph_uncensored_norm_x[, i])
  sd = sd(lymph_uncensored_norm_x[, i])
  lymph_uncensored_norm_x[, i] <- (lymph_uncensored_norm_x[, i] - mean)/sd
}
lymph_uncensored_norm_y <- log(lymph_uncensored[, 1])
lymph_uncensored_norm <- cbind(lymph_uncensored_norm_x, lymph_uncensored_norm_y)

library(caret)
holdout_set <- createFolds(lymph_uncensored_norm[, 501], k = 4, list = FALSE, 
                           returnTrain = FALSE)
test <- lymph_uncensored_norm[holdout_set == 1, ]
train <- lymph_uncensored_norm[holdout_set != 1, ]

scode_3 <- "
data {
  int<lower=0> N;
  int<lower=0> P;
  int<lower=0> K;
  matrix[N, P] X;
  vector[N] Y;
  real<lower=0, upper=1> w;
}
parameters {
  matrix[P, K] lambda;
  matrix[N, K] F;
  vector<lower=0>[P] tao_f;
  real<lower=0> a_f;
  real<lower=0> b_f;
  matrix<lower=0>[P, K] phi;
  real<lower=0> nu;
  real<lower=0> delta1;
  real<lower=0> delta2;
  real<lower=0> a1;
  real<lower=0> a2;
  
  real alpha;
  vector[K] theta;
  real<lower=0> tao_y;
  real<lower=0> a_y;
  real<lower=0> b_y;
}
transformed parameters {
  vector<lower=0>[K] tao_lambda;
  for (k in 1:K)
    tao_lambda[k] = delta1 + (k - 1) * delta2;
}
model {
  delta1 ~ gamma(a1, 1);
  delta2 ~ gamma(a2, 1);
  tao_y ~ gamma(a_y, b_y);

  for (p in 1:P) {
    tao_f[p] ~ gamma(a_f, b_f);
    for (k in 1:K) {
      phi[p, k] ~ gamma(nu/2, nu/2);
      lambda[p, k] ~ normal(0, 1/sqrt(phi[p, k] * tao_lambda[k]));}}
  
  for (n in 1:N) {
    target += normal_lpdf(Y[n]| alpha + F[n, ] * theta, 1/sqrt(tao_y));
    for (p in 1:P) {
      target += w * normal_lpdf(X[n, p]| F[n, ] * (lambda[p, ])', 1/sqrt(tao_f[p]));}}
}"


data_i <- list(
  N = nrow(train), 
  P = 500, 
  K = 20, 
  X = train[, 1:500], 
  Y = train[, 501], 
  w = 1
)

fit_5 <- stan(model_name = "latent_weighted", model_code = scode_3, data = data_i, 
            chains = 1, iter = 500)
fit_5_summary <- as.matrix(summary(fit_5)$summary)
lambda_est_5 <- matrix(fit_5_summary[1:10000, 1], nrow = 500, ncol = 20, byrow = TRUE)
alpha_est_5 <- fit_5_summary[nrow(train) * 20 + 20508, 1]
theta_est_5 <- as.vector(fit_5_summary[(nrow(train) * 20 + 20509):(nrow(train) * 20 + 20528), 1])

pred_5 <- exp(as.matrix(test[, -501]) %*% lambda_est_5 %*% solve(t(lambda_est_5) %*% lambda_est_5) %*% theta_est_5 + alpha_est_5)
err_5 <- pred_5 - exp(test[, 501])

library(glmnet)
lasso_1 <- glmnet(as.matrix(train[, 1:500]), train[, 501], family = "gaussian", 
                  alpha = 1)
lasso_1_coef <- coef(lasso_1)
lasso_1_coef <- lasso_1_coef[, 100]
pred_lasso <- exp(predict(lasso_1, as.matrix(test[, -501]), s = 0.01))
err_lasso <- pred_lasso - exp(test[, 501])

library(caret)
fit_control <- trainControl(method = "boot")
train_data <- train_set[-1]
test_data <- test_set[-1]
library(bartMachine)
bart_model <- train(lymph_imputed_norm_y_new ~ ., data = train_data, 
                    method = "bartMachine", metric = "RMSE", 
                    trControl = fit_control, verbose = FALSE, tuneLength = 5)
library(rpart)
cart_model <- train(lymph_imputed_norm_y_new ~ ., data = train_data, 
                    method = "rpart", metric = "RMSE", trControl = fit_control, 
                    verbose = FALSE, tuneLength = 5)