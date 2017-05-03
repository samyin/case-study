Coord_Brain<-matrix(,68,2)

Coord_Brain[1,]<-c(-7.1,8.9)
Coord_Brain[2,]<-c(-0.8,16.4)
Coord_Brain[3,]<-c(-4.7,14.2)
Coord_Brain[4,]<-c(-1.3,12.1)
Coord_Brain[5,]<-c(-1.1,1.9)
Coord_Brain[6,]<-c(-3.9,15.3)
Coord_Brain[7,]<-c(-4.3,8)
Coord_Brain[8,]<-c(-5,4.2)
Coord_Brain[9,]<-c(-6.4,10.8)
Coord_Brain[10,]<-c(-1.1,7.6)
Coord_Brain[11,]<-c(-3.3,1.5)
Coord_Brain[12,]<-c(-3,19.7)
Coord_Brain[13,]<-c(-1.7,5.5)
Coord_Brain[14,]<-c(-0.9,20.2)
Coord_Brain[15,]<-c(-7.6,11.3)
Coord_Brain[16,]<-c(-3.8,11.4)
Coord_Brain[17,]<-c(-1.5,13.5)
Coord_Brain[18,]<-c(-6.3,16)
Coord_Brain[19,]<-c(-5.2,21)
Coord_Brain[20,]<-c(-6.2,18.3)
Coord_Brain[21,]<-c(-1.2,3.3)
Coord_Brain[22,]<-c(-5.7,9.4)
Coord_Brain[23,]<-c(-0.7,10.8)
Coord_Brain[24,]<-c(-5.3,12.1)
Coord_Brain[25,]<-c(-1.3,4.8)
Coord_Brain[26,]<-c(-0.9,19.1)
Coord_Brain[27,]<-c(-4.5,19.5)
Coord_Brain[28,]<-c(-1.5,16.8)
Coord_Brain[29,]<-c(-2.6,4.2)
Coord_Brain[30,]<-c(-6.8,12.7)
Coord_Brain[31,]<-c(-6.8,8.2)
Coord_Brain[32,]<-c(-0.7,23)
Coord_Brain[33,]<-c(-3.5,17.3)
Coord_Brain[34,]<-c(-5.4,10.6)
Coord_Brain[35:68,]<-Coord_Brain[1:34,]
Coord_Brain[35:68,1]<--Coord_Brain[35:68,1]

# Convert weight matrix to adjacency matrix
W_bin <- W
W_bin[W_bin > 0] = 1

# Compute distance matrix
D <- matrix(nrow = 68, ncol = 68)
for (i in 1:68) {
  for (j in 1:68) {
    D[i, j] <- sqrt((ROI[i, 1] - ROI[j, 1])^2 + (ROI[i, 2] - ROI[j, 2])^2 + 
                      (ROI[i, 3] - ROI[j, 3])^2)/100
  }
}

# Compute indicator for same-Hemisphere
H <- matrix(0, nrow = 68, ncol = 68)
H[1:34, 1:34] = H[35:68, 35:68] = 1


W_per <- W_bin
order <- shuffle(114)
for (i in 1:114) {
  W_per[ , , i] <- W_bin[ , , order[i]]
}
par <- matrix(nrow = 6, ncol = 3)
mis <- matrix(nrow = 6, ncol = 19)

mle <- function(x) {
  f_1 <- 0
  f_2 <- 0
  f_3 <- 0
  for (i in 1:76) {
    for (j in 2:68) {
      for (k in 1:(j - 1)) {
        # three partial derivatives
        r <- x[1] + x[2]*D[j, k] + x[3]*H[j, k]
        f_1 <- f_1 + W_per[j, k, i]*dnorm(r)/pnorm(r) - 
          (1 - W_per[j, k, i])*dnorm(r)/pnorm(-r)
        f_2 <- f_2 + W_per[j, k, i]*dnorm(r)*D[j, k]/pnorm(r) - 
          (1 - W_per[j, k, i])*dnorm(r)*D[j, k]/pnorm(-r)
        f_3 <- f_3 + W_per[j, k, i]*dnorm(r)*H[j, k]/pnorm(r) - 
          (1 - W_per[j, k, i])*dnorm(r)*H[j, k]/pnorm(-r)
      }
    }
  }
  for (i in 77:95) {
    for (j in 2:68) {
      for (k in 1:(j - 1)) {
        # three partial derivatives
        r <- x[1] + x[2]*D[j, k] + x[3]*H[j, k]
        f_1 <- f_1 + W_per[j, k, i]*dnorm(r)/pnorm(r) - 
          (1 - W_per[j, k, i])*dnorm(r)/pnorm(-r)
        f_2 <- f_2 + W_per[j, k, i]*dnorm(r)*D[j, k]/pnorm(r) - 
          (1 - W_per[j, k, i])*dnorm(r)*D[j, k]/pnorm(-r)
        f_3 <- f_3 + W_per[j, k, i]*dnorm(r)*H[j, k]/pnorm(r) - 
          (1 - W_per[j, k, i])*dnorm(r)*H[j, k]/pnorm(-r)
      }
    }
  }
  c(f_1, f_2, f_3)
}

# Solve for the MLE estimate
library("rootSolve")
par[6, ] <- multiroot(f = mle, start = c(1, 1, 1))$root

for (p in 96:114) {
# Simulation
W_est <- matrix(nrow = 68, ncol = 68)
for (i in 2:68) {
  for (j in 1:(i - 1)) {
    W_est[i, j] = W_est[j, i] = rbinom(1, 1, pnorm(par[6, 1] + par[6, 2]*D[i, j] + 
                                                     par[6, 3]*H[i, j]))
  }
}
diag(W_est) <- 0

mis[6, p - 95] <- sum(W_est != W_per[ , , p])/(68*67)
}

# Evaluate with ROC
#library(pROC)
#roc(W_bin[ , , p], W_est, plot = TRUE)
# p is a subject indicator in the test set

pdf("error.pdf")
plot(as.vector(mis), xlab = "", ylab = "", 
     main = "Misclassification Error Rate in Cross Validation", xaxt = "n", 
     cex.main = 1)
dev.off()
summary(as.vector(mis))

apply(par, 2, mean)
summary(auc)
