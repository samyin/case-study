attach(covariate)
covariate <- cbind(INDICATOR, covariate[, 4:5], CCI, covariate[, 7:13])

pdf("plot_1.pdf")
hist(Age, freq = FALSE, main = "Distribution of Subject Age (n = 114)")
dev.off()

pdf("plot_2.pdf")
boxplot(covariate[, 4:6], cex.main = 1, main = "Boxplot of CCI, FSIQ, and CAQ")
dev.off()

pdf("plot_3.pdf")
boxplot(covariate[, 7:11], cex.axis = 0.57, cex.main = 1, 
        main = "Boxplot of the Big Five Personality Traits score")
dev.off()

pdf("plot_4.pdf")
library(corrplot)
CM <- cor(covariate[, 2:11], use = "complete.obs")
corrplot(CM, method = "number")
title("Correlation Matrix of Cognitive Traits")
dev.off()

pdf("plot_5.pdf")
# Adjacency Matrix
pdf("plot_binary.pdf")
W_bin <- W
W_bin[W_bin > 0] = 1
W_bin_mean <- apply(W_bin, c(1, 2), mean)
qgraph(W_bin_mean, layout = Coord_Brain, vsize = c(2.5, 4))
title(main = "Average Brain Connectivity Diagram (Binary, n = 114)")
dev.off()

# Connectedness
c <- rep(0, 114)
for (i in 1:114) {
  s <- W_bin[ , , i]
  c[i] <- (sum(s > 0))/(68*67)
}

# The ratio of inter-hemispheric connection
lr_ratio <- rep(0, 114)
for (i in 1:114) {
  lr_ratio[i] <- (sum(W_bin[1:34, 35:68, i]) + 
                    sum(W_bin[35:68, 1:34, i]))/(sum(W_bin[1:68, 1:68, i]))
}

hist(c, col = "lightblue", xlab = "", xlim = c(0, 0.6), freq = FALSE, 
     main = "Distribution of Connectedness and Inter-hemispheric Connection Ratio", 
     cex.main = 0.8)
hist(lr_ratio, freq = FALSE, col = "lightgreen", add = TRUE)
l <- list(c, lr_ratio)
pdf("2hist.pdf")
multhist(l, freq = FALSE, breaks = 10, ylab = "Density", col = c("lightblue", "lightgreen"))
title(main = "Distribution of Connectedness and Inter-hemispheric Connection Ratio", cex.main = 1)
legend("topleft", cex = 0.75, legend = c("Connectedness", "Inter-hemispheric Connection Ratio"), 
       col = c("lightblue", "lightgreen"), lty = 1)
dev.off()

pdf("plot_6.pdf")
# Frequncy of Connection across brain
average <- apply(W_bin, c(1,2), sum)/114
heatmap(indicators, Rowv = NA, Colv = NA)
dev.off()

library(qgraph)
library(igraph)
library(sna)

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

pdf("plot_7.pdf")
w <- apply(W, c(1, 2), median)
qgraph(w, layout = Coord_Brain, vsize = c(2.5, 5))
title(main = "Brain Connectivity based on Median Weight Matrix (n = 114)")
dev.off()