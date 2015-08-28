# Chris Hamm NGS-2015 sesion on sex-linked dosage compensation

set.seed(13342143)
sessionInfo()
setwd("~/Desktop/Projects/NGS-course/SLDC")

#### First we will make up, errrr, simulate data
# simulate read counts under a log normal distribution
go.1 <- rlnorm(n = 1e4, meanlog = 3, sdlog = 1)
head(go.1)
hist(go.1, col = "grey", main = "", las = 1, breaks = 50)

#### Now we will simulate data for males and females under equal expression
#Male data
M.exp <- go.1 + rnorm(1e4, mean = 0, sd = 10)
M.exp[M.exp < 0.1] <- 0.1 # we will be using log2, and don't want the world to implode.

#Female data
F.exp <- go.1 + rnorm(1e4, mean = 0, sd = 10)
F.exp[F.exp < 0.1] <- 0.1

# Combine into data.frame
df.2 <- data.frame(exp = go.1, male = M.exp, female = F.exp, chr = as.factor(rep(1:10)))
head(df.2)
tail(df.2)
str(df.2)

#### We can explore the data a bit:
plot(M.exp, F.exp, pch = 19, las = 1, col = rgb(0, 0, 0, 0.2))

#### Make your own M:A plot
plot(log2(go.1), log2(M.exp) - log2(F.exp), pch = 19, col = rgb(0, 0, 0, 0.2), ylim = c(-10, 10), xlim = c(0, 10), las = 1, xlab = expression(paste("log Expression")), ylab = expression(paste("log"[2], " Male - log"[2], " Female")))

# plot males count density
plot(density(log2(M.exp)), lwd = 3, main = "") 
# overlay female count density
lines(density(log2(F.exp)), lwd = 3, main = "", col = "red", lty = 2) 
# overlay the base expression
lines(density(log2(df.2$exp)), lwd = 3, col = "dodgerblue", lty = 4) 

#### Note the bimodality? That is a consequence of not wanting the world to implode. We'll get rid of that in short order.

#### Let's plot the log2 expression for males and females by chromosome:
par(mfrow = c(1, 2))
# Chromosome 1 is the sex chromosome
boxplot(log2(df.2$male) ~ df.2$chr, las =1, pch = 19, col = c("red", rep("dodgerblue", 9)), ylim = c(-1, 9), main = "male", outline = FALSE) 
boxplot(log2(df.2$female) ~ df.2$chr, las =1, pch = 19, col = c("red", rep("dodgerblue", 9)), ylim = c(-1, 9), main = "female", outline = FALSE)
par(mfrow = c(1, 1))

#### Cool, Males and Females show the same expression by chromosome. We can look at this other ways:
# This is chromosome by chromosome
boxplot((log2(df.2$male) - log2(df.2$female)) ~ df.2$chr, pch = 19, col = c("red", rep("dodgerblue", 9)), las = 1, ylim = c(-3, 3), main = expression(paste("log"[2], "Male - log"[2], "Female")), outline = FALSE) 

# Autosomes vs. sex-chromosome
plot(density(log2(df.2$male[df.2$chr != 1]) - log2(df.2$female[df.2$chr != 1])), lwd = 3, las = 1, main = "", xlab = "")
lines(density(log2(df.2$male[df.2$chr == 1]) - log2(df.2$female[df.2$chr == 1])), lwd = 3, col = "red")

#### Now let's add a dosage effect that is **_equal_** in both sexes (70% reduction)
df.2$M.z <- df.2$male
df.2$M.z[df.2$chr == 1] <- df.2$male[df.2$chr == 1] * 0.7

df.2$F.z <- df.2$female
df.2$F.z[df.2$chr == 1] <- df.2$female[df.2$chr == 1] * 0.7
head(df.2)

par(mfrow = c(1, 2))
boxplot(log2(df.2$M.z) ~ df.2$chr, las =1, pch = 19, col = c("red", rep("dodgerblue", 9)), ylim = c(-1, 9), main = "male", outline = FALSE)
boxplot(log2(df.2$F.z) ~ df.2$chr, las =1, pch = 19, col = c("red", rep("dodgerblue", 9)), ylim = c(-1, 9), main = "female", outline = FALSE)
par(mfrow = c(1, 1))

#### The effect is pretty visible in the chromosome plots
#### Now let's look at what happens if we double male expression
df.2$trt <- df.2$male
df.2$trt[df.2$chr == 1] <- df.2$male[df.2$chr == 1] * 2
head(df.2)

par(mfrow = c(1, 2))
boxplot(log2(df.2$trt) ~ df.2$chr, las =1, pch = 19, col = c("red", rep("dodgerblue", 9)), ylim = c(-1, 10), main = "male", outline = FALSE)
boxplot(log2(df.2$F.z) ~ df.2$chr, las =1, pch = 19, col = c("red", rep("dodgerblue", 9)), ylim = c(-1, 10), main = "female", outline = FALSE)
par(mfrow = c(1, 1))

#### You can really see male expression go up:
# with density plots for male vs. female
plot(density(log2(df.2$trt[df.2$chr != 1]) - log2(df.2$F.z[df.2$chr != 1])), lwd = 3, las = 1, main = "", xlab = "", ylim = c(0, 0.65))
lines(density(log2(df.2$trt[df.2$chr == 1]) - log2(df.2$F.z[df.2$chr == 1])), lwd = 3, col = "red")

# and with M:F boxplots
boxplot((log2(df.2$trt) - log2(df.2$female)) ~ df.2$chr, pch = 19, col = c("red", rep("dodgerblue", 9)), las = 1, ylim = c(-3, 4), main = expression(paste("log"[2], "Male - log"[2], "Female")), outline = FALSE)

#### Now let's play with some imaginary data from the silkmoth Bombyx mori
T <- FALSE
T
Bmori.dat <- read.csv("Data/Bmori-data.csv", header = T)
head(Bmori.dat) 

Bmori.dat <- read.csv("Data/Bmori-data.csv", header = TRUE)
head(Bmori.dat)

TRUE <- FALSE
TRUE <- F
TRUE <- "anything"
str(Bmori.dat)

#### The data are:
# f68 = Female control (Fcont)
# m68 = Male control (Mcont)
# f69 = Female experimental (Fexp)
# m70 = Male experimental (Mexp)
#### Now we'll look at our data a bit:

# Make your own M:A plots
M.Mcont_Fcont<- log2(Bmori.dat$rpkm.m68) - log2(Bmori.dat$rpkm.f67) 
A.Mcont_Fcont <- log2((Bmori.dat$rpkm.m68 + Bmori.dat$rpkm.f67) / 2)

M.Mcont_Fexp <- log2(Bmori.dat$rpkm.m68) - log2(Bmori.dat$rpkm.f69)
A.Mcont_Fexp <- log2((Bmori.dat$rpkm.m68 + Bmori.dat$rpkm.f69) / 2)

M.Mexp_Fcont <- log2(Bmori.dat$rpkm.m70) - log2(Bmori.dat$rpkm.f67)
A.Mexp_Fcont <- log2((Bmori.dat$rpkm.m70 + Bmori.dat$rpkm.f67) / 2)

M.Mexp_Fexp <- log2(Bmori.dat$rpkm.m70) - log2(Bmori.dat$rpkm.f69)
A.Mexp_Fexp <- log2((Bmori.dat$rpkm.m70 + Bmori.dat$rpkm.f69) / 2)

# make the four pairwise plots
par(mfrow = c(2, 2))
plot(A.Mcont_Fcont, M.Mcont_Fcont, main = "Mcont - Fcont", pch = 19, col = rgb(0, 0, 0, 0.1), cex = 0.8, ylim = c(-4, 6), xlim = c(-5, 17), ylab = expression(paste("log"[2], " expression")), xlab = "")
abline(h = 0, lwd = 3, lty = 1, col = "red")

plot(A.Mcont_Fcont, M.Mcont_Fexp, main = "Mcont - Fexp", pch = 19, col = rgb(0, 0, 0, 0.1), cex = 0.8, ylim = c(-8, 8), xlim = c(-5, 17), xlab = "")
abline(h = 0, lwd = 3, lty = 1, col = "red")

plot(A.Mexp_Fcont, M.Mexp_Fcont, main = "Mexp - Fcont", pch = 19, col = rgb(0, 0, 0, 0.1), cex = 0.8, ylim = c(-6, 8), xlim = c(-5, 17), ylab = expression(paste("log"[2], " expression")), xlab = "Expected count")
abline(h = 0, lwd = 3, lty = 1, col = "red")

plot(A.Mexp_Fexp, M.Mexp_Fexp, main = "Mexp - Fexp", pch = 19, col = rgb(0, 0, 0, 0.1), cex = 0.8, ylim = c(-5, 7), xlim = c(-5, 17), xlab = "Expected count")
abline(h = 0, lwd = 3, lty = 1, col = "red")
par(mfrow = c(1, 1))

#### Now we plot Male and female control groups by chromosome

par(mfrow = c(1, 2))
boxplot(log2(Bmori.dat$rpkm.m68[Bmori.dat$rpkm.m68 >= 1]) ~ Bmori.dat$chr[Bmori.dat$rpkm.m68 >= 1], pch = 19, notch = TRUE, col = c("red", rep("grey", 27)), na.rm = TRUE, las = 1, ylim = c(-1, 15), main = "Control (male)", ylab = expression(paste("log"[2]," expression (RPKM)")), outline = FALSE)

boxplot(log2(Bmori.dat$rpkm.m70[Bmori.dat$rpkm.m70 >= 1]) ~ Bmori.dat$chr[Bmori.dat$rpkm.m70 >= 1],  pch = 19, notch = TRUE, col = c("red", rep("grey", 27)), na.rm = TRUE, las = 1, ylim = c(-1, 15), main = expression(paste(italic("Expr"), " knockdown (male)")), ylab = expression(paste("log"[2]," expression (RPKM)")), outline = FALSE)
par(mfrow = c(1, 1))

#### We'll write a function to speed things up a bit:
plot.ratios <- function(chr = Bmori.dat$chr, minRPKM = 1, sampA, sampB, ...) {
  rpkm <- data.frame(chr, sampA, sampB)
  filter <- (sampA >= minRPKM & sampB >= minRPKM)
  rpkm <- rpkm[filter,]
  rpkm$lratio <- log2(rpkm$sampA) - log2(rpkm$sampB)
  boxplot(rpkm$lratio ~ rpkm$chr, ...)
}
plot.ratios(sampA = Bmori.dat$rpkm.m70, sampB=Bmori.dat$rpkm.f69, minRPKM = 0.1,  main = expression(paste(italic(Expr), " Male vs. Female")), outline = FALSE,  col = c("red", rep("grey", 27)), notch = TRUE, ylab = expression(paste("log" [2], "(Male) - log" [2], "(Female)")), las = 1)

#### Clearly the experimental treatment has resulted in the upregulation of the sex chromosome

#### We can also plot sex-chromosome vs. autosomes:
plot.density.za <- function(chr = Bmori.dat$chr, minRPKM = 1, sampA, sampB, ...) {
  rpkm <- data.frame(chr, sampA, sampB)
  filter <- (sampA > minRPKM & sampB > minRPKM)
  rpkm <- rpkm[filter,]
  rpkm$lratio <- log2(rpkm$sampA) - log2(rpkm$sampB)
  rpkm$za <- ifelse(rpkm$chr == 1, "Z", "A")
  Zvals <- rpkm$lratio[rpkm$za == "Z"]
  Avals <- rpkm$lratio[rpkm$za == "A"]
  plot(density(Avals), las = 1, ...)
  abline(v = median(Avals), lty = 2, lwd = 3)
  lines(density(Zvals), col = "red", lwd = 3)
  abline(v = median(Zvals), lty = 2, col = "red", lwd = 3)
}

plot.density.za(sampA = Bmori.dat$rpkm.m70, sampB = Bmori.dat$rpkm.f69, minRPKM = 1,  main = expression(paste(italic("Masc"), ": Male vs. Female ")), xlim = c(-2,2), lwd = 3, xlab = expression(paste("Log"[2], " (Male:Female)")), ylab = "")
legend("topright", legend = c("Z Chromosome", "Autosomes"), lty = 1, col = c("red", "black"), bty = "n", lwd = 3)


get.za.stats <- function(chr = Bmori.dat$chr, minRPKM = 1, samp, plot.it = FALSE, ...) {
  rpkm <- data.frame(chr, samp)
  rpkm <- rpkm[samp > minRPKM,]
  rpkm$za <- ifelse(rpkm$chr == 1, "Z", "A")
  if (plot.it == TRUE) {
    boxplot(log2(rpkm$samp) ~ rpkm$za, ...)
  }
  Zvals <- rpkm$samp[rpkm$za == "Z"]
  Avals <- rpkm$samp[rpkm$za == "A"]
  Zmean <- round(mean(Zvals), digits = 3)
  Zmedian <- round(median(Zvals), digits = 3)
  Amean <- round(mean(Avals), digits = 3)
  Amedian <- round(median(Avals), digits = 3)
  za.mean <- round(Zmean/Amean, digits = 4)
  za.median <- round(Zmedian/Amedian, digits = 4)
  wilcox <- wilcox.test(Zvals, Avals)
  out.stats <- list(	"Zmean" = Zmean, "Amean" = Amean, "Z:A_mean"= za.mean, "Zmedian" = Zmedian,  "Amedian" = Amedian, "Z:A_median" = za.median, "ZvA_MWU-p" = signif(wilcox$p.value, digits = 4), "n_Zgenes" = length(Zvals), "n_Agenes" = length(Avals) 
  )				
  return(out.stats)
}
sapply(Bmori.dat[3:6], FUN = function(x) get.za.stats(samp = x))
