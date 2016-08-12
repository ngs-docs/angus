library(ggplot2)

w_rep1 <- 'ORE_w_r1_ATCACG_L001'
w_rep2 <- 'ORE_w_r2_GTTTCG_L002'
sde3_rep1 <- 'ORE_sdE3_r1_GTGGCC_L004'
sde3_rep2 <- 'ORE_sdE3_r2_TGACCA_L005'

quants <- list()
quants$w_rep1 <- read.csv(file.path(paste0('quants/',w_rep1,'/quant.sf')), sep='\t')
quants$w_rep2 <- read.csv(file.path(paste0('quants/',w_rep2,'/quant.sf')), sep='\t')
quants$sde3_rep1 <- read.csv(file.path(paste0('quants/',sde3_rep1,'/quant.sf')), sep='\t')
quants$sde3_rep2 <- read.csv(file.path(paste0('quants/',sde3_rep2,'/quant.sf')), sep='\t')

w <- data.frame(x=quants$w_rep1$TPM, y=quants$w_rep2$TPM)
sde3 <- data.frame(x=quants$sde3_rep1$TPM, y=quants$sde3_rep2$TPM)
z <- data.frame(x=quants$w_rep1$TPM, y=quants$sde3_rep2$TPM)

plotw <- ggplot(w, aes(x=x+0.1,y=y+0.1)) + geom_point() + scale_x_log10() + scale_y_log10()
plotsde3 <- ggplot(sde3, aes(x=x+0.1,y=y+0.1)) + geom_point() + scale_x_log10() + scale_y_log10()
plotz <- ggplot(z, aes(x=x+0.1,y=y+0.1)) + geom_point() + scale_x_log10() + scale_y_log10()

plot_density_w <- ggplot(melt(w), aes(x=log(value + 0.1), fill=variable)) + geom_density(alpha = 0.5)
plot_density_sde3 <- ggplot(melt(sde3), aes(x=log(value + 0.1), fill=variable)) + geom_density(alpha = 0.5)
plot_density_z <- ggplot(melt(z), aes(x=log(value + 0.1), fill=variable)) + geom_density(alpha = 0.5)

cor(quants$w_rep1$TPM + 0.1, quants$w_rep2$TPM + 0.1)
plotw
plot_density_w

cor(quants$sde3_rep2$TPM + 0.1, quants$sde3_rep1$TPM + 0.1)
plotsde3
plot_density_sde3

cor(quants$w_rep1$TPM + 0.1, quants$sde3_rep1$TPM + 0.1)
plotz
plot_density_z