dataset <- "Cytobank_43324_4FI"

#stuff <- readRDS(file=paste0(dataset, "_res.rds"))
#qvals <- stuff$results$FDR
#coords <- stuff$coords    
#is.sig <- qvals <= 0.05
#sig.coords <- coords[is.sig,]
#require(Rtsne)
#set.seed(100)
#tsne.out <- Rtsne(sig.coords, perplexity=10)

tsne.out <- list(Y=read.table(paste0(dataset, "_coords.txt"), row.names=1))
sig.coords <- tsne.out$Y

##############################################################
# Plotting the behaviour of the new signalling population.

## This code block can be used to select the population from the t-SNE plot.
## For simplicity's sake, I've hard-coded the relevant hyperspheres in 'of.interest'.
#require(scran)
#app <- selectorPlot(tsne.out$Y[,1], tsne.out$Y[,2])
#saved <- shiny::runApp(app)
#which(saved[[1]])

of.interest <- c(444, 1007, 1080, 1236, 2324, 2743, 2828, 3408, 4825, 4974, 5281, 5807, 5985, 6187)
of.interest <- rownames(sig.coords)[of.interest]

# plot(tsne.out$Y[,1], tsne.out$Y[,2])
# points(tsne.out$Y[of.interest,1], tsne.out$Y[of.interest,2], col="red")

raw <- readRDS(file.path("../../refdata", paste0(dataset, "_counts.rds")))
props <- t(t(raw$counts[of.interest,])/raw$totals)*100
timings <- as.integer(sub(".*_([0-9]+).fcs", "\\1", colnames(raw$counts)))

pdf(file.path("pics", "novel_timecourse.pdf"))
plot(0, 0, ylab="Cell abundance (%)", xlab="Time (days)", xlim=range(timings), ylim=range(props), 
     type="n", cex.lab=1.4, cex.axis=1.2)
all.cols <- rainbow(nrow(props))
for (i in seq_len(nrow(props))) {
    lines(timings, props[i,], col="grey80", lwd=2)
}
lines(timings, colMeans(props), col="black", lwd=3)
dev.off()

##############################################################
# Making a plot for the specific time points.

library(edgeR)
ti <- factor(timings)
fit <- glmFit(raw$counts[rownames(sig.coords),], design=model.matrix(~0 + ti), offset=log(raw$totals), dispersion=1, prior.count=3)
induction <- (fit$coefficients[,levels(ti)=="1"] - fit$coefficients[,levels(ti)=="0"])/log(2)
withdrawal <- (fit$coefficients[,levels(ti)=="17"] - fit$coefficients[,levels(ti)=="16"])/log(2)

library(cydar)
mfc <- 5

png(file.path("pics", paste0(dataset, "_induction.png")), width=6, height=6, units="in", res=300)
plotCellLogFC(tsne.out$Y[,1], tsne.out$Y[,2], induction, max.logFC=mfc, 
              main="After induction", xlab="t-SNE1", ylab="t-SNE2", 
              cex.axis=1.2, cex.lab=1.4, cex.main=1.4)
dev.off()

png(file.path("pics", paste0(dataset, "_withdrawal.png")), width=6, height=6, units="in", res=300)
out <- plotCellLogFC(tsne.out$Y[,1], tsne.out$Y[,2], withdrawal, max.logFC=mfc,
              main="After withdrawal", xlab="t-SNE1", ylab="t-SNE2", 
              cex.axis=1.2, cex.lab=1.4, cex.main=1.4)
dev.off()

png(file.path("pics", "extra_bar.png"), width=0.9, height=6, units="in", res=300)
par(mar=c(0,0,0,0))
plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
start.loc <- seq(-0.5, 0.5, length.out=length(out))
interval <- diff(start.loc)[1]
rect(-0.5, start.loc, 0.5, start.loc+interval, col=out, border=NA)
text(0, -0.5, pos=1, -mfc, cex=1.5)
text(0, 0.5+interval,  pos=3, mfc, cex=1.5)
text(-0.9, 0, pos=1, srt=90, "Log-fold change", cex=1.5)
dev.off()
