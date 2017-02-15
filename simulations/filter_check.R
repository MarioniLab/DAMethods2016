# Constructing cell count data for 10 samples in a one-way layout.
# Exploring the effect of the mean on the p-value, so setting up counts with the maximum possible difference at each mean count.

groups <- factor(rep(LETTERS[1:2], each=5))
counts <- matrix(0, ncol=10, nrow=2000)
second <- groups=="B"
counts[,second] <- 1:nrow(counts)/100 

# Running edgeR with offsets of zero (i.e., no normalization necessary, assuming all samples are of equal size).
# Using a NB dispersion of 1 based on real data.
# Using the LRT to compute the p-value to be liberal, i.e., see the smallest p-value we can expect.

design <- model.matrix(~groups)
library(edgeR)
fit <- glmFit(counts, design, dispersion=1, offset=0)
res <- glmLRT(fit)

# Computing the adjusted p-value threshold with various correction factors.
# We see a max correction of >10000 in the BMMC data set, and >200 in the iPSC data set.
# > blah <- readRDS("Cytobank_43324_4FI_res.rds")
# > plot(blah$results$PValue, blah$results$FDR/blah$results$PValue, log="xy")

pdf("plot_filter.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))

mu <- rowMeans(counts)
logp <- log10(res$table$PValue)
plot(mu, logp, xlab="Mean count per hypersphere", ylab=expression("Minimum"~log[10]~"p-value"), cex.axis=1.2, cex.lab=1.4, pch=16, cex=0.5)

correction <- c(2, 4, 6)
colors <- c("orange", "sienna", "red")

for (cx in seq_along(correction)) {
    raw.threshold <- 5*10^-correction[cx]
    threshold <- log10(raw.threshold)
    x.pt <- max(mu[logp > threshold])
    segments(-10, threshold, x.pt, threshold, lwd=2, lty=2, col=colors[cx])
    segments(x.pt, threshold, x.pt, threshold-10, lwd=2, lty=2, col=colors[cx])
    text(-0.3, threshold + 0.2, format(scientific=TRUE, raw.threshold), col=colors[cx], pos=4, offset=0)
}

dev.off()

# The aim is to ask; what is the largest mean where the smallest possible p-value is still above the rejection threshold?
# The exact p-value threshold depends on the effect of multiple testing and the desired FDR threshold.
# With strong correction, the largest mean is quite high (~5) and requires filtering to get rid of it.
# Even with no correction, though, the mean filter is non-zero.

