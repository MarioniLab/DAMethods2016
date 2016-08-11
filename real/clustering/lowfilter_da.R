#####################################################################
# This checks the correspondence when the filter threshold is reduced.

# Loading data and counts.
require(cydar)
dataset <- "Cytobank_43324_4FI"
out <- readRDS(file=file.path("../../refdata", paste0(dataset, "_counts.rds")))

# Setting up the design matrix.
timings <- as.integer(sub(".*_([0-9]+).fcs", "\\1", colnames(out$counts)))
design <- model.matrix(~splines::ns(timings, 3))

# Testing for differential proportions across time, using the spline.
require(edgeR)
y <- DGEList(out$counts, lib.size=out$total, genes=out$coordinates)
keep <- aveLogCPM(y) >= aveLogCPM(1, mean(out$total)) # REDUCED FILTER!
y <- y[keep,]

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
res <- glmQLFTest(fit, coef=2:ncol(design))

qvals <- spatialFDR(out$coordinates[keep,], res$table$PValue)

#####################################################################
# Checking correspondence.

it <- 1L
collected.failed <- list()
all.hypers <- t(out$coordinates)#[keep,][qvals <= 0.05,])

for (k in c(50, 100, 200)) { 
    incoming <- read.table(sprintf("da_clusters_%i.tsv", k), check.names=FALSE)
    m <- match(colnames(out$coordinates), colnames(incoming))
    is.sig <- which(incoming$FDR <= 0.05)
    
    failed <- 0L
    it <- 1L
    max.dist <- 0
    for (i in is.sig) {
        distances <- sqrt(colSums((all.hypers - as.numeric(incoming[i,m]))^2))
        closest <- which.min(distances)
        if (distances[closest] > sqrt(ncol(out$coordinates))*0.5) { 
            failed <- failed + 1L
        }
    }
    collected.failed[[as.character(k)]] <- failed
}

collected.failed

#####################################################################
