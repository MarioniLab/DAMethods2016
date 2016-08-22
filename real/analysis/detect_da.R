# Testing for differential abundance and saving the results.

require(cydar)
for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    cd <- readRDS(file.path("../../refdata", paste0(dataset, "_raw.rds")))
    out <- readRDS(file.path("../../refdata", paste0(dataset, "_counts.rds")))

    # Setting up the design matrix.
    timings <- as.integer(sub(".*_([0-9]+).fcs", "\\1", colnames(out$counts)))
    design <- model.matrix(~splines::ns(timings, 3))

    # Testing for differential proportions across time, using the spline.
    require(edgeR)
    y <- DGEList(out$counts, lib.size=out$total, genes=out$coordinates)
    keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$total))
    y <- y[keep,]

    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit, coef=2:ncol(design))

    # Computing fold changes.
    redesign <- model.matrix(~timings)
    refit <- glmFit(y, redesign, dispersion=1)
    reres <- glmLRT(refit)

    qvals <- spatialFDR(out$coordinates[keep,], res$table$PValue)
    saveRDS(list(coords=y$genes, 
                 ranges=intensityRanges(cd, p=0.01),
                 results=data.frame(logFC=reres$table$logFC, AveLogCPM=y$AveLogCPM, Dispersion=fit$dispersion, PValue=res$table$PValue, FDR=qvals)),
            file=paste0(dataset, "_res.rds"))
}

