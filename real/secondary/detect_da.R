# Testing for differential abundance and saving the results.

require(cydar)
out <- readRDS(file.path("../../refdata", "Cytobank_44185.rds"))

# Setting up the design matrix.
patient <- sub(".*_H([0-9]+)\\..*", "\\1", colnames(out))
treatment <- sub(".*NoDrug_([^_]+)_.*", "\\1", colnames(out))
design <- model.matrix(~0 + treatment + patient)
colnames(design) <- sub("-", "", colnames(design))

# Testing for differential proportions across time, using the spline.
require(edgeR)
y <- DGEList(assay(out), lib.size=out$totals)
keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$totals))
out <- out[keep,]
y <- y[keep,]
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

res <- glmQLFTest(fit, contrast=makeContrasts(treatmentIL10 - treatmentBasal1, levels=design))
qvals <- spatialFDR(intensities(out), res$table$PValue)
saveRDS(list(coords=intensities(out), ranges=intensityRanges(out, p=0.01),
    results=data.frame(logFC=res$table$logFC, AveLogCPM=y$AveLogCPM, Dispersion=fit$dispersion, PValue=res$table$PValue, FDR=qvals)),
file="IL10_res.rds")

#res <- glmQLFTest(fit, contrast=makeContrasts(treatmentIL3 - treatmentBasal1, levels=design))
#qvals <- spatialFDR(intensities(out), res$table$PValue)
#saveRDS(list(coords=intensities(out), ranges=intensityRanges(cd, p=0.01),
#    results=data.frame(logFC=res$table$logFC, AveLogCPM=y$AveLogCPM, Dispersion=fit$dispersion, PValue=res$table$PValue, FDR=qvals)),
#file="IL3_res.rds")

