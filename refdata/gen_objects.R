# Converts the FCS data into a matrix for further processing through the pipeline.
# Also generates hypersphere counts using default parameters for countCells.

require(cyder)
require(ncdfFlow)
for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    x <- read.ncdfFlowSet(list.files(file.path(dataset, "fcs"), full=TRUE))

    # Removing some of the unnecessary channels.
    toignore <- c(match(c("Cell_length", "Time", "beadDist", "barcode"), colnames(x)), grep("^BC-[0-9]", colnames(x)), grep("DNA", colnames(x)))
    x <- x[,-toignore]
    cd <- prepareCellData(x)

    out <- countCells(cd, BPPARAM=MulticoreParam(3), downsample=10, tol=0.5)
    saveRDS(cd, file=paste0(dataset, "_raw.rds"))
    saveRDS(out, file=paste0(dataset, ".rds"))
}


