# Converts the FCS data into a matrix for further processing through the pipeline.
# Also generates hypersphere counts using default parameters for countCells.

require(cydar)
require(ncdfFlow)
for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    x <- read.ncdfFlowSet(list.files(file.path(dataset, "fcs"), full=TRUE))

    # Removing some of the unnecessary channels.
    toignore <- c(match(c("Cell_length", "Time", "beadDist", "barcode"), colnames(x)), grep("^BC-[0-9]", colnames(x)), grep("DNA", colnames(x)))
    x <- x[,-toignore]

    # Counting the data into hyperspheres.
    cd <- prepareCellData(x)
    out <- countCells(cd, BPPARAM=SerialParam(), downsample=10, tol=0.5)
    saveRDS(out, file=paste0(dataset, ".rds"))
}


