# Converts the FCS data into a matrix for further processing through the pipeline.
# Also generates hypersphere counts using default parameters for countCells.

require(cydar)
require(ncdfFlow)

all.collected <- list()
for (dataset in c("Cytobank_44185_H1", "Cytobank_44185_H2", "Cytobank_44185_H3", "Cytobank_44185_H4", "Cytobank_44185_H5")) {
    x <- read.ncdfFlowSet(list.files(file.path(dataset, "fcs"), full=TRUE))

    # Removing some of the unnecessary channels.
    descriptions <- as.character(parameters(x[[1]])$desc)
    toignore <- c(match(c("Time", "Cell_length", "PhenoGraph"), descriptions), grep("^BC[0-9]", descriptions),
                  match(c("DNA1", "DNA2"), descriptions), match(c("Viability"), descriptions))
    x <- x[,-toignore]

    all.collected[[dataset]] <- x
}

######################################################################
# Creating some pictures showing the expected shift between barcoding batches (using the healthy sample).

dir.create("pics", showWarning=FALSE)
all.markers <- seq_len(ncol(x[[1]]))
s1 <- grep("NoDrug_Basal1", sampleNames(all.collected[[1]]))
s2 <- grep("NoDrug_Basal1", sampleNames(all.collected[[2]]))
s3 <- grep("NoDrug_Basal1", sampleNames(all.collected[[3]]))
s4 <- grep("NoDrug_Basal1", sampleNames(all.collected[[4]]))
s5 <- grep("NoDrug_Basal1", sampleNames(all.collected[[5]]))
descriptions <- as.character(parameters(all.collected[[1]][[1]])$desc)

pdf("pics/ecdf_raw.pdf")
for (m in all.markers) {
    plot(ecdf(exprs(all.collected[[1]][[s1]])[,m]), col="red", main=descriptions[m], xlab="Intensity", cex.axis=1.3, cex.lab=1.5)
    plot(ecdf(exprs(all.collected[[2]][[s2]])[,m]), add=TRUE, col="blue")
    plot(ecdf(exprs(all.collected[[3]][[s3]])[,m]), add=TRUE, col="grey50")
    plot(ecdf(exprs(all.collected[[4]][[s4]])[,m]), add=TRUE, col="forestgreen")
    plot(ecdf(exprs(all.collected[[5]][[s5]])[,m]), add=TRUE, col="black")
}
dev.off()

#####################################################################
# No need to specify experimental design, all batches have the same design.

out <- normalizeBatch(all.collected, batch.comp=NULL, mode="range")   

pdf("pics/ecdf_rnorm.pdf")
for (m in all.markers) {
    plot(ecdf(out[[1]][[s1]][,m]), col="red", main=descriptions[m], xlab="Intensity", cex.axis=1.3, cex.lab=1.5)
    plot(ecdf(out[[2]][[s2]][,m]), add=TRUE, col="blue")
    plot(ecdf(out[[3]][[s3]][,m]), add=TRUE, col="grey50")
    plot(ecdf(out[[4]][[s4]][,m]), add=TRUE, col="forestgreen")
    plot(ecdf(out[[5]][[s5]][,m]), add=TRUE, col="black")
}
dev.off()

#####################################################################
# Repeating with quantile normalization as a demonstration.

out.alt <- normalizeBatch(all.collected, batch.comp=NULL, mode="quantile")   
pdf("pics/ecdf_qnorm.pdf")
for (m in all.markers) {
    plot(ecdf(out.alt[[1]][[s1]][,m]), col="red", main=descriptions[m], xlab="Intensity", cex.axis=1.3, cex.lab=1.5)
    plot(ecdf(out.alt[[2]][[s2]][,m]), add=TRUE, col="blue")
    plot(ecdf(out.alt[[3]][[s3]][,m]), add=TRUE, col="grey50")
    plot(ecdf(out.alt[[4]][[s4]][,m]), add=TRUE, col="forestgreen")
    plot(ecdf(out.alt[[5]][[s5]][,m]), add=TRUE, col="black")
}
dev.off()

#####################################################################
# Cleaning up the memory before preparing the cell data.

x <- unlist(out, recursive=FALSE)
rm(out, all.collected)
gc()

for (sample in seq_along(x)) {
    colnames(x[[sample]]) <- descriptions 
}

# Preparing the data for counting.
cd <- prepareCellData(x)
saveRDS(cd, file="Cytobank_44185_raw.rds")

# Actually counting cells into hyperspheres.
out <- countCells(cd, BPPARAM=SerialParam(), downsample=10, tol=0.5)
saveRDS(out, file="Cytobank_44185_counts.rds")

#####################################################################
# End.
