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
cols <- c("red", "blue", "grey50", "forestgreen", "black")

pdf("pics/ecdf_raw2.pdf")
for (m in all.markers) {
    collected <- list(exprs(all.collected[[1]][[s1]])[,m],
                      exprs(all.collected[[2]][[s2]])[,m], 
                      exprs(all.collected[[3]][[s3]])[,m],
                      exprs(all.collected[[4]][[s4]])[,m],
                      exprs(all.collected[[5]][[s5]])[,m])
    cydar:::multiIntHist(collected, cols=cols, main=descriptions[m], cex.axis=1.3, cex.lab=1.5, cex.main=1.5, lwd=2)
}
dev.off()

collated <- do.call(diffIntDistr, all.collected)
pdf("pics/diff_raw.pdf", width=10, height=20)
par(mfrow=c(8,4), mar=c(2.1, 2.1, 2.1, 1.1))
for (m in names(collated$difference)) { 
    current <- collated$difference[[m]] * 100
    is.inter <- collated$index==0L
    is.intra <- !is.inter & lower.tri(collated$index)
    boxplot(list(Inter=current[is.inter], Intra=current[is.intra]), main=m)
}
dev.off()

#####################################################################
# No need to specify experimental design, all batches have the same design.

out <- normalizeBatch(all.collected, batch.comp=NULL, mode="range")   
pdf("pics/ecdf_rnorm2.pdf")
for (m in all.markers) {
    collected <- list(out[[1]][[s1]][,m],
                      out[[2]][[s2]][,m], 
                      out[[3]][[s3]][,m],
                      out[[4]][[s4]][,m],
                      out[[5]][[s5]][,m])
    cydar:::multiIntHist(collected, cols=cols, main=descriptions[m], cex.axis=1.3, cex.lab=1.5, cex.main=1.5, lwd=2)
}
dev.off()

rcollated <- do.call(diffIntDistr, out)
pdf("pics/diff_rnorm.pdf", width=10, height=20)
par(mfrow=c(8,4), mar=c(2.1, 2.1, 2.1, 1.1))
for (m in names(rcollated$difference)) { 
    current <- rcollated$difference[[m]] * 100
    is.inter <- rcollated$index==0L
    is.intra <- !is.inter & lower.tri(rcollated$index)
    boxplot(list(Inter=current[is.inter], Intra=current[is.intra]), main=m)
}
dev.off()

#####################################################################
# Repeating with quantile normalization as a demonstration.

out.alt <- normalizeBatch(all.collected, batch.comp=NULL, mode="warp")   

pdf("pics/ecdf_wnorm2.pdf")
for (m in all.markers) {
    collected <- list(out.alt[[1]][[s1]][,m],
                      out.alt[[2]][[s2]][,m], 
                      out.alt[[3]][[s3]][,m],
                      out.alt[[4]][[s4]][,m],
                      out.alt[[5]][[s5]][,m])
    cydar:::multiIntHist(collected, cols=cols, main=descriptions[m], cex.axis=1.3, cex.lab=1.5, cex.main=1.5, lwd=2)
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

# Counting cells into hyperspheres.
cd <- prepareCellData(x)

diag <- neighborDistances(cd)
pdf("pics/distances.pdf")
boxplot(diag, ylab=expression(r[0]*sqrt(2)), ylim=c(0, 0.9), xlab="Neighbours", cex.axis=1.2, cex.lab=1.4, outline=FALSE)
abline(h=0.55, col="red", lwd=2, lty=2)
dev.off()

out <- countCells(cd, BPPARAM=SerialParam(), downsample=10, tol=0.55)
saveRDS(out, file="Cytobank_44185.rds")

#####################################################################
# End.
