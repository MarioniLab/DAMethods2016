#####################################
# Simulating a data set with one subpopulation containing a shift in one dimension.

source("functions.R")
set.seed(500)
samples <- c(1,1,2,2)
design <- model.matrix(~factor(samples))

ncells <- 20000
nmarkers <- 30
signaller <- 30

###################################
# Running the simulation.

require(cydar)
library(edgeR)

collected.res <- list()
for (mode in c("large", "small")) {
    if (mode=="large") {
        lfc <- log10(40)
    } else {
        lfc <- log10(2)
    }

    for (type in c("major", "minor")) { 
        for (it in 1:50) { 
            combined <- matrix(rnorm(ncells*nmarkers, 1, sd=0.5), ncol=nmarkers)
            sample.ids <- sample(length(samples), ncells, replace=TRUE)

            affected <- which(samples[sample.ids]==2)
            if (type=="minor") {
                # Only taking 10%
                affected <- sample(affected, length(affected)/10)
            }
            combined[affected,signaller] <- combined[affected,signaller] + lfc 

            # Splitting by sample in preparation for counting.
            colnames(combined) <- paste0("X", seq_len(nmarkers))
            by.sample <- list()
            for (x in seq_along(samples)) { 
                by.sample[[x]] <- combined[sample.ids==x,,drop=FALSE]
            }
            names(by.sample) <- paste0("Y", seq_along(by.sample))
           
            for (anal in c("standard", "shifting")) {
                if (anal=="standard") {
                    cd <- prepareCellData(by.sample)
                    out <- countCells(cd, downsample=10, tol=0.5, BPPARAM=SerialParam())
                    
                    # Testing for DA.
                    y <- DGEList(assay(out), lib.size=out$totals)
                    keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$totals))
                    out <- out[keep,]
                    y <- y[keep,]
                    y <- estimateDisp(y, design)
                    fit <- glmQLFit(y, design, robust=TRUE)
                    
                    res <- glmQLFTest(fit)
                    coords <- intensities(out)
                    qval <- spatialFDR(coords, res$table$PValue)

                } else {
                    cd <- prepareCellData(by.sample, markers=colnames(combined)[-signaller])
                    out <- countCells(cd, downsample=10, tol=0.5, BPPARAM=SerialParam())
                    out2 <- medIntensities(out, markers=signaller)
                    med.signal <- assay(out2, "med.X30")
                    cell.count <- assay(out2, "counts") 
                    
                    el <- new("EList", list(E=med.signal, weights=cell.count))
                    fit <- lmFit(el, design)
                    fit <- eBayes(fit, robust=TRUE)
                    res <- topTable(fit, coef=ncol(design), n=Inf, sort.by="none")
                    qval <- spatialFDR(intensities(out), res$P.Value)
                }
                
                # Checking if any of the positions contain the affected cells.
                # For DA, we need to pick up hyperspheres at the shifted position.
                # For shifting, we need to pick up hyperspheres containing affected cells.
                da.hypersphere <- qval <= 0.05 & !is.na(qval)
                if (anal=="standard") {
                    detected <- any(abs(coords[da.hypersphere,signaller] - (1+lfc)) < 0.5 & res$table$logFC[da.hypersphere] > 0)
                } else {
                    candidates <- da.hypersphere & res$logFC > 0
                    for (x in which(candidates)) {
                        cell.ids <- cellData(out)$cell.id[unpackIndices(cellAssignments(out)[x])[[1]]]
                        if (!any(cell.ids %in% affected)) candidates[x] <- FALSE                            
                    }
                    detected <- any(candidates)
                }
                
                cur.sim <- paste(type, mode, anal, sep=".")
                if (is.null(collected.res[[cur.sim]])) {
                    collected.res[[cur.sim]] <- 0
                } 
                if (detected) {
                    collected.res[[cur.sim]] <- collected.res[[cur.sim]] + 1
                }
            }
        }
    }
}

stuff <- strsplit(names(collected.res), split="\\.")
stuff <- do.call(rbind, stuff)
colnames(stuff) <- c("Size", "Shift", "Analysis")
write.table(file="results_signal.txt", data.frame(stuff, Detected=unlist(collected.res)),
            row.names=FALSE, sep="\t", quote=FALSE)

###################################
# Making a pretty plot.

blah <- read.table("results_signal.txt", header=TRUE, stringsAsFactors=FALSE)
blah.standard <- blah[blah$Analysis=="standard",]
blah.shifting <- blah[blah$Analysis=="shifting",]

stopifnot(identical(blah.standard$Size, blah.shifting$Size))
stopifnot(identical(blah.standard$Shift, blah.shifting$Shift))

p <- rbind(blah.standard$Detected, blah.shifting$Detected)/50

replacement <- c(major="All", minor="10%", large="40-fold", small="2-fold")
colnames(p) <- paste0(replacement[blah.standard$Size], ", ", replacement[blah.standard$Shift])

pdf("plot_signal.pdf")
cols <- c("grey30", "grey80")
loc <- barplot(p, beside=TRUE, ylab="Detection frequency", cex.axis=1.2, cex.lab=1.4, col=cols, cex.names=1.2)
stder <- sqrt(p*(1-p)/50)
upper <- stder + p
segments(loc, p, loc, upper)
segments(loc-0.1, upper, loc+0.1, upper)

legend("topright", fill=cols, c("DA", "Direct"), cex=1.2)
dev.off()
