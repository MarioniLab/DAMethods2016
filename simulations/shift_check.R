# This checks how edgeR performs when intensities are systematically shifted between samples.
# We use Poisson sampling as any loss of control is more evident in that simulation.

require(cydar)
require(edgeR)
source("functions.R")

ofile <- "results_shift_alpha.txt"
existing <- FALSE
saved.FC <- FALSE

############################################ 
# Setting up

for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    x <- readRDS(file.path("../refdata", paste0(dataset, ".rds")))
    set.seed(12321)

    for (it in seq_len(10)) {
        for (shift in c(0, 0.1, 0.2, 0.3)) { 
            current.exprs <- resampleCells(x, setting=1L)

            # Adding a random, sample-specific intensity shift to all markers.
            for (lib in seq_along(current.exprs)) {
                whee <- current.exprs[[lib]]
                for (m in seq_len(ncol(whee))) {
                    whee[,m] <- whee[,m] + rnorm(1, mean=0, sd=shift) # Common intensity shift to all cells.
                }
                current.exprs[[lib]] <- whee
            }

            cd <- prepareCellData(current.exprs)
            for (expansion in c(TRUE, FALSE)) { 
                if (expansion) { 
                    tol <- expandRadius(cd, tol=0.5)
                } else {
                    tol <- 0.5
                }

                # Setting up the experimental design.
                out <- countCells(cd, BPPARAM=SerialParam(), tol=tol, downsample=10)
                groupings <- rep(1:2, length.out=ncol(out))
                design <- model.matrix(~factor(groupings))
        
                # Fully fledged edgeR analysis (just to get to p-values).
                y <- DGEList(assay(out), lib.size=out$totals)
                keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$totals))
                out <- out[keep,]
                y <- y[keep,]
                
                y <- estimateDisp(y, design)
                fit <- glmQLFit(y, design, robust=TRUE)
                res <- glmQLFTest(fit)

                write.table(file=ofile, data.frame(Dataset=dataset, Shift=shift, Expansion=expansion,
                                                   edgeR.error1=sum(res$table$PValue <= 0.01)/nrow(res),
                                                   edgeR.error5=sum(res$table$PValue <= 0.05)/nrow(res),
                                                   edgeR.Dispersion=y$common),
                            quote=FALSE, sep="\t", col.names=!existing, append=existing, row.names=FALSE)
                existing <- TRUE
                gc()    
            }
        }
    }
}

############################################ 
# Making bar plots of type I error rates and dispersions.

stuff <- read.table("results_shift_alpha.txt", header=TRUE)
bydataset <- split(stuff[,-1], stuff$Dataset)

all.mean1 <- all.se1 <- all.disp <- all.dispse <- list()
datasets <- mode <- list()
for (x in names(bydataset)) {
    current.d <- bydataset[[x]]
    byscheme <- split(current.d[,-1], current.d$Shift)
    
    collected.mean1 <- collected.se1 <- collected.disp <- collected.dispse <- list()
    for (s in rev(names(byscheme))) {
        current.s <- byscheme[[s]]
        is.expanded <- current.s$Expansion
        collected.mean1[[s]] <- c(mean(current.s$edgeR.error1[!is.expanded]),
                                  mean(current.s$edgeR.error1[is.expanded]))
        collected.se1[[s]] <- sqrt(c(var(current.s$edgeR.error1[!is.expanded]),
                                     var(current.s$edgeR.error1[is.expanded])))/nrow(current.s)

#        current.s$edgeR.Dispersion <- (current.s$edgeR.Dispersion)^0.25
        collected.disp[[s]] <- c(mean(current.s$edgeR.Dispersion[!is.expanded]),
                                 mean(current.s$edgeR.Dispersion[is.expanded]))
        collected.dispse[[s]] <- sqrt(c(var(current.s$edgeR.Dispersion[!is.expanded]),
                                        var(current.s$edgeR.Dispersion[is.expanded])))/nrow(current.s)
    }

    if (x=="Cytobank_43324_4FI") {
        transfect <- "Oct4-GFP"
    } else if (x=="Cytobank_43324_NN") {
        transfect <- "Nanog-Neo"
    } else {
        transfect <- "Nanog-GFP"
    }

    all.mean1[[transfect]] <- do.call(cbind, collected.mean1)
    all.se1[[transfect]] <- do.call(cbind, collected.se1)
    all.disp[[transfect]] <- do.call(cbind, collected.disp)
    all.dispse[[transfect]] <- do.call(cbind, collected.dispse)
    datasets[[transfect]] <- rep(transfect, nrow(all.mean1[[1]]))
    mode[[transfect]] <- c("default", "expanded")
}

datasets <- unlist(datasets)
mode <- unlist(mode)
all.mean1 <- do.call(rbind, all.mean1)
all.se1 <- do.call(rbind, all.se1)
all.disp <- do.call(rbind, all.disp)
all.dispse <- do.call(rbind, all.dispse)

o <- order(as.numeric(colnames(all.mean1)))
all.mean1 <- all.mean1[,o]
all.se1 <- all.se1[,o]
all.disp <- all.disp[,o]
all.dispse <- all.dispse[,o]
               
colours <- c("plum", "purple", "lightblue", "blue", "palegreen", "forestgreen")

pdf("plot_shift.pdf")
all.upper <- all.mean1 + all.se1
out <- barplot(all.mean1, beside=TRUE, col=colours, xlab="Intensity shift per marker",
               ylab="Observed type I error rate", cex.axis=1.2, cex.lab=1.4, cex.names=1.2, 
               ylim=c(0, max(all.upper)))
segments(out, all.mean1, out, all.upper)
segments(out-0.1, all.upper, out+0.1, all.upper)
abline(h=0.01, col="red", lwd=2, lty=2)

legend(min(out), max(all.mean1), fill=colours, sprintf("%s (%s)", datasets, mode), cex=1.2)

all.upper <- all.disp + all.dispse
out <- barplot(all.disp, beside=TRUE, col=colours, xlab="Intensity shift per marker",
               ylab="Common NB dispersion", cex.axis=1.2, cex.lab=1.4, cex.names=1.2,
               ylim=c(0, max(all.upper)))
segments(out, all.disp, out, all.upper)
segments(out-0.1, all.upper, out+0.1, all.upper)
dev.off()

############################################ 
# End.

