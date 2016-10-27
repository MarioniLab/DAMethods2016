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
    x <- readRDS(file.path("../refdata", paste0(dataset, "_raw.rds")))
    set.seed(12321)

    for (it in seq_len(10)) {
        for (shift in c(0, 0.1, 0.25, 0.5)) { 
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
                    tol <- 0.5 + shift
                } else {
                    tol <- 0.5
                }

                # Setting up the experimental design.
                out <- countCells(cd, BPPARAM=SerialParam(), tol=tol, downsample=10)
                groupings <- rep(1:2, length.out=ncol(out$counts))
                design <- model.matrix(~factor(groupings))
        
                # Fully fledged edgeR analysis (just to get to p-values).
                y <- DGEList(out$counts, lib.size=out$total)
                keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$total))
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
# Plotting the power results.

#to.save <- read.table("results_edgeR_power.txt", header=TRUE)
#lr <- log10(to.save$edgeR) - log10(to.save$MW)
#col <- rep("grey", length(lr))
#col[lr < 0 & abs(to.save$logFC) > 2] <- "purple"
#col[lr > 0 & abs(to.save$logFC) > 2] <- "darkorange"
#sizes <- -log10(to.save$edgeR)#, to.save$MW))
#
#pdf("plot_power.pdf")
#par(mar=c(5.1, 5.1, 4.1, 2.1))
#plot(to.save$logFC, lr, col=col, pch=16, 
#     xlab=expression(Log[2]*"-fold change in abundance"), 
#     ylab=expression(Log[10]*"-ratio of p-values (edgeR/MW)"), 
#                     cex.axis=1.2, cex.lab=1.4, cex=sizes)
#text.loc <- round(min(to.save$logFC))
#text(text.loc, 0, sum(lr < 0 & abs(to.save$logFC) > 2), col="purple", pos=1, cex=1.5)
#text(text.loc, 0, sum(lr > 0 & abs(to.save$logFC) > 2), col="darkorange", pos=3, cex=1.5)
#legend(text.loc, max(lr), pch=1, pt.cex=c(1, 2, 3), legend=c("0.1", "0.01", "0.001"), cex=1.2)
#dev.off()
#
## Making bar plots of type I error rates.
#
#stuff <- read.table("results_edgeR_alpha.txt", header=TRUE)
#bydataset <- split(stuff[,-1], stuff$Dataset)
#
#all.mean1 <- all.se1 <- all.mean5 <- all.se5 <- list()
#all.mean1.mw <- all.se1.mw <- all.mean5.mw <- all.se5.mw <- list()
#for (x in names(bydataset)) {
#    current.d <- bydataset[[x]]
#    byscheme <- split(current.d[,-1], current.d$Setting)
#    
#    collected.mean1 <- collected.se1 <- collected.mean5 <- collected.se5 <- list()
#    collected.mean1.mw <- collected.se1.mw <- collected.mean5.mw <- collected.se5.mw <- list()
#    for (s in rev(names(byscheme))) { 
#        current.s <- byscheme[[s]]
#        collected.mean1[[s]] <- mean(current.s$edgeR.error1)
#        collected.se1[[s]] <- sqrt(var(current.s$edgeR.error1)/nrow(current.s))
#        collected.mean5[[s]] <- mean(current.s$edgeR.error5)
#        collected.se5[[s]] <- sqrt(var(current.s$edgeR.error5)/nrow(current.s))
#        collected.mean1.mw[[s]] <- mean(current.s$MW.error1)
#        collected.se1.mw[[s]] <- sqrt(var(current.s$MW.error1)/nrow(current.s))
#        collected.mean5.mw[[s]] <- mean(current.s$MW.error5)
#        collected.se5.mw[[s]] <- sqrt(var(current.s$MW.error5)/nrow(current.s))
#    }
#
#    if (x=="Cytobank_43324_4FI") {
#        transfect <- "Oct4-GFP"
#    } else if (x=="Cytobank_43324_NN") {
#        transfect <- "Nanog-Neo"
#    } else {
#        transfect <- "Nanog-GFP"
#    }
#
#    all.mean1[[transfect]] <- unlist(collected.mean1)
#    all.se1[[transfect]] <- unlist(collected.se1)
#    all.mean5[[transfect]] <- unlist(collected.mean5)
#    all.se5[[transfect]] <- unlist(collected.se5)
#    all.mean1.mw[[transfect]] <- unlist(collected.mean1.mw)
#    all.se1.mw[[transfect]] <- unlist(collected.se1.mw)
#    all.mean5.mw[[transfect]] <- unlist(collected.mean5.mw)
#    all.se5.mw[[transfect]] <- unlist(collected.se5.mw)
#}
#
#all.mean1 <- do.call(rbind, all.mean1)
#all.se1 <- do.call(rbind, all.se1)
#all.mean5 <- do.call(rbind, all.mean5)
#all.se5 <- do.call(rbind, all.se5)
#
#all.means <- t(rbind(all.mean1, all.mean5))
#all.se <- all.means + t(rbind(all.se1, all.se5))

############################################ 
# End.

