# This checks that edgeR runs properly on CyToF data. 
# We do so by pooling all cells from all real samples into a single matrix, and then sample from that matrix to constitute each sample.

require(cyder)
require(edgeR)
source("functions.R")

ofile <- "results_edgeR_alpha.txt"
existing <- FALSE
saved.FC <- FALSE

############################################ 
# Setting up

for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    x <- readRDS(file.path("../refdata", paste0(dataset, "_raw.rds")))
    set.seed(12321)

    for (it in seq_len(10)) {
        for (setting in 1:2) {
            setname <- ifelse(setting==1L, "Random", "Noisy")
            current.exprs <- resampleCells(x, setting=setting)
    
            # Setting up the experimental design.
            cd <- prepareCellData(current.exprs)
            out <- countCells(cd, BPPARAM=SerialParam(), downsample=10)
            groupings <- rep(1:2, length.out=ncol(out$counts))
            design <- model.matrix(~factor(groupings))

            # Fully fledged edgeR analysis (just to get to p-values).
            y <- DGEList(out$counts, lib.size=out$total)
            keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$total))
            y <- y[keep,]

            y <- estimateDisp(y, design)
            fit <- glmQLFit(y, design, robust=TRUE)
            res <- glmQLFTest(fit)

            # Trying a Mann-Whitney analysis.
            props <- t(t(y$counts)/out$total)
            propA <- props[,groupings==1L]
            propB <- props[,groupings==2L]
            mw.out <- sapply(seq_len(nrow(y)), FUN=function(i) {
                   suppressWarnings(wilcox.test(propA[i,], propB[i,])$p.value)
            })

            write.table(file=ofile, data.frame(Dataset=dataset, Setting=setname, 
                                               edgeR.error1=sum(res$table$PValue <= 0.01)/nrow(res),
                                               edgeR.error5=sum(res$table$PValue <= 0.05)/nrow(res),
                                               edgeR.Dispersion=y$common,
                                               MW.error1=sum(mw.out <= 0.01)/length(mw.out),
                                               MW.error5=sum(mw.out <= 0.05)/length(mw.out)),
                        quote=FALSE, sep="\t", col.names=!existing, append=existing, row.names=FALSE)
            existing <- TRUE

            if (setting==2L && !saved.FC) {
                # Saving p-values and log-fold changes for visualization later.
                to.save <- data.frame(logFC=res$table$logFC, edgeR=res$table$PValue, MW=mw.out)
                write.table(file="results_edgeR_power.txt", to.save, row.names=FALSE, sep="\t", quote=FALSE)
                saved.FC <- TRUE
            }

            gc()    
        }
    }
}

############################################ 
# Plotting the power results.

to.save <- read.table("results_edgeR_power.txt", header=TRUE)
lr <- log10(to.save$edgeR) - log10(to.save$MW)
col <- rep("grey", length(lr))
col[lr < 0 & abs(to.save$logFC) > 2] <- "purple"
col[lr > 0 & abs(to.save$logFC) > 2] <- "darkorange"
sizes <- -log10(to.save$edgeR)#, to.save$MW))

pdf("plot_power.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
plot(to.save$logFC, lr, col=col, pch=16, 
     xlab=expression(Log[2]*"-fold change in abundance"), 
     ylab=expression(Log[10]*"-ratio of p-values (edgeR/MW)"), 
                     cex.axis=1.2, cex.lab=1.4, cex=sizes)
text.loc <- round(min(to.save$logFC))
text(text.loc, 0, sum(lr < 0 & abs(to.save$logFC) > 2), col="purple", pos=1, cex=1.5)
text(text.loc, 0, sum(lr > 0 & abs(to.save$logFC) > 2), col="darkorange", pos=3, cex=1.5)
legend(text.loc, max(lr), pch=1, pt.cex=c(1, 2, 3), legend=c("0.1", "0.01", "0.001"), cex=1.2)
dev.off()

# Making bar plots of type I error rates.

stuff <- read.table("results_edgeR_alpha.txt", header=TRUE)
bydataset <- split(stuff[,-1], stuff$Dataset)

all.mean1 <- all.se1 <- all.mean5 <- all.se5 <- list()
for (x in names(bydataset)) {
    current.d <- bydataset[[x]]
    byscheme <- split(current.d[,-1], current.d$Setting)
    
    collected.mean1 <- collected.se1 <- collected.mean5 <- collected.se5 <- list()
    for (s in rev(names(byscheme))) { 
        current.s <- byscheme[[s]]
        collected.mean1[[s]] <- mean(current.s$edgeR.error1)
        collected.se1[[s]] <- sqrt(var(current.s$edgeR.error1)/nrow(current.s))
        collected.mean5[[s]] <- mean(current.s$edgeR.error5)
        collected.se5[[s]] <- sqrt(var(current.s$edgeR.error5)/nrow(current.s))
    }

    if (x=="Cytobank_43324_4FI") {
        transfect <- "Oct4-GFP"
    } else if (x=="Cytobank_43324_NN") {
        transfect <- "Nanog-Neo"
    } else {
        transfect <- "Nanog-GFP"
    }

    all.mean1[[transfect]] <- unlist(collected.mean1)
    all.se1[[transfect]] <- unlist(collected.se1)
    all.mean5[[transfect]] <- unlist(collected.mean5)
    all.se5[[transfect]] <- unlist(collected.se5)
}

all.mean1 <- do.call(rbind, all.mean1)
all.se1 <- do.call(rbind, all.se1)
all.mean5 <- do.call(rbind, all.mean5)
all.se5 <- do.call(rbind, all.se5)

all.means <- t(rbind(all.mean1, all.mean5))
all.se <- all.means + t(rbind(all.se1, all.se5))

pdf("plot_alpha.pdf")
par(mar=c(6.5, 5.1, 4.1, 2.1))
spacing <- matrix(c(1, 0), 2, ncol(all.means))
spacing[1,4] <- 2
out <- barplot(all.means, beside=TRUE, space=spacing, las=2, ylim=c(0, 0.06), ylab="Observed type I error rate", 
               cex.lab=1.4, cex.axis=1.1, cex.names=1.2)
segments(out, all.means, out, all.se)
segments(out+0.1, all.se, out-0.1, all.se)
segments(out[1,1]-0.5, 0.01, out[2,3]+0.5, 0.01, col="red", lwd=2, lty=2)
segments(out[1,4]-0.5, 0.05, out[2,6]+0.5, 0.05, col="red", lwd=2, lty=2)

legend(out[1,1]-0.5, 0.06, legend=c("Technical", "Biological"), fill=grey.colors(2), cex=1.2)

dev.off()

############################################ 
# End.

