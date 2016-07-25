# This checks that edgeR runs properly on CyToF data. 
# We do so by pooling all cells from all real samples into a single matrix, and then sample from that matrix to constitute each sample.
# We repeat this, once with random sample and again as a Dirichlet process (where each cell is a table) to mimic correlations.

require(cyder)
require(edgeR)

ofile <- "results_FDR.txt"
existing <- FALSE

#####################################
# Assessment 

assessFDR <- function(coords, true.diff, boundary=0.5, plot=FALSE) {
    xbin.id <- ceiling(coords[,1]/boundary)
    ybin.id <- ceiling(coords[,2]/boundary)
    all.ids <- paste0(xbin.id, ".", ybin.id)
    all.tests <- table(all.ids)
    false.pos <- table(all.ids[!true.diff])

    if (plot) {
        xrange <- range(coords[,1])
        yrange <- range(coords[,2])
        xspan <- diff(xrange)
        yspan <- diff(yrange)
        
        newspan <- max(xspan, yspan) * 1.01 # A slight increase to avoid problems with numerical precision.
        ymean <- mean(yrange)
        yrange <- (yrange - ymean)/yspan * newspan + ymean
        xmean <- mean(xrange)
        xrange <- (xrange - xmean)/xspan * newspan + xmean

        plot(0, 0, xlim=xrange, ylim=yrange, type="n", xlab="PC1", ylab="PC2", cex.axis=1.2, cex.lab=1.4)
        all.colors <- rev(grey.colors(max(false.pos)+1))[-(1)]
        for (i in unique(all.ids)) {
            xo <- as.numeric(sub("\\..*", "", i)) * boundary
            yo <- as.numeric(sub(".*\\.", "", i)) * boundary
            if (!is.na(false.pos[i]) && all.tests[i]==false.pos[i]) { 
                col <- all.colors[false.pos[i]]
            } else {
                col <- "red"
            }
            rect(xo - boundary, yo - boundary, xo, yo, col=col, border=NA)
        } 
    }

    m <- match(names(false.pos), names(all.tests))
    return(sum(false.pos/all.tests[m])/length(all.tests))
}

############################################ 
# Setting up

for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    x <- readRDS(file.path("../refdata", paste0(dataset, "_raw.rds")))
    groupings <- rep(1:2, length.out=nsamples)
    nDA <- 0.1
    set.seed(12321)

    for (it in seq_len(50)) {
        current.exprs <- resampleCells(x, setting=2L)
        
        # Adding a large DA subpopulation to both groups.
        for (i in seq_len(nsamples)) {
            if (groupings[i]==1L) {
                loc <- 1
            } else {
                loc <- 0
            }
            extras <- matrix(loc, round(to.sample*nDA), ncol(pooled))
            current.exprs[[i]] <- rbind(current.exprs[[i]], extras)
        }
    
        # edgeR.
        cd <- prepareCellData(current.exprs)
        out <- countCells(cd, BPPARAM=MulticoreParam(2), downsample=10)
 
        # Testing for differential proportions
        y <- DGEList(out$counts, lib.size=out$total)
        y$genes <- out$coordinates
        keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$total))
        y <- y[keep,]
        
        design <- model.matrix(~factor(groupings))
        y <- estimateDisp(y, design, robust=TRUE)
        fit <- glmQLFit(y, design)
        res <- glmQLFTest(fit)

        # Figuring out which hyperspheres contain the DA spot(s).
        # This should be easy, as it's quite a clear distinction.
        is.DA <- rowMeans(y$counts) > 500 & abs(res$table$logFC) > 1

        # Controlling the FDR, spatially or naively.
        all.results <- list()
        for (con in c("naive", "spatial")) {
            if (con=="naive") {
                qval <- p.adjust(res$table$PValue, method="BH")
            } else {
                qval <- spatialFDR(y$genes, res$table$PValue)
            }
            is.sig <- qval <= 0.05
            coords <- prcomp(y$genes[is.sig,])$x[,1:2]

            for (width in c(0.5, 1, 2)) { 
                if (existing || con=="spatial" || width!=0.5) {
                    fdr <- assessFDR(coords, is.DA[is.sig], boundary=width)
                } else {
                    # Making an image of what the analysis might look like.
                    pdf("FDR_setup.pdf")
                    fdr <- assessFDR(coords, is.DA[is.sig], boundary=width, plot=TRUE)
                    dev.off()
                }
                all.results[[paste0(con, "_", width)]] <- fdr
            }
        }     
     
        write.table(file=ofile, data.frame(Dataset=dataset, rbind(unlist(all.results))),
                    quote=FALSE, sep="\t", col.names=!existing, append=existing, row.names=FALSE)
        existing <- TRUE
        gc()
    }
}

############################################ 
# Plotting

stuff <- read.table("results_FDR.txt", header=TRUE)
bydataset <- split(stuff[,-1], stuff$Dataset)

all.naive.means <- list()
all.naive.se <- list()
all.weight.means <- list()
all.weight.se <- list()
for (x in names(bydataset)) {
    current.d <- bydataset[[x]]
    
    is.naive <- grep("naive", colnames(current.d))
    is.weight <- grep("spatial", colnames(current.d))
    cur.naive.means <- colMeans(current.d[,is.naive])
    cur.naive.se <- sqrt(apply(current.d[,is.naive], 2, var)/nrow(current.d))
    cur.weight.means <- colMeans(current.d[,is.weight])
    cur.weight.se <- sqrt(apply(current.d[,is.weight], 2, var)/nrow(current.d))

    if (x=="Cytobank_43324_4FI") {
        transfect <- "Oct4-GFP"
    } else if (x=="Cytobank_43324_NN") {
        transfect <- "Nanog-Neo"
    } else {
        transfect <- "Nanog-GFP"
    }

    all.naive.means[[transfect]] <- cur.naive.means    
    all.naive.se[[transfect]] <- cur.naive.se
    all.weight.means[[transfect]] <- cur.weight.means    
    all.weight.se[[transfect]] <- cur.weight.se
}

all.naive.means <- do.call(cbind, all.naive.means)
all.naive.se <- do.call(cbind, all.naive.se)
all.weight.means <- do.call(cbind, all.weight.means)
all.weight.se <- do.call(cbind, all.weight.se)

all.means <- cbind(all.naive.means, all.weight.means)
all.se <- all.means + cbind(all.naive.se, all.weight.se)

pdf("plot_FDR.pdf", width=10, height=6)
par(mar=c(6.5, 5.1, 4.1, 2.1))
out <- barplot(all.means, beside=TRUE, ylab="Observed spatial FDR", ylim=c(0, 1), cex.axis=1.1, cex.lab=1.4, cex.name=1.2)

segments(out, all.means, out, all.se)
segments(out-0.1, all.se, out+0.1, all.se)
abline(h=0.05, col="red", lwd=2, lty=2)

legend("topright", legend=paste("Width of", sub(".*_", "", rownames(all.means))), fill=grey.colors(3), cex=1.2)

par(xpd=TRUE)
height <- -0.15
segments(out[1,1], height, out[3,3], height, lwd=2)
text(out[2,2], height, "Na\u00EFve", cex=1.4, pos=1)
segments(out[1,4], height, out[3,6], height, lwd=2)
text(out[2,5], height, "Weighted", cex=1.4, pos=1)

dev.off()

############################################ 
# End.
