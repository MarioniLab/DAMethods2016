# This checks that the spatial FDR is being properly controlled.
# We use the same simulations as in edgeR_check.R, with biological variability.
# (Note that this results in some liberalness at low p-values as the counts aren't exactly NB-distributed.)

require(cydar)
require(edgeR)
source("functions.R")
ofile <- "results_FDR.txt"
existing <- FALSE
plotgen <- TRUE 

############################################ 
# Setting up

for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    curdata <- readRDS(file.path("../refdata", paste0(dataset, ".rds")))
    nsamples <- ncol(curdata)
    groupings <- rep(1:2, length.out=nsamples)
    set.seed(12321)

    for (it in seq_len(50)) {
        current.exprs <- resampleCells(curdata, setting=2L)

        # Adding a large DA subpopulation to both groups.
        current.exprs <- addPointDifference(current.exprs, which(groupings==1L), loc=1, prop.DA=0.1)
        current.exprs <- addPointDifference(current.exprs, which(groupings==2L), loc=0, prop.DA=0.1)
    
        # edgeR.
        cd <- prepareCellData(current.exprs)
        out <- countCells(cd, BPPARAM=SerialParam(), downsample=10)
 
        # Testing for differential proportions
        y <- DGEList(assay(out), lib.size=out$totals)
        keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$totals))
        out <- out[keep,]
        y <- y[keep,]
        
        design <- model.matrix(~factor(groupings))
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design, robust=TRUE)
        res <- glmQLFTest(fit)

        # Figuring out which hyperspheres contain the DA spot(s).
        tolerance <- 0.5*sqrt(ncol(current.exprs[[1]]))
        index <- rowData(out)$center.cell
        true.centres <- t(cellIntensities(out)[,index])
        is.up <- sqrt(rowSums((true.centres - 1)^2)) < tolerance  
        is.down <- sqrt(rowSums((true.centres - 0)^2)) < tolerance 
        is.DA <- is.up | is.down

        # Controlling the FDR, spatially or naively.
        all.results <- list()
        hypercoords <- intensities(out)
        for (con in c("naive", "spatial")) {
            if (con=="naive") {
                qval <- p.adjust(res$table$PValue, method="BH")
            } else {
                qval <- spatialFDR(hypercoords, res$table$PValue)
            }
            is.sig <- qval <= 0.05
            all.results[[paste0(con, "_up")]] <- sum(is.sig & is.up)
            all.results[[paste0(con, "_down")]] <- sum(is.sig & is.down)

            # Assessing the FDR using partitions of varying size.
            for (width in c(0.2, 0.4, 0.6, 0.8, 1)) { 
                partitions <- apply(hypercoords[is.sig,,drop=FALSE], 1, function(x) { paste(floor(x/width), collapse=".") })
                all.tests <- table(partitions)
                false.pos <- table(partitions[!is.DA[is.sig]])
                m <- match(names(false.pos), names(all.tests))
                fdr <- sum(false.pos/all.tests[m])/length(all.tests)
                all.results[[paste0(con, "_", width)]] <- fdr
            }

            if (plotgen) {
                # Plotting an example.
                plotgen <- FALSE
                coords <- prcomp(hypercoords[is.sig,])$x[,1:2]
                boundary <- 0.5
                partitions <- apply(coords, 1, function(x) { paste(ceiling(x/boundary), collapse=".") })
                all.tests <- table(partitions)
                false.pos <- table(partitions[!is.DA[is.sig]])

                xbin.id <- ceiling(coords[,1]/boundary)
                ybin.id <- ceiling(coords[,2]/boundary)
                all.ids <- paste0(xbin.id, ".", ybin.id)

                xrange <- range(coords[,1])
                yrange <- range(coords[,2])
                xspan <- diff(xrange)
                yspan <- diff(yrange)
                
                newspan <- max(xspan, yspan) * 1.01 # A slight increase to avoid problems with numerical precision.
                ymean <- mean(yrange)
                yrange <- (yrange - ymean)/yspan * newspan + ymean
                xmean <- mean(xrange)
                xrange <- (xrange - xmean)/xspan * newspan + xmean
                
                pdf("FDR_setup.pdf")
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
                dev.off()
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
stuff <- stuff[,!grepl("up|down", colnames(stuff))]
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

legend("topright", legend=paste("Width of", sub(".*_", "", rownames(all.means))), fill=grey.colors(nrow(all.means)), cex=1.2)

par(xpd=TRUE)
height <- -0.15
segments(out[1,1], height, out[5,3], height, lwd=2)
text(out[3,2], height, "Na\u00EFve", cex=1.4, pos=1)
segments(out[1,4], height, out[5,6], height, lwd=2)
text(out[3,5], height, "Weighted", cex=1.4, pos=1)

dev.off()

############################################ 
# End.
