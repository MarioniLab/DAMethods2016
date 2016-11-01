# This checks that citrus runs properly on CyToF data. 
# We do so by pooling all cells from all real samples into a single matrix, and then sample from that matrix to constitute each sample.
# We also add some DA subpopulations that should be fairly obvious (comprising 10% of the toal sample).

require(citrus)
require(flowCore)
require(Biobase)
source("functions.R")

odir <- "temp_citrus"
dir.create(odir)
ofile <- "results_citrus.txt"
existing <- FALSE

############################################ 
# Setting up

for (dataset in c("Cytobank_43324_4FI", "Cytobank_43324_NG", "Cytobank_43324_NN")) {
    rawdata <- readRDS(file.path("../refdata", paste0(dataset, "_raw.rds")))
    nsamples <- length(attributes(rawdata)$samples)
    groupings <- rep(1:2, length.out=nsamples)
    set.seed(12321)

    for (it in seq_len(50)) {
        current.exprs <- resampleCells(rawdata, setting=2L)

        # Adding a large DA subpopulation to both groups.
        current.exprs <- addPointDifference(current.exprs, which(groupings==1L), loc=1, prop.DA=0.1)
        current.exprs <- addPointDifference(current.exprs, which(groupings==2L), loc=0, prop.DA=0.1)

        ##########################################################################
        # Running CITRUS
        ##########################################################################

        # Writing to FCS files.
        out.files <- dumpToFile(odir, current.exprs)

        # Running CITRUS.
        cit.out <- file.path(odir, "citrusOutput")
        dir.create(cit.out, showWarning=FALSE)
        all.files <- data.frame(default=basename(unlist(out.files)))
        all.markers <- colnames(current.exprs[[1]])
        min.pct <- 0.05 # Default

        results <- citrus.full(
            fileList=all.files,
            labels=groupings,
            clusteringColumns=all.markers,
            dataDirectory=odir,
            outputDirectory=cit.out,
            family="classification",
            modelTypes="sam",
            nFolds=1,
            
            fileSampleSize=1000, # Default
            featureType="abundances",
            minimumClusterSizePercent=min.pct,
            transformColumns=NULL, # Already transformed, no scaling.
            transformCofactor=NULL,
            scaleColumns=NULL
        )

        # Get clusters at a nominal FDR of 5%.
        sig.clusters <- as.integer(results$conditionRegressionResults$default$sam$differentialFeatures$fdr_0.05[["clusters"]])

        ##########################################################################
        # Relaxed evaluation:
        ##########################################################################

        # Identifying all clusters containing cells with all-zero or all-unity intensity.
        found.up <- found.down <- NULL
        for (x in seq_along(results$citrus.foldClustering$allClustering$clusterMembership)) {
            cur.clust <- results$citrus.foldClustering$allClustering$clusterMembership[[x]]
            selected <- results$citrus.combinedFCSSet$data[cur.clust,all.markers,drop=FALSE]
            if (any(rowSums(abs(selected-1)<1e-6)==length(all.markers))) {
                found.up <- append(found.up, x)
            } 
            if (any(rowSums(abs(selected-0)<1e-6)==length(all.markers))) {
                found.down <- append(found.down, x)
            }
        }

        up <- sig.clusters %in% found.up
        down <- sig.clusters %in% found.down
        relaxed.fdr <- 1-sum(up|down)/length(up)
        relaxed.up <- any(up)
        relaxed.down <- any(down)
    
        ##########################################################################
        # Strict evaluation:
        ##########################################################################

        threshold <- 0.5*sqrt(length(all.markers)) 
        strict.up <- strict.down <- FALSE
        failed <- 0

        for (x in sig.clusters) {
            cur.clust <- results$citrus.foldClustering$allClustering$clusterMembership[[x]]
            selected <- results$citrus.combinedFCSSet$data[cur.clust,]
            central <- apply(selected[,all.markers], 2, median)

            pass <- FALSE
            if (sqrt(sum((central - 1)^2)) <= threshold) {
                pass <- TRUE
                strict.up <- TRUE 
            }
            if (sqrt(sum((central - 0)^2)) <= threshold) {
                pass <- TRUE
                strict.down <- TRUE
            }
            if (!pass) { failed <- failed + 1 }
        }

        strict.fdr <- failed/length(sig.clusters)

        ##########################################################################
        # Saving results 
        ##########################################################################

        write.table(data.frame(Dataset=dataset, RelaxedFDR=relaxed.fdr, RelaxedUp=relaxed.up, RelaxedDown=relaxed.down,
                               StrictFDR=strict.fdr, StrictUp=strict.up, StrictDown=strict.down),
                    file=ofile, append=existing, row.names=FALSE, col.names=!existing, sep="\t", quote=FALSE)
        existing <- TRUE
    }
}

# Cleaning up.
unlink(odir, recursive=TRUE)

############################################ 

stuff <- read.table("results_citrus.txt", header=TRUE)
bydataset <- split(stuff[,-1], stuff$Dataset)

all.relaxed.means <- list()
all.relaxed.se <- list()
all.strict.means <- list()
all.strict.se <- list()
for (x in names(bydataset)) {
    current.d <- bydataset[[x]]
    
    is.relaxed <- "RelaxedFDR"==colnames(current.d)
    is.strict <- "StrictFDR"==colnames(current.d)
    cur.relaxed.means <- mean(current.d[,is.relaxed])
    cur.relaxed.se <- sqrt(var(current.d[,is.relaxed])/nrow(current.d))
    cur.strict.means <- mean(current.d[,is.strict])
    cur.strict.se <- sqrt(var(current.d[,is.strict])/nrow(current.d))
    
    if (x=="Cytobank_43324_4FI") {
        transfect <- "Oct4-GFP"
    } else if (x=="Cytobank_43324_NN") {
        transfect <- "Nanog-Neo"
    } else {
        transfect <- "Nanog-GFP"
    }
    
    all.relaxed.means[[transfect]] <- cur.relaxed.means    
    all.relaxed.se[[transfect]] <- cur.relaxed.se
    all.strict.means[[transfect]] <- cur.strict.means    
    all.strict.se[[transfect]] <- cur.strict.se
}

all.means <- cbind(unlist(all.relaxed.means), unlist(all.strict.means))
all.se <- cbind(unlist(all.relaxed.se), unlist(all.strict.se))

pdf("plot_citrus.pdf", width=7, height=6)
par(mar=c(5.1, 4.1, 2.1, 8.1))
out <- barplot(all.means, beside=TRUE, ylim=c(0, 0.3), ylab="Observed FDR", cex.axis=1.2, cex.lab=1.4)
segments(out, all.means, out, all.means+all.se)
segments(out-0.1, all.means+all.se, out+0.1)

mtext(c("Relaxed", "Strict"), at=colMeans(out), side=1, line=1, cex=1.4)
abline(h=0.05, col="red", lwd=2, lty=2)

par(xpd=TRUE)
legend(max(out)+0.5, 0.3, fill=grey.colors(3), rownames(all.means), cex=1.2)
dev.off()
