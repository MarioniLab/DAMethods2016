# This checks that edgeR runs properly on CyToF data. 
# We do so by pooling all cells from all real samples into a single matrix, and then sample from that matrix to constitute each sample.
# We repeat this, once with random sample and again as a Dirichlet process (where each cell is a table) to mimic correlations.

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
    nDA <- 0.1
    set.seed(12321)

    for (it in seq_len(50)) {
        current.exprs <- resampleCells(rawdata, setting=2L)

        # Adding a large DA subpopulation to both groups.
        for (i in seq_len(nsamples)) {
            if (groupings[i]==1L) {
                loc <- 1
            } else {
                loc <- 0
            }
            to.sample <- nrow(current.exprs[[i]])
            extras <- matrix(loc, round(to.sample*nDA), ncol(current.exprs[[i]]))
            current.exprs[[i]] <- rbind(current.exprs[[i]], extras)
        }

        ##########################################################################
        # Running CITRUS
        ##########################################################################

        # Writing to FCS files.
        out.files <- list()
        for (f in seq_along(current.exprs)) {
            curexp <- current.exprs[[f]]
            p <- AnnotatedDataFrame(data.frame(name=colnames(curexp), 
                                               desc=colnames(curexp),
                                               range=apply(curexp, 2, function(x) { diff(range(x)) }),
                                               minRange=apply(curexp, 2, min),
                                               maxRange=apply(curexp, 2, max)))
            ff <- flowFrame(curexp, p)
            fname <- file.path(odir, paste0(f, ".fcs"))
            suppressWarnings(write.FCS(file=fname, ff))
            out.files[[f]] <- fname
        }

        # Running CITRUS.
        cit.out <- file.path(odir, "citrusOutput")
        all.files <- data.frame(default=basename(unlist(out.files)))
        all.markers <- colnames(current.exprs[[1]])
        min.pct <- 0.05

        results <- citrus.full(
            fileList=all.files,
            labels=groupings,
            clusteringColumns=all.markers,
            dataDirectory=odir,
            outputDirectory=cit.out,
            family="classification",
            modelTypes="sam",
            nFolds=1,
            
            fileSampleSize=1000,
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
            selected <- results$citrus.combinedFCSSet$data[cur.clust,]
            if (any(rowSums(selected==1)==length(all.markers))) {
                found.up <- append(found.up, x)
            } 
            if (any(rowSums(selected==0)==length(all.markers))) {
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
            if (sum((central - 1)^2) <= threshold) {
                pass <- TRUE
                strict.up <- TRUE 
            }
            if (sum((central - 0)^2) <= threshold) {
                pass <- TRUE
                strict.down <- TRUE
            }
            if (!pass) { failed <- failed + 1 }
        }

        strict.fdr <- failed/length(sig.clusters)

        ##########################################################################
        # Saving results 
        ##########################################################################

        write.table(data.frame(RelaxedFDR=relaxed.fdr, RelaxedUp=relaxed.up, RelaxedDown=relaxed.down,
                               StrictFDR=strict.fdr, StrictUp=strict.up, StrictDown=strict.down),
                    file=ofile, append=existing, row.names=FALSE, col.names=!existing, sep="\t", quote=FALSE)
        existing <- TRUE
    }
}

# Cleaning up.
unlink(odir, recursive=TRUE)
