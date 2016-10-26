#####################################
# Simulating a data set with two DA subpopulations merged into a larger subpopulation. 

source("functions.R")
odir <- "temp_cluster"
dir.create(odir)

set.seed(500)
samples <- c(1,1,2,2)
design <- model.matrix(~factor(samples))

ncells <- 20000
nda <- 20

nmarkers <- 30
threshold <- 0.5*sqrt(nmarkers) 

down.pos <- 2
up.pos <- 1

###################################
# Running the simulation.

detected.clust <- detected.hyper <- detected.citrus <- list()
for (it in 1:50) { 
    combined <- rbind(matrix(rnorm(ncells*nmarkers, 1.5, 0.6), ncol=nmarkers),
                      matrix(rnorm(nda*nmarkers, down.pos, 0.3), ncol=nmarkers),
                      matrix(rnorm(nda*nmarkers, up.pos, 0.3), ncol=nmarkers))
    sample.id <- c(sample(length(samples), ncells, replace=TRUE), 
                   sample(which(samples==1), nda, replace=TRUE), 
                   sample(which(samples==2), nda, replace=TRUE))
    true.direction <- rep(1:3, c(ncells, nda, nda))

    ##########################################################################
    #### With clustering, at verying cut depths.
    ##########################################################################

    adist <- dist(combined)
    mytree <- hclust(adist)
    collected <- list()
    for (strat in c("hierachical", "kmeans")) { 
    for (k in c(5, 50, 500)) { 
        if (strat=="hierarchical") {
            clusters <- cutree(mytree, k=k)
        } else {
            clusters <- kmeans(combined, centers=k, iter.max=100)$cluster
        }

        # Assembling counts per cluster.
        all.counts <- list()
        by.cluster <- split(sample.id, clusters)
        for (i in by.cluster) {
            all.counts <- c(all.counts, list(tabulate(i, nbins=length(samples))))
        }
        all.counts <- do.call(rbind, all.counts)
        rownames(all.counts) <- seq_len(k)
        
        # Testing each cluster for DA.
        require(edgeR)
        y <- DGEList(all.counts)
        keep <- aveLogCPM(y) >= aveLogCPM(5, mean(colSums(all.counts)))
        y <- y[keep,]
        y <- estimateDisp(y, design)
        fit <- glmQLFit(y, design, robust=TRUE)
        res <- glmQLFTest(fit)
        
        # Checking if the cluster center is close to the center of the altered populations.
        current.collected <- logical(2) 
        for (p in 1:2) { 
            ref <- ifelse(p==1, down.pos, up.pos)
            is.sig <- which(p.adjust(res$table$PValue, method="BH")<=0.05 & (res$table$logFC < 0)==(p==1))

            affirmed <- FALSE
            for (is in is.sig) {
                current.cluster <- as.integer(rownames(res)[is]) # Getting the cluster ID.
                clustered <- clusters==current.cluster
                center <- apply(combined[clustered,,drop=FALSE], 2, median) # Getting the cluster centre.
                if (sqrt(sum((center - ref)^2)) <= threshold) {
                    affirmed <- TRUE
                }
            }
            if (affirmed) { current.collected[p] <- TRUE }
        }
        collected[[paste0(strat, k, "_down")]] <- current.collected[1]
        collected[[paste0(strat, k, "_up")]] <- current.collected[2]
    }
    }
    detected.clust[[it]] <- unlist(collected)
   
    ##########################################################################
    #### With hypersphere counts.
    ##########################################################################

    # Splitting by sample in preparation for counting.
    colnames(combined) <- paste0("X", seq_len(nmarkers))
    by.sample <- list()
    for (x in seq_along(samples)) { 
        by.sample[[x]] <- combined[sample.id==x,,drop=FALSE]
    }
    names(by.sample) <- paste0("Y", seq_along(by.sample))
    
    require(cydar)
    cd <- prepareCellData(by.sample)
    out <- countCells(cd, downsample=10, tol=0.5, BPPARAM=SerialParam())
    
    # Testing for DA.
    y <- DGEList(out$counts, lib.size=out$totals, genes=out$coordinates)
    keep <- aveLogCPM(y) >= aveLogCPM(5, mean(out$total))
    y <- y[keep,]
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit)
    qval <- spatialFDR(y$genes, res$table$PValue)
   
    # Checking if the position is within range of the center of the altered population.
    da.hypersphere <- qval <= 0.05
    centers <- t(y$genes[da.hypersphere,]) 
    detected.hyper[[it]] <- c(cydar_down=any(sqrt(colSums((centers - down.pos)^2)) <= threshold & res$table$logFC[da.hypersphere] < 0), 
                              cydar_up=any(sqrt(colSums((centers - up.pos)^2)) <= threshold & res$table$logFC[da.hypersphere] > 0))
    
    ##########################################################################
    #### With CITRUS.
    ##########################################################################

    library(citrus)
    
    # Writing to FCS files.
    out.files <- dumpToFile(odir, by.sample)
    
    # Running CITRUS.
    cit.out <- file.path(odir, "citrusOutput")
    dir.create(cit.out, showWarning=FALSE)
    all.files <- data.frame(default=basename(unlist(out.files)))
    all.markers <- colnames(by.sample[[1]])
    min.pct <- 0.05
    
    results <- citrus.full(
                           fileList=all.files,
                           labels=samples,
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

    strict.up <- strict.down <- FALSE
    failed <- 0
    for (x in sig.clusters) {
        cur.clust <- results$citrus.foldClustering$allClustering$clusterMembership[[x]]
        selected <- results$citrus.combinedFCSSet$data[cur.clust,]
        central <- apply(selected[,all.markers], 2, median)
        
        pass <- FALSE
        if (sum((central - up.pos)^2) <= threshold) {
            pass <- TRUE
            strict.up <- TRUE 
        }
        if (sum((central - down.pos)^2) <= threshold) {
            pass <- TRUE
            strict.down <- TRUE
        }
        if (!pass) { failed <- failed + 1 }
    }
    detected.citrus[[it]] <- c(CITRUS_down=strict.down, CITRUS_up=strict.up)

    ##########################################################################
    #### Other stuff.
    ##########################################################################

    # Plotting the setup.
    if (it==1L) { 
        pdf("cluster_setup.pdf")
        pca <- prcomp(combined)
        plot(pca$x[,1], pca$x[,2], pch=16, col="grey80", xlab="PC1", ylab="PC2")
        points(pca$x[true.direction==2L,1], pca$x[true.direction==2L,2], pch=16, col="blue")
        points(pca$x[true.direction==3L,1], pca$x[true.direction==3L,2], pch=16, col="red")
        dev.off()
    }
}

collected.results <- data.frame(do.call(rbind, detected.clust), do.call(rbind, detected.hyper), do.call(rbind, detected.citrus))
write.table(file="results_cluster.txt", collected.results, sep="\t", quote=FALSE, row.names=FALSE)
unlink(odir, recursive=TRUE)

###############################################
# Making a barplot.

all.sums <- colMeans(collected.results)
all.sums <- all.sums[!grepl("kmeans", names(all.sums))]
every.second <- seq(from=1, to=length(all.sums), by=2)
all.stats <- rbind(all.sums[every.second], all.sums[every.second+1]) 
colnames(all.stats) <- c("k = 5", "k = 50", "k = 500", " ", "CITRUS")

pdf("plot_cluster.pdf")
out <- barplot(all.stats, beside=TRUE, ylab="Detection frequency", ylim=c(0, 1), col=c("blue", "red"),
        cex.axis=1.2, cex.lab=1.4, cex.names=1.4)
legend(1, 1, fill=c("blue", "red"), legend=c("First", "Second"), cex=1.2)
mtext(at=mean(out[,4]), side=1, line=2, text="Hyper-\nspheres", cex=1.2)

# Adding binomial standard errors.
se <- sqrt(all.stats * (1-all.stats)/nrow(collected.results))
combo <- all.stats + se
segments(out, all.stats, out, combo)
segments(out - 0.1, combo, out + 0.1, combo)

par(xpd=TRUE)
yline <- -0.1
segments(min(out[,1]), yline, max(out[,3]), yline, lwd=2)
mtext(at=mean(out[,1:3]), side=1, line=3, text="Hierarachical clustering", cex=1.2)
dev.off()
