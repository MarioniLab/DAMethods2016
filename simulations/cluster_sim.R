#####################################
# Simulating a data set with two DA subpopulations merged into a larger subpopulation. 

set.seed(500)
samples <- c(1,1,2,2)
design <- model.matrix(~factor(samples))

ncells <- 20000
nda <- 50

nmarkers <- 30
threshold <- 0.5*sqrt(nmarkers) 

down.pos <- 1.8
up.pos <- 1.2

###################################
# Running the simulation.

detected.clust <- detected.hyper <- list()
for (it in 1:20) { 
    combined <- rbind(matrix(rnorm(ncells*nmarkers, 1.5, 0.6), ncol=nmarkers),
                      matrix(rnorm(nda*nmarkers, down.pos, 0.3), ncol=nmarkers),
                      matrix(rnorm(nda*nmarkers, up.pos, 0.3), ncol=nmarkers))
    sample.id <- c(sample(length(samples), ncells, replace=TRUE), 
                   sample(which(samples==1), nda, replace=TRUE), 
                   sample(which(samples==2), nda, replace=TRUE))
    true.direction <- rep(1:3, c(ncells, nda, nda))

    #### With clustering, at verying cut depths.
    adist <- dist(combined)
    mytree <- hclust(adist)
    collected <- list()
    for (k in c(20, 50, 100)) { 
        clusters <- cutree(mytree, k=k)

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
        collected[[as.character(k)]] <- current.collected
    }
    detected.clust[[it]] <- unlist(collected)
    
    #### With hypersphere counts.

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
    fit <- glmQLFit(y, design)
    res <- glmQLFTest(fit)
    qval <- spatialFDR(y$genes, res$table$PValue)
   
    # Checking if the position is within range of the center of the altered population.
    da.hypersphere <- qval <= 0.05
    centers <- t(y$genes[da.hypersphere,]) 
    detected.hyper[[it]] <- c(any(sqrt(colSums((centers - down.pos)^2)) <= threshold & res$table$logFC[da.hypersphere] < 0), 
                              any(sqrt(colSums((centers - up.pos)^2)) <= threshold & res$table$logFC[da.hypersphere] > 0))
                           
    # Plotting the results.
    if (it==1L) { 
        pdf("cluster_setup.pdf")
        pca <- prcomp(combined)
        plot(pca$x[,1], pca$x[,2], pch=16, col="grey80", xlab="PC1", ylab="PC2")
        points(pca$x[true.direction==2L,1], pca$x[true.direction==2L,2], pch=16, col="blue")
        points(pca$x[true.direction==3L,1], pca$x[true.direction==3L,2], pch=16, col="red")
        dev.off()
    }
}

collected.results <- cbind(do.call(rbind, detected.clust), do.call(rbind, detected.hyper))
colnames(collected.results) <- paste(rep(c("D", "U"), 4), rep(c("20", "50", "100", "Hyper"), each=2), sep=".")
write.table(file="results_cluster.txt", collected.results, sep="\t", quote=FALSE, row.names=FALSE)

###############################################
# Making a barplot.

all.sums <- colSums(collected.results)/nrow(collected.results)
every.second <- seq(from=1, to=ncol(collected.results), by=2)
all.stats <- rbind(all.sums[every.second], all.sums[every.second+1]) 
colnames(all.stats) <- c("k = 20", "k = 50", "k = 100", "Hyperspheres")

pdf("plot_cluster.pdf")
barplot(all.stats, beside=TRUE, ylab="Detection frequency", ylim=c(0, 1), col=c("blue", "red"),
        cex.axis=up.pos, cex.lab=1.4, cex.names=1.4)
legend(1, 1, fill=c("blue", "red"), legend=c("First", "Second"), cex=up.pos)
dev.off()
