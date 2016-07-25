library(ncdfFlow)
dataset <- "Cytobank_43324_4FI"
new.fcs <- list.files(file.path("../../refdata", dataset, "fcs"), full=TRUE)
x <- read.ncdfFlowSet(new.fcs)
toignore <- c(match(c("Cell_length", "Time", "beadDist", "barcode"), colnames(x)), grep("^BC-[0-9]", colnames(x)), grep("DNA", colnames(x)))
x <- x[,-toignore]

timings <- as.integer(sub(".*_([0-9]+).fcs", "\\1", basename(new.fcs)))
design <- model.matrix(~splines::ns(timings, 3))

# Downsampling; taking 1000 from each sample prior to clustering, a la CITRUS.
set.seed(10)
collected <- list()
samples <- list()
for (i in seq_along(sampleNames(x))) {
    current <- exprs(x[[i]])
    collected[[i]] <- current[sample(nrow(current), 1000),]
    samples[[i]] <- rep(i, nrow(collected[[i]]))
}
collected <- do.call(rbind, collected)
samples <- unlist(samples)

# Clusting using hierarchical methods.
mytree <- hclust(dist(collected))
for (k in c(50, 100, 200)) {
    clusters <- cutree(mytree, k=k)
    
    # Assembling counts per cluster.
    all.counts <- all.coords <- list()
    by.cluster <- split(seq_along(clusters), clusters)
    for (i in seq_along(by.cluster)) {
        chosen <- by.cluster[[i]]
        all.counts[[i]] <- tabulate(samples[chosen], nbins=nrow(design))
        all.coords[[i]] <- apply(collected[chosen,,drop=FALSE], 2, median)
    }
    all.counts <- do.call(rbind, all.counts)
    all.coords <- do.call(rbind, all.coords)
    
    # Testing each cluster for DA.
    require(edgeR)
    y <- DGEList(all.counts, genes=all.coords)
    keep <- aveLogCPM(y) >= aveLogCPM(5, mean(colSums(all.counts)))
    y <- y[keep,]
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design)
    res <- glmQLFTest(fit, coef=2:ncol(design))

    # Computing fold changes.
    redesign <- model.matrix(~timings)
    refit <- glmFit(y, redesign, dispersion=1)
    reres <- glmLRT(refit)
    
    # Storing the result.
    ofile <- sprintf("da_clusters_%i.tsv", k)
    output <- cbind(y$genes, logFC=reres$table$logFC, PValue=res$table$PValue, 
                    FDR=p.adjust(res$table$PValue, method="BH"))
    write.table(output, file=ofile, sep="\t", quote=FALSE)
}

