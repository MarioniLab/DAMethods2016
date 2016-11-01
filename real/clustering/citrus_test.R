library(ncdfFlow)
dataset <- "Cytobank_43324_4FI"
ref.dir <- file.path("../../refdata", dataset, "fcs")
new.fcs <- list.files(ref.dir, full=TRUE)
x <- read.ncdfFlowSet(new.fcs)
toignore <- c(match(c("Cell_length", "Time", "beadDist", "barcode"), colnames(x)), grep("^BC-[0-9]", colnames(x)), grep("DNA", colnames(x)))

markers <- colnames(x)[-toignore]
timings <- as.integer(sub(".*_([0-9]+).fcs", "\\1", basename(new.fcs)))

# Running CITRUS.
library(citrus)
cit.out <- "citrusOutput"
dir.create(cit.out, showWarning=FALSE)
all.files <- data.frame(default=basename(new.fcs))
min.pct <- 0.05

set.seed(1000)
results <- citrus.full(
                       fileList=all.files,
                       labels=timings,
                       clusteringColumns=markers,
                       dataDirectory=ref.dir,
                       outputDirectory=cit.out,
                       family="continuous",
                       modelTypes="sam",
                       nFolds=1,
                       
                       fileSampleSize=1000,
                       featureType="abundances",
                       minimumClusterSizePercent=min.pct,
                       transformColumns=NULL, # Already transformed, no scaling.
                       transformCofactor=NULL,
                       scaleColumns=NULL
                       )

# Saving the results.
is.sig <- results$conditionRegressionResults$default$sam$differentialFeatures$fdr_0.05$clusters
collected <- list()
for (clust in is.sig) {
    cur.clust <- results$citrus.foldClustering$allClustering$clusterMembership[[clust]]
    cur.exprs <- results$citrus.combinedFCSSet$data[cur.clust,]
    center <- apply(cur.exprs[,markers], 2, median)
    collected[[length(collected)+1]] <- center
}
collected <- do.call(rbind, collected)
write.table(file="da_citrus.tsv", collected, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

unlink(cit.out, recursive=TRUE)
