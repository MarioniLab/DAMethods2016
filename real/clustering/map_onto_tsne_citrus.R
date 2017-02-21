# This maps the cluster centers onto the t-SNE plot.

ori.dir <- "../analysis"
ref <- readRDS(file.path(ori.dir, "Cytobank_43324_4FI_res.rds"))
tsne.coords <- read.table(file.path(ori.dir, "Cytobank_43324_4FI_coords.txt"))
all.hypers <- t(ref$coords[ref$results$FDR <= 0.05,])

incoming <- read.table("da_citrus.tsv", check.names=FALSE, header=TRUE)
m <- match(colnames(ref$coords), colnames(incoming))

current.x <- current.y <- list()
failed <- 0L
it <- 1L
max.dist <- 0
for (i in seq_len(nrow(incoming))) {
    distances <- sqrt(colSums((all.hypers - as.numeric(incoming[i,m]))^2))
    closest <- which.min(distances)
    if (distances[closest] <= sqrt(ncol(ref$coords))*0.5) { 
        current.x[[it]] <- tsne.coords[closest,1]
        current.y[[it]] <- tsne.coords[closest,2]
        it <- it + 1L
    } else {
        failed <- failed + 1L
    }
    max.dist <- max(max.dist, distances[closest])
}
failed

png("Cytobank_43324_4FI_citrus.png", width=7, height=6, units="in", res=300)
plot(tsne.coords[,1], tsne.coords[,2], xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, pch=1, col="grey50")
points(tsne.coords[,1], tsne.coords[,2], xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, pch=16, col="white", cex=0.9)
points(unlist(current.x), unlist(current.y), pch=16)
dev.off()
