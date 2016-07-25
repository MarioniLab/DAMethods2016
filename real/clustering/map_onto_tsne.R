# This maps the cluster centers onto the t-SNE plot.

ori.dir <- "../analysis"
ref <- readRDS(file.path(ori.dir, "Cytobank_43324_4FI_res.rds"))
tsne.coords <- read.table(file.path(ori.dir, "Cytobank_43324_4FI_coords.txt"))
all.hypers <- t(ref$coords[as.character(tsne.coords[,1]),])

collected.x <- collected.y <- collected.col <- collected.pch <- list()
it <- 1L
max.logFC <- 1
zero.col <- 0.8
collected.failed <- list()

for (k in c(50, 100, 200)) { 
    if (k==50) { 
        pch <- 22
    } else if (k==100) { 
        pch <- 23
    } else if (k==200) { 
        pch <- 24
    }
    incoming <- read.table(sprintf("da_clusters_%i.tsv", k), check.names=FALSE)
    m <- match(colnames(ref$coords), colnames(incoming))
    is.sig <- which(incoming$FDR <= 0.05)

    current.x <- current.y <- list()
    failed <- 0L
    it <- 1L
    max.dist <- 0
    for (i in is.sig) {
        distances <- sqrt(colSums((all.hypers - as.numeric(incoming[i,m]))^2))
        closest <- which.min(distances)
        if (distances[closest] <= sqrt(ncol(ref$coords))*0.5) { 
            current.x[[it]] <- tsne.coords[closest,2]
            current.y[[it]] <- tsne.coords[closest,3]
            it <- it + 1L
        } else {
            failed <- failed + 1L
        }
        max.dist <- max(max.dist, distances[closest])
    }

    kc <- as.character(k)
    collected.x[[kc]] <- unlist(current.x)
    collected.y[[kc]] <- unlist(current.y)
    all.logFC <- pmin(pmax(incoming$logFC[is.sig], -max.logFC), max.logFC)
    collected.col[[kc]] <- plotrix::color.scale(all.logFC, c(0, zero.col, 1), c(0, zero.col, 0), 
                                                c(1, zero.col, 0), xrange = c(-max.logFC, max.logFC))
    collected.pch[[kc]] <- rep(pch, length(is.sig))
    collected.failed[[kc]] <- c(failed, max.dist)
}

collected.failed

png("Cytobank_43324_4FI_clusters.png", width=7, height=6, units="in", res=300)
layout(cbind(1,2), c(10, 2))
plot(tsne.coords[,2], tsne.coords[,3], xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, pch=1, col="grey50")
points(tsne.coords[,2], tsne.coords[,3], xlab="t-SNE1", ylab="t-SNE2", cex.axis=1.2, cex.lab=1.4, cex.main=1.4, pch=16, col="white", cex=0.9)
points(unlist(collected.x), unlist(collected.y), bg=unlist(collected.col), pch=unlist(collected.pch))

# Adding a legend.
par(mar=c(5.1, 0, 4.1, 0.1))
plot(0,0,type="n", axes=FALSE, xlab="", ylab="")
legend("topleft", pch=22:24, legend=paste0("k = ", c(50, 100, 200)), cex=1.2)

# Adding a colourbar
length.out <- 100
values <- seq(from = -max.logFC, to = max.logFC, length.out = length.out)
out <- plotrix::color.scale(values, c(0, zero.col, 1), c(0, zero.col, 0), c(1, zero.col, 0), xrange = c(-max.logFC, max.logFC))
start.loc <- seq(-0.5, 0.5, length.out=length(out)) - 0.3
rect(-0.5, start.loc, 0.5, start.loc+diff(start.loc)[1], col=out, border=NA)
text(0, min(start.loc), pos=1, -max.logFC, cex=1.2)
text(0, max(start.loc),  pos=3, max.logFC, cex=1.2)
text(-0.9, mean(start.loc), pos=1, srt=90, "Log-FC per day", cex=1)
dev.off()
