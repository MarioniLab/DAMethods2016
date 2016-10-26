####################################################################################
# This checks the effect of the hypersphere radius on the log-fold change.

generateSphere <- function(npts, radius, centre) { 
    ndim <- length(centre)
    coords <- matrix(rnorm(npts*ndim), ncol=ndim)
    mult <- radius * runif(npts)^(1/ndim) / sqrt(rowSums(coords^2))
    t(t(coords * mult) + centre)
}

pdf("plot_radius.pdf")
par(mar=c(5.1, 5.1, 4.1, 2.1))
for (scenario in 1:4) {
    ncells <- 10000
    alt.cells <- 10000
    alt.rad <- 0.5
    nmarkers <- 30
    samples <- c(1,2)
    loc <- 0.5
    xlim <- c(0, 4)

    if (scenario==1L) {
        ;
    } else if (scenario==2L) {
        alt.cells <- 5000
    } else if (scenario==3L) {
        loc <- 0.8
        xlim <- c(0, 5)
    } else if (scenario==4L) {
        alt.rad <- 0.2
    }

    set.seed(10000)
    nonda <- generateSphere(ncells, 0.5*sqrt(nmarkers), numeric(nmarkers))
    da <- generateSphere(alt.cells, alt.rad*sqrt(nmarkers), rep(loc, nmarkers))
    combined <- rbind(nonda, da)
    sample.id <- c(sample(2, ncells, replace=TRUE), rep(1, alt.cells))
    
    colnames(combined) <- paste0("X", seq_len(nmarkers))
    by.sample <- list()
    true.diff <- list()
    for (x in seq_along(samples)) { 
        selected <- which(sample.id==x)
        by.sample[[x]] <- combined[selected,,drop=FALSE]
        suppressWarnings(true.diff[[x]] <- min(which(selected > ncells)))
    }
    names(by.sample) <- paste0("Y", seq_along(by.sample))
    true.diff <- unlist(true.diff)
   
    require(cydar)
    cd <- prepareCellData(by.sample)
    for (tol in c(0.5, 0.6, 0.7, 0.8)) { 
        out <- countCells(cd, downsample=10, tol=tol, BPPARAM=SerialParam())
        
        distances <- sqrt(rowSums((out$coordinates - loc)^2))
        logfc <- log2((out$counts[,1]+1)/(out$counts[,2]+1))
        columnar <- as.integer(sub("c", "", rownames(out$counts)))
        s.origin <- attributes(cd)$sample.id[columnar]
        c.origin <- attributes(cd)$cell.id[columnar]
        suppressWarnings(col <- ifelse(true.diff[s.origin+1] <= c.origin, "red", "grey"))
        plot(distances, logfc, col=col, pch=16, main=sprintf("Radius = %.2f, Scenario = %i", tol, scenario),
             xlim=xlim, ylim=c(-2, 10), xlab="Distance from DA subpopulation centre", ylab=expression("Log"[2]~"fold change in abundance"),
             cex.lab=1.6)

        threshold <- alt.rad * sqrt(nmarkers)
        abline(v=threshold, col="blue", lwd=2, lty=2) 
        abline(h=1, col="black", lwd=2, lty=2) 
    }
}
dev.off()

####################################################################################
# Adding a plot describing what the simulation is about.

pdf("radius_setup.pdf")
plot(0.5, 0, type="n", axes=FALSE, xlab="", ylab="", xlim=c(-1,1), ylim=c(-1,1))
library(plotrix)
first.col <- rgb(0.6, 0.6, 0.6, 0.2)
nda.centre <- c(-0.4, 0)
draw.circle(nda.centre[1], nda.centre[2], 0.5, border=NA, col=first.col)
text(nda.centre[1], nda.centre[2]-0.55, "Non-DA subpopulation", pos=1, cex=0.7)

second.col <- rgb(1, 0, 0, 0.2)
da.centre <- c(0.4, 0)
draw.circle(da.centre[1], da.centre[2], 0.5, border=NA, col=second.col)
text(da.centre[1], da.centre[2]-0.55, "DA subpopulation", pos=1, cex=0.7)

set.seed(100)
cur.pts <- generateSphere(100, 0.5, nda.centre)
points(cur.pts[,1], cur.pts[,2], col="grey80", pch=16)
cur.pts <- generateSphere(100, 0.5, da.centre)
points(cur.pts[,1], cur.pts[,2], col=rgb(1, 0.5, 0.5), pch=16)

# Adding the centres.
points(nda.centre[1], nda.centre[2], pch=4, cex=1.5) 
points(da.centre[1], da.centre[2], pch=4, cex=1.5) 
arrows(nda.centre[1] + 0.02, nda.centre[2], da.centre[1] - 0.02, da.centre[2], code=3, length=0.05)
text((nda.centre[1]+da.centre[1])/2, (nda.centre[2]+da.centre[2])/2, "Distance between subpopulations", pos=1, cex=0.7)

# Adding the radius.
arrows(da.centre[1]+0.02, da.centre[2], da.centre[1] + 0.5*cos(0),  da.centre[2] + 0.5*sin(0), code=3, length=0.05)
text(da.centre[1]+ 0.5/2, 0.0, "Subpopulation size", cex=0.7, pos=1)

# Introducing the hypersphere centre.
hyper.centre <- c(-0.05, 0.3)
points(hyper.centre[1], hyper.centre[2], pch=16, col="grey80")
points(hyper.centre[1], hyper.centre[2])
draw.circle(hyper.centre[1], hyper.centre[2], 0.2, border="black")
text(hyper.centre[1], hyper.centre[2]+0.2, "Hypersphere", pos=3, cex=0.7)
text(hyper.centre[1], hyper.centre[2], "Centre", pos=3, cex=0.4, offset=0.25)

hyper.loc <- c(-0.00, 0.26)
points(hyper.loc[1], hyper.loc[2], pch=16, col="dodgerblue")
arrows(hyper.centre[1], hyper.centre[2], hyper.loc[1], hyper.loc[2], length=0.05)
text(hyper.loc[1], hyper.loc[2], "Median\nlocation", pos=1, cex=0.4, offset=0.35)

arrows(hyper.centre[1] - 0.02, hyper.centre[2], hyper.centre[1] - 0.2 + 0.001, hyper.centre[2], length=0.05)
text(hyper.centre[1] - 0.1, hyper.centre[2], "Radius", pos=1, cex=0.4, offset=0.35)


dev.off()
