require(cydar)
dir.create("pics")

for (il in c("IL-10")) { #, "IL-3")) {
    true.name <- sub("-", "", il)
    stuff <- readRDS(paste0(true.name, "_res.rds"))
    qvals <- stuff$results$FDR
    coords <- stuff$coords    

    # Visualizing with t-SNE.
    require(Rtsne)
    is.sig <- qvals <= 0.05
    sig.coords <- coords[is.sig,]
    set.seed(100)
    change.up <- (stuff$results$logFC > 0)[is.sig]
    tsne.out.1 <- Rtsne(sig.coords[change.up,], perplexity=10)
    change.down <- (stuff$results$logFC < 0)[is.sig]
    tsne.out.2 <- Rtsne(sig.coords[change.down,], perplexity=10)

    tsne.out <- list(Y=matrix(0, ncol=2, nrow=nrow(sig.coords)))
    tsne.out$Y[change.up,] <- tsne.out.1$Y 
    tsne.out$Y[change.up,1] <- tsne.out$Y[change.up,1] - max(tsne.out$Y[change.up,1]) - 10
    tsne.out$Y[change.down,] <- tsne.out.2$Y 
    tsne.out$Y[change.down,1] <- tsne.out$Y[change.down,1] - min(tsne.out$Y[change.down,1]) + 10
    write.table(tsne.out$Y, file=paste0(true.name, "_coords.txt"), sep="\t", quote=FALSE, col.names=FALSE)

    for (plotmode in c(TRUE, FALSE)) { 
        # Plotting log-FCs
        if (plotmode) { 
            png(file.path("pics", paste0(true.name, "_logFC.png")), width=6.2, height=6, units="in", res=300)
            plottype <- "p"
        } else {
            pdf(file.path("pics", paste0(true.name, "_logFC.pdf")), width=6.2, height=6)
            plottype <- "n"
        }

        layout(cbind(1,2), widths=c(10, 1))
        par(mar=c(5.1,4.1,4.1,1.1))
        out <- plotCellLogFC(tsne.out$Y[,1], tsne.out$Y[,2], stuff$results$logFC[is.sig], 
                      main=paste("Effect of", il, "treatment"), xlab="t-SNE1", ylab="t-SNE2", 
                      cex.axis=1.2, cex.lab=1.4, cex.main=1.4, max.logFC=3, cex=1, type=plottype)
        par(mar=c(0,0,0,0))
        plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
        start.loc <- seq(-0.5, 0.5, length.out=length(out))
        rect(-0.5, start.loc, 0.5, start.loc+diff(start.loc)[1], col=out, border=NA)
        text(0, -0.5, pos=1, names(out)[1], cex=1.2)
        text(0, 0.5,  pos=3, names(out)[length(out)], cex=1.2)
        text(-0.9, 0, pos=1, srt=90, "Log-FC", cex=1)
        dev.off()

        # Plotting intensities
        reranges <- stuff$ranges
        if (plotmode) { 
            png(file.path("pics", paste0(true.name, "_markers.png")), width=20, height=20, units="in", res=300)
        } else {
            pdf(file.path("pics", paste0(true.name, "_markers.pdf")), width=20, height=20)
        }
        
        lmat <- cbind(matrix(seq_len(6*6), ncol=6, nrow=6), 37)
        layout(lmat, widths=c(rep(1, 6), 0.2))
        for (i in order(colnames(sig.coords))) {
            par(mar=c(2.1, 2.1, 2.1, 2.1))
            out <- plotCellIntensity(tsne.out$Y[,1], tsne.out$Y[,2], sig.coords[,i], 
                                     irange=reranges[,i], main=colnames(sig.coords)[i], 
                                     xlab="t-SNE1", ylab="t-SNE2", cex=1, type=plottype)
        }
        for (i in (1+ncol(sig.coords)):36) {
            # Using up the remaining panels before we get to the colour bar.
            plot.new()
        }
        par(mar=c(0,0,0,0))
        plot(0,0, type="n", axes=FALSE, ylab="", xlab="", ylim=c(-1, 1), xlim=c(-1, 0.5))
        start.loc <- seq(-0.5, 0.5, length.out=length(out))
        interval <- diff(start.loc)[1]
        rect(-0.5, start.loc, 0.5, start.loc+interval, col=out, border=NA)
        text(0, -0.5, pos=1, "Low", cex=1.5)
        text(0, 0.5+interval,  pos=3, "High", cex=1.5)
        text(-0.9, 0, pos=1, srt=90, "Marker intensity", cex=1.5)
        
        dev.off()
    }
}
